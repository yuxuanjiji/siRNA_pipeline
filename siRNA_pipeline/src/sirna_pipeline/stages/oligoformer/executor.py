#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""OligoFormer 效率预测执行器（内置，external 契约实现）。

契约：
    <python> executor.py --input <统一记录CSV> --output <输出CSV>
        [--repo <OligoFormer部分目录>] [--model <best_model.pth>]
        [--cds <CDS FASTA>] [--seed N]

对输入记录逐条：由 (guide_checked, cds_start, target_mRNA_19) 重建 19+19+19 mRNA
上下文（边界以 X 补足 57nt，与官方 infer.py 的 cRNA 构造一致），并按官方
infer 流程给出 efficiency 分：RNA-FM 嵌入(scripts/RNA-FM.sh) → loader 数据集 →
模型前向 → pred[:,1]×1.341（越大越好，官方口径）。

输出 CSV：variant_id, oligo_efficacy, oligo_available（每输入行一行）。

依赖与前提（由 stages/oligoformer/adapter.py 先做 status 探测，不满足则不调用本脚本）：
  - 运行环境 = OligoFormer 的 .venv（torch/pandas/numpy/sklearn）；
  - repo 下存在 scripts/RNA-FM.sh 与 RNA-FM 资产（RNA-FM/redevelop +
    pretrained/extract_embedding.yml），本机无资产时该步会失败——属于前置探测项。
本文件不改动 OligoFormer 仓库任何文件。
"""
from __future__ import annotations

import argparse
import csv
import os
import random
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

RNA_BASES = {"A", "C", "G", "U", "X"}


def _rc(seq: str) -> str:
    return seq.translate(str.maketrans("AUCGT", "UAGCA"))[::-1]


def _norm(seq: str) -> str:
    return "".join(seq.split()).upper().replace("T", "U")


def _resolve_bash() -> str:
    """优先 Git Bash（Windows 下 System32\\bash.exe 是 WSL 启动器，不可用于本脚本）。"""
    for cand in (r"C:\Program Files\Git\bin\bash.exe",
                 r"C:\Program Files (x86)\Git\bin\bash.exe"):
        if Path(cand).exists():
            return cand
    return shutil.which("bash") or "bash"


def _load_cds(fasta: str) -> str:
    name, seqs = None, []
    with open(fasta, encoding="utf-8") as fh:
        for ln in fh:
            ln = ln.strip()
            if not ln:
                continue
            if ln.startswith(">"):
                continue
            seqs.append(ln)
    return _norm("".join(seqs))


def _mrna57(cds: str, w19: str | None, cds_start_1based: int) -> str:
    """重建 57nt = 19 上游 + 靶窗19 + 19 下游；边界以 X 补足（对齐官方 cRNA）。"""
    n = len(cds)
    i0 = max(0, cds_start_1based - 1)
    win = w19 if w19 and len(w19) == 19 else cds[i0:i0 + 19]
    win = _norm(win)
    win = (win + "X" * 19)[:19]                # 靶窗不足 19 以 X 补足
    # 上游
    l = cds[max(0, i0 - 19):i0]
    l = "X" * (19 - len(l)) + l if len(l) < 19 else l[-19:]
    # 下游
    r = cds[i0 + 19:i0 + 38]
    r = r + "X" * (19 - len(r)) if len(r) < 19 else r[:19]
    return l + win + r


def _read_input(path: Path) -> list[dict]:
    with open(path, encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh))


def _find_ckpt(repo_root: Path, given: str | None) -> Path | None:
    """定位 RNA-FM 权重：优先显式路径，其次仓库常见位置。"""
    if given:
        p = Path(str(given))
        return p if p.exists() else None
    for rel in ("RNA-FM/redevelop/pretrained/RNA-FM_pretrained.pth",
                "RNA-FM/pretrained/RNA-FM_pretrained.pth"):
        p = repo_root / rel
        if p.exists():
            return p
    return None


def _md5seq(seq: str) -> str:
    import hashlib
    return hashlib.md5(seq.encode("utf-8")).hexdigest()[:16]


def _cache_root(repo_root: Path) -> Path:
    """跨 run 复用的嵌入缓存（按序列 md5 命名，之后重跑几乎瞬时）。"""
    return repo_root / "data" / "RNAFM_cache"


def _embed_repo_fm(uniq: dict, cache: Path, repo_root: Path, ckpt: Path,
                   seed: int, batch_size: int = 32, log=None,
                   device: str = "cpu") -> None:
    """用 RNA-FM 仓库自带 `fm` 包（纯 PyTorch，无需 fairseq）产出 (L,640) 嵌入。

    uniq: {"siRNA": {md5: seq}, "mRNA": {...}}（已按序列去重）
    写入 cache/<kind>/representations/<md5>.npy；已存在的跳过（缓存命中）。
    device: 'cpu' 或 'cuda'（GPU 机器上把模型与张量搬到显卡，见 README 环境清单）。
    """
    import argparse as _argparse
    import numpy as np
    import torch

    log = log or (lambda *a: None)
    sys.path.insert(0, str(repo_root / "RNA-FM"))
    torch.serialization.add_safe_globals([_argparse.Namespace])
    torch.manual_seed(seed)
    random.seed(seed)
    import fm  # type: ignore

    dev = torch.device(device)
    t0 = __import__("time").time()
    model, alphabet = fm.pretrained.rna_fm_t12(str(ckpt))
    model.to(dev)
    model.eval()
    conv = alphabet.get_batch_converter()
    log("RNA-FM 权重加载完成 %.1fs（%s，device=%s）"
        % (__import__("time").time() - t0, ckpt.name, dev))

    for kind, mapping in uniq.items():
        rep = cache / kind / "representations"
        rep.mkdir(parents=True, exist_ok=True)
        todo = [(n, s) for n, s in mapping.items() if not (rep / (n + ".npy")).exists()]
        log("RNA-FM[%s]: 唯一序列 %d 条，需计算 %d 条（其余命中缓存）"
            % (kind, len(mapping), len(todo)))
        t0 = __import__("time").time()
        for i in range(0, len(todo), batch_size):
            chunk = todo[i:i + batch_size]
            _, _, toks = conv([(n, s) for n, s in chunk])
            with torch.no_grad():
                out = model(toks.to(dev), repr_layers=[12])
            reps = out["representations"][12]
            for j, (name, seq) in enumerate(chunk):
                np.save(rep / (name + ".npy"),
                        reps[j, 1:len(seq) + 1].cpu().numpy().astype("float32"))
        if todo:
            log("RNA-FM[%s]: %d 条完成，用时 %.1fs（%.3fs/条）"
                % (kind, len(todo), __import__("time").time() - t0,
                   (__import__("time").time() - t0) / max(1, len(todo))))


def _materialize(token_dir: Path, cache: Path, uniq: dict, log=None) -> None:
    """在本次 run 的 token 目录下生成 fasta 与 representations（缓存硬链接/拷贝）。"""
    log = log or (lambda *a: None)
    for kind, mapping in uniq.items():
        rep = token_dir / kind / "representations"
        rep.mkdir(parents=True, exist_ok=True)
        fa = token_dir / ("siRNA.fa" if kind == "siRNA" else "mRNA.fa")
        with open(fa, "w", encoding="utf-8") as fh:
            for name, seq in mapping.items():
                fh.write(">%s\n%s\n" % (name, seq))
        n_link = n_copy = 0
        for name in mapping:
            src = cache / kind / "representations" / (name + ".npy")
            dst = rep / (name + ".npy")
            if dst.exists():
                continue
            try:
                os.link(src, dst)
                n_link += 1
            except OSError:
                shutil.copyfile(src, dst)
                n_copy += 1
        log("materialize[%s]: %d 条（link %d / copy %d）"
            % (kind, len(mapping), n_link, n_copy))


def _build_dataset(rows: list[dict], cds: str, work: Path, repo_root: Path,
                   rnafm_python: str | None = None, fm_mode: str = "auto",
                   rnafm_model: str | None = None, fm_batch: int = 32,
                   seed: int = 42, log=None, device: str = "cpu"):
    """生成 (mRNA/siRNA df, token 名, 数据目录) —— 按官方 loader 相对路径布局。

    嵌入路线：
      * repo-fm（默认优先）：RNA-FM 仓库自带 `fm` 包，纯 PyTorch CPU，无需 fairseq；
      * bash：官方 scripts/RNA-FM.sh（fairseq 环境，如 py3.8/3.10 conda env）。
    两者都按**序列去重**，并复用 `data/RNAFM_cache`（重跑近瞬时）。
    """
    import pandas as pd
    log = log or (lambda *a: None)
    token = "exec_%d" % os.getpid()
    recs = []
    for r in rows:
        g = _norm(str(r.get("guide_checked", "")))
        w19 = _norm(str(r.get("target_mRNA_19", "")))
        try:
            start = int(float(str(r.get("cds_start", ""))))
        except (TypeError, ValueError):
            start = 1
        if len(g) != 19:
            g = (g + "X" * 19)[:19]
        m57 = _mrna57(cds, w19, start)
        recs.append({"siRNA": g, "mRNA": m57, "cds_start": start})

    df = pd.DataFrame({"siRNA": [x["siRNA"] for x in recs],
                       "mRNA": [x["mRNA"] for x in recs]})
    # 新版 pandas 会把全字符串列推断为 StringDtype，官方 calculate_td 写入 float 会
    # TypeError；转 object 保持与官方旧 pandas 行为一致
    df = df.astype(object)

    # td 特征：直接复用官方 infer.calculate_td（同表同列序，保证与训练一致）
    sys.path.insert(0, str(repo_root / "scripts"))
    import infer as infer_mod  # type: ignore
    df = infer_mod.calculate_td(df)
    df = df.reset_index(drop=True)
    for col in ("siRNA", "mRNA"):
        df[col] = df[col].astype(str)

    # 按序列去重（同一窗口 16 条共享 mRNA；guide 亦可能重复）
    uniq = {"siRNA": {}, "mRNA": {}}
    for seq in df["siRNA"].tolist():
        uniq["siRNA"].setdefault(_md5seq(seq), seq)
    for seq in df["mRNA"].tolist():
        uniq["mRNA"].setdefault(_md5seq(seq), seq)

    data_dir = repo_root / "data" / "infer" / token
    if data_dir.exists():
        shutil.rmtree(data_dir)
    data_dir.mkdir(parents=True, exist_ok=True)

    mode = fm_mode or "auto"
    ckpt = _find_ckpt(repo_root, rnafm_model)
    if mode == "auto":
        mode = "repo-fm" if (ckpt and (repo_root / "RNA-FM" / "fm").exists()) else "bash"
    log("FM 路线=%s（RNA-FM 权重=%s）" % (mode, ckpt.name if ckpt else "未找到"))

    if mode == "repo-fm":
        if ckpt is None:
            raise FileNotFoundError("未找到 RNA-FM 权重（可用 --rnafm-model 指定）")
        cache = _cache_root(repo_root)
        _embed_repo_fm(uniq, cache, repo_root, ckpt, seed=seed,
                       batch_size=fm_batch, log=log, device=device)
        _materialize(data_dir, cache, uniq, log=log)
    else:
        # 官方 bash 脚本路线（fairseq 环境）；先落去重 fasta，再调用脚本
        cache = data_dir
        _materialize_fasta_only(data_dir, uniq)
        script = repo_root / "scripts" / "RNA-FM.sh"
        if not script.exists():
            raise FileNotFoundError(f"缺少 RNA-FM.sh：{script}")
        env = dict(os.environ)
        py = Path(str(rnafm_python or sys.executable)).resolve()
        if not py.exists():
            raise FileNotFoundError(f"RNA-FM 解释器不存在：{py}")
        env["PATH"] = str(py.parent) + os.pathsep + env.get("PATH", "")
        run = subprocess.run([_resolve_bash(), str(script), str(data_dir)],
                             cwd=str(repo_root), capture_output=True, text=True,
                             encoding="utf-8", errors="replace", timeout=7200,
                             env=env)
        if run.returncode != 0:
            raise RuntimeError("RNA-FM.sh 失败（请检查 repo/RNA-FM 资产）："
                               + (run.stdout + run.stderr)[-1200:])
    return df, token, data_dir


def _materialize_fasta_only(token_dir: Path, uniq: dict) -> None:
    """bash 路线：只落去重 fasta（脚本会自行产出 representations）。"""
    for kind, mapping in uniq.items():
        fa = token_dir / ("siRNA.fa" if kind == "siRNA" else "mRNA.fa")
        with open(fa, "w", encoding="utf-8") as fh:
            for name, seq in mapping.items():
                fh.write(">%s\n%s\n" % (name, seq))


def _infer(df, token: str, repo_root: Path, model_path: Path, seed: int,
           device: str = "cpu"):
    """官方前向口径：loader 数据集 → Oligo 模型 → pred[:,1]×1.341。"""
    import numpy as np  # noqa: F401
    import torch
    from torch.utils.data import DataLoader
    torch.manual_seed(seed)
    random.seed(seed)
    dev = torch.device(device)

    sys.path.insert(0, str(repo_root / "scripts"))
    from loader import data_process_loader_infer  # type: ignore
    import model as model_mod  # type: ignore

    params = {"batch_size": 1, "shuffle": False, "num_workers": 0,
              "drop_last": False}
    # 官方 loader 用相对路径 './data/infer/<token>/...' 读 fasta 与嵌入 → 需 cwd=仓库根
    cwd = os.getcwd()
    os.chdir(str(repo_root))
    try:
        ds = DataLoader(data_process_loader_infer(df.index.values, df, token),
                        **params)

        m = model_mod.Oligo()      # 类默认参数即训练口径（vocab 26 / dim 128 / head 8 / L1）
        state = torch.load(str(model_path), map_location=dev)
        if isinstance(state, dict) and "state_dict" in state:
            state = state["state_dict"]
        m.load_state_dict(state)
        m.to(dev)
        m.eval()

        outs = []
        with torch.no_grad():
            for i, data in enumerate(ds):
                if len(data) != 5:  # 防御：loader 版本兼容
                    continue
                out = m(*[d.to(dev) for d in data])
                outs.append(out[0][:, 1].item())
    finally:
        os.chdir(cwd)
    return [v * 1.341 for v in outs]       # 官方 ×1.341 对齐 PITA 分数量纲


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--repo", required=True, help="OligoFormer部分 目录")
    ap.add_argument("--model", required=True, help="best_model.pth 路径")
    ap.add_argument("--cds", required=True, help="靶基因 CDS FASTA（重建 57nt）")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--work-dir", default=None, help="中间目录（默认 repo/data/infer）")
    ap.add_argument("--rnafm-python", default=None,
                    help="跑 RNA-FM 的独立解释器（如 py3.8/3.10 conda 环境），"
                         "留空则用当前 PATH 中的 python")
    ap.add_argument("--fm-mode", default="auto", choices=("auto", "repo-fm", "bash"),
                    help="嵌入路线：auto=有仓库 fm 包+权重则用 repo-fm，否则 bash 官方脚本")
    ap.add_argument("--rnafm-model", default=None,
                    help="RNA-FM 权重 .pth（默认自动定位 RNA-FM/**/RNA-FM_pretrained.pth）")
    ap.add_argument("--fm-batch", type=int, default=32, help="RNA-FM 批大小")
    ap.add_argument("--device", default="auto",
                    help="auto|cpu|cuda[：N] —— GPU 机器上装 CUDA 版 torch 后自动用卡")
    args = ap.parse_args(argv)

    def _log(*a):
        print("[executor]", *a, flush=True)

    import torch as _torch
    device = args.device
    if str(device).lower() == "auto":
        device = "cuda" if _torch.cuda.is_available() else "cpu"
    if str(device).startswith("cuda") and not _torch.cuda.is_available():
        _log("请求 device=%s 但 CUDA 不可用 → 回退 CPU" % device)
        device = "cpu"
    _log("device = %s（torch %s）" % (device, _torch.__version__))
    _log("CPU threads = %d" % _torch.get_num_threads())

    repo_root = Path(args.repo).resolve()
    cds = _load_cds(args.cds)
    rows = _read_input(Path(args.input))
    if not rows:
        raise SystemExit("无输入行")
    # 路径绝对化：_infer() 内部会 chdir 到 repo 根以适配官方 loader 的相对路径，
    # 若 --model/--cds 为相对路径会在 chdir 后失效（FileNotFoundError）。
    model_path = Path(args.model).resolve()
    work_arg = Path(args.work_dir).resolve() if args.work_dir else None

    tmp_root = work_arg if work_arg else repo_root
    # 沙箱/受限环境下 tempfile 的 0o700 属性操作会被拒，故自建目录并忽略错误清理
    work_dir = tmp_root / ("exec_oligo_%d" % os.getpid())
    shutil.rmtree(work_dir, ignore_errors=True)
    work_dir.mkdir(parents=True, exist_ok=True)
    df, token, data_dir = (None, None, None)
    try:
        df, token, data_dir = _build_dataset(
            rows, cds, work_dir, repo_root,
            rnafm_python=args.rnafm_python, fm_mode=args.fm_mode,
            rnafm_model=args.rnafm_model, fm_batch=args.fm_batch,
            seed=args.seed, log=_log, device=device)
        scores = _infer(df, token, repo_root, model_path, args.seed,
                        device=device)
    finally:
        if data_dir is not None:
            shutil.rmtree(data_dir, ignore_errors=True)
        shutil.rmtree(work_dir, ignore_errors=True)
    if len(scores) != len(rows):
        raise SystemExit("推断行数与输入不一致：%d != %d" % (len(scores), len(rows)))

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["variant_id", "oligo_efficacy", "oligo_available"])
        for r, s in zip(rows, scores):
            w.writerow([r.get("variant_id", ""), "%.10g" % s, "1"])
    print("executor done: rows=%d" % len(rows), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
