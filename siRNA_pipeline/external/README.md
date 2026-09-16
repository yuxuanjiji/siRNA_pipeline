# external/ —— 外部工具与数据（不随本仓库分发，按需获取）

## 1. ViennaRNA（结构检测与 siRNADiscovery 特征真值）

结构检测 `stages/structure/legacy/structure_detector.py` 自动探测顺序：
`import RNA`(Python 绑定) → `RNAfold/RNAcofold/RNAplfold`(PATH 或 `VIENNARNA_BIN`) →
否则纯 Python 最近邻近似（`structure_reliable=0`，不静默）。

安装（三种任选）：
1. **Windows 安装包（本机已采用，验证通过）**：`ViennaRNA Package` 2.7.2 64bit →
   默认装在 `C:\Program Files (x86)\ViennaRNA Package`（含 `RNAfold.exe/RNAcofold.exe/RNAplfold.exe`，
   无需 Python 绑定）；随后设 `VIENNARNA_BIN` 指向该目录或加入 PATH。
2. conda（含 python 绑定）：`conda install -c bioconda viennarna`
3. 源码编译（`项目搭建/结构检测/ViennaRNA-2.7.2/` 为源码副本，需 configure/make/make install）。
运行前确认 `import RNA` 可用，或把可执行文件目录加入 `VIENNARNA_BIN`（CLI 层 `level=1` 即可，
`structure_reliable=1`）。

## 2. 脱靶 PITA + TargetScan（任务9）

原方案是 OligoFormer 仓库内 bash/perl 脚本（`参考模型/OligoFormer-main/off-target/`），
需要 bash + perl + ViennaRNA，属于强外部环境。
本工程以 **外部执行器契约** 接入（`stages/offtarget/adapter.py`）：

```
<cmd> --input <统一记录CSV> --output <输出CSV>
输出 CSV 必须含 variant_id；可含 offtarget_pita_score/pita_filter/
offtarget_targetscan_score/offtarget_targetscan_filter/offtarget_pass 子集
```

就绪后启用：`configs/pipeline.yaml` 中 `offtarget.enabled: true, external_cmd: <cmd>`。
未就绪保持默认关闭（管道 graceful skip，不伪造预测）。

## 3. 参考模型（效率预测/双模型一致性）

- `参考模型/OligoFormer-main`：与 `项目搭建/Oligoformer部分/OligoFormer-main` 为同源副本（已 MD5 验证）；
  本工程以其为引用基准，**不复制**；**运行副本在 `工程根/OligoFormer部分/`**（model/scripts/.venv 齐全）。
- `参考模型/siRNADiscovery-2`：GNN（tensorflow+stellargraph+RNA-FM/RPISeq）。
  ⚠️ **该资产目前不在本仓库内，且团队已放弃"双模型一致性检验"计划（2026-09-10 决定，任务18 取消）**；
  主管道只用 OligoFormer 作效率辅助分。`task14` 的 `Simone.csv`（Sciabola 2013）仍来自该仓库，可作纯数据基准使用。
- OligoFormer 已以**效率辅助分**接入主管道（`stages/oligoformer/`，configs 默认 enabled=true）：
  内置 executor 复用其 model/loader/infer 口径输出 `oligo_efficacy`；缺模型权重或
  **RNA-FM 资产**（`<repo>/RNA-FM` + `pretrained/extract_embedding.yml`，需另外获取）时
  自动 graceful skip（rank α=1），取舍见 `docs/design/oligoformer_keep_delete.md`。

## 4. 数据集单源

`数据集/`（SFRP1-mRNA.txt、mismatch 校准集、Hu/Mix/Taka/Simone 效率基准）只读引用
（`configs/paths.yaml`），不在本工程内复制。

## 5. BLAST 近全长同源脱靶层（脱靶模块第三层）

- 源码：`BLAST/blast_offtarget.py`（**零修改、只读引用**，由 `stages/offtarget/blast_adapter.py`
  以子进程按 CLI 契约调用）；口径/判据/参数详见 `BLAST/README.txt`；
- 环境（enabled 时才需要）：运行脚本的 python 装有 pandas；NCBI BLAST+（blastn）；
  BLAST 库（makeblastdb 产物，如 GENCODE v46 转录本库，注意其数据许可）；
- 集成说明与配置样例：`docs/design/blast_offtarget_integration.md`。
