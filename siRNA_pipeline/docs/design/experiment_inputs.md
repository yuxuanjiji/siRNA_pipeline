# 任务14-17实验输入契约

统一入口是 `python scripts/run_experiments.py`。脚本只读取本地数据，不下载或补造论文数据；所有结果写入 `outputs/analysis/`。

## 任务14

`--input-dir` 指向 OligoFormer 目录，目录下需要 `data/{Hu,Taka,Mix,Simone}.csv` 与根目录对应的 `{name}_predictions.csv`。数据集必须有 `label`；预测文件必须有 `pred` 或 `pred_efficiency`。若还要计算纯热力主分与主分+辅助，数据集还必须包含 `dG_total`、`dG_seed`、`delta_deltaG_ends`。缺列时报告只会给纯 DL，并明确标记另外两种方法不可计算。

## 任务15

Holen 输入必须正好 50 行，Birmingham 输入必须正好 362 行；两者都需要 `label` 和预测列。预测可以直接在输入中，也可以用 `--prediction-file` 按行补入。按 guide 侧/mRNA 侧分层需要提供 `mismatch_side`（值建议为 `guide` 或 `mRNA`）。

## 任务16

消融输入至少需要 `label`、`dl_score`、`thermo_score`。脚本遍历 `$\\alpha_{thermo} \\in \\{0,0.2,0.4,0.6,0.8,1\\}$`，输出 Spearman/Pearson/R² 和推荐点；若存在 `mutation_site`，同时输出五突变位点分组摘要。

## 任务17

阳性对照 CSV 至少需要 `control_id`、`guide_seq`；`--rank` 指向已经跑完结构和排序的 `rank_final.csv`。输出包含是否匹配、`rules_pass`、`structure_pass`、`final_rank` 和排名分位。匹配不到的序列不会被当作“不通过”，而会单独列为 `unmatched`。

## 最终冻结

`freeze` 要求 `configs/paths.yaml` 配置的 SFRP1 FASTA 存在，并要求设置 `VIENNARNA_BIN`。满足后运行 `predict.py --run-name final_frozen`，由主管道发布 `outputs/results/results.csv`；缺任一条件只写 `freeze_status.json`，不会覆盖已有结果。
