BLAST 脱靶检测模块（脱靶验证第三层：近全长序列同源，独立模块）
================================================================
定位：脱靶验证的 BLAST 层，与 offtarget_module\（off-target/PITA 种子区打分 +
     热力学可及性两层）平行，共同构成项目脱靶验证三件套。
     本层负责"近全长切割型"脱靶——用序列同源搜索回答：
     人类转录组里还有哪些基因含与 guide 反向互补的近全长位点。

来源
  算法    : Altschul et al. 1990, J Mol Biol 215:403-410（BLAST 原理）
  实现    : NCBI BLAST+ 2.14.1（Camacho et al. 2009, BMC Bioinformatics 10:421）
  参数    : NCBI BLAST Command Line Applications User Manual
           （blastn-short 任务、word_size 7——19nt 短 query 专用，
             Camacho 2009 论文不含此 task，勿引错）
  数据库  : GENCODE v46 蛋白编码全转录本（约 11 万 ENST，含 5'/3' UTR）

正式判据（tolerant 容错口径，2026-09-06 第三轮修正定稿）
  一条候选判"脱靶" <= 其任一 BLAST 命中同时满足：
    1) 命中基因非靶基因（默认 SFRP1）
    2) 反向互补方向命中（sstart > send：转录本含 rc(guide)，guide 可结合；
       plus 命中=与 guide 同向相同的同义片段，不能配对，不计）
    3) 正确配对数 (length - mismatch) >= 17
    4) 错配数 mismatch <= 1
  数值依据：判据为项目自定（对称性：候选自身对靶基因即 18/19+1 错配；
  机制见 Becker 2019 Mol Cell 75:741-755 错配位置高通量实测、
  Ui-Tei 2008 NAR 36:2136-2151、Haley & Zamore 2004 NSMB 11:599-606、
  Ameres 2007 Cell 130:101-112）；敏感性 = 三轮口径对比
  （strict16 剔 1,015 / strict17 剔 275 / +minus 剔 137 / tolerant 剔 823，
  后者抓到 SFRP2/SFRP5 同源家族 51 条，100% 全互补口径仅 2 条）。

环境要求（Windows 已验证；Linux 同理）
  Python 3.8+ 且装有 pandas；NCBI BLAST+ 可执行文件（blastn/makeblastdb）
  BLAST 库：E:\blastdb\gencode_v46_pc_simple（makeblastdb 产物，
  库文件不在本包内，路径通过 --db 传入）

输入
  CSV：须含 siRNA 列（19nt guide，U/T 均可，模块内部统一 U->T）
  或 FASTA：>名字 + 19nt 序列

运行（在 BLAST\ 目录下）
  python blast_offtarget.py --input <candidates.csv> --out-dir <输出目录> ^
      --prefix sfrp1
  全量复跑正式归档（3,605 条）：
  python blast_offtarget.py --input ..\..\08_结果档案\benchmark_results\sfrp1_pipeline\sfrp1_structure_pass.csv --out-dir <输出目录> --prefix sfrp1_full
  可选参数：
    --target-gene <名字>     排除的靶基因（默认 SFRP1，按子串匹配）
    --min-paired / --max-mismatch  判据阈值（默认 17 / 1）
    --require-no-gaps        额外要求 gapopen==0（正式归档口径未启用）
    --make-db <转录本.fa> --db-prefix <前缀>   从 FASTA 重建库
    --threads N              BLAST 线程数（默认 8）

输出（--out-dir 下）
  <prefix>_queries.fa     实际送检 FASTA（U->T 后）
  <prefix>_blast_raw.tsv  BLAST 原始 12 列输出
  <prefix>_hits_all.tsv   全部 >=16nt 命中（审计留痕，89,395 行口径）
  <prefix>_flagged.csv    判脱靶候选
  <prefix>_clean.csv      通过候选
  <prefix>_stats.json     参数 + 计数 + 脱靶基因 top

与 SFRP1 项目连接
  正式归档：08_结果档案\benchmark_results\sfrp1_pipeline\
    offtarget_hits_all_len16.tsv        3,605 条全量命中留痕
    offtarget_audit_full_tolerant.csv   每候选 tolerant_offtarget 标志
    sfrp1_FINAL_tolerant.csv            终榜（2,782 通过 -> 0.95 门槛 2,385）
  历史脚本（本包的来源，保留于 benchmark\ 供追溯）：
    filter_offtarget.py -> offtarget_audit.py -> offtarget_minus_correction.py
    -> offtarget_tolerant_final.py（口径三轮演进的完整记录见
    offtarget_audit_report.md）

复现性
  同 BLAST+ 2.14.1 + 同库 + 同参数 -> 命中逐行一致；
  30 条候选冒烟验证见 10_中间与临时\tmp_blast_smoke\（命中逐行、
  脱靶标志均与归档一致）。
