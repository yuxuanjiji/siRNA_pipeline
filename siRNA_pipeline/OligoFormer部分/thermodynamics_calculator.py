# -*- coding: utf-8 -*-
"""
siRNA热力学特征计算工具
========================
功能：
  1. 双链整体ΔG（最近邻模型）
  2. 双链末端ΔΔG（5'端 vs 3'端前2nt ΔG差，用于链选择）
  3. seed区结合能（引导链位置2-8，共7nt的ΔG）
  4. 引导链MFE（最小自由能，需ViennaRNA/RNAfold；未安装时输出NaN）
  5. OligoFormer兼容的24维td特征

参考：
  - SantaLucia (1998) RNA最近邻热力学参数
  - OligoFormer infer.py 中的 calculate_td 函数

用法：
  python thermodynamics_calculator.py --input candidates.csv --output candidates_td.csv
  python thermodynamics_calculator.py --sirna GCGCUCAUCAUUGUGCUGC --mrna GCAGCACAAUGAUGAGUGC
"""

import argparse
import csv
import os
import sys
import subprocess
import numpy as np

# ============================================================
# RNA最近邻热力学参数（SantaLucia 1998，单位 kcal/mol）
# ============================================================
DeltaG = {
    'AA': -0.93, 'UU': -0.93, 'AU': -1.10, 'UA': -1.33,
    'CU': -2.08, 'AG': -2.08, 'CA': -2.11, 'UG': -2.11,
    'GU': -2.24, 'AC': -2.24, 'GA': -2.35, 'UC': -2.35,
    'CG': -2.36, 'GG': -3.26, 'CC': -3.26, 'GC': -3.42,
    'init': 4.09, 'endAU': 0.45, 'sym': 0.43
}

DeltaH = {
    'AA': -6.82, 'UU': -6.82, 'AU': -9.38, 'UA': -7.69,
    'CU': -10.48, 'AG': -10.48, 'CA': -10.44, 'UG': -10.44,
    'GU': -11.40, 'AC': -11.40, 'GA': -12.44, 'UC': -12.44,
    'CG': -10.64, 'GG': -13.39, 'CC': -13.39, 'GC': -14.88,
    'init': 3.61, 'endAU': 3.72, 'sym': 0
}


def antiRNA(RNA):
    """返回RNA的反向互补链"""
    comp = {'A': 'U', 'U': 'A', 'C': 'G', 'G': 'C', 'T': 'A', 'a': 'u', 'u': 'a', 'c': 'g', 'g': 'c'}
    return ''.join([comp.get(c, 'X') for c in RNA[::-1]])


def Calculate_DGH(seq):
    """计算一条RNA序列的双链ΔG和ΔH（假设与互补链配对）"""
    seq = seq.upper().replace('T', 'U')
    DG = DeltaG['init']
    # 末端AU惩罚
    ends = seq[0] + seq[-1]
    DG += (ends.count('A') + ends.count('U')) * DeltaG['endAU']
    # 对称性校正
    if antiRNA(seq) == seq:
        DG += DeltaG['sym']
    # 最近邻累加
    for i in range(len(seq) - 1):
        dimer = seq[i] + seq[i+1]
        if dimer in DeltaG:
            DG += DeltaG[dimer]

    DH = DeltaH['init']
    DH += (ends.count('A') + ends.count('U')) * DeltaH['endAU']
    if antiRNA(seq) == seq:
        DH += DeltaH['sym']
    for i in range(len(seq) - 1):
        dimer = seq[i] + seq[i+1]
        if dimer in DeltaH:
            DH += DeltaH[dimer]
    return DG, DH


def Calculate_end_diff(siRNA):
    """计算5'端与3'端前2nt的ΔG差值（链选择指标）
    正值表示5'端更不稳定（有利于引导链选择），负值表示3'端更不稳定
    """
    siRNA = siRNA.upper().replace('T', 'U')
    _5 = siRNA[:2]   # 5'端前2nt
    _3 = siRNA[-2:]  # 3'端后2nt
    count = 0
    if _5 in ['AC', 'AG', 'UC', 'UG']:
        count += 1
    elif _5 in ['GA', 'GU', 'CA', 'CU']:
        count -= 1
    if _3 in ['AC', 'AG', 'UC', 'UG']:
        count += 1
    elif _3 in ['GA', 'GU', 'CA', 'CU']:
        count -= 1
    dg5 = DeltaG.get(_5, 0)
    dg3 = DeltaG.get(_3, 0)
    return round(dg5 - dg3 + count * 0.45, 2)


def Calculate_seed_dG(siRNA):
    """计算seed区（引导链位置2-8，共7nt）的双链ΔG
    seed区是siRNA与mRNA结合的关键区域，结合能越强（ΔG越负）脱靶风险越高
    """
    siRNA = siRNA.upper().replace('T', 'U')
    seed = siRNA[1:8]  # 位置2-8（0-based index 1-7）
    dg, _ = Calculate_DGH(seed)
    return round(dg, 2)


def Calculate_total_dG(siRNA):
    """计算整条19nt引导链的双链ΔG"""
    dg, _ = Calculate_DGH(siRNA)
    return round(dg, 2)


_mfe_available = None  # 缓存：None=未检测, True=可用, False=不可用

def Calculate_MFE(siRNA):
    """计算引导链的MFE（最小自由能二级结构）
    需要ViennaRNA/RNAfold；未安装时返回NaN
    """
    global _mfe_available
    siRNA = siRNA.upper().replace('T', 'U')

    # 如果已经确认不可用，直接返回NaN
    if _mfe_available is False:
        return float('nan')

    try:
        # 尝试用Python的RNA模块（ViennaRNA Python绑定）
        import RNA
        _mfe_available = True
        fc = RNA.fold_compound(siRNA)
        (ss, mfe) = fc.mfe()
        return round(mfe, 2)
    except ImportError:
        pass

    try:
        # 尝试调用命令行RNAfold（缩短超时到3秒）
        result = subprocess.run(
            ['RNAfold', '--noPS', '-T', '37'],
            input=siRNA + '\n',
            capture_output=True, text=True, timeout=3
        )
        _mfe_available = True
        lines = result.stdout.strip().split('\n')
        if len(lines) >= 2:
            import re
            m = re.search(r'\(([-\d.]+)\)', lines[1])
            if m:
                return round(float(m.group(1)), 2)
    except (FileNotFoundError, subprocess.TimeoutExpired):
        _mfe_available = False
        return float('nan')

    return float('nan')


def calculate_oligoformer_td(siRNA):
    """计算OligoFormer兼容的24维td特征
    完全复现 infer.py 中的 calculate_td 函数
    """
    siRNA = siRNA.upper().replace('T', 'U')
    td = {}

    # 1. ends: 5'端vs3'端ΔG差异
    td['ends'] = Calculate_end_diff(siRNA)

    # 2. DG_1: 前2nt ΔG
    td['DG_1'] = DeltaG.get(siRNA[0:2], 0)

    # 3. DH_1: 前2nt ΔH
    td['DH_1'] = DeltaH.get(siRNA[0:2], 0)

    # 4. U_1: 第1位是否为U
    td['U_1'] = int(siRNA[0] == 'U')

    # 5. G_1: 第1位是否为G
    td['G_1'] = int(siRNA[0] == 'G')

    # 6. DH_all: 整条ΔH
    _, td['DH_all'] = Calculate_DGH(siRNA)

    # 7. U_all: U含量比例
    td['U_all'] = siRNA.count('U') / 19

    # 8. UU_1: 前2nt是否为UU
    td['UU_1'] = int(siRNA[0:2] == 'UU')

    # 9. G_all: G含量比例
    td['G_all'] = siRNA.count('G') / 19

    # 10. GG_1: 前2nt是否为GG
    td['GG_1'] = int(siRNA[0:2] == 'GG')

    # 11. GC_1: 前2nt是否为GC
    td['GC_1'] = int(siRNA[0:2] == 'GC')

    # 12. GG_all: GG二核苷酸频率
    dimers = [siRNA[j] + siRNA[j+1] for j in range(18)]
    td['GG_all'] = dimers.count('GG') / 18

    # 13. DG_2: 第2-3nt ΔG
    td['DG_2'] = DeltaG.get(siRNA[1:3], 0)

    # 14. UA_all: UA二核苷酸频率
    td['UA_all'] = dimers.count('UA') / 18

    # 15. U_2: 第2位是否为U
    td['U_2'] = int(siRNA[1] == 'U')

    # 16. C_1: 第1位是否为C
    td['C_1'] = int(siRNA[0] == 'C')

    # 17. CC_all: CC二核苷酸频率
    td['CC_all'] = dimers.count('CC') / 18

    # 18. DG_18: 最后2nt ΔG
    td['DG_18'] = DeltaG.get(siRNA[17:19], 0)

    # 19. CC_1: 前2nt是否为CC
    td['CC_1'] = int(siRNA[0:2] == 'CC')

    # 20. GC_all: GC二核苷酸频率
    td['GC_all'] = dimers.count('GC') / 18

    # 21. CG_1: 前2nt是否为CG
    td['CG_1'] = int(siRNA[0:2] == 'CG')

    # 22. DG_13: 第13-14nt ΔG
    td['DG_13'] = DeltaG.get(siRNA[12:14], 0)

    # 23. UU_all: UU二核苷酸频率
    td['UU_all'] = dimers.count('UU') / 18

    # 24. A_19: 最后1位是否为A
    td['A_19'] = int(siRNA[18] == 'A')

    return td


def compute_all_features(siRNA, mrna=None):
    """计算全部热力学特征，返回字典"""
    siRNA = siRNA.upper().replace('T', 'U')
    result = {
        'siRNA': siRNA,
        'dG_total': Calculate_total_dG(siRNA),
        'dG_seed': Calculate_seed_dG(siRNA),
        'delta_deltaG_ends': Calculate_end_diff(siRNA),
        'MFE_guide': Calculate_MFE(siRNA),
        'GC_content': round((siRNA.count('G') + siRNA.count('C')) / 19 * 100, 1),
    }
    # OligoFormer 24维td —— 必须严格按照 infer.py 中列创建顺序
    OLIGO_TD_ORDER = [
        'ends', 'DG_1', 'DH_1', 'U_1', 'G_1', 'DH_all', 'U_all', 'UU_1',
        'G_all', 'GG_1', 'GC_1', 'GG_all', 'DG_2', 'UA_all', 'U_2', 'C_1',
        'CC_all', 'DG_18', 'CC_1', 'GC_all', 'CG_1', 'DG_13', 'UU_all', 'A_19'
    ]
    td = calculate_oligoformer_td(siRNA)
    result['td'] = ','.join([str(td[k]) for k in OLIGO_TD_ORDER])
    result.update(td)
    return result


def process_csv(input_path, output_path, sirna_col='siRNA', mrna_col='mRNA'):
    """批量处理CSV文件，添加热力学特征列"""
    rows = []
    with open(input_path, encoding='utf-8-sig') as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames
        for row in reader:
            rows.append(row)

    print(f"Loaded {len(rows)} records from {input_path}")

    # 计算特征
    new_fields = ['dG_total', 'dG_seed', 'delta_deltaG_ends', 'MFE_guide', 'GC_content', 'td']
    out_fieldnames = list(fieldnames) + [f for f in new_fields if f not in fieldnames]

    out_rows = []
    for i, row in enumerate(rows):
        sirna = row.get(sirna_col, '').strip()
        if not sirna or len(sirna) != 19:
            print(f"  Warning: row {i} siRNA length != 19, skipping td")
            row['dG_total'] = ''
            row['dG_seed'] = ''
            row['delta_deltaG_ends'] = ''
            row['MFE_guide'] = ''
            row['GC_content'] = ''
            row['td'] = ''
        else:
            feats = compute_all_features(sirna)
            row['dG_total'] = feats['dG_total']
            row['dG_seed'] = feats['dG_seed']
            row['delta_deltaG_ends'] = feats['delta_deltaG_ends']
            row['MFE_guide'] = feats['MFE_guide']
            row['GC_content'] = feats['GC_content']
            row['td'] = feats['td']
        out_rows.append(row)
        if (i + 1) % 10 == 0:
            print(f"  Processed {i+1}/{len(rows)}")

    with open(output_path, 'w', newline='', encoding='utf-8-sig') as f:
        writer = csv.DictWriter(f, fieldnames=out_fieldnames)
        writer.writeheader()
        writer.writerows(out_rows)

    print(f"\nSaved to {output_path}")
    print(f"Columns added: {', '.join(new_fields)}")


def main():
    parser = argparse.ArgumentParser(description='siRNA热力学特征计算工具')
    parser.add_argument('--input', type=str, help='输入CSV文件路径')
    parser.add_argument('--output', type=str, help='输出CSV文件路径')
    parser.add_argument('--sirna', type=str, help='单条siRNA序列（19nt）')
    parser.add_argument('--mrna', type=str, help='对应mRNA靶标序列')
    parser.add_argument('--sirna_col', type=str, default='siRNA', help='CSV中siRNA列名')
    parser.add_argument('--mrna_col', type=str, default='mRNA', help='CSV中mRNA列名')
    args = parser.parse_args()

    if args.sirna:
        # 单条模式
        feats = compute_all_features(args.sirna, args.mrna)
        print("\n" + "=" * 60)
        print("siRNA热力学特征")
        print("=" * 60)
        print(f"siRNA序列:       {feats['siRNA']}")
        print(f"GC含量:          {feats['GC_content']}%")
        print(f"整体双链ΔG:      {feats['dG_total']} kcal/mol")
        print(f"seed区ΔG(2-8):  {feats['dG_seed']} kcal/mol")
        print(f"末端ΔΔG(5'-3'): {feats['delta_deltaG_ends']} kcal/mol")
        print(f"引导链MFE:       {feats['MFE_guide']} kcal/mol")
        print(f"\nOligoFormer td (24维):")
        print(f"  {feats['td']}")
    elif args.input:
        # 批量模式
        output = args.output or args.input.replace('.csv', '_td.csv')
        process_csv(args.input, output, args.sirna_col, args.mrna_col)
    else:
        parser.print_help()
        print("\n示例:")
        print("  单条: python thermodynamics_calculator.py --sirna GCGCUCAUCAUUGUGCUGC")
        print("  批量: python thermodynamics_calculator.py --input mismatch.csv --output mismatch_td.csv")


if __name__ == '__main__':
    main()
