# -*- coding: utf-8 -*-
"""task15: 综合排名 vs 实验效率 Spearman（按数据真实相关性取向）"""
import pandas as pd, numpy as np
from scipy.stats import spearmanr, pearsonr
base = r'D:\HuaweiMoveData\Users\及宇轩\Desktop\生科挑战赛\siRNA_pipeline\OligoFormer部分\data'
td = pd.read_csv(base+r'\mismatch_td.csv', encoding='utf-8-sig')
pred = pd.read_csv(r'mismatch_predictions.csv', encoding='utf-8-sig')

df = td.copy()
df['oligo_pred'] = pred['pred_efficiency'].values
df['label'] = df['label'].astype(float)

# 剔除MFE_guide(全NaN)，用三特征+oligo
THERM_RAW = {'dG_total': None, 'dG_seed': None, 'delta_deltaG_ends': None}
def mm(s, hb):
    s = s.astype(float); r = (s-s.min())/(s.max()-s.min()+1e-9)
    return r if hb else 1-r

# 按 raw Pearson 符号定方向
for c in THERM_RAW:
    cc = np.sign(df[c].astype(float).corr(df['label'])) or 1
    df['s_'+c] = mm(df[c], (cc>0))
THERM = ['s_dG_total','s_dG_seed','s_delta_deltaG_ends']
df['s_oligo'] = mm(df['oligo_pred'], True)
df['S_thermo'] = df[THERM].mean(axis=1)
df['S_combo']  = 0.7*df['S_thermo'] + 0.3*df['s_oligo']

def rep(col,name):
    r,p = spearmanr(df[col], df['label'])
    rp,pp = pearsonr(df[col], df['label'])
    print('  %-30s Spearman=%.3f (p=%.3f)  Pearson=%.3f' % (name, r, p, rp))

print('=== task15 Spearman vs 实验label (n=%d, Holen错配集) ===' % len(df))
rep('S_thermo','纯热力学主排序(3特征等权)')
rep('s_oligo','纯OligoFormer')
rep('S_combo','综合(0.7热+0.3奥)')
df.to_csv('mismatch_ranked.csv', index=False, encoding='utf-8-sig')
print('saved mismatch_ranked.csv')
