# -*- coding: utf-8 -*-
import csv, os, shutil
base = r'D:\HuaweiMoveData\Users\及宇轩\Desktop\生科挑战赛\siRNA_pipeline\OligoFormer部分\data'
comp = str.maketrans('ACGU','UGCA')
def rc(s): return s.translate(comp)[::-1]

# WT行验证: pos67=20 -> 57nt窗[1:58] -> pos57=19(中心)
START, END = 1, 58
rows = list(csv.DictReader(open(base+r'\mismatch.csv', encoding='utf-8-sig')))
fa = [l.strip() for l in open(base+r'\fasta\mismatch_mRNA.fa') if l.strip()]
fa_mrnas = [fa[i] for i in range(1, len(fa), 2)]

# 备份
shutil.copy(base+r'\mismatch.csv', base+r'\mismatch_67nt_backup.csv')
shutil.copy(base+r'\fasta\mismatch_mRNA.fa', base+r'\fasta\mismatch_mRNA_67nt_backup.fa')

# CSV: mRNA列统一截[1:58]
for r in rows:
    assert len(r['mRNA'])==67
    r['mRNA'] = r['mRNA'][START:END]
cols = list(rows[0].keys())
with open(base+r'\mismatch.csv','w',encoding='utf-8-sig',newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=cols); w.writeheader(); w.writerows(rows)

# fasta: 同截法，保持RNA0..17顺序
with open(base+r'\fasta\mismatch_mRNA.fa','w') as fh:
    for i, s in enumerate(fa_mrnas):
        assert len(s)==67
        w57 = s[START:END]
        fh.write('>RNA%d\n%s\n' % (i, w57))

# 验证WT行中心
r18 = rows[18]
t = rc(r18['siRNA']); p = r18['mRNA'].find(t)
print('CSV row18 WT pos in 57nt:', p, '(expect 19)')
print('mRNA lens:', len(rows[0]['mRNA']))
print('done')
