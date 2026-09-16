import pandas as pd
import numpy as np
import re
import copy
import os

path = './'
_FLANK = [19]

def antiRNA(RNA):
    antiRNA = []
    for i in RNA:
        if i == 'A' or i == 'a':
            antiRNA.append('U')
        elif i == 'U' or i == 'u' or i == 'T' or i == 't':
            antiRNA.append('A')
        elif i == 'C' or i == 'c':
            antiRNA.append('G')
        elif i == 'G' or i == 'g':
            antiRNA.append('C')
    return ''.join(antiRNA[::-1])


Hu_mRNA = pd.read_csv(path + 'data/fasta/Hu_mRNA.fa', header=None)[1::2].reset_index(drop=True)
for i in range(Hu_mRNA.shape[0]):
    Hu_mRNA.iloc[i,0] = 'X' * 300 + Hu_mRNA.iloc[i,0].replace('T','U') + 'X' * 300


Hu_o = pd.read_csv(path + 'data/Hu.csv')


Hu = Hu_o.copy(deep=True)

Hu['flanking'] = Hu['y']



for FLANK in _FLANK:
    print('FLANK:',FLANK)
    os.makedirs(path + 'flanking/' + str(FLANK) + '/fasta/', exist_ok=True)

    for i in range(Hu.shape[0]):
        count = 0
        for j in range(Hu_mRNA.shape[0]):
            res = re.search(antiRNA(Hu.iloc[i,0]),Hu_mRNA.iloc[j,0])
            if res is not None:
                Hu.iloc[i,-1] = Hu_mRNA.iloc[j,0][(res.span(0)[0] - FLANK) : (res.span(0)[1] + FLANK)]
                count += 1
        if count == 0:
            print(i,j,'here')

        if count > 1:
            print(i,j)

    with open(path + 'flanking/' + str(FLANK) + '/fasta/' + 'Hu_mRNA.fa','w') as f:
        for i in range(Hu.shape[0]):
            f.write('>RNA' + str(i) + '\n')
            f.write(Hu.iloc[i,-1] + '\n')
    Hu_o['mRNA'] = Hu['flanking']
    Hu_o.to_csv(path + 'flanking/' + str(FLANK) + '/Hu.csv',index=False)