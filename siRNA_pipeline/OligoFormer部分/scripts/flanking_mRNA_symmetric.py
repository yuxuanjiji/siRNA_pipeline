import pandas as pd
import numpy as np
import re 
import copy
import os  # 添加 os 模块用于创建目录

path = './'#path 是当前目录的路径
# _FLANK = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,50,60,70,80,90,100,200,300]
_FLANK = [19]#可以选取侧翼长度,最终应该选择为19,尝试更换为其他数值,观察是否与论文结论一致
def antiRNA(RNA):#函数功能:形成siRNA序列互补链并且反转
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
    return ''.join(antiRNA[::-1])#''.join函数的作用是将列表中的元素用空字符串连接起来，得到一个字符串


Hu_mRNA = pd.read_csv(path + 'data/fasta/Hu_mRNA.fa',header=None)[1::2].reset_index(drop=True)#需要mRNA序列，所以需要读取fasta文件
for i in range(Hu_mRNA.shape[0]):
    Hu_mRNA.iloc[i,0] = 'X' * 300 + Hu_mRNA.iloc[i,0].replace('T','U') + 'X' * 300 #在mRNA序列的左右两侧添加300个X，是为了后面进行侧翼序列的提取
#第 i 行，第 0 列 的那个 mRNA 序列
new_mRNA = pd.read_csv(path + 'data/fasta/new_mRNA.fa',header=None)[1::2].reset_index(drop=True)
for i in range(new_mRNA.shape[0]):
    new_mRNA.iloc[i,0] = 'X' * 300 + new_mRNA.iloc[i,0].replace('T','U') + 'X' * 300 

Taka_mRNA = pd.read_csv(path + 'data/fasta/Taka_mRNA.fa',header=None)[1::2].reset_index(drop=True)
for i in range(Taka_mRNA.shape[0]):
    Taka_mRNA.iloc[i,0] = 'X' * 300 + Taka_mRNA.iloc[i,0].replace('T','U') + 'X' * 300 

Hu_o = pd.read_csv(path + 'data/Hu.csv')
new_o = pd.read_csv(path + 'data/new.csv')
Taka_o = pd.read_csv(path + 'data/Taka.csv')#csv文件中有siRNA序列,沉默效率和y值

Hu = Hu_o.copy(deep=True)#深拷贝一个副本
new = new_o.copy(deep=True)
Taka = Taka_o.copy(deep=True)
Hu['flanking'] = Hu['y']# 创建 flanking 列，用 y 先填充，占位置，后面再改成真正的 mRNA 序列
new['flanking'] = new['y']
Taka['flanking'] = Taka['y']


for FLANK in _FLANK:
    print('FLANK:',FLANK)
    
    # 创建输出目录
    fasta_dir = path + 'flanking/' + str(FLANK) + '/fasta/'
    os.makedirs(fasta_dir, exist_ok=True)
    
    for i in range(Hu.shape[0]):
        count = 0
        for j in range(Hu_mRNA.shape[0]):
            res = re.search(antiRNA(Hu.iloc[i,0]),Hu_mRNA.iloc[j,0])
            if res is not None:
                Hu.iloc[i,-1] = Hu_mRNA.iloc[j,0][(res.span(0)[0] - FLANK) : (res.span(0)[1] + FLANK)]
                count += 1
        if count == 0:
            print(i,j,'here') 
    for i in range(new.shape[0]):
        count = 0
        for j in range(new_mRNA.shape[0]):
            res = re.search(antiRNA(new.iloc[i,0]),new_mRNA.iloc[j,0])
            if res is not None:
                new.iloc[i,-1] = new_mRNA.iloc[j,0][(res.span(0)[0] - FLANK) : (res.span(0)[1] + FLANK)]
                count += 1
                break
        if count == 0:
            print(i,j,'here')
    for i in range(Taka.shape[0]):
        count = 0
        for j in range(Taka_mRNA.shape[0]):
            res = re.search(antiRNA(Taka.iloc[i,0]),Taka_mRNA.iloc[j,0])
            if res is not None:
                Taka.iloc[i,-1] = Taka_mRNA.iloc[j,0][(res.span(0)[0] - FLANK) : (res.span(0)[1] + FLANK)]
                count += 1
                #span = 返回一个元组，元组的第一个元素是匹配到的字符串的起始位置，第二个元素是匹配到的字符串的结束位置
        if count > 1:
            print(i,j)
    # for i in `ls ./`; do mv $i/*.fa $i/fasta/ ; done
    with open(path + 'flanking/' + str(FLANK) + '/fasta/' + 'Hu_mRNA.fa','w') as f:
        for i in range(Hu.shape[0]):
            f.write('>RNA' + str(i) + '\n')
            f.write(Hu.iloc[i,-1] + '\n')
    Hu_o['mRNA'] = Hu['flanking']
    Hu_o.to_csv(path + 'flanking/' + str(FLANK) + '/Hu.csv',index=False)
    with open(path + 'flanking/' + str(FLANK) + '/fasta/' + 'new_mRNA.fa','w') as f:
        for i in range(new.shape[0]):
            f.write('>RNA' + str(i) + '\n')
            f.write(new.iloc[i,-1] + '\n')
    new_o['mRNA'] = new['flanking']
    new_o.to_csv(path + 'flanking/' + str(FLANK) + '/new.csv',index=False)
    with open(path + 'flanking/' + str(FLANK) + '/fasta/' + 'Taka_mRNA.fa','w') as f:
        for i in range(Taka.shape[0]):
            f.write('>RNA' + str(i) + '\n')
            f.write(Taka.iloc[i,-1] + '\n')
    Taka_o['mRNA'] = Taka['flanking']
    Taka_o.to_csv(path + 'flanking/' + str(FLANK) + '/Taka.csv',index=False)