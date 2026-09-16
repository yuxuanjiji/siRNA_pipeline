import re
import os

path = './'
_FLANK = [19]

def antiRNA(RNA):
    result = []
    for i in RNA:
        if i == 'A' or i == 'a':
            result.append('U')
        elif i == 'U' or i == 'u' or i == 'T' or i == 't':
            result.append('A')
        elif i == 'C' or i == 'c':
            result.append('G')
        elif i == 'G' or i == 'g':
            result.append('C')
    return ''.join(result[::-1])

def read_fasta(filename):
    sequences = []
    with open(filename, 'r') as f:
        lines = f.readlines()
        for i in range(1, len(lines), 2):
            sequences.append(lines[i].strip())
    return sequences

def read_csv(filename):
    data = []
    with open(filename, 'r') as f:
        header = f.readline().strip().split(',')
        for line in f:
            parts = line.strip().split(',')
            row = dict(zip(header, parts))
            data.append(row)
    return header, data

# 读取 mRNA 序列
Hu_mRNA = read_fasta(path + 'data/fasta/Hu_mRNA.fa')
for i in range(len(Hu_mRNA)):
    Hu_mRNA[i] = 'X' * 300 + Hu_mRNA[i].replace('T', 'U') + 'X' * 300

# 读取 raw CSV 数据
raw_header, raw_data = read_csv(path + 'data/process_demo/raw/Hu.csv')

# 读取 data CSV 数据
data_header, data_data = read_csv(path + 'data/Hu.csv')

# 比较同一 siRNA 的 label 值
print("=" * 80)
print("比较 raw 和 data 中同一 siRNA 的 label 值")
print("=" * 80)

for i in range(min(20, len(raw_data))):
    raw_seq = raw_data[i]['siRNA']
    raw_label = raw_data[i]['label']

    # 在 data 中找相同的 siRNA
    for j in range(len(data_data)):
        if data_data[j]['siRNA'] == raw_seq:
            data_label = data_data[j]['label']
            ratio = float(raw_label) / float(data_label) if float(data_label) != 0 else 0
            print(f"siRNA: {raw_seq[:20]}...")
            print(f"  raw label:  {raw_label}")
            print(f"  data label: {data_label}")
            print(f"  比值 (raw/data): {ratio:.6f}")
            print()
            break

print("=" * 80)
print("分析归一化系数")
print("=" * 80)

# 计算平均比值
ratios = []
for i in range(min(100, len(raw_data))):
    raw_seq = raw_data[i]['siRNA']
    raw_label = raw_data[i]['label']

    for j in range(len(data_data)):
        if data_data[j]['siRNA'] == raw_seq:
            data_label = data_data[j]['label']
            if float(data_label) != 0:
                ratio = float(raw_label) / float(data_label)
                ratios.append(ratio)
            break

if ratios:
    import statistics
    print(f"样本数量: {len(ratios)}")
    print(f"平均比值: {statistics.mean(ratios):.6f}")
    print(f"标准差:   {statistics.stdev(ratios):.6f}")
    print(f"最小值:   {min(ratios):.6f}")
    print(f"最大值:   {max(ratios):.6f}")