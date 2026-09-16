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

def write_csv(filename, header, data):
    with open(filename, 'w', newline='') as f:
        f.write(','.join(header) + '\n')
        for row in data:
            row_values = []
            for col in header:
                if col in row:
                    val = str(row[col])
                    # 如果值包含逗号，需要用引号包围
                    if ',' in val:
                        val = f'"{val}"'
                    row_values.append(val)
                else:
                    row_values.append('')
            f.write(','.join(row_values) + '\n')

# 读取 mRNA 序列
Hu_mRNA = read_fasta(path + 'data/fasta/Hu_mRNA.fa')
for i in range(len(Hu_mRNA)):
    Hu_mRNA[i] = 'X' * 300 + Hu_mRNA[i].replace('T', 'U') + 'X' * 300

Mix_mRNA = read_fasta(path + 'data/fasta/Mix_mRNA.fa')
for i in range(len(Mix_mRNA)):
    Mix_mRNA[i] = 'X' * 300 + Mix_mRNA[i].replace('T', 'U') + 'X' * 300

Taka_mRNA = read_fasta(path + 'data/fasta/Taka_mRNA.fa')
for i in range(len(Taka_mRNA)):
    Taka_mRNA[i] = 'X' * 300 + Taka_mRNA[i].replace('T', 'U') + 'X' * 300

# 读取 CSV 数据
hu_header, Hu = read_csv(path + 'data/Hu.csv')
mix_header, new = read_csv(path + 'data/Mix.csv')
taka_header, Taka = read_csv(path + 'data/Taka.csv')

for FLANK in _FLANK:
    print('FLANK:', FLANK)
    os.makedirs(path + 'flanking/' + str(FLANK) + '/fasta/', exist_ok=True)
    
    # 处理 Hu
    for i, row in enumerate(Hu):
        count = 0
        for j, mrna in enumerate(Hu_mRNA):
            res = re.search(antiRNA(row['siRNA']), mrna)
            if res is not None:
                start = res.span(0)[0] - FLANK
                end = res.span(0)[1] + FLANK
                row['mRNA'] = mrna[start:end]
                count += 1
        if count == 0:
            print(f"Hu 第 {i} 条 siRNA {row['siRNA']} 未找到匹配")
    
    # 处理 Mix
    for i, row in enumerate(new):
        count = 0
        for j, mrna in enumerate(Mix_mRNA):
            res = re.search(antiRNA(row['siRNA']), mrna)
            if res is not None:
                start = res.span(0)[0] - FLANK
                end = res.span(0)[1] + FLANK
                row['mRNA'] = mrna[start:end]
                count += 1
                break
        if count == 0:
            print(f"Mix 第 {i} 条 siRNA {row['siRNA']} 未找到匹配")
    
    # 处理 Taka
    for i, row in enumerate(Taka):
        count = 0
        for j, mrna in enumerate(Taka_mRNA):
            res = re.search(antiRNA(row['siRNA']), mrna)
            if res is not None:
                start = res.span(0)[0] - FLANK
                end = res.span(0)[1] + FLANK
                row['mRNA'] = mrna[start:end]
                count += 1
        if count == 0:
            print(f"Taka 第 {i} 条 siRNA {row['siRNA']} 未找到匹配")
        elif count > 1:
            print(f"Taka 第 {i} 条 siRNA 找到 {count} 个匹配")
    
    # 写入 Hu
    with open(path + 'flanking/' + str(FLANK) + '/fasta/Hu_mRNA.fa', 'w') as f:
        for i, row in enumerate(Hu):
            f.write('>RNA' + str(i) + '\n')
            f.write(row['mRNA'] + '\n')
    write_csv(path + 'flanking/' + str(FLANK) + '/Hu.csv', hu_header, Hu)
    
    # 写入 Mix
    with open(path + 'flanking/' + str(FLANK) + '/fasta/Mix_mRNA.fa', 'w') as f:
        for i, row in enumerate(new):
            f.write('>RNA' + str(i) + '\n')
            f.write(row['mRNA'] + '\n')
    write_csv(path + 'flanking/' + str(FLANK) + '/Mix.csv', mix_header, new)
    
    # 写入 Taka
    with open(path + 'flanking/' + str(FLANK) + '/fasta/Taka_mRNA.fa', 'w') as f:
        for i, row in enumerate(Taka):
            f.write('>RNA' + str(i) + '\n')
            f.write(row['mRNA'] + '\n')
    write_csv(path + 'flanking/' + str(FLANK) + '/Taka.csv', taka_header, Taka)

print("处理完成！")