import re

def read_csv(filename):
    data = []
    with open(filename, 'r') as f:
        header = f.readline().strip().split(',')
        for line in f:
            parts = line.strip().split(',')
            row = dict(zip(header, parts))
            data.append(row)
    return data

# 读取 Taka 数据
raw_data = read_csv('./data/process_demo/raw/Taka.csv')
data_data = read_csv('./data/Taka.csv')

print("=" * 80)
print("raw/Taka.csv 前5条 siRNA:")
print("=" * 80)
for i in range(min(5, len(raw_data))):
    print(f"[{i}] {raw_data[i]['siRNA']}")

print("\n" + "=" * 80)
print("data/Taka.csv 前5条 siRNA:")
print("=" * 80)
for i in range(min(5, len(data_data))):
    print(f"[{i}] {data_data[i]['siRNA']}")

print("\n" + "=" * 80)
print("检查 siRNA 长度:")
print("=" * 80)
raw_lens = set(len(row['siRNA']) for row in raw_data)
data_lens = set(len(row['siRNA']) for row in data_data)
print(f"raw siRNA 长度: {raw_lens}")
print(f"data siRNA 长度: {data_lens}")

print("\n" + "=" * 80)
print("检查 data/Taka.csv 是否有 mRNA 列:")
print("=" * 80)
print(f"data/Taka.csv 列名: {list(data_data[0].keys())}")