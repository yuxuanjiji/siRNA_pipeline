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

# 读取 raw 和 data 数据
raw_data = read_csv('./data/process_demo/raw/Hu.csv')
data_data = read_csv('./data/Hu.csv')

# 找 raw 中 label > 1 的记录
print("=" * 80)
print("raw/Hu.csv 中 label > 1 的记录：")
print("=" * 80)

count = 0
for row in raw_data:
    try:
        label = float(row['label'])
        if label > 1:
            count += 1
            if count <= 10:
                print(f"siRNA: {row['siRNA'][:30]}... | label: {label}")
    except:
        pass

print(f"\n总数: {count} 条记录 label > 1")

# 找 raw 和 data 中 label 的最大值和最小值
print("\n" + "=" * 80)
print("raw/Hu.csv label 统计：")
print("=" * 80)

raw_labels = []
for row in raw_data:
    try:
        raw_labels.append(float(row['label']))
    except:
        pass

print(f"最大值: {max(raw_labels):.6f}")
print(f"最小值: {min(raw_labels):.6f}")
print(f"样本数: {len(raw_labels)}")

print("\n" + "=" * 80)
print("data/Hu.csv label 统计：")
print("=" * 80)

data_labels = []
for row in data_data:
    try:
        data_labels.append(float(row['label']))
    except:
        pass

print(f"最大值: {max(data_labels):.6f}")
print(f"最小值: {min(data_labels):.6f}")
print(f"样本数: {len(data_labels)}")

print("\n" + "=" * 80)
print("验证归一化系数：")
print("=" * 80)
print(f"raw 最大值 / data 最大值 = {max(raw_labels):.6f} / {max(data_labels):.6f} = {max(raw_labels)/max(data_labels):.6f}")