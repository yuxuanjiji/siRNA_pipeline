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
print("Taka.csv 详细分析")
print("=" * 80)

# 样本数量
print(f"\nraw/Taka.csv 样本数: {len(raw_data)}")
print(f"data/Taka.csv 样本数: {len(data_data)}")

# raw label 统计
raw_labels = [float(row['label']) for row in raw_data if row['label']]
print(f"\nraw/Taka.csv label 统计：")
print(f"  最大值: {max(raw_labels):.6f}")
print(f"  最小值: {min(raw_labels):.6f}")
print(f"  平均值: {sum(raw_labels)/len(raw_labels):.6f}")

# data label 统计
data_labels = [float(row['label']) for row in data_data if row['label']]
print(f"\ndata/Taka.csv label 统计：")
print(f"  最大值: {max(data_labels):.6f}")
print(f"  最小值: {min(data_labels):.6f}")
print(f"  平均值: {sum(data_labels)/len(data_labels):.6f}")

# 计算归一化系数
ratio = max(raw_labels) / max(data_labels)
print(f"\n归一化系数 (raw最大值 / data最大值): {ratio:.6f}")

# 尝试通过 siRNA 匹配来找对应关系
print("\n" + "=" * 80)
print("通过 siRNA 匹配 raw 和 data 的对应关系")
print("=" * 80)

matched = 0
unmatched_raw = []
for i, raw_row in enumerate(raw_data):
    raw_seq = raw_row['siRNA']
    raw_label = raw_label_val = float(raw_row['label'])
    for j, data_row in enumerate(data_data):
        if data_row['siRNA'] == raw_seq:
            data_label_val = float(data_row['label'])
            ratio_i = raw_label_val / data_label_val if data_label_val != 0 else 0
            matched += 1
            if matched <= 10:
                print(f"matched[{i}] -> data[{j}]:")
                print(f"  siRNA: {raw_seq[:30]}...")
                print(f"  raw label:  {raw_label_val:.6f}")
                print(f"  data label: {data_label_val:.6f}")
                print(f"  比值: {ratio_i:.6f}")
            break
    else:
        unmatched_raw.append(i)

print(f"\n匹配的记录数: {matched}")
print(f"未匹配的 raw 记录数: {len(unmatched_raw)}")

# 如果匹配，检查比值是否一致
if matched > 0:
    ratios = []
    for raw_row in raw_data:
        raw_seq = raw_row['siRNA']
        raw_label_val = float(raw_row['label'])
        for data_row in data_data:
            if data_row['siRNA'] == raw_seq:
                data_label_val = float(data_row['label'])
                if data_label_val != 0:
                    ratios.append(raw_label_val / data_label_val)
                break
    if ratios:
        print(f"\n比值统计 ({len(ratios)} 个样本):")
        print(f"  平均值: {sum(ratios)/len(ratios):.6f}")
        print(f"  最小值: {min(ratios):.6f}")
        print(f"  最大值: {max(ratios):.6f}")