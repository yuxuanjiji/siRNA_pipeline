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
unnorm_data = read_csv('./data/unnorm/Taka.csv')
data_data = read_csv('./data/Taka.csv')

print("=" * 80)
print("unnorm/Taka.csv vs data/Taka.csv 分析")
print("=" * 80)

# 样本数量
print(f"\nunnorm/Taka.csv 样本数: {len(unnorm_data)}")
print(f"data/Taka.csv 样本数: {len(data_data)}")

# unnorm label 统计
unnorm_labels = [float(row['label']) for row in unnorm_data if row['label']]
print(f"\nunnorm/Taka.csv label 统计：")
print(f"  最大值: {max(unnorm_labels):.6f}")
print(f"  最小值: {min(unnorm_labels):.6f}")
print(f"  平均值: {sum(unnorm_labels)/len(unnorm_labels):.6f}")

# data label 统计
data_labels = [float(row['label']) for row in data_data if row['label']]
print(f"\ndata/Taka.csv label 统计：")
print(f"  最大值: {max(data_labels):.6f}")
print(f"  最小值: {min(data_labels):.6f}")
print(f"  平均值: {sum(data_labels)/len(data_labels):.6f}")

# 通过 siRNA 匹配来找对应关系
print("\n" + "=" * 80)
print("通过 siRNA 匹配 unnorm 和 data 的对应关系")
print("=" * 80)

matched = 0
ratios = []
for i, unnorm_row in enumerate(unnorm_data):
    unnorm_seq = unnorm_row['siRNA']
    unnorm_label = float(unnorm_row['label'])
    for j, data_row in enumerate(data_data):
        if data_row['siRNA'] == unnorm_seq:
            data_label = float(data_row['label'])
            ratio = unnorm_label / data_label if data_label != 0 else 0
            ratios.append(ratio)
            matched += 1
            if matched <= 10:
                print(f"matched[{i}] -> data[{j}]:")
                print(f"  siRNA: {unnorm_seq[:30]}...")
                print(f"  unnorm label: {unnorm_label:.6f}")
                print(f"  data label:   {data_label:.6f}")
                print(f"  比值: {ratio:.6f}")
            break

print(f"\n匹配的记录数: {matched}")

if ratios:
    print(f"\n比值统计 ({len(ratios)} 个样本):")
    print(f"  平均值: {sum(ratios)/len(ratios):.6f}")
    print(f"  最小值: {min(ratios):.6f}")
    print(f"  最大值: {max(ratios):.6f}")

    # 检查是否所有比值都相同
    unique_ratios = set(ratios)
    print(f"  不同比值数量: {len(unique_ratios)}")
    if len(unique_ratios) <= 5:
        print(f"  比值种类: {sorted(unique_ratios)}")