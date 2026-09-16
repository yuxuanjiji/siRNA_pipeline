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

ratios = []
for i, unnorm_row in enumerate(unnorm_data):
    unnorm_seq = unnorm_row['siRNA']
    unnorm_label = float(unnorm_row['label'])
    for j, data_row in enumerate(data_data):
        if data_row['siRNA'] == unnorm_seq:
            data_label = float(data_row['label'])
            if data_label != 0:
                ratio = unnorm_label / data_label
                ratios.append(ratio)
            break

print(f"有效样本数（data label != 0）: {len(ratios)}")

if ratios:
    print(f"比值统计:")
    print(f"  平均值: {sum(ratios)/len(ratios):.6f}")
    print(f"  最小值: {min(ratios):.6f}")
    print(f"  最大值: {max(ratios):.6f}")

    unique_ratios = set(ratios)
    print(f"  不同比值数量: {len(unique_ratios)}")

    if len(unique_ratios) <= 20:
        for r in sorted(unique_ratios):
            count = ratios.count(r)
            print(f"    {r:.6f}: {count} 次")

print(f"\n可能的归一化系数:")
print(f"  unnorm最大值 / data最大值 = 0.9799 / 1.0 = 0.9799")
print(f"  1 / 0.9799 = {1/0.9799:.6f}")

factor = 1 / 0.9799
print(f"\nfactor = 1 / 0.9799 = {factor:.6f}")

correct = 0
wrong = 0
for i, unnorm_row in enumerate(unnorm_data):
    unnorm_seq = unnorm_row['siRNA']
    unnorm_label = float(unnorm_row['label'])
    for j, data_row in enumerate(data_data):
        if data_row['siRNA'] == unnorm_seq:
            data_label = float(data_row['label'])
            if data_label != 0:
                expected = unnorm_label / 0.9799
                if abs(expected - data_label) < 0.0001:
                    correct += 1
                else:
                    wrong += 1
                    if wrong <= 5:
                        print(f"  不匹配: unnorm={unnorm_label}, expected={expected:.6f}, data={data_label}")
            break

print(f"\n验证结果: 正确={correct}, 错误={wrong}")