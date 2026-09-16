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

datasets = ['Mix', 'Taka']

for dataset in datasets:
    print("=" * 80)
    print(f"{dataset}.csv 分析")
    print("=" * 80)

    raw_data = read_csv(f'./data/process_demo/raw/{dataset}.csv')
    data_data = read_csv(f'./data/{dataset}.csv')

    # raw 统计
    raw_labels = []
    for row in raw_data:
        try:
            raw_labels.append(float(row['label']))
        except:
            pass

    print(f"\nraw/{dataset}.csv label 统计：")
    print(f"  最大值: {max(raw_labels):.6f}")
    print(f"  最小值: {min(raw_labels):.6f}")
    print(f"  label > 1 的记录数: {sum(1 for x in raw_labels if x > 1)}")

    # data 统计
    data_labels = []
    for row in data_data:
        try:
            data_labels.append(float(row['label']))
        except:
            pass

    print(f"\ndata/{dataset}.csv label 统计：")
    print(f"  最大值: {max(data_labels):.6f}")
    print(f"  最小值: {min(data_labels):.6f}")
    print(f"  label > 1 的记录数: {sum(1 for x in data_labels if x > 1)}")

    # 计算归一化系数
    if max(data_labels) > 0:
        ratio = max(raw_labels) / max(data_labels)
        print(f"\n归一化系数 (raw最大值 / data最大值): {ratio:.6f}")

    # 验证归一化
    print("\n验证 (比较同一 siRNA):")
    for i in range(min(3, len(raw_data))):
        raw_seq = raw_data[i]['siRNA']
        raw_label = raw_data[i]['label']
        for j in range(len(data_data)):
            if data_data[j]['siRNA'] == raw_seq:
                data_label = data_data[j]['label']
                print(f"  siRNA: {raw_seq[:20]}...")
                print(f"    raw:  {raw_label}")
                print(f"    data: {data_label}")
                print(f"    比值: {float(raw_label)/float(data_label):.6f}")
                break

    print()