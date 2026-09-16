import csv
import os

# 输入文件路径
input_dir = 'c:/Users/jiyuxuan/Downloads/OligoFormer-main/OligoFormer-main/data/process_demo/raw/'
output_dir = 'c:/Users/jiyuxuan/Downloads/OligoFormer-main/OligoFormer-main/data/process_demo/raw_no_mRNA/'

# 确保输出目录存在
os.makedirs(output_dir, exist_ok=True)

# 处理的文件名列表
datasets = ['Hu.csv', 'Mix.csv', 'Taka.csv']

for dataset in datasets:
    input_path = input_dir + dataset
    output_path = output_dir + dataset
    
    with open(input_path, 'r', newline='', encoding='utf-8') as infile:
        reader = csv.reader(infile)
        header = next(reader)
        
        # 找到 mRNA 列的索引
        if 'mRNA' in header:
            mrna_idx = header.index('mRNA')
            # 保留除 mRNA 外的列
            new_header = [col for col in header if col != 'mRNA']
            
            with open(output_path, 'w', newline='', encoding='utf-8') as outfile:
                writer = csv.writer(outfile)
                writer.writerow(new_header)
                
                for row in reader:
                    # 删除 mRNA 列
                    new_row = [cell for i, cell in enumerate(row) if i != mrna_idx]
                    writer.writerow(new_row)
            
            print(f"已处理 {dataset}，删除了 mRNA 列")
        else:
            # 如果没有 mRNA 列，直接复制文件
            with open(output_path, 'w', newline='', encoding='utf-8') as outfile:
                writer = csv.writer(outfile)
                writer.writerow(header)
                for row in reader:
                    writer.writerow(row)
            print(f"{dataset} 中不存在 mRNA 列，直接复制")

print(f"\n处理完成！结果保存在: {output_dir}")