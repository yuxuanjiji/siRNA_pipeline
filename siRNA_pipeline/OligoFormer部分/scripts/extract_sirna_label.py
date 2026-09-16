import re

with open('flanking/19/Taka.csv', 'r') as f:
    content = f.read()

lines = content.split('\n')
header = lines[0].split(',')
siRNA_idx = header.index('siRNA')
label_idx = header.index('label')

# 使用正则表达式匹配CSV行（处理td列跨行的情况）
pattern = r'^([^,]+),([^,]+),([^,]+),([^,]+),(".+?"|[^,]+)(?:,|$)'

output_rows = ['siRNA,label']
count = 0

for line in lines[1:]:
    if not line.strip():
        continue
    # 找到siRNA和label的位置
    parts = line.split(',')
    if len(parts) >= 5:
        siRNA = parts[siRNA_idx]
        label = parts[label_idx]
        output_rows.append(f'{siRNA},{label}')
        count += 1

with open('flanking/19/Taka_simplified.csv', 'w') as f:
    f.write('\n'.join(output_rows))

print(f"已完成，共 {count} 行数据")
print(f"输出文件: flanking/19/Taka_simplified.csv")