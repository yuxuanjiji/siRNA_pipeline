import os

data_dir = r'c:\Users\jiyuxuan\Downloads\OligoFormer-main\OligoFormer-main\data'

print("="*60)
print("对比用户标准答案 vs 当前文件")
print("="*60)

# 用户提供的标准答案（从用户选择的前26行）
standard_sequences = [
    "AAAUUGAAAGGAAUUGUAUAAAUCAAUUAACAUAUUAGCUGAGUUGXXXXXXXXXXX",  # RNA0
    "UAAAAUUGAAAGGAAUUGUAUAAAUCAAUUAACAUAUUAGCUGAGUUGXXXXXXXXX",  # RNA1
    "CUUAUUUUUAGAUAAAAUUGAAAGGAAUUGUAUAAAUCAAUUAACAUAUUAGCUGAG",  # RNA2
    "UUGGUCUGCUUAUUUUUAGAUAAAAUUGAAAGGAAUUGUAUAAAUCAAUUAACAUAU",  # RNA3
    "UGGAAUCUUAUGUAACUUUCUUAUUUAAUUUUGGUCUGCUUAUUUUUAGAUAAAAUU",  # RNA4
    "AUGGAAUCUUAUGUAACUUUCUUAUUUAAUUUUGGUCUGCUUAUUUUUAGAUAAAAU",  # RNA5
    "AAAUGGAAUCUUAUGUAACUUUCUUAUUUAAUUUUGGUCUGCUUAUUUUUAGAUAAA",  # RNA6
    "UCAAAUGGAAUCUUAUGUAACUUUCUUAUUUAAUUUUGGUCUGCUUAUUUUUAGAUA",  # RNA7
    "UUAUGUAUGCUCAAAUGGAAUCUUAUGUAACUUUCUUAUUUAAUUUUGGUCUGCUUA",  # RNA8
    "UGGCCUUAUGUAUGCUCAAAUGGAAUCUUAUGUAACUUUCUUAUUUAAUUUUGGUCU",  # RNA9
    "UAUCAUGGCCUUAUGUAUGCUCAAAUGGAAUCUUAUGUAACUUUCUUAUUUAAUUUU",  # RNA10
    "AAAGUAUCAUGGCCUUAUGUAUGCUCAAAUGGAAUCUUAUGUAACUUUCUUAUUUAA",  # RNA11
    "CACAUUAAAGUAUCAUGGCCUUAUGUAUGCUCAAAUGGAAUCUUAUGUAACUUUCUU",  # RNA12
]

# 读取当前FASTA文件的前13条序列
with open(os.path.join(data_dir, 'fasta', 'Hu_mRNA.fa'), 'r') as f:
    content = f.read()
    current_sequences = []
    for line in content.split('\n'):
        if not line.startswith('>') and line.strip():
            current_sequences.append(line.strip())

print(f"\n标准答案序列数量: {len(standard_sequences)}")
print(f"当前FASTA序列数量: {len(current_sequences)}")

print("\n" + "="*60)
print("逐条对比前13条序列")
print("="*60)
for i in range(min(13, len(standard_sequences), len(current_sequences))):
    std = standard_sequences[i]
    curr = current_sequences[i]
    match = "✓" if std == curr else "✗"
    print(f"\nRNA{i}:")
    print(f"  标准: {std}")
    print(f"  当前: {curr}")
    print(f"  一致: {match}")
    if std != curr:
        print(f"  差异位置: ", end="")
        for j in range(min(len(std), len(curr))):
            if std[j] != curr[j]:
                print(f"{j}(标准={std[j]},当前={curr[j]})", end=" ")
        print()

print("\n" + "="*60)
print("用户选择的是哪一行？")
print("="*60)
print("用户选择的行内容:")
print("CACAUUAAAGUAUCAUGGCCUUAUGUAUGCUCAAAUGGAAUCUUAUGUAACUUUCUU")
print("\n这是RNA12的序列，长度:", len("CACAUUAAAGUAUCAUGGCCUUAUGUAUGCUCAAAUGGAAUCUUAUGUAACUUUCUU"))