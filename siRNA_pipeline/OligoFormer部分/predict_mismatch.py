# -*- coding: utf-8 -*-
"""
OligoFormer错配siRNA预测脚本
用法: python predict_mismatch.py
输出: mismatch_predictions.csv (含预测效率和真实效率)
"""
import os
import sys
import csv
import numpy as np
import pandas as pd
import torch
import torch.nn.functional as F
from torch.utils.data import DataLoader

# 添加scripts目录到路径
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(SCRIPT_DIR, "scripts"))

from loader import data_process_loader
from model import Oligo

# 配置
DATASET = "mismatch"
DATA_PATH = os.path.join(SCRIPT_DIR, "data")
BEST_MODEL = os.path.join(SCRIPT_DIR, "model", "best_model.pth")
OUTPUT_FILE = os.path.join(SCRIPT_DIR, "mismatch_predictions.csv")
BATCH_SIZE = 16

device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
print(f"Device: {device}")

# 加载数据
df = pd.read_csv(os.path.join(DATA_PATH, f"{DATASET}.csv"), dtype=str)
print(f"Loaded {len(df)} records")

params = {"batch_size": BATCH_SIZE, "shuffle": False, "num_workers": 0, "drop_last": False}
ds = DataLoader(data_process_loader(df.index.values, df.label.values, df.y.values, df, DATASET, DATA_PATH + "/"), **params)

# 加载模型
model = Oligo(vocab_size=26, embedding_dim=128, lstm_dim=32, n_head=8, n_layers=1, lm1=19, lm2=19).to(device)
model.load_state_dict(torch.load(BEST_MODEL, map_location=device))
model.eval()
print(f"Model loaded: {BEST_MODEL}")

# 预测
all_pred_efficacy = []
all_pred_class = []
all_true_label = []
all_true_y = []

with torch.no_grad():
    for i, data in enumerate(ds):
        siRNA = data[0].to(device)
        mRNA = data[1].to(device)
        siRNA_FM = data[2].to(device)
        mRNA_FM = data[3].to(device)
        label = data[4].to(device)
        y = data[5]
        td = data[6].to(device)

        pred, _, _ = model(siRNA, mRNA, siRNA_FM, mRNA_FM, td)
        efficacy = pred[:, 1].cpu().numpy()  # 正类概率 = 预测效率
        pred_cls = torch.argmax(pred, dim=-1).cpu().numpy()

        all_pred_efficacy.extend(efficacy)
        all_pred_class.extend(pred_cls)
        all_true_label.extend(label.cpu().numpy())
        all_true_y.extend(y.numpy())

# 计算指标
from sklearn.metrics import accuracy_score, f1_score, precision_score, recall_score, roc_auc_score, matthews_corrcoef
from scipy.stats import pearsonr, spearmanr

true_eff = np.array(all_true_label)
pred_eff = np.array(all_pred_efficacy)
true_cls = np.array(all_true_y)
pred_cls = np.array(all_pred_class)

# PCC
pcc, pcc_p = pearsonr(true_eff, pred_eff)
# Spearman
spr, spr_p = spearmanr(true_eff, pred_eff)
# RMSE
rmse = np.sqrt(np.mean((true_eff - pred_eff) ** 2))
# MAE
mae = np.mean(np.abs(true_eff - pred_eff))
# 分类指标
acc = accuracy_score(true_cls, pred_cls)
f1 = f1_score(true_cls, pred_cls)
auc = roc_auc_score(true_cls, pred_eff) if len(np.unique(true_cls)) > 1 else 0
mcc = matthews_corrcoef(true_cls, pred_cls)

print("\n" + "=" * 60)
print("预测结果评估指标")
print("=" * 60)
print(f"PCC (Pearson):  {pcc:.4f} (p={pcc_p:.4f})")
print(f"Spearman:       {spr:.4f} (p={spr_p:.4f})")
print(f"RMSE:           {rmse:.4f}")
print(f"MAE:            {mae:.4f}")
print(f"Accuracy:       {acc:.4f}")
print(f"F1 Score:       {f1:.4f}")
print(f"AUC:            {auc:.4f}")
print(f"MCC:            {mcc:.4f}")
print(f"y=1 (高效):     {np.sum(true_cls==1)} / {len(true_cls)}")
print(f"预测y=1:        {np.sum(pred_cls==1)} / {len(pred_cls)}")

# 保存预测结果
with open(OUTPUT_FILE, "w", newline="", encoding="utf-8-sig") as f:
    writer = csv.writer(f)
    writer.writerow(["siRNA", "mRNA", "true_efficiency", "pred_efficiency", "true_y", "pred_y", "abs_error"])
    for i in range(len(df)):
        writer.writerow([
            df.iloc[i]["siRNA"],
            df.iloc[i]["mRNA"],
            f"{true_eff[i]:.4f}",
            f"{pred_eff[i]:.4f}",
            int(true_cls[i]),
            int(pred_cls[i]),
            f"{abs(true_eff[i] - pred_eff[i]):.4f}"
        ])

print(f"\n预测结果已保存: {OUTPUT_FILE}")
