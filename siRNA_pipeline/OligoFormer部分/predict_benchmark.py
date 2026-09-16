# -*- coding: utf-8 -*-
"""task14: OligoFormer 公开基准预测 (Hu/Taka/Mix)
用法: python predict_benchmark.py Hu
"""
import os, sys, csv, argparse
import numpy as np, pandas as pd, torch
from torch.utils.data import DataLoader
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(SCRIPT_DIR, "scripts"))
from loader import data_process_loader
from model import Oligo

ap = argparse.ArgumentParser()
ap.add_argument("dataset")
ap.add_argument("--batch", type=int, default=16)
args = ap.parse_args()
DS = args.dataset
DATA_PATH = os.path.join(SCRIPT_DIR, "data")
device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
df = pd.read_csv(os.path.join(DATA_PATH, f"{DS}.csv"), dtype=str)
print(f"[{DS}] Loaded {len(df)} records")
params = {"batch_size": args.batch, "shuffle": False, "num_workers": 0, "drop_last": False}
ds = DataLoader(data_process_loader(df.index.values, df.label.values, df.y.values, df, DS, DATA_PATH + "/"), **params)
model = Oligo(vocab_size=26, embedding_dim=128, lstm_dim=32, n_head=8, n_layers=1, lm1=19, lm2=19).to(device)
model.load_state_dict(torch.load(os.path.join(SCRIPT_DIR,"model","best_model.pth"), map_location=device))
model.eval()
preds=[]; trues=[]
with torch.no_grad():
    for data in ds:
        siRNA=data[0].to(device); mRNA=data[1].to(device)
        siFM=data[2].to(device); mFM=data[3].to(device)
        label=data[4]; td=data[6].to(device)
        out,_,_=model(siRNA,mRNA,siFM,mFM,td)
        preds.extend(out[:,1].cpu().numpy())
        trues.extend(label.numpy())
from scipy.stats import pearsonr, spearmanr
t=np.array(trues); p=np.array(preds)
pcc,_=pearsonr(t,p); spr,_=spearmanr(t,p)
rmse=np.sqrt(((t-p)**2).mean())
print(f"[{DS}] n={len(t)}  PCC={pcc:.4f}  Spearman={spr:.4f}  RMSE={rmse:.4f}")
output = df.loc[:len(preds) - 1, ["siRNA", "mRNA", "label"]].copy()
output["pred"] = p
output.to_csv(os.path.join(SCRIPT_DIR, f"{DS}_predictions.csv"), index=False)
