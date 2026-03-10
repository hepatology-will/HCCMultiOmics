import os
import glob
import pickle
import pandas as pd
import numpy as np
import loompy as lp # 必须要安装 loompy
from dask.diagnostics import ProgressBar
from arboreto.algo import grnboost2

# --- 修正点 1: 新版 API 路径 ---
from ctxcore.rnkdb import FeatherRankingDatabase 
from pyscenic.utils import modules_from_adjacencies, load_motifs
from pyscenic.prune import prune2df, df2regulons
from pyscenic.aucell import aucell
# 移除了 pyscenic.export，因为新版经常报错

# ================= 配置区域 (保持不变) =================
F_EXPRESSION = "hep_mal_p1_counts.csv"
F_TFS = "hs_hgnc_tfs.txt"
F_MOTIFS = "motifs-v9-nr.hgnc-m0.001-o0.0.tbl"
# 请确认这是您刚才下载的那个文件
F_DB = "hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather"

F_OUTPUT_LOOM = "hep_mal_p1_SCENIC.loom"
F_OUTPUT_CSV = "hep_mal_p1_AUC.csv" # 新增：CSV 备份

N_CPU = 40
# =======================================================

def run():
    print("Step 1: Loading Expression Matrix...")
    ex_matrix = pd.read_csv(F_EXPRESSION, index_col=0)
    print(f"Matrix shape: {ex_matrix.shape} (Cells x Genes)")
    
    tf_names = load_tf_names(F_TFS)
    print(f"Loaded {len(tf_names)} TFs.")

    print("Step 2: Inferring Co-expression Network (GRNBoost2)...")
    adjacencies = grnboost2(ex_matrix, tf_names=tf_names, verbose=True)

    print("Step 3: Pruning Network with Motifs (CTX)...")
    modules = modules_from_adjacencies(adjacencies, ex_matrix)
    
    # 修正点: 这里的调用方式适配新版
    dbs = [FeatherRankingDatabase(fname=F_DB, name="hg38_10kb")]
    
    df = prune2df(dbs, modules, F_MOTIFS)
    regulons = df2regulons(df)
    print(f"Found {len(regulons)} regulons.")
    
    # 备份 Regulons
    with open("regulons_final.p", "wb") as f:
        pickle.dump(regulons, f)

    print("Step 4: Calculating Regulon Activity (AUCell)...")
    # 计算 AUC 矩阵
    auc_mtx = aucell(ex_matrix, regulons, num_workers=N_CPU)

    print("Step 5: Exporting Results...")
    
    # --- 方案 A: 保存为 CSV (最推荐，R 读取最稳) ---
    auc_mtx.to_csv(F_OUTPUT_CSV)
    print(f"AUC Matrix saved to {F_OUTPUT_CSV} (Use this for R visualization!)")

    # --- 方案 B: 手动保存为 LOOM (适配新版 loompy) ---
    # 创建基本的 loom 文件
    row_attrs = {"Gene": np.array(ex_matrix.columns)}
    col_attrs = {"CellID": np.array(ex_matrix.index)}
    
    lp.create(F_OUTPUT_LOOM, ex_matrix.T.values, row_attrs, col_attrs)
    
    # 将 AUC 结果写入 Loom
    with lp.connect(F_OUTPUT_LOOM, mode="r+") as ds:
        # 将 regulons AUC 添加到 loom 的全局属性中或作为 layer
        # 为了简单起见，我们这里主要依赖 CSV，但 loom 也会生成
        ds.ra['RegulonsAUC'] = auc_mtx.T.values 
        
    print(f"Loom file saved to {F_OUTPUT_LOOM}")

def load_tf_names(fname):
    with open(fname, "r") as f:
        tfs = [line.strip() for line in f.readlines()]
    return tfs

if __name__ == "__main__":
    run()
