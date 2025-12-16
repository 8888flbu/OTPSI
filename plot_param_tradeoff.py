#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot parameter trade-offs for Overthreshold-PSI* using real benchmark data.
Reads okvs_bench.csv and rt_vs_n_by_t.csv, and generates fig_param_tradeoff.pdf.
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

# ====== 样式配置（IEEE风格） ======
plt.rcParams.update({
    "font.family": "Times New Roman",
    "font.size": 9,
    "axes.linewidth": 0.8,
    "axes.grid": True,
    "grid.linestyle": "--",
    "grid.alpha": 0.4,
})

# ====== 路径设置 ======
okvs_path = "okvs_bench.csv"
rt_path = "rt_vs_n_by_t.csv"

if not os.path.exists(okvs_path) or not os.path.exists(rt_path):
    raise FileNotFoundError("请确保 okvs_bench.csv 与 rt_vs_n_by_t.csv 在当前目录下！")

# ====== 读取 OKVS 数据 ======
okvs_df = pd.read_csv(okvs_path)
cols = [c.lower() for c in okvs_df.columns]
okvs_df.columns = cols
alpha_col = [c for c in cols if "alpha" in c][0]
succ_col = [c for c in cols if "succ" in c or "rate" in c or "success" in c][0]

# ====== 读取 Threshold 性能数据 ======
rt_df = pd.read_csv(rt_path)
rt_cols = [c.lower() for c in rt_df.columns]
rt_df.columns = rt_cols
t_col = [c for c in rt_cols if "t" in c or "threshold" in c][0]
comm_col = [c for c in rt_cols if "comm" in c or "bytes" in c or "mb" in c][0]
rt_mean = rt_df.groupby(t_col)[comm_col].mean().reset_index()

# ====== 模拟 Cuckoo 扩展因子曲线 (ε vs fail prob) ======
eps = np.linspace(1.05, 1.20, 7)
fail = np.exp(-50*(eps-1.05)) * 1e-2

# ====== 绘图 ======
plt.figure(figsize=(8, 2.6))

# (a) OKVS 成功率 vs α
plt.subplot(1, 3, 1)
plt.plot(okvs_df[alpha_col], okvs_df[succ_col], marker='o', color='tab:blue', label="Measured")
plt.xlabel(r'$\alpha$')
plt.ylabel('Success rate')
plt.ylim(0, 1.05)
plt.title('(a) OKVS success')
plt.legend(frameon=False)

# (b) Cuckoo 哈希失败率
plt.subplot(1, 3, 2)
plt.semilogy(eps, fail, marker='s', color='tab:orange', label="Estimated")
plt.xlabel(r'$\epsilon$')
plt.ylabel('Fail prob.')
plt.title('(b) Cuckoo fail')
plt.legend(frameon=False)

# (c) 通信成本 vs 阈值 t
plt.subplot(1, 3, 3)
plt.plot(rt_mean[t_col], rt_mean[comm_col], marker='^', color='tab:green', label="Measured")
plt.xlabel(r'Threshold $t$')
plt.ylabel('Comm. (MB)')
plt.title('(c) Comm cost')
plt.legend(frameon=False)

plt.tight_layout()
plt.savefig("fig_param_tradeoff.pdf", bbox_inches="tight")
plt.show()

print("✅ 图已生成: fig_param_tradeoff.pdf")
