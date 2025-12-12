# plot_okvs.py
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("okvs_bench.csv")
enc = df[df['phase']=="encode"]

plt.plot(enc['n'], enc['time_ms'], marker='o')
plt.xlabel("n (number of key-value pairs)")
plt.ylabel("Encode time (ms)")
plt.title("RB-OKVS Encode Benchmark")
plt.grid(True)
plt.tight_layout()
plt.savefig("okvs_bench.png", dpi=150)
print("Saved okvs_bench.png")

