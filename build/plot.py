import csv
import matplotlib.pyplot as plt

xs_hi, xs_lo = [], []
with open('result.csv') as f:
    r = csv.DictReader(f)
    for row in r:
        xs_hi.append(int(row['hi']))
        xs_lo.append(int(row['lo']))

plt.figure()
plt.scatter(xs_hi, xs_lo)
plt.title('Estimated Intersection (by tags)')
plt.xlabel('hi'); plt.ylabel('lo')
plt.grid(True)
plt.tight_layout()
plt.savefig("result.png", dpi=160)
print("Saved to result.png")

