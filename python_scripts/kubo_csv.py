import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
from os.path import join as pjoin

# ------------------------- Usage Check ------------------------- #
if len(sys.argv) != 3:
    print("Usage: python plot_kubo_from_txt.py <path_to_txt_folder> <species_name>")
    sys.exit(1)

path = sys.argv[1]
species_name = sys.argv[2]

txt_file = pjoin(path, f"per_particle_kubo_{species_name}.txt")

# ------------------------- Load TXT ---------------------------- #
df = pd.read_csv(txt_file, comment='#', header=None)

df.columns = [
    "particle_index",
    "tau_int",
    "tau_exp",
    "kubo",
    "velocity_mean",
    "omega_b_max",
]

# Convert numerics safely
for col in df.columns:
    df[col] = pd.to_numeric(df[col], errors="coerce")

# --------------------- Basic Stats ----------------------------- #
print("Min Kubo:", df["kubo"].min())
print("Max Kubo:", df["kubo"].max())

# -------------------- Particle with Max Kubo -------------------- #
idx_max = df["kubo"].idxmax()
row_max = df.loc[idx_max]

print("\n=== Particle with Maximum Kubo Number ===")
print(f"Particle Index : {row_max['particle_index']}")
print(f"Kubo Value     : {row_max['kubo']}")
print(f"tau_int        : {row_max['tau_int']}")
print(f"tau_exp        : {row_max['tau_exp']}")

# -------------------- Kubo Threshold Counts -------------------- #
count_ge_1 = (df["kubo"] >= 1).sum()
count_gt_1 = (df["kubo"] > 1).sum()
count_lt_1 = (df["kubo"] < 1).sum()

print("\n=== Kubo Threshold Counts ===")
print(f"Kubo >= 1 : {count_ge_1}")
print(f"Kubo > 1  : {count_gt_1}")
print(f"Kubo < 1  : {count_lt_1}")
print(f"Total     : {len(df)}")

# -------------------- Particles with Kubo > 1 -------------------- #
df_kubo_gt1 = df[df["kubo"] > 1]
df_kubo_gt1 = df[(df["kubo"] > 0.1) & (df["kubo"] < 2)]


print("\n=== Particles with Kubo > 1 (First 10) ===")
if len(df_kubo_gt1) == 0:
    print("No particles have Kubo > 1.")
else:
    print(df_kubo_gt1[["particle_index", "kubo", "tau_int", "tau_exp"]].head(10).to_string(index=False))

# -------------------- Histogram (Kubo < 10) ---------------------- #
kubo_clipped = df["kubo"][df["kubo"] < 10]
#kubo_clipped = df["kubo"]

plt.figure()
plt.hist(kubo_clipped, bins = 200, edgecolor='black')
plt.yscale('log')
plt.xlabel(r"$K$")
plt.ylabel(r"$f(K)$")
plt.title(f"Kubo Number Distribution — {species_name}")
plt.tight_layout()
plt.show()

# -------------------- Velocity vs Kubo Scatter Plot -------------- #
plt.figure()
plt.scatter( df["velocity_mean"], (df["kubo"]), s=10, alpha=0.5)
plt.ylabel(r"$K$")
plt.xlabel(r"$\langle v \rangle$")
plt.title(f"Velocity vs Kubo — {species_name}")
#plt.ticklabel_format(style='sci', axis='y')
# Source - https://stackoverflow.com/a
# Posted by kklocker
# Retrieved 2025-11-24, License - CC BY-SA 4.0
plt.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))


plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()
