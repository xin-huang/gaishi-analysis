# Copyright 2025 Xin Huang
#
# GNU General Public License v3.0
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, please see
#
#    https://www.gnu.org/licenses/gpl-3.0.en.html


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


df = pd.read_csv(snakemake.input.perf, sep="\t")
model_name = f"Model_{snakemake.wildcards.model_id}"
cutoff = float(snakemake.wildcards.cutoff)

cols = [
    "Model",
    "Replicate",
    "Cutoff",
    "Precision",
    "Recall",
    "L_TT_sample",
    "L_IT_sample",
    "L_TP_sample",
    "L_FP_sample",
    "L_TN_sample",
    "L_FN_sample",
]
df = df[cols]

group = df.groupby(["Model", "Cutoff"], as_index=False).mean()
row = group.query("Model == @model_name and Cutoff == @cutoff").iloc[0]

tp = row["L_TP_sample"]
fp = row["L_FP_sample"]
tn = row["L_TN_sample"]
fn = row["L_FN_sample"]

cm = np.array([
    [tn, fp],
    [fn, tp],
])

cm_pct = cm / cm.sum(axis=1, keepdims=True) * 100

fig, ax = plt.subplots()
im = ax.imshow(cm_pct, cmap="Blues")

for i in range(cm_pct.shape[0]):
    for j in range(cm_pct.shape[1]):
        ax.text(j, i, f"{cm_pct[i, j]:.1f}%", ha="center", va="center")

ax.set_xticks([0, 1])
ax.set_yticks([0, 1])
ax.set_xticklabels(["Non-introgressed", "Introgressed"])
ax.set_yticklabels(["Non-introgressed", "Introgressed"])

ax.set_xlabel("Predicted")
ax.set_ylabel("True")
ax.set_title(f"Model={model_name} | Cutoff={cutoff}")

plt.tight_layout()
plt.savefig(snakemake.output.plot, bbox_inches="tight")
plt.close()
