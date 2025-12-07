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


import pandas as pd
import matplotlib.pyplot as plt

perf_files = list(snakemake.input.performance)

dfs = [pd.read_csv(f, sep="\t") for f in perf_files]
full_df = pd.concat(dfs, ignore_index=True)

cols = [
    "Replicate",
    "Cutoff",
    "Precision",
    "Recall",
    "True_tracts_length",
    "Inferred_tracts_length",
    "Overlapped_length",
]
full_df = full_df[cols]

pr_df = (
    full_df.groupby("Cutoff", as_index=False)
    .agg(
        Precision_mean=("Precision", "mean"),
        Precision_std=("Precision", "std"),
        Recall_mean=("Recall", "mean"),
        Recall_std=("Recall", "std"),
    )
    .sort_values("Cutoff")
)

plt.figure(figsize=(5, 5), dpi=300)

plt.errorbar(
    pr_df["Recall_mean"],
    pr_df["Precision_mean"],
    xerr=pr_df["Recall_std"],
    yerr=pr_df["Precision_std"],
    fmt="o-",
)

plt.xlim([0,1])
plt.ylim([0,1])
plt.xlabel("Recall")
plt.ylabel("Precision")
plt.title("Precision–Recall Curve")
plt.grid(True)
plt.tight_layout()
plt.savefig(snakemake.output.plot, bbox_inches="tight")
plt.close()
