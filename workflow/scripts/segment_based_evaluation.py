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
import pyranges as pr


def calc_pr(
    ntrue_tracts: int, ninferred_tracts: int, ntrue_positives: int
) -> tuple[float, float]:
    """
    Calculate segment-based precision and recall.

    Parameters
    ----------
    ntrue_tracts : int
        Length of true introgressed fragments.
    ninferred_tracts : int
        Length of inferred introgressed fragments.
    ntrue_positives : int
        Length of fragments belonging to true positives.

    Returns
    -------
    precision : float
        Estimated precision.
    recall : float
        Estimated recall.
    """
    if float(ninferred_tracts) == 0:
        precision = np.nan
    else:
        precision = ntrue_positives / float(ninferred_tracts)
    if float(ntrue_tracts) == 0:
        recall = np.nan
    else:
        recall = ntrue_positives / float(ntrue_tracts)

    return precision, recall


def evaluate(
    true_tract_file: str, inferred_tract_file: str, cutoff: float, output: str
) -> None:
    """
    Evaluate model performance using segment-based precision and recall.

    Parameters
    ----------
    true_tract_file : str
        Path to the file containing true introgressed fragments.
    inferred_tract_file : str
        Path to the file containing inferred introgressed fragments.
    cutoff : float
        Probability threshold used to classify a fragment as introgressed.
    output : str
        Path to the output file where model performance metrics
        (e.g., precision and recall) will be written.
    """
    try:
        true_tracts = pd.read_csv(
            true_tract_file,
            sep="\t",
            header=None,
            names=["Chromosome", "Start", "End", "Sample"],
        )
    except pd.errors.EmptyDataError:
        true_tracts_samples = []
    else:
        true_tracts_samples = true_tracts["Sample"].unique()
        true_tracts = pr.PyRanges(true_tracts).merge(by="Sample")

    try:
        inferred_tracts = pd.read_csv(
            inferred_tract_file,
            sep="\t",
            header=None,
            names=["Chromosome", "Start", "End", "Sample"],
        )
    except pd.errors.EmptyDataError:
        inferred_tracts_samples = []
    else:
        inferred_tracts_samples = inferred_tracts["Sample"].unique()
        inferred_tracts = pr.PyRanges(inferred_tracts).merge(by="Sample")

    res = pd.DataFrame(
        columns=[
            "Cutoff",
            "Precision",
            "Recall",
            "True_tracts_length",
            "Inferred_tracts_length",
            "Overlapped_length",
        ]
    )

    sum_ntrue_tracts = 0
    sum_ninferred_tracts = 0
    sum_ntrue_positives = 0

    for s in np.intersect1d(true_tracts_samples, inferred_tracts_samples):
        ind_true_tracts = true_tracts[true_tracts.Sample == s].merge(by="Sample")
        ind_inferred_tracts = inferred_tracts[inferred_tracts.Sample == s].merge(
            by="Sample"
        )
        ind_overlaps = ind_true_tracts.intersect(ind_inferred_tracts)

        ntrue_tracts = np.sum([x[1].End - x[1].Start for x in ind_true_tracts])
        ninferred_tracts = np.sum([x[1].End - x[1].Start for x in ind_inferred_tracts])
        ntrue_positives = np.sum([x[1].End - x[1].Start for x in ind_overlaps])

        #precision, recall = cal_pr(ntrue_tracts, ninferred_tracts, ntrue_positives)
        sum_ntrue_tracts += ntrue_tracts
        sum_ninferred_tracts += ninferred_tracts
        sum_ntrue_positives += ntrue_positives

    for s in np.setdiff1d(true_tracts_samples, inferred_tracts_samples):
        # ninferred_tracts = 0
        ind_true_tracts = true_tracts[true_tracts.Sample == s]

        ntrue_tracts = np.sum([x[1].End - x[1].Start for x in ind_true_tracts])
        sum_ntrue_tracts += ntrue_tracts

    for s in np.setdiff1d(inferred_tracts_samples, true_tracts_samples):
        # ntrue_tracts = 0
        ind_inferred_tracts = inferred_tracts[inferred_tracts.Sample == s]

        ninferred_tracts = np.sum([x[1].End - x[1].Start for x in ind_inferred_tracts])
        sum_ninferred_tracts += ninferred_tracts

    total_precision, total_recall = calc_pr(
        sum_ntrue_tracts, sum_ninferred_tracts, sum_ntrue_positives
    )
    res.loc[len(res.index)] = [
        cutoff,
        total_precision,
        total_recall,
        sum_ntrue_tracts,
        sum_ninferred_tracts,
        sum_ntrue_positives,
    ]

    res.fillna("NaN").to_csv(output, sep="\t", index=False)


evaluate(
    true_tract_file=snakemake.input.true_tracts,
    inferred_tract_file=snakemake.input.inferred_tracts,
    cutoff=snakemake.params.cutoff,
    output=snakemake.output.tsv,
)
