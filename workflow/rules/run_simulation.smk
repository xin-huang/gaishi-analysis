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


def get_model_params(wildcards, input):
    df = pd.read_csv(input.tsv, sep="\t")
    mid = int(wildcards.model_id)
    model = df.iloc[mid]

    return model.to_dict()


rule run_msprime_simulation:
    input:
        tsv="config/msprime_simulation_params.tsv",
        demes="config/ArchIE_3D19.yaml",
    output:
        ts="results/simulation/model_{model_id}/model_{model_id}.ts",
        vcf="results/simulation/model_{model_id}/model_{model_id}.vcf",
        bed="results/simulation/model_{model_id}/model_{model_id}.true.tracts.bed",
        ref_list="results/simulation/model_{model_id}/model_{model_id}.ref.list",
        tgt_list="results/simulation/model_{model_id}/model_{model_id}.tgt.list",
        src_list="results/simulation/model_{model_id}/model_{model_id}.src.list",
        #seed_file="results/simulation/model_{model_id}/model_{model_id}.seedmsprime",
    params:
        model=get_model_params,
        ref_id="Reference",
        tgt_id="Target",
        src_id="Source",
        reps=1,
        seed=3821,
    script:
        "../scripts/msprime_simulation.py"
