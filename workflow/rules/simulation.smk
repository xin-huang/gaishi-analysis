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


np.random.seed(3821)
seed_list = np.random.randint(1, 2**31, n_rep)


rule run_msprime_simulation:
    input:
        tsv="config/msprime_simulation_params.tsv",
        demes="config/ArchIE_3D19.yaml",
    output:
        ts="results/simulation/model_{model_id}/rep_{rep}/model_{model_id}.rep_{rep}.ts",
        vcf="results/simulation/model_{model_id}/rep_{rep}/model_{model_id}.rep_{rep}.vcf",
        bed="results/simulation/model_{model_id}/rep_{rep}/model_{model_id}.rep_{rep}.true.tracts.bed",
        ref_list="results/simulation/model_{model_id}/rep_{rep}/model_{model_id}.rep_{rep}.ref.list",
        tgt_list="results/simulation/model_{model_id}/rep_{rep}/model_{model_id}.rep_{rep}.tgt.list",
        src_list="results/simulation/model_{model_id}/rep_{rep}/model_{model_id}.rep_{rep}.src.list",
        seed_file="results/simulation/model_{model_id}/rep_{rep}/model_{model_id}.rep_{rep}.seedmsprime",
    params:
        model=get_model_params,
        ref_id="Reference",
        tgt_id="Target",
        src_id="Source",
        seed=lambda wildcards: seed_list[int(wildcards.rep)],
    resources:
        mem_gb=16,
    script:
        "../scripts/msprime_simulation.py"


rule extract_biallelic_snps:
    input:
        vcf=rules.run_msprime_simulation.output.vcf,
    output:
        vcf="results/simulation/model_{model_id}/rep_{rep}/model_{model_id}.rep_{rep}.biallelic.snps.vcf.gz",
    shell:
        """
        bcftools view {input.vcf} -v snps -m 2 -M 2 -g ^miss | bgzip -c > {output.vcf}
        tabix -p vcf {output.vcf}
        """
