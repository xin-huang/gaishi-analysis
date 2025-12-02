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


rule render_gaishi_config_template:
    input:
        tsv="config/msprime_simulation_params.tsv",
        template="config/gaishi.config.template.yaml",
    output:
        config="results/gaishi/model_{model_id}/gaishi.model_{model_id}.config.yaml",
    params:
        model=get_model_params,
        ref_id="Reference",
        tgt_id="Target",
        src_id="Source",
        output_prefix="gaishi.lr",
        output_dir="results/gaishi/model_{model_id}",
        nfeature=10000,
        nprocess=16,
        seedmsprime=4836,
        vcf_file=rules.extract_biallelic_snps.output.vcf,
        ref_ind_file=rules.run_msprime_simulation.output.ref_list,
        tgt_ind_file=rules.run_msprime_simulation.output.tgt_list,
    template_engine:
        "yte"


rule run_gaishi_train:
    input:
        demes="config/ArchIE_3D19.yaml",
        config=rules.render_gaishi_config_template.output.config,   
    output:
        model="results/gaishi/model_{model_id}/gaishi.trained.model",
    resources:
        cpus=16,
    conda:
        "../envs/gaishi.yaml",
    shell:
        """
        gaishi train --demes {input.demes} --config {input.config} --output {output.model}
        """


rule run_gaishi_infer:
    input:
        model=rules.run_gaishi_train.output.model,
        config=rules.render_gaishi_config_template.output.config,
        vcf=rules.extract_biallelic_snps.output.vcf,
    output:
        pred="results/gaishi/model_{model_id}/gaishi.pred.tsv",
    resources:
        mem_gb=16, cpus=2,
    conda:
        "../envs/gaishi.yaml",       
    shell:
        """
        gaishi infer --model {input.model} --config {input.config} --output {output.pred}
        """
