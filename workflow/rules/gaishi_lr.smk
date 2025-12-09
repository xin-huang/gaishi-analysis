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


np.random.seed(4836)
seed_list = np.random.randint(1, 2**31, n_rep)

cutoffs = np.arange(0.1, 1.0, 0.1).round(1)
cutoffs = np.concatenate([np.array([0.001, 0.01]), cutoffs, np.array([0.99, 0.999])])


rule render_lr_config_template:
    input:
        tsv="config/msprime_simulation_params.tsv",
        template="config/gaishi.config.template.yaml",
    output:
        config="results/gaishi/model_{model_id}/lr/rep_{rep}/gaishi.model_{model_id}.rep_{rep}.config.yaml",
    params:
        model=get_model_params,
        ref_id="Reference",
        tgt_id="Target",
        src_id="Source",
        train_output_prefix="gaishi.lr.train",
        infer_output_prefix="gaishi.lr.infer",
        output_dir="results/gaishi/model_{model_id}/lr/rep_{rep}",
        nfeature=1000000,
        nprocess=16,
        seedmsprime=lambda wildcards: int(seed_list[int(wildcards.rep)]),
        vcf_file=rules.extract_biallelic_snps.output.vcf,
        ref_ind_file=rules.run_msprime_simulation.output.ref_list,
        tgt_ind_file=rules.run_msprime_simulation.output.tgt_list,
    template_engine:
        "yte"


rule run_lr_train:
    input:
        demes="config/ArchIE_3D19.yaml",
        config=rules.render_lr_config_template.output.config,   
    output:
        model="results/gaishi/model_{model_id}/lr/rep_{rep}/gaishi.trained.model",
    resources:
        mem_gb=64, cpus=16,
    conda:
        "../envs/gaishi.yaml",
    shell:
        """
        gaishi train --demes {input.demes} --config {input.config} --output {output.model}
        """


rule run_lr_infer:
    input:
        model=rules.run_lr_train.output.model,
        config=rules.render_lr_config_template.output.config,
        vcf=rules.extract_biallelic_snps.output.vcf,
    output:
        pred="results/gaishi/model_{model_id}/lr/rep_{rep}/gaishi.pred.tsv",
    resources:
        mem_gb=256, cpus=16,
    conda:
        "../envs/gaishi.yaml",       
    shell:
        """
        gaishi infer --model {input.model} --config {input.config} --output {output.pred}
        """


rule get_inferred_tracts:
    input:
        pred=rules.run_lr_infer.output.pred,
    output:
        inferred_tracts=temp("results/gaishi/model_{model_id}/lr/rep_{rep}/gaishi.pred.cutoff_{cutoff}.inferred.tracts.bed"),
    params:
        cutoff="{cutoff}",
    shell:
        r"""
        awk -v cutoff={params.cutoff} -v OFS="\t" '$6>cutoff{{print $1,$2,$3,$4}}' {input.pred} | sed '1d' > {output.inferred_tracts}
        """


rule evaluate_lr_segment_based:
    input:
        true_tracts=rules.run_msprime_simulation.output.bed,
        inferred_tracts=rules.get_inferred_tracts.output.inferred_tracts,
    output:
        tsv=temp("results/gaishi/model_{model_id}/lr/rep_{rep}/gaishi.pred.cutoff_{cutoff}.performance.tsv"),
    params:
        cutoff="{cutoff}",
        seq_len=200_000_000,
    script:
        "../scripts/segment_based_evaluation.py"


rule summarize_performance:
    input:
        perf=expand(
            "results/gaishi/model_{model_id}/lr/rep_{rep}/gaishi.pred.cutoff_{cutoff}.performance.tsv",
            cutoff=cutoffs,
            allow_missing=True,
        ),
    output:
        perf=temp("results/gaishi/model_{model_id}/lr/rep_{rep}/gaishi.pred.performance.tsv"),
    shell:
        r"""
        cat {input.perf} | grep -v Cutoff | awk -v rep={wildcards.rep} '{{print rep"\t"$0}}' > {output.perf}
        sed -i '1iReplicate\tCutoff\tPrecision\tRecall\tL_TT_sample\tL_IT_sample\tL_TP_sample\tL_FP_sample\tL_TN_sample\tL_FN_sample' {output.perf}
        """


rule collect_replicate_performance:
    input:
        perf=expand(
            "results/gaishi/model_{model_id}/lr/rep_{rep}/gaishi.pred.performance.tsv",
            rep=range(n_rep),
            allow_missing=True,
        ),
    output:
        perf="results/gaishi/model_{model_id}/lr/gaishi.model_{model_id}.lr.pred.performance.tsv",
    shell:
        r"""
        cat {input.perf} | grep -v Cutoff | awk -v id={wildcards.model_id} '{{print "Model_"id"\t"$0}}' > {output.perf}
        sed -i '1iModel\tReplicate\tCutoff\tPrecision\tRecall\tL_TT_sample\tL_IT_sample\tL_TP_sample\tL_FP_sample\tL_TN_sample\tL_FN_sample' {output.perf}
        """


rule plot_confusion_matrix:
    input:
        perf=rules.collect_replicate_performance.output.perf,
    output:
        plot="results/gaishi/model_{model_id}/plots/gaishi.model_{model_id}.lr.pred.cutoff_{cutoff}.confusion.matrix.png",
    script:
       "../scripts/plot_confusion_matrix.py"


rule plot_pr_curve:
    input:
        perf=rules.collect_replicate_performance.output.perf,
    output:
        plot="results/gaishi/model_{model_id}/plots/gaishi.model_{model_id}.lr.pred.pr.curve.png",
    script:
        "../scripts/plot_pr_curve.py"
