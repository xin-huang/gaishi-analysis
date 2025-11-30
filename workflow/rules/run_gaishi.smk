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


rule run_gaishi_train:
    input:
        demes="config/ArchIE_3D19.yaml",
        config="config/gaishi.train.config.yaml",   
    output:
        model="results/gaishi/gaishi.trained.model",
    resources:
        cpus=2,
    conda:
        "../envs/gaishi-env.yaml",
    shell:
        """
        gaishi train --demes {input.demes} --config {input.config} --output {output.model}
        """


rule run_gaishi_infer:
    input:
        model=rules.run_gaishi_train.output.model,
        config="config/gaishi.infer.config.yaml",
    output:
        pred="results/gaishi/gaishi.pred.tsv",
    resources:
        cpus=2,
    conda:
        "../envs/gaishi-env.yaml",       
    shell:
        """
        gaishi infer --model {input.model} --config {input.config} --output {output.pred}
        """
