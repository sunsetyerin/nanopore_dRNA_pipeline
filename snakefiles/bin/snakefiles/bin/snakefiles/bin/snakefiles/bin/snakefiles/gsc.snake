#!/usr/bin python3
# coding=utf-8

__author__ = "Jean-Michel Garant"
__copyright__ = "Copyright (C) 2020, " + __author__
__email__ = "jmgarant@bcgsc.ca"
__license__ = "GPLv3"

include: "snakefiles/comparative.snake"

configfile: "configs/sample_config.yaml"

rule slurmgpu:
    """
    Jobs submission on dlhosts (guppy).
    """
    input:
        expand(
            os.path.join(
                "data",
                "{sample}",
                "guppy_basecall",
                "sequencing_summary.txt"),
            sample=list(config["samples"].keys())),

rule numbers:
    """
    Jobs submission on numbers (tailfindr).
    """
    input:
        expand(
            os.path.join(
                "data",
                "{sample}",
                "tailfindr",
                "polya_tails.csv"),
            sample=list(config["samples"].keys())),

rule general:
    """
    Set of jobs to be run on any cpu ressources.
    """
    input:
        expand(os.path.join(
            "data",
            "{sample}",
            "tailfindr",
            "plots",
            "polya_tails_len_bins.html"
            ), sample=list(config["samples"].keys())),
        expand(os.path.join(
            "data",
            "{sample}",
            "flair",
            "pass_catenated.isoform.read.map.txt",
            ), sample=list(config["samples"].keys())),
#        expand(os.path.join(
#            "data",
#            "{sample}",
#            "epinano",
#            "pass_catenated.tsv.per.site.var.csv",
#            ), sample=list(config["samples"].keys())),
        expand(os.path.join(
            "data",
            "{sample}",
            "epinano",
            "epinano_summary.txt"
            ), sample=list(config["samples"].keys())),
