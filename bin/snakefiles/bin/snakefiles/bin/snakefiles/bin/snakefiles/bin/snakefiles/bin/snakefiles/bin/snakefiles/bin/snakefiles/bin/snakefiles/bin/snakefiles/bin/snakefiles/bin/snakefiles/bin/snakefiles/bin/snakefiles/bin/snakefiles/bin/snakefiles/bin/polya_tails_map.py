#!/usr/bin/env python3
# coding=utf-8

"""
Match the ONT reads and their polyA tail length to the gene expression
data. The result is written in a table and displayed in a plotly
scatter plot.

Outputs are a .tsv and an .html files.

The script is called by the SnakeMake rule "match_read_polya_tpm" and is
only meant to be used in this context. A custom error is raised upon
inappropriate usage of this script.
"""

__author__ = "Jean-Michel Garant"
__copyright__ = "Copyright (C) 2020, " + __author__
__email__ = "jmgarant@bcgsc.ca"
__license__ = "GPLv3"

import collections          as col

import pandas               as pd
import plotly.graph_objects as go
import plotly.express       as px

def join_exp_tails_on_map(
        df_exp,
        df_tails,
        df_map
        ):
    """
    Merge the dataframes containing gene expression data, polyA tail
    lengths of reads by joining using an intermediate dataframe mapping
    reads to genes.
    :return: a pandas dataframe
    """
    df_map["transcript"], df_map["gene"] = df_map["transcript_gene"].str.split(
        "_", 1).str
    df_map["reads"] = df_map["reads"].str.split(",")
    df_map = df_map.drop("transcript_gene", axis=1).explode("reads")
    merged_df = df_tails.merge(
        df_map,
        left_on="read_id",
        right_on="reads",
        ).merge(
            df_exp,
            left_on="gene",
            right_on="Gene"
            )
    group_gene = merged_df.groupby([
        "gene",
        "Gene",
        "TPM",
        "Count",
        "IlluminaTPM",
        "TPMGroup",
        ], sort=True)
    return group_gene.agg({
        "samples_per_nt": "describe",
        "tail_length": "describe",
        "transcript": lambda x: list(col.Counter(list(x)).items()),
        }).reset_index().set_index("gene")

def main():
    """
    Top-level function.
    Manages input, write to files.
    """
    df_map = pd.read_csv(
        snakemake.input["read_map"], # pylint: disable=undefined-variable
        sep="\t",
        names=["transcript_gene", "reads"])
    df_tails = pd.read_csv(
        snakemake.input["polya_tails"], # pylint: disable=undefined-variable
        sep=","
        ).drop(["tail_start", "tail_end", "file_path"], 1)
    df_exp = pd.read_csv(
        snakemake.input["gene_expression"], # pylint: disable=undefined-variable
        sep="\t",
        ).drop(columns=[
            #"Count",
            "RPKM",
            "logCount",
            "logRPKM",
            "logTPM",
            "Group",
            "log2IlluminaTPM",
            "log2NanoporeTPM",
            ])

    group_gene = join_exp_tails_on_map(df_exp, df_tails, df_map)

#    group_gene.columns = [
#        '_'.join(col).strip("_") for col in group_gene.columns.values
#        ]
    group_gene.to_csv(
        snakemake.output["gene_to_polya"], # pylint: disable=undefined-variable
        sep="\t",
        )

if __name__ == "__main__" and 'snakemake' in globals():
    main()
else:
    raise Exception('Custom Error: This software is meant to be part of '
                    'a SnakeMake workflow.')
