#!/usr/bin/env python3
# coding=utf-8

"""
Summarize the output of TailfindR jobs from multiple samples using
plotly with a combination of binned histogram and boxplot.

Output is an .html file.

The script is called by the SnakeMake rule "bins_tailfindr" and is only
meant to be used in this context. A custom error is raised upon
inappropriate usage of this script.
"""

__author__ = "Jean-Michel Garant"
__copyright__ = "Copyright (C) 2020, " + __author__
__email__ = "jmgarant@bcgsc.ca"
__license__ = "GPLv3"

import subprocess
import pandas         as pd
import plotly.express as px
import plotly.io      as pio
from scipy.stats import ks_2samp, ttest_ind

def add_group_count_legend(figure_object):
    """
    Modifies the name of traces in a plotly Figure object so the
    content of a legend shows the count of points in each colored
    traces.

    In plotly express scatter, the different traces are defined by the
    groups generated using the 'color' argument on a column of the
    dataframe.
    :return: a plotly figure object
    """
    for data_group in figure_object.data:
        data_group.name = "{} | {}".format(data_group.name.split("=")[1],
                                           str(len(data_group.customdata))
                                           )
    return figure_object

def histo_bin_polya_len(dataf,
                        marginal=False,
                        ):
    """
    Generate a plotly express binned histogram with/without marginal box
    plot on x axis.
    :return: a plotly figure object
    """
    out_fig = px.histogram(dataf,
                           x="tail_length",
                           color="sample",
                           barmode="overlay",
                           histnorm="percent",
                           nbins=454,
                           marginal=marginal,
                           template="presentation",
                           range_x=[0, 400],
                           )
    out_fig.update_layout(xaxis_title="PolyA Tail Length",
                          yaxis_title="Proportion of Counts (%) in the Sample",
                          )
    return out_fig

def full_scattergl(dataf,
                   x_col="tail_length",
                   y_col="read_length",
                   marginal=False,
                   range_x=None,
                   range_y=None,
                   ):
    """
    Generate a plotly express scatterplotgl (accelerated for large sets
    html) with/without marginal box plot on both axis.
    :return: a plotly figure object
    """
    out_fig = px.scatter(dataf,
                         x=x_col,
                         y=y_col,
                         color="sample",
                         marginal_x=marginal,
                         marginal_y=marginal,
                         template="presentation",
                         range_x=range_x,
                         range_y=range_y,
                         render_mode='webgl',
                         )
    return out_fig

def main():
    """
    Top-level function.
    Manages input, write to files.
    """
    # fetch and merge dataframes
    # TODO changed the limit of two samples in comparative
    df_tailfindr = []
    df_guppy = []
    dframes = []
    for (sample_id,
         sample_tailfindr,
         sample_guppy
         ) in zip(snakemake.params["samples_id"][0:2],  # pylint: disable=undefined-variable
               snakemake.input["samples_tailfindr"][0:2],  # pylint: disable=undefined-variable
               snakemake.input["samples_guppy"][0:2]  # pylint: disable=undefined-variable
               ):
        df_tailfindr = pd.read_csv(sample_tailfindr,
                                   sep=",",
                                   )
        df_guppy = pd.read_csv(sample_guppy,
                               sep="\t",
                               usecols=["read_id", "sequence_length_template"],
                               )
        df_sample = pd.merge(left=df_tailfindr,
                             right=df_guppy,
                             how="left",
                             on="read_id",
                             )
        df_sample["sample"] = sample_id
        dframes.append(df_sample)
    # Compare distributions before concat
    print(ks_2samp(dframes[0]["tail_length"],
                   dframes[1]["tail_length"]))
    print(ttest_ind(dframes[0]["tail_length"],
                    dframes[1]["tail_length"],
                    equal_var=False,
                    nan_policy='omit',
                    ))
    dataf = pd.concat(dframes)  # [df_sample_a, df_sample_b])

    # Binned hisogram of polya lengths
    polya_dist_fig = histo_bin_polya_len(dataf,
                                         marginal="box")
    # add_group_count_legend(out_fig)
    polya_dist_fig.write_html(snakemake.output["fig_html"]) # pylint: disable=undefined-variable
    polya_dist_fig = histo_bin_polya_len(dataf,
                                         marginal=False)
    polya_dist_fig.write_image(file=snakemake.output["fig_svg"], # pylint: disable=undefined-variable
                               format="svg")
#    pio._orca.shutdown_server()

    # Comparative of variable TODO recomment after
    fig = full_scattergl(dataf,
                         x_col="tail_length",
                         y_col="samples_per_nt",
                         marginal=False,
                         range_x=[0, 400],
                         range_y=[0, 200],
                         )
    fig.write_html(
        snakemake.output["fig_html"].replace(".html",# pylint: disable=undefined-variable
                                             "_polya_samplesnt.html"))
    fig.write_image(
        file=snakemake.output["fig_svg"].replace(".svg", # pylint: disable=undefined-variable
                                                 "_polya_samplesnt.svg"),
        format="svg")
#    pio._orca.shutdown_server()

    fig = full_scattergl(dataf,
                         x_col="tail_length",
                         y_col="sequence_length_template",
                         marginal=False,
                         range_x=[0, 400],
                         range_y=[0, 10000],
                         )
    fig.write_html(
        snakemake.output["fig_html"].replace(".html",# pylint: disable=undefined-variable
                                             "_polya_readlen.html"))
    fig.write_image(
        file=snakemake.output["fig_svg"].replace(".svg", # pylint: disable=undefined-variable
                                                 "_polya_readlen.svg"),
        format="svg")
#    pio._orca.shutdown_server()

    fig = full_scattergl(dataf,
                         x_col="sequence_length_template",
                         y_col="samples_per_nt",
                         marginal=False,
                         range_x=[0, 10000],
                         range_y=[0, 200],
                         )
    fig.write_html(
        snakemake.output["fig_html"].replace(".html",# pylint: disable=undefined-variable
                                             "_readlen_samplesnt.html"))
    fig.write_image(
        file=snakemake.output["fig_svg"].replace(".svg", # pylint: disable=undefined-variable
                                                 "_readlen_samplesnt.svg"),
        format="svg")

if __name__ == "__main__" and "snakemake" in globals():
    main()
else:
    raise Exception('Custom Error: This software is meant to be part of '
                    'a SnakeMake workflow.')
