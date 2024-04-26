#!/usr/bin/env python3
# coding=utf-8

"""
Summarize the output of TailfindR jobs using plotly with a combination
of binned histogram and boxplot.

Output is an .html file.

The script is called by the SnakeMake rule "summary_tailfindr" and is
only meant to be used in this context. A custom error is raised upon
inappropriate usage of this script.
"""

__author__ = "Jean-Michel Garant"
__copyright__ = "Copyright (C) 2020, " + __author__
__email__ = "jmgarant@bcgsc.ca"
__license__ = "GPLv3"

import pandas         as pd
import plotly.express as px

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

def custom_plot(dataf):
                # header,
    """
    Generate a plotly express binned histogram with marginal box plot on
    x axis.
    :return: a plotly figure object
    """
    # df.columns = header
    out_fig = px.histogram(dataf,
                           x="tail_length",
                           nbins=908,
                           marginal="box",
                           template="none",
                           )
    out_fig.update_layout(xaxis_title="PolyA_tail_length")
    return out_fig

def main():
    """
    Top-level function.
    Manages input, write to files.
    """
    dataf = pd.read_csv(snakemake.input["polya_csv"], sep=",") # pylint: disable=undefined-variable
    out_fig = custom_plot(dataf)
        # snakemake.params["polya_csv_header"],
    out_fig.write_html(snakemake.output["fig_html"]) # pylint: disable=undefined-variable

if __name__ == "__main__" and "snakemake" in globals():
    main()
else:
    raise Exception('Custom Error: This software is meant to be part of'
                    'a SnakeMake workflow.')
