#!/usr/bin/env python3
# coding=utf-8

"""
"""

__author__ = "Jean-Michel Garant"
__copyright__ = "Copyright (C) 2020, " + __author__
__email__ = "jmgarant@bcgsc.ca"
__license__ = "GPLv3"

from collections import OrderedDict

import pandas               as pd
import plotly.graph_objects as go
import plotly.express       as px

def add_group_count_legend(figure_object):
    """
    Modifies the name of traces in a plotly Figure object so the
    content of a legend shows the count of points in each colored
    traces.

    In plotly express scatter, the different traces are defined by the
    groups generated using the 'color' argument on a column of the
    dataframe.
    """
    for data_group in figure_object.data:
        if len(data_group.name.split("=")) >= 2:
            data_group.name = " | ".join(data_group.name.split("=")[1],
                                         len(data_group.customdata))
    return figure_object

def drop_tail_len_box(dataframe,
                  x,
                  y,
                  groups_column, #  unused
                  title,
                  ):
    figs = {}
    for num, key in enumerate(y.keys()):
        print("let's get to it {}, {}".format(num,key))
        figs[key] = px.box(
            dataframe.loc[dataframe[key].notnull()],
            x=x,
            y=key,
            #log_y=True,
#            marginal_x="violin",
#            marginal_y="violin",
#            color=groups_column,
            hover_data=dataframe.columns,
            title=title,
            points=False,
#            render_mode='webgl',
            )
        print("One of the traces was generated")
        if num != 0:
            figs[key].update_traces(visible=False)
    dat = tuple()
    for key in y.keys():
        dat = dat + figs[key]["data"]
    fig = go.Figure(data=dat)
    fig.layout.update(
        figs[list(y.keys())[0]].layout,
#        xaxis_title="Mean Tails Length of Reads in Gene (nts)",
#        yaxis_title="Reads per Millions (RPM)",
        hovermode="closest",
        )
    fig.update_traces(marker_line=dict(width=1, color='Gray'))
    print("color layout changed")
    fig.update_layout(
        template="none",
        updatemenus=[
            go.layout.Updatemenu(
                ##type="buttons",
                buttons=[
                    dict(
                        args=[{"visible": [False]*num*len(figs[key]["data"])
                                          + [True]*len(figs[key]["data"])
                                          + [False]*(len(y)-num-1
                                                     )*len(figs[key]["data"])},
                              {"yaxis.title.text": y[key]},
                              ],
                        label=y[key],
                        method="update",
                        ) for  num, key in enumerate(y)
                    ],
                direction="down",
                pad={"r": 0, "t": 0},
                showactive=True,
                x=0,
                xanchor="right",
                y=1,
                yanchor="bottom"
                ),
            go.layout.Updatemenu(
                type="buttons",
                buttons=[
                    dict(
                        args=[{"yaxis.type": axis_type,
                               "yaxis2.type": axis_type},
                              ],
                        label=axis_type,
                        method="relayout",
                        ) for axis_type in ["log", "linear"]
                    ],
                direction="right",
                pad={"r": 0, "t": 0.02},
                showactive=True,
                x=0,
                xanchor="right",
                y=1,
                yanchor="top"
                ),
            ]
        )
    print("Major layout redefined")
    return fig

def one_violin(dataframe,
                  x,
                  y,
                  title,
                  ):
    fig = px.violin(
        dataframe,
        x=x,
        y="tail_length",
#        log_y=True,
#        marginal_x="violin",
#        marginal_y="violin",
#        color=groups_column,
        title=title,
#        points=True,
#        render_mode='webgl',
        )
    fig.update_traces(marker_line=dict(width=1, color='Gray'))
    fig.update_layout(template="presentation")
    return fig

def drop_tail_len_bar(dataframe,
                      x,
                      groups_column, #  unused
                      title,
                      ):
    figs = {}
    for num, key in enumerate(x.keys()):
        print("let's get to it {}, {}".format(num,key))
        figs[key] = px.histogram(
            dataframe.loc[dataframe[key].notnull()],
            x=key,
#            marginal_x="violin",
#            marginal_y="violin",
            color=groups_column,
#            hover_data=dataframe.columns,
            title=title,
            barmode="overlay",
            )
        print("One of the traces was generated")
        if num != 0:
            figs[key].update_traces(visible=False)
    dat = tuple()
    for key in x.keys():
        dat = dat + figs[key]["data"]
    fig = go.Figure(data=dat)
    fig.layout.update(
        figs[list(x.keys())[0]].layout,
#        xaxis_title="Mean Tails Length of Reads in Gene (nts)",
#        yaxis_title="Reads per Millions (RPM)",
        hovermode="closest",
        )
    print("color layout changed")
    fig.update_layout(
        template="none",
        updatemenus=[
            go.layout.Updatemenu(
                ##type="buttons",
                buttons=[
                    dict(
                        args=[{"visible": [False]*num*len(figs[key]["data"])
                                          + [True]*len(figs[key]["data"])
                                          + [False]*(len(x)-num-1
                                                     )*len(figs[key]["data"])},
                              {"xaxis.title.text": x[key]},
                              ],
                        label=x[key],
                        method="update",
                        ) for  num, key in enumerate(x)
                    ],
                direction="down",
                pad={"r": 0, "t": 0},
                showactive=True,
                x=0,
                xanchor="right",
                y=1,
                yanchor="bottom"
                ),
            go.layout.Updatemenu(
                type="buttons",
                buttons=[
                    dict(
                        args=[{"yaxis.type": axis_type,
                               "yaxis2.type": axis_type},
                              ],
                        label=axis_type,
                        method="relayout",
                        ) for axis_type in ["log", "linear"]
                    ],
                direction="right",
                pad={"r": 0, "t": 0.02},
                showactive=True,
                x=0,
                xanchor="right",
                y=1,
                yanchor="top"
                ),
            ]
        )
    print("Major layout redefined")
    return fig

def one_histogram(dataframe,
                  x,
#                  y,
                  groups_column,
                  title,
                  ):
    fig = px.histogram(
        dataframe,
        x=x,
        barmode="overlay",
#        y="tail_length",
#        log_y=True,
#        marginal_x="violin",
#        marginal_y="violin",
        color=groups_column,
#        title=title,
        range_x=[0, 100],
#        nbins=100,
#        points=True,
#        render_mode='webgl',
        )
#    fig.update_traces(marker_line=dict(width=1, color='Gray'))
    fig.update_layout(template="presentation")
    return fig

def main():
    """
    Top-level function.
    Manages input, write to files.
    """
    main_df = pd.read_csv(
        snakemake.input["mapped"],
        sep="\t",
        index_col=0,
        header=0, # 2 rows columns headers
        )
    # Somehow pandas does not natively support the multirow headers it
    # produces when writing grouped dataframe...

    # Dynamic figure generated using a dataframe and an ordered
    # dictionnary of the axis to use
    all_dict = OrderedDict({col:col for col in main_df.select_dtypes(
        include="number").columns})
    print(main_df)
    print(all_dict)
    if snakemake.output["plotly_fig"].endswith(".html"):
        fig = drop_tail_len_bar(main_df,
                                #snakemake.params["x_axis"],
                                all_dict,
                                snakemake.params["group_column"],
                                "PolyA tails length for "
                                "{} reads".format(len(main_df)),
                                )
        print("fig object generate")
        fig = add_group_count_legend(fig)
        fig.write_html(snakemake.output["plotly_fig"])
    if snakemake.output["plotly_fig"].endswith((".png", ".jpg", ".svg")):
        fig = one_histogram(main_df,
                        snakemake.params["x_axis"],
                        snakemake.params["group_column"],
                        "PolyA tails length"
                        #"for {} reads".format(len(main_df)),
                        )
        fig.write_image(snakemake.output["plotly_fig"])

if __name__ == "__main__" and 'snakemake' in globals():
    main()
else:
    raise Exception('Custom Error: This software is meant to be part of '
                    'a SnakeMake workflow.')
