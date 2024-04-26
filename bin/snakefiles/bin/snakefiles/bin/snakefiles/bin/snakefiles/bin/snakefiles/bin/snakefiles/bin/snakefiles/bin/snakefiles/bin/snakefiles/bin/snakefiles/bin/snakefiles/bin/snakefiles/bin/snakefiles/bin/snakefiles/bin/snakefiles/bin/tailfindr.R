#!/usr/bin/env R
# coding=utf-8

#' TailfindR script for SnakeMake rule "tailfindr".
#' 
#' Runs TailfindR, using the function find_tails which writes a
#' polya_tails.csv file for each .fast5 directories provided with an
#' option to plot the raw data traces for each reads.
#'
#' The script is called by the SnakeMake rule "tailfindr" and is only
#' meant to be used in this context.

library(tailfindr)

find_tails(fast5_dir = dirname(snakemake@input[["fast5_file"]]),
           save_dir = dirname(snakemake@output[["tails_out"]]),
           csv_filename = 'polya_tails.csv',
           basecall_group = "Basecall_1D_001",
           num_cores = snakemake@params[["num_cores_per_jobs"]],
           save_plots = snakemake@params[["tail_plots"]],
           plot_debug_traces = snakemake@params[["debug_traces"]],
           plotting_library = 'rbokeh')
