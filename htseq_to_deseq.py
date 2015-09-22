#!/usr/bin/env python3
"""This script will combine multiple htseq-count output files to produce a table
for DESeq

Copyright 2015 Oleksandr Moskalenko <om@rc.ufl.edu>

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

    Unless required by applicable law or agreed to in writing, software
    distributed under the License is distributed on an "AS IS" BASIS,
    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
    See the License for the specific language governing permissions and
    limitations under the License.

Version 1.0.0
"""

import logging, subprocess, re, pickle, os, argparse, sys, csv, pandas

# CONSTANTS:

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--verbose", default=False, action='store_true', help="Verbose)")
    parser.add_argument("-d", "--debug", default=False, action='store_true', help="Debug)")
    parser.add_argument("-p", "--pattern", default="_htseq_count", dest='pattern', help="File name pattern to split on")
    parser.add_argument("-o", "--output", default="output.tsv", dest='outfile', help="Verbose)")
    parser.add_argument("htseqcfiles", type=str, nargs='+', help='htseq-count files to combine')
    args = parser.parse_args()
    return args

def setup_logger(args):
     log = logging.getLogger(__name__)
     ch = logging.StreamHandler()
     formatter = logging.Formatter('%(asctime)s:%(levelname)s - %(message)s')
     ch.setFormatter(formatter)
     if args.debug:
         ch.setLevel(logging.DEBUG)
         log.setLevel('DEBUG')
     else:
         ch.setLevel(logging.INFO)
         log.setLevel('INFO')
     log.addHandler(ch)
     args.log = log


def process_htseq_count_data(args, infile_list):
    combined_dict = {}
    gene_id_list = []
    for infile in infile_list:
        file_dict = {}
        if args.pattern:
            pattern = args.pattern
        else:
            pattern = '_htseq_counts'
        sample_name = os.path.split(infile)[1].split(pattern)[0]
        with open(infile, 'r') as in_fh:
            reader = csv.reader(in_fh, delimiter='\t')
            file_dict = dict(filter(None, reader))
            gene_ids = list(file_dict.keys())
            gene_id_list.extend(gene_ids)
        combined_dict[sample_name] = file_dict
    if args.verbose:
        print("Total number of gene_ids: {0}".format(len(gene_ids)))
    return combined_dict, gene_ids


def combine_counts(args, data, gene_ids):
    """Combine multiple data dictionaries with per-sample htseq-counts"""
    log = args.log
    if args.verbose:
        log.info("Combining samples")
    main_df = pandas.DataFrame(index=gene_ids)
    main_df.index.name = 'gene_id'
    processed_samples = []
    for sample in data:
        if args.debug:
            if processed_samples:
                print("All processed samples: '{0}'".format(", ".join(processed_samples)))
            print("Processing sample '{0}'".format(sample))
        if sample not in processed_samples:
            sample_df = pandas.DataFrame.from_dict(data[sample], orient='index')
            sample_df.index.name = 'gene_id'
            sample_df.columns = [sample]
            processed_samples.append(sample)
        main_df = main_df.merge(sample_df, how='left', left_index=True, right_index=True)
    if args.debug:
        log.debug("Processed samples: {0}".format(", ".join(processed_samples)))
#    print(main_df)
    return main_df


def write_outfile(args):
    out_fh = open(args.outfile, 'w')


def main():
    args = parse_args()
    setup_logger(args)
    log = args.log
    if not args.outfile:
        print("Using the default output filename 'output.tsv'")
    htseq_count_file_list = [x for x in args.htseqcfiles]
    if args.verbose:
        log.info("Reading in the input files: {0}".format(", ".join(htseq_count_file_list)))
    data_dict, gene_ids = process_htseq_count_data(args, htseq_count_file_list)
    combined_df = combine_counts(args, data_dict, gene_ids)
    combined_df.to_csv(args.outfile, sep='\t', na_rep='NA')


if __name__=='__main__':
    main()
