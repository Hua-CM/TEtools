# -*- coding: utf-8 -*-
# @Time    : 2022/7/26 10:26
# @Author  : Zhongyi Hua
# @File    : Distribution.py
# @Usage   :
# @Note    : 
# @E-mail  : njbxhzy@hotmail.com

from BCBio import GFF
from collections import defaultdict
from copy import copy
from Bio.SeqFeature import FeatureLocation
import pandas as pd


def gen_region(gff_path):
    gff = [_ for _ in GFF.parse(gff_path)]
    scaf_region_dict = defaultdict(defaultdict)
    for scaffold in gff:
        region_dict = defaultdict()
        region_dict[0] = 'inter'
        # construct a dummy head and a dummy_end
        dummy_head = copy(scaffold.features[0])
        dummy_head.type = 'gene'
        dummy_head.location = FeatureLocation(0, 0)
        dummy_end = copy(scaffold.features[0])
        dummy_end.type = 'gene'
        dummy_end.location = FeatureLocation(len(scaffold), len(scaffold))

        tmp_lst = [dummy_head] + scaffold.features + [dummy_end]
        for idx in range(1, len(tmp_lst) - 1):
            gene = tmp_lst[idx]
            if gene.type == 'gene':
                gene_start = int(gene.location.start)
                gene_end = int(gene.location.end)
                next_gene_start = int(tmp_lst[idx+1].location.start)
                pre_gene_end = int(tmp_lst[idx-1].location.end)

                region_dict[gene_start] = 'inner'
                region_dict[gene_end] = 'ud2k'
                if (gene_end + 2000) < next_gene_start:
                    region_dict[gene_end + 2000] = 'ud5k'
                    if (gene_end + 5000) < next_gene_start:
                        region_dict[gene_end + 5000] = 'inter'
                if (gene_start - 2000) > pre_gene_end:
                    region_dict[gene_start - 2000] = 'ud2k'
                    if (gene_start - 5000) > pre_gene_end:
                        region_dict[gene_start - 5000] = 'ud5k'
        # add 0 position
        min_pos = min(list(region_dict.keys()))
        if min_pos > 5000:
            region_dict[0] = 'inter'
        elif min_pos > 2000:
            region_dict[0] = 'ud5k'
        else:
            region_dict[0] = 'ud2k'

        scaf_region_dict[scaffold.id] = region_dict
    return scaf_region_dict


def determine_type(scaffold, start: int, end: int, scaf_region_dict):
    start, end = int(start), int(end)
    region_dict = scaf_region_dict.get(scaffold)
    type_count = defaultdict(int)
    if not region_dict:
        type_count['inter'] = end - start + 1
        return type_count
    pre_start = max(_start for _start in sorted(list(region_dict.keys())) if start - _start > 0)

    _type = region_dict.get(pre_start)
    l_p = start
    r_p = start
    while r_p < end:
        r_p += 1
        if r_p in region_dict:
            type_count[_type] += (r_p - l_p) + 1
            _type = region_dict.get(r_p)
            l_p = r_p
    # The last one
    type_count[_type] += (end - l_p)
    return type_count


def parse_args():
    import argparse
    cmdparser = argparse.ArgumentParser(prog="TE distribution determination")
    cmdparser.add_argument('-i', '--input', required=True,
                           help='RM_TRIPs.R output')
    cmdparser.add_argument('-g', '--gff', required=True,
                           help='Genome gff file')
    cmdparser.add_argument('-o', '--output', required=True,
                           help='Determination result')
    args = cmdparser.parse_args()
    return args


def main(args):
    gff_dict = gen_region(args.gff)
    te_df = pd.read_table(args.input)
    distribution = te_df.\
        apply(lambda x: determine_type(x['qry_id'], x['merged_qrystart'], x['merged_qryend'], gff_dict), axis=1)
    te_types = [max(_, key=_.get) for _ in distribution.to_list()]
    te_df2 = pd.DataFrame(distribution.to_list()).fillna(0)
    te_df2['region'] = te_types
    te_df = pd.concat([te_df, te_df2], axis=1)
    te_df.to_csv(args.output, sep='\t', index=False)


if __name__ == '__main__':
    main(parse_args())
