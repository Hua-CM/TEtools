# -*- coding: utf-8 -*-
# @Time    : 2022/7/29 12:43
# @Author  : Zhongyi Hua
# @File    : solo_detect.py
# @Usage   :
# @Note    : 
# @E-mail  : njbxhzy@hotmail.com

import pandas as pd
from Bio.Blast.Applications import NcbiblastnCommandline, NcbimakeblastdbCommandline, NcbitblastnCommandline
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from tempfile import mkdtemp
import os
from collections import defaultdict
import re


def del_directory(_dir):
    for r_, _d_, f_ in os.walk(_dir):
        for files in f_:
            os.remove(os.path.join(r_, files))
        os.removedirs(r_)


def get_seq(_seq, _start, _end):
    """
    :param _seq:
    :param _start: !!! should be LTR start position on genome (begin from 1, not 0)
    :param _end: !!! should be LTR end position on genome (begin from 1, not 0)
    :return:
    """
    if len(_seq.seq) < 6000:
    # filter too short scaffolds
        return None
    elif _start <= 3000:
        return _seq.seq[0: _end+3000]
    elif _end >= len(_seq.seq)-3000:
        return _seq.seq[_start-3001:]
    else:
        return _seq.seq[_start-3001: _end+3000]


def parse_args():
    import argparse
    cmdparser = argparse.ArgumentParser(prog="Solo LTR determination")
    cmdparser.add_argument('-g', '--genome', required=True,
                           help='<file path> Genome sequence fasta')
    cmdparser.add_argument('-i', '--intact', required=True,
                           help='<file path> intact LTR fasta')
    cmdparser.add_argument('--gag', required=True,
                           help='<file path> Gag-pal protein sequence fasta')
    cmdparser.add_argument('-t', '--threads', default=4, type=int,
                           help='<int> BLAST threads number. Default: 4')
    cmdparser.add_argument('-o', '--output', required=True,
                           help='<file path> output file name. (Generate Automately)')
    args = cmdparser.parse_args()
    return args


def main(args):
    tmp_dir = mkdtemp(dir='.')
    args.genome = os.path.abspath(args.genome)
    args.intact = os.path.abspath(args.intact)
    args.gag = os.path.abspath(args.gag)
    args.output = os.path.abspath(args.output)

    os.symlink(args.genome, os.path.join(tmp_dir, 'genome.fasta'))
    mkdb1_cmd1 = NcbimakeblastdbCommandline(input_file=os.path.join(tmp_dir, 'genome.fasta'),
                                            dbtype='nucl')
    mkdb1_cmd1()
    blast_cmd1 = NcbiblastnCommandline(query=args.intact,
                                       db=os.path.join(tmp_dir, 'genome.fasta'),
                                       evalue=1e-10,
                                       out=os.path.join(tmp_dir, 'blastn.out'),
                                       outfmt='6 qaccver saccver pident length qlen qstart qend sstart send sstrand evalue',
                                       max_hsps=1,
                                       num_threads=args.threads
                                       )
    blast_cmd1()
    # tidy1
    blastn_res = pd.read_table(os.path.join(tmp_dir, 'blastn.out'), header=None,
                  names=['qaccver', 'saccver', 'pident', 'length', 'qlen', 'qstart',
                         'qend', 'sstart', 'send', 'sstrand', 'evalue'])
    blastn_res = blastn_res[(blastn_res['pident'] > 90) & (blastn_res['pident'] < 100) & (blastn_res['length']/blastn_res['qlen'] > 0.9)]
    # remove duplicated sequences
    blastn_res = blastn_res[['saccver', 'sstart', 'send']].drop_duplicates()
    blastn_res = blastn_res.sort_values(by=['saccver', 'sstart', 'send'])
    putative_lst = [['dummy',0,0]]
    for _idx, _row in blastn_res.iterrows():
        _start, _end = (_row['sstart'], _row['send']) if _row['sstart'] < _row['send'] else (_row['send'], _row['sstart'])
        pre_ele = putative_lst[-1]
        if not pre_ele[0] == _row['saccver']:
            putative_lst.append([_row['saccver'], _start, _end])
        else:
            if _start < pre_ele[2]:
                if _end > pre_ele[2]:
                    pre_ele[2] = _row['send']
            else:
                putative_lst.append([_row['saccver'], _start, _end])
    putative_lst = putative_lst[1:]
    # get sequence
    out_seqs = []
    seq_dicts = SeqIO.to_dict(SeqIO.parse(os.path.join(tmp_dir, 'genome.fasta'), 'fasta'))
    for ltr in putative_lst:
        _seq = seq_dicts.get(ltr[0])
        _short_seq = get_seq(_seq, ltr[1], ltr[2])
        if _short_seq:
            short_seq_record = SeqRecord(_short_seq, id='_'.join([ltr[0], str(ltr[1]), str(ltr[2])]), description='')
            out_seqs.append(short_seq_record)
    SeqIO.write(out_seqs, os.path.join(tmp_dir, 'putative_LTR.fasta'), 'fasta')
    # Blast Gag-pol
    mkdb1_cmd2 = NcbimakeblastdbCommandline(input_file=os.path.join(tmp_dir, 'putative_LTR.fasta'),
                                            dbtype='nucl')
    mkdb1_cmd2()
    blast_cmd2 = NcbitblastnCommandline(query=args.gag,
                                        db=os.path.join(tmp_dir, 'putative_LTR.fasta'),
                                        evalue=1e-8,
                                        out=os.path.join(tmp_dir, 'tblastn.out'),
                                        outfmt='6 qaccver saccver pident length qlen qstart qend sstart send sstrand evalue',
                                        max_hsps=1,
                                        num_threads=args.threads
                                        )
    blast_cmd2()
    # tidy2
    tblastn_res = pd.read_table(os.path.join(tmp_dir, 'tblastn.out'),
                                header=None,
                                names=['qaccver', 'saccver', 'pident', 'length', 'qlen', 'qstart',
                                       'qend', 'sstart', 'send', 'sstrand', 'evalue'])
    tblastn_res = tblastn_res[(tblastn_res['pident'] > 30) & (tblastn_res['length'] / tblastn_res['qlen'] > 0.5)]
    truncted = list(set(tblastn_res['saccver'].to_list()))
    putative_lst = ['_'.join([ltr[0], str(ltr[1]), str(ltr[2])]) for ltr in putative_lst]
    solo = [ltr for ltr in putative_lst if ltr not in truncted]
    outs = os.path.splitext(args.output)
    with open(outs[0] + '_truncted' + outs[1], 'w') as f_out:
        f_out.write('\n'.join(truncted))
    with open(outs[0] + '_solo' + outs[1], 'w') as f_out:
        f_out.write('\n'.join(solo))
    # cluster into family
    os.system(f'cat {os.path.join(tmp_dir, "putative_LTR.fasta")} {args.intact} > {os.path.join(tmp_dir, "all_LTR.fasta")}')
    mkdb1_cmd3 = NcbimakeblastdbCommandline(input_file=os.path.join(tmp_dir, "all_LTR.fasta"),
                                            dbtype='nucl')
    mkdb1_cmd3()
    blast_cmd3 = NcbiblastnCommandline(query=os.path.join(tmp_dir, "all_LTR.fasta"),
                                       db=os.path.join(tmp_dir, "all_LTR.fasta"),
                                       evalue=1e-3,
                                       out=os.path.join(tmp_dir, 'all_LTR.blastout'),
                                       outfmt=6,
                                       num_threads=args.threads
                                       )
    blast_cmd3()
    os.system(f'silix -i 0.6 -r 0.7  {os.path.join(tmp_dir, "all_LTR.fasta")} {os.path.join(tmp_dir, "all_LTR.blastout")} > {os.path.join(tmp_dir, "LTR.clstr")}')
    # output
    ## seqid type cluster
    tmp_dict = defaultdict()
    for _ in solo:
        tmp_dict[_] = 'solo'
    for _ in truncted:
        tmp_dict[_] = 'truncted'

    tmp_lst = []
    with open(os.path.join(tmp_dir, "LTR.clstr")) as f_in:
        for line in f_in.read().strip().split('\n'):
            _clstr, _seq_id = line.split('\t')
            tmp_lst.append({'cluster': _clstr, 'seq_id': _seq_id, 'type': tmp_dict.get(_seq_id, 'intact')})
    family_df = pd.DataFrame(tmp_lst)
    family_df.to_csv(args.output, sep='\t', index=False)
    del_directory(tmp_dir)


if __name__ == '__main__':
    main(parse_args())
