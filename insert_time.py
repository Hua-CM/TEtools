# -*- coding: utf-8 -*-
# @Time    : 2022/7/21 18:40
# @Author  : Zhongyi Hua
# @File    : insert_time.py
# @Usage   :
# @Note    : 
# @E-mail  : njbxhzy@hotmail.com
from math import log


def cal_insert(i, miu):
    """
    calculate insert time
    :param i: 0~1, the identity of L-LTR and R-LTR
    :param miu: neutral mutation rate
    :return: insert time
    """
    return (-3/4 * log(-1/3+4/3*i)) / (2*miu)


def parse_edta_gff(gff_path):
    iden_lst = []
    type_lst = []
    with open(gff_path) as f_in:
        lines = f_in.read().strip().split('\n')
        for line in lines:
            if not line.startswith("#"):
                eles = line.split('\t')
                if eles[2] == 'repeat_region':
                    for _ in eles[8].split(';'):
                        if _.startswith('ltr_identity'):
                            iden_lst.append(float(_.split('=')[1]))
                        elif _.startswith('Classification'):
                            type_lst.append(_.split('=')[1])
    return iden_lst, type_lst


def parse_args():
    import argparse
    cmdparser = argparse.ArgumentParser(prog="TE insert time determination")
    cmdparser.add_argument('-i', '--input', required=True,
                           help='<file path> *.mod.EDTA.raw/*.mod.LTR.intact.gff3')
    cmdparser.add_argument('-n', type=float, required=True,
                           help='<float> The neutral mutation rate, e.g. 1e-10')
    cmdparser.add_argument('-o', '--output', required=True,
                           help='<file path> output list')
    args = cmdparser.parse_args()
    return args


def main(args):
    iden_lst, type_lst = parse_edta_gff(args.input)
    time_list = [cal_insert(_, 3e-9) for _ in iden_lst]
    with open(args.output, 'w') as f_out:
        f_out.write('\n'.join([_type + '\t' + str(_time) for _type, _time in zip(type_lst, time_list)]))


if __name__ == '__main__':
    main(parse_args())
