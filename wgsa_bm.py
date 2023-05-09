#!/gs/gsfs0/users/yizhao/miniconda2/bin/python
# coding: utf-8
from __future__ import print_function

import sys

import scipy

# print(sys.version)
from scipy import stats
import logging
import re
from collections import Counter
from itertools import compress, product, chain
import copy
from subprocess import Popen, PIPE
import os
import sqlite3
import itertools
import math
from scipy.stats import wilcoxon
import multiprocessing
from multiprocessing.managers import BaseManager
import time
import socket

if sys.version_info[0] == 2:
    import Queue
elif sys.version_info[0] == 3:
    import queue
from numpy import mean, median, ptp, var, std, array, cov, corrcoef
import numpy as np
import itertools
import traceback
from bs4 import BeautifulSoup
from selenium import webdriver
from liftover import get_lifter

log_format = "%(asctime)s %(levelname)-8s %(process)d [%(funcName)s]%(message)s(%(filename)s line%(lineno)d)"


# logging.basicConfig(filename="_wgsa.log", level=logging.DEBUG, format=log_format, filemode="w")
def filter1(rate, bystro_tsv, output):
    """
    According to the result file tsv of bystro, establish and filter the original vcf file. 
    The filter condition is the frequency of occurrence in the population<=rate
    :param output: str
    :param bystro_tsv: str
    :type rate: float
    """
    logging.basicConfig(filename="filter1.log", level=logging.DEBUG, format=log_format, filemode="w")

    def adjust_f(f, ref):
        if ref > 0.5:
            if f > 0.5:
                return f
            else:
                return 1 - f
        else:
            if f < 0.5:
                return f
            else:
                return 1 - f

    assert type(rate) == float
    assert type(bystro_tsv) == str
    assert type(output) == str
    assert rate < 0.5
    iCounter = 0
    with open(bystro_tsv, "r") as fp_in, open(output, "w") as fp_out:
        fp_out.write(fp_in.readline())  # head
        while True:
            line = fp_in.readline()
            if iCounter % 10000 == 0 and iCounter > 0:
                print("handled {} lines".format(iCounter))
            if not line:
                break
            iCounter += 1
            line_list = line.strip("\n").split("\t")
            if line_list[68] == "!":
                fp_out.write(line)
                iCounter += 1
                continue
            f1 = float(line_list[68])
            f2 = float(line_list[14])
            if adjust_f(f1, f2) >= rate:
                iCounter += 1
                continue
            fp_out.write(line)
            iCounter += 1
    print("handled {} lines".format(iCounter))
    print("all done!")


def isnumber(aString):
    try:
        float(aString)
        return True
    except:
        return False


def data_reduce(file_in, col, output):
    logging.basicConfig(filename="data_reduce.log", level=logging.DEBUG, format=log_format, filemode="w")
    assert type(file_in) == str
    assert type(col) == str
    assert type(output) == str
    col_list = [int(i) for i in col.split(",")]
    iCounter = 0
    with open(file_in, "r") as fp_in, open(output, "w") as fp_out:
        fp_out.write(fp_in.readline())
        iCounter += 1
        while True:
            line_list = fp_in.readline().strip("\n").split("\t")
            iCounter += 1

            if not line_list or line_list == ['']:
                break
            for i in col_list:
                if len(line_list) < i:
                    print("row={0} {1}".format(iCounter, line_list))
                    continue

                tmp = filter(lambda x: isnumber(x), re.split(r";|\|", line_list[i - 1]))
                if not tmp:
                    print("row={2} col={0} [{1}]".format(i, line_list[i - 1], iCounter))
                else:
                    line_list[i - 1] = str(max([float(j) for j in tmp]))
            fp_out.write("{}\n".format("\t".join(line_list)))


def get_words(file_in, col_list):
    """
    For a given tsv file, analyze which words appear in the specified column (words are separated 
    by semicolons)）
    :type col_list: str
    :param col_list: 以逗号分隔的待分析列序号（从1开始）
    :param file_in: tsv文件
    :return:
    """
    logging.basicConfig(filename="get_words.log", level=logging.DEBUG, format=log_format, filemode="w")
    assert type(file_in) == str
    assert type(col_list) == str

    def analyze_column(name, all_data, col_num):
        """

        :param name:
        :param col_num: 列序号（从1开始）
        :return:
        """
        my_data = [i[col_num - 1] for i in all_data]
        # (";".join(my_data)).split(";")
        print("column name = {}".format(name))
        print([i[0] for i in Counter(re.split(";|\|", ";".join(my_data))).most_common()])
        print("\n")

    col_list = [int(i) for i in col_list.split("/")]
    # print col_list
    with open(file_in, "r") as fp:
        head = fp.readline().strip("\n").split("\t")
        all_data = [i.strip("\n").split("\t") for i in fp.readlines()]
    for i in col_list:
        analyze_column(head[i - 1], all_data, i)


def reduce_words(input_file, config_file, output):
    """
    According to the word_order, the specified column in the tsv is deleted, and only the most pre-
    examined word in the word list is kept.Which columns and word lists are analyzed are recorded 
    in a config file in the following format:
    col_num1    word1_1 word1_2 word1_3
    col_num2    word2_1 word2_2
    :param config_file: word_order文件
    :param input_file:
    :param output:
    :return:
    """
    logging.basicConfig(filename="reduce_words.log", level=logging.DEBUG, format=log_format, filemode="w")
    logging.debug("input_file=[{0}] confit=[{1}] output=[{2}]".format(input_file, config_file, output))
    col_list = []
    word_order_list = []
    with open(config_file, "r") as fp:
        while True:
            data_list = fp.readline()
            if not data_list:
                break
            if data_list.startswith("#"):
                continue
            data_list = data_list.strip("\n").strip("\r").split("\t")
            if len(data_list) <= 1:
                continue
            col_list.append(int(data_list[0]))
            word_order_list.append(dict(zip(data_list[1:], xrange(len(data_list) - 1))))  # list[dict{word:order}]

    def select_word(word_list, order_dict):
        selected = word_list[0]
        for word in word_list:
            if word not in order_dict:
                logging.debug("word = {0} not in order_dict = {1}".format(word, order_dict))
            if selected not in order_dict:
                logging.debug("word = {0} not in order_dict = {1}".format(selected, order_dict))
            if order_dict[word] < order_dict[selected]:
                selected = word
        return selected

    def handle_column(col_num, word_order_dict, data_pointer):
        """
        :param word_order_dict:
        :param col_num: 从1开始
        :param data_pointer:
        :return:
        """
        logging.debug("handling " + str(col_num))
        for line_data in data_pointer:
            word_list = [i[0] for i in Counter(re.split(";|\|", line_data[col_num - 1])).most_common()]
            if len(word_list) > 1:
                line_data[col_num - 1] = select_word(word_list, word_order_dict)
                continue
            line_data[col_num - 1] = word_list[0]

    with open(input_file, "r") as fp:
        all_data = fp.readlines()
    data = [i.strip("\n").split("\t") for i in all_data[1:]]
    head_data = all_data[0]
    del all_data
    for i in xrange(len(col_list)):
        handle_column(col_list[i], word_order_list[i], data)

    with open(output, "w") as fp:
        fp.write(head_data)
        fp.write("\n".join(["\t".join(i) for i in data]))


def filter_by_word(file_in, word_file, mode, output):
    """
    Filter the tsv file according to word, filter condition: as long as one of the selected 
    columns is hit, keep the row

    :param mode:    0 Keep the row as long as one of the selected columns hits
                    1 All selected columns must be hit to keep the row
    :param file_in:
    :param word_file:
    :param output:
    :return:
    """
    logging.basicConfig(filename="filter_by_word.log", level=logging.DEBUG, format=log_format, filemode="w")
    assert type(file_in) == str
    assert type(word_file) == str
    assert type(mode) == str
    assert type(output) == str

    def should_keep(mode, line_list, col_list, word_set_list):
        for i in xrange(len(col_list)):
            if len(line_list) < col_list[i]:
                return False
            if mode == "0":
                if line_list[col_list[i] - 1] in word_set_list[i]:
                    return True
            else:
                if line_list[col_list[i] - 1] not in word_set_list[i]:
                    return False
        return mode == "1"

    col_list = []
    word_set_list = []
    with open(word_file, "r") as fp:
        while True:
            data_list = fp.readline()
            if not data_list:
                break
            if data_list.startswith("#"):
                continue
            data_list = data_list.strip("\n").split("\t")
            if len(data_list) <= 1:
                continue
            col_list.append(int(data_list[0]))
            word_set_list.append(set(data_list[1:]))
    with open(file_in, "r") as fp_in, open(output, "w") as fp_out:
        fp_out.write(fp_in.readline())
        while True:
            line = fp_in.readline()
            if not line:
                break
            line_list = line.strip("\n").split("\t")
            if should_keep(mode, line_list, col_list, word_set_list):
                fp_out.write(line)


def reclassify_ALT(lof_file, output):
    logging.basicConfig(filename="reclassify_ALT.log", level=logging.DEBUG, format=log_format, filemode="w")
    assert type(lof_file) == str
    assert type(output) == str
    name_col = 14
    alt_col = 16
    lof_alt_col = 8
    sample_begin_col = 21
    with open(lof_file, "r") as fp_in:
        lof_lines = fp_in.readlines()
    lof_list = filter(lambda x: len(x) > 20, [i.strip("\n").split("\t") for i in lof_lines])
    name_list = [i[name_col - 1] for i in lof_list]
    with open(output, "w") as fp_out, open(output + ".unusual", "w") as fp_out2:
        for i in lof_list:
            alt_list = i[alt_col - 1].split(",")
            if name_list.count(i[name_col - 1]) < len(alt_list):
                if len(alt_list) > 2:
                    fp_out2.write("{}\n".format("\t".join(i)))
                    continue
                assert len(alt_list) == 2
                # two alt, 1 lof
                lof_alt_index = alt_list.index(i[lof_alt_col - 1])
                none_lof_value = str(2 - lof_alt_index)
                # start modifying data
                i[alt_col - 1] = alt_list[lof_alt_index]
                for index in xrange(len(i) - sample_begin_col + 1):
                    GT_list = i[sample_begin_col - 1 + index].split(":")
                    GT_list[0] = GT_list[0].replace(none_lof_value, "0")
                    i[sample_begin_col - 1 + index] = ":".join(GT_list)
                fp_out.write("{}\n".format("\t".join(i)))
            else:

                fp_out.write("{}\n".format("\t".join(i)))


def filter_alt(vcf_file, output, output2):
    logging.basicConfig(filename="filter_alt.log", level=logging.DEBUG, format=log_format, filemode="w")
    iCounter = 0
    with open(vcf_file, "r") as fp_in, open(output, "w") as fp_out, open(output2, "w") as fp_out2:
        while True:
            vcf_line = fp_in.readline()
            if vcf_line == "":
                break
            if vcf_line.startswith("#"):
                fp_out.write(vcf_line)
                fp_out2.write(vcf_line)
                continue
            vcf_list = vcf_line.strip("\r").strip("\n").strip("\r").split("\t")
            if len(vcf_list) < 9:
                continue

            alt_list = vcf_list[4].split(",")
            GT_element_list = []
            for GT in vcf_list[9:]:
                GT_element_list.extend(GT.split(":")[0].split("/"))
            # all the numbers in samples
            GT_set = set(GT_element_list)

            alt_num_list = [i + 1 for i in list(xrange(len(alt_list)))]
            my_filter = [str(i) in GT_set for i in alt_num_list]
            new_alt_list = list(compress(alt_list, my_filter))
            alt_num_list = list(compress(alt_num_list, my_filter))

            new_alt_num_list = [i + 1 for i in list(xrange(len(alt_num_list)))]
            org2new_num_dict = dict(zip(alt_num_list, new_alt_num_list))

            # change the data
            vcf_list[4] = ",".join(new_alt_list)
            for GT_index in xrange(len(vcf_list) - 9):
                GT_list = vcf_list[9 + GT_index].split(":")
                GT_list2 = GT_list[0].split("/")
                for index2 in xrange(len(GT_list2)):
                    if GT_list2[index2] == ".":
                        continue
                    if int(GT_list2[index2]) in org2new_num_dict:
                        GT_list2[index2] = org2new_num_dict[int(GT_list2[index2])]
                    else:
                        GT_list2[index2] = 0
                GT_list[0] = "/".join([str(i) for i in GT_list2])
                vcf_list[9 + GT_index] = ":".join(GT_list)
            result_str = "{}\n".format("\t".join(vcf_list))
            fp_out.write(result_str)
            if new_alt_list != alt_list:
                fp_out2.write("{0}\t{1}".format(",".join(alt_list), result_str))
            iCounter += 1
            if iCounter % 1000 == 0:
                print("handled {} lines".format(iCounter))
    print("all done!")


def line_left_normalization(vcf_line):
    """
    Perform left_normalization on a row of vcf data and return
    :param vcf_line:
    :return:
    """

    if vcf_line.startswith("#"):
        return vcf_line
    if not vcf_line.strip():
        return ""
    vcf_list = vcf_line.strip("\r").strip("\n").strip("\r").split("\t")
    if len(vcf_list) < 2:
        print('vcf line is [{}]'.format(vcf_line))
        sys.stdout.flush()
    pos = int(vcf_list[1])
    ref = vcf_list[3]
    alt_list = vcf_list[4].split(",")
    alt_len_list = [len(i) for i in alt_list]
    ileft_delete = 0
    for index in xrange(min(min(alt_len_list), len(ref)) - 1):
        ch_list = list(set(map(lambda x: x[index], alt_list)))
        if len(ch_list) == 1 and ref[index] == ch_list[0]:
            ileft_delete += 1
        else:
            break
    pos = pos + ileft_delete
    ref = ref[ileft_delete:]
    alt_list = map(lambda x: x[ileft_delete:], alt_list)
    alt_len_list = [len(i) for i in alt_list]
    if len(ref) == min(alt_len_list) == max(alt_len_list):
        iright_delete = 0
        for index in xrange(len(ref) - 1):
            ch_list = list(set(map(lambda x: x[-1 * index - 1], alt_list)))
            if len(ch_list) == 1 and ref[-1 * index - 1] == ch_list[0]:
                iright_delete += 1
            else:
                break
        if iright_delete > 0:
            ref = ref[:-1 * iright_delete]
            alt_list = map(lambda x: x[:-1 * iright_delete], alt_list)

    vcf_list[1] = str(pos)
    # if not vcf_list[2].startswith("rs"):
    #     vcf_list[2] = "{0}_{1}".format(vcf_list[0], vcf_list[1])
    vcf_list[3] = ref
    vcf_list[4] = ",".join(alt_list)
    return "\t".join(vcf_list)


def left_normalization(vcf_file, output, mode=0):
    logging.basicConfig(filename="left_normalization.log", level=logging.DEBUG, format=log_format, filemode="w")
    if mode == 0:
        icounter = 0
        with open(vcf_file, "r") as fp_in, open(output, "w") as fp_out:
            while True:
                vcf_line = fp_in.readline()
                if not vcf_line:
                    break
                if vcf_line.startswith("#"):
                    fp_out.write(vcf_line)
                    continue
                if not vcf_line.strip():
                    continue
                icounter += 1
                if icounter > 1:
                    fp_out.write("\n")
                fp_out.write("{}".format(line_left_normalization(vcf_line)))
    else:
        while True:
            vcf_line = sys.stdin.readline()
            if not vcf_line:
                break
            if vcf_line.startswith("#"):
                sys.stdout.write(vcf_line)
                continue
            if not vcf_line.strip():
                continue
            sys.stdout.write("{}\n".format(line_left_normalization(vcf_line)))


def reduce_duplicate(annovar_result, col_ref, col_alt, col_chr, col_pos, list_file, output):
    logging.basicConfig(filename="reduce_duplicate.log", level=logging.DEBUG, format=log_format, filemode="w")
    logging.debug("1")
    assert type(annovar_result) == str
    assert type(col_ref) == int
    assert type(col_alt) == int
    assert type(col_chr) == int
    assert type(col_pos) == int
    assert type(list_file) == str
    assert type(output) == str
    with open(list_file, "r") as fp:
        order_list = fp.readline().strip("\r").strip("\n").strip("\r").split("\t")
    col_key = int(order_list[0])
    order_list = order_list[1:]
    # print order_list
    order_dict = dict(zip(order_list, xrange(len(order_list))))
    # logging.debug("order_dict=[{}]".format(order_dict))
    data_dict = {}
    with open(annovar_result, "r") as fp:
        while True:
            data_line = fp.readline()
            if not data_line:
                break
            data_list = data_line.split("\t")
            if len(data_list) < 12:
                continue
            data_key = "{0}_{1}_{2}_{3}".format(data_list[col_ref - 1], data_list[col_alt - 1],
                                                data_list[col_chr - 1], data_list[col_pos - 1])
            # logging.debug("handling data_key=[{}]".format(data_key))
            if data_list[col_key - 1] not in order_list:
                logging.debug("[{0}] is not in order list. org data = [{1}]".format(data_list[col_key - 1], data_line))
                continue
            if data_key not in data_dict:
                data_dict[data_key] = [data_line, order_dict[data_list[col_key - 1]]]
            else:
                if order_dict[data_list[col_key - 1]] < data_dict[data_key][1]:
                    data_dict[data_key] = [data_line, order_dict[data_list[col_key - 1]]]

    with open(output, "w") as fp:
        for value in data_dict.itervalues():
            fp.write(value[0])


class vcf_line_info:
    def __init__(self, vcf_list, head_list):
        # type: (list[str], list[str]) -> None
        self.vcf_list = vcf_list  # type: list[str]
        self.head_list = []  # type: list[list[str]]
        self.head_list.append(head_list)
        self.alt_list = self.vcf_list[4].split(",")
        self.alt_num = len(self.alt_list)
        self.head_ref_alt_key_set = set(["{0}_{1}".format(head_list[5], head_list[6])])

    def add_head(self, head_list):
        self.head_list.append(head_list)
        self.head_ref_alt_key_set.add("{0}_{1}".format(head_list[5], head_list[6]))


def left_normalization_one_alt_annovar(ref, alt):
    alt_len = len(alt)
    ileft_delete = 0
    for i in xrange(min(alt_len, len(ref))):
        if ref[i] == alt[i]:
            ileft_delete += 1
        else:
            break
    ref = ref[ileft_delete:]
    alt = alt[ileft_delete:]
    iright_delete = 0
    if len(ref) == len(alt):
        for i in xrange(len(ref)):
            if ref[-1 * i - 1] == alt[-1 * i - 1]:
                iright_delete += 1
            else:
                break
        if iright_delete > 0:
            ref = ref[:-1 * iright_delete]
            alt = alt[:-1 * iright_delete]
    ref = "-" if ref == "" else ref
    alt = "-" if alt == "" else alt
    return [ref, alt, ileft_delete]


def regroup_alt(line_info):
    # type: (vcf_line_info) -> str
    left_normalized_ref_alt_list = [[i[5], i[6]] for i in line_info.head_list]  # [[ref, alt], [ref, alt]]
    filter_list = []
    for alt in line_info.alt_list:
        [ref_new, alt_new, ileft_delete] = left_normalization_one_alt_annovar(line_info.vcf_list[3], alt)

        if [ref_new, alt_new] not in left_normalized_ref_alt_list:
            filter_list.append(False)
        else:
            filter_list.append(True)
    org_alt_num_list = [i + 1 for i in xrange(len(line_info.alt_list))]
    new_alt_list = list(compress(line_info.alt_list, filter_list))
    filtered_org_alt_num_list = list(compress(org_alt_num_list, filter_list))
    line_info.vcf_list[4] = ",".join(new_alt_list)

    # build new alt num list
    new_alt_num_list = []
    for index in xrange(len(line_info.alt_list)):
        if not filter_list[index]:
            new_alt_num_list.append("0")
        else:
            new_alt_num_list.append(str(filtered_org_alt_num_list.index(index + 1) + 1))
    # dict org alt num --> new alt num
    dict_org2new_alt_num = dict(zip([str(i) for i in org_alt_num_list], new_alt_num_list))

    for i in xrange(len(line_info.vcf_list) - 9):
        gt_list = line_info.vcf_list[i + 9].split(":")
        gt_num_list = gt_list[0].split("/")
        for index in xrange(len(gt_num_list)):
            if gt_num_list[index] in dict_org2new_alt_num:
                gt_num_list[index] = dict_org2new_alt_num[gt_num_list[index]]

        gt_list[0] = "/".join(gt_num_list)
        line_info.vcf_list[i + 9] = ":".join(gt_list)

    return "{}\n".format("\t".join(line_info.vcf_list))


def rebuild_vcf(file_in, output):
    """
    According to the filtered annovar results, restore vcf. The alt that does not appear is treated as ref in the genotype.（group non-damaging alt alleles back to reference）
    """
    logging.basicConfig(filename="rebuild_vcf.log", level=logging.DEBUG, format=log_format, filemode="w")
    assert type(file_in) == str
    assert type(output) == str
    with open(file_in, "r") as fp:
        data_list = filter(lambda x: len(x) > 19, [i.strip().split("\t") for i in fp.readlines()])
    multiple_dict = {}
    with open(output, "w") as fp:
        for i in data_list:
            if "," not in i[14]:
                fp.write("{}\n".format("\t".join(i[10:])))
            else:
                vcf_key = "{0}_{1}".format(i[10], i[11])  # orgvcf:chr_pos
                if vcf_key not in multiple_dict:
                    multiple_dict[vcf_key] = vcf_line_info(i[10:], i[:10])  # A line of vcf may have multiple annovar comments (head)
                else:
                    multiple_dict[vcf_key].add_head(i[:10])
        for snp in multiple_dict.itervalues():
            if snp.alt_num == len(snp.head_ref_alt_key_set):
                fp.write("{}\n".format("\t".join(snp.vcf_list)))
            else:
                fp.write(regroup_alt(snp))


class bystro_anno_info:
    @staticmethod
    def build_annovar_ref_alt_key_from_bystro(bystro_ref, bystro_alt):
        """
        bystro style ref and alt --> annovar style ref and alt
        :type bystro_ref: str
        :type bystro_alt: str
        """
        if bystro_alt.startswith("-"):
            ref = bystro_ref.replace("|", "")
            ret = [ref, "-"]
        elif bystro_alt.startswith("+"):
            alt = bystro_alt.replace("+", "")
            ret = ["-", alt]
        else:
            ret = [bystro_ref, bystro_alt]
        return ret

    def __init__(self, vcf_pos, bystro_ref, bystro_alt):
        assert type(vcf_pos) == str
        assert type(bystro_ref) == str
        assert type(bystro_alt) == str
        self.vcf_pos = vcf_pos
        self.anno_list = []
        annovar_format_ref_alt_list = self.build_annovar_ref_alt_key_from_bystro(bystro_ref, bystro_alt)
        self.anno_list.append(annovar_format_ref_alt_list)
        self.anno_ref_alt_set = set(["_".join(annovar_format_ref_alt_list)])

    def add_anno(self, bystro_ref, bystro_alt):
        annovar_format_ref_alt_list = self.build_annovar_ref_alt_key_from_bystro(bystro_ref, bystro_alt)
        self.anno_list.append(annovar_format_ref_alt_list)
        self.anno_ref_alt_set.add("_".join(annovar_format_ref_alt_list))


def rebuild_vcf_bystro(tsv_in, vcf_in, output, need_sample=False):
    logging.basicConfig(filename="rebuild_vcf_bystro_{}.log".format(os.path.basename(tsv_in)), level=logging.DEBUG,
                        format=log_format, filemode="w")

    def get_vcf_variance_list(vcf_list):
        ret = []
        alt_list = vcf_list[4].split(",")
        for alt in alt_list:
            [new_ref, new_alt, ileft_delete] = left_normalization_one_alt_annovar(vcf_list[3], alt)
            ret.append("{0}_{1}".format(new_ref, new_alt))
        return ret

    # tsv --> tsv_dict
    with open(tsv_in, "r") as fp:
        tsv_data = [i.split("\t") for i in fp.readlines()]
    tsv_dict = {}
    for tsv_list in tsv_data[1:]:
        vcf_key = "{0}_{1}".format(tsv_list[0], tsv_list[15])  # chr_vcfpos
        if vcf_key not in tsv_dict:
            tsv_dict[vcf_key] = bystro_anno_info(tsv_list[15], tsv_list[16], tsv_list[4])  # vcf_pos ref alt
        else:
            tsv_dict[vcf_key].add_anno(tsv_list[16], tsv_list[4])  # ref alt

    with open(vcf_in, "r") as fp_in, open(output, "w") as fp_out:
        icounter = 0
        while True:
            vcf_line = fp_in.readline()
            if not vcf_line:
                break
            if vcf_line.startswith("#") or not vcf_line.strip():
                continue
            icounter += 1
            if icounter > 0 and icounter % 10000 == 0:
                logging.debug("handled {} vcf lines".format(icounter))
            vcf_list = vcf_line.strip().split("\t")
            vcf_key = "{0}_{1}".format(vcf_list[0], vcf_list[1])  # chr_vcfpos

            if vcf_key not in tsv_dict:
                # logging.debug("{0} not in tsv_dict".format(vcf_key))
                continue
            anno = tsv_dict[vcf_key]  # type: bystro_anno_info
            vcf_ref_alt_list = get_vcf_variance_list(vcf_list)  # annovar style [ref_alt, ref_alt]
            vcf_ref_alt_set = set(vcf_ref_alt_list)

            if vcf_ref_alt_set & anno.anno_ref_alt_set == vcf_ref_alt_set:
                if need_sample:
                    fp_out.write(vcf_line)
                else:
                    fp_out.write("{}\n".format("\t".join(vcf_list[:7])))
                continue
            if vcf_ref_alt_set & anno.anno_ref_alt_set == set([]):
                continue
            # # if len(alt_list) < len(anno.anno_ref_alt_set):
            # #     logging.debug("{0} unusual".format(vcf_key))
            # #     continue
            # alt_list = vcf_list[4].split(",")
            # ref = vcf_list[3]
            # filter_list = []
            # for alt in alt_list:
            #     # left normalization
            #     [new_ref, new_alt, ileft_delete] = left_normalization_one_alt_annovar(ref, alt)
            #
            #     if "{0}_{1}".format(new_ref, new_alt) not in anno.anno_ref_alt_set:
            #         filter_list.append(False)
            #     else:
            #         filter_list.append(True)
            filter_list = [i in anno.anno_ref_alt_set for i in vcf_ref_alt_list]

            new_alt_list = list(compress(vcf_list[4].split(","), filter_list))
            vcf_list[4] = ",".join(new_alt_list)

            if need_sample:
                org_alt_num_list = [i + 1 for i in xrange(len(vcf_ref_alt_list))]
                filtered_org_alt_num_list = list(compress(org_alt_num_list, filter_list))
                # build new alt num list
                new_alt_num_list = []
                for index in xrange(len(vcf_ref_alt_list)):
                    if not filter_list[index]:
                        new_alt_num_list.append("0")
                    else:
                        new_alt_num_list.append(str(filtered_org_alt_num_list.index(index + 1) + 1))
                # dict org alt num --> new alt num
                dict_org2new_alt_num = dict(zip([str(i) for i in org_alt_num_list], new_alt_num_list))
                for i in xrange(len(vcf_list) - 9):
                    gt_list = vcf_list[i + 9].split(":")
                    gt_num_list = gt_list[0].split("/")
                    for index in xrange(len(gt_num_list)):
                        if gt_num_list[index] in dict_org2new_alt_num:
                            gt_num_list[index] = dict_org2new_alt_num[gt_num_list[index]]

                    gt_list[0] = "/".join(gt_num_list)
                    vcf_list[i + 9] = ":".join(gt_list)
                fp_out.write("{}\n".format("\t".join(vcf_list)))
            else:
                fp_out.write("{}\n".format("\t".join(vcf_list[:7])))
        logging.debug("all done")


class VepAnnoInfo:
    @staticmethod
    def __vep_key(vep_pos, vep_alt):
        # return "{0}:{1}".format(vep_pos, vep_alt)
        ret = vep_pos.split(":")
        ret[1] = ret[1].split("-")[0]
        ret.append(vep_alt)
        return "_".join(ret)  # "chr_pos_alt"

    def __init__(self, org_id, vep_pos, vep_alt):
        assert type(org_id) == str
        assert type(vep_pos) == str
        assert type(vep_alt) == str
        self.org_id = org_id
        self.anno_set = set([])
        self.anno_set.add(self.__vep_key(vep_pos, vep_alt))

    def add_anno(self, vep_pos, vep_alt):
        self.anno_set.add(self.__vep_key(vep_pos, vep_alt))


def rebuild_vcf_vep(lof_in, vcf_in, need_head, output):
    """
    The first step is to preliminarily screen the vcf according to the id in the lof file and save
    it in the form of a file
    :param lof_in:
    :param vcf_in:
    :param need_head:
    :param output:
    :return:
    """
    logging.basicConfig(filename="rebuild_vcf_vep.log", level=logging.DEBUG, format=log_format, filemode="w")
    assert type(lof_in) == str
    assert type(vcf_in) == str
    assert type(output) == str
    assert type(need_head) == bool

    # lof_in --> lof_dict
    with open(lof_in, "r") as fp:
        lof_data = [i.split("\t") for i in fp.readlines()]
    # lof_data = [[i[0], i[1], i[2]] for i in lof_data]
    lof_dict = {}
    for lof_list in lof_data:
        if lof_list[0] not in lof_dict:
            lof_dict[lof_list[0]] = VepAnnoInfo(lof_list[0], lof_list[1], lof_list[2])
        else:
            lof_dict[lof_list[0]].add_anno(lof_list[1], lof_list[2])
    icounter = 0
    with open(vcf_in, "r") as fp_vcf, open(output, "w") as fp_out:
        while True:
            if icounter > 0 and icounter % 10000 == 0:
                logging.debug("handled {} vcf lines".format(icounter))
            vcf_line = fp_vcf.readline()
            icounter += 1
            if not vcf_line:
                break
            if vcf_line.startswith("#"):
                if need_head:
                    fp_out.write(vcf_line)
                continue
            vcf_list = vcf_line.split("\t")

            # not lof
            if vcf_list[2] not in lof_dict:
                continue

            # is lof
            fp_out.write(vcf_line)


def rebuild_vcf_vep2(lof_in, vcf_in, need_head, output):
    """
    In the second step, the vcf filtered by id is reconstructed, and the alt that 
    does not appear is treated as ref. Delete the same id but not lof.
    :param lof_in:
    :param vcf_in:
    :param need_head:
    :param output:
    :return:
    """
    assert type(lof_in) == str
    assert type(vcf_in) == str
    assert type(output) == str
    assert type(need_head) == bool
    logging.basicConfig(filename="rebuild_vcf_vep2.log", level=logging.DEBUG, format=log_format, filemode="w")
    # lof_in --> lof_dict
    with open(lof_in, "r") as fp:
        lof_data = [i.split("\t") for i in fp.readlines()]
    # lof_data = [[i[0], i[1], i[2]] for i in lof_data]
    lof_dict = {}
    for lof_list in lof_data:
        if lof_list[0] not in lof_dict:
            lof_dict[lof_list[0]] = VepAnnoInfo(lof_list[0], lof_list[1], lof_list[2])
        else:
            lof_dict[lof_list[0]].add_anno(lof_list[1], lof_list[2])
    icounter = 0
    with open(vcf_in, "r") as fp_vcf, open(output, "w") as fp_out:
        should_continue = False
        while True:
            should_continue = False
            if icounter > 0 and icounter % 10000 == 0:
                logging.debug("handled {} vcf lines".format(icounter))
            vcf_line = fp_vcf.readline()
            icounter += 1
            if not vcf_line:
                break
            if vcf_line.startswith("#"):
                if need_head:
                    fp_out.write(vcf_line)
                continue
            vcf_list = vcf_line.split("\t")

            # not lof
            if vcf_list[2] not in lof_dict:
                continue

            # is lof
            anno_info = lof_dict[vcf_list[2]]
            anno_list = list(anno_info.anno_set)  # type: list[str]
            anno_alt_list = [i.split("_")[2] for i in anno_list]
            vcf_pos = int(vcf_list[1])
            vcf_alt_list = vcf_list[4].split(",")
            anno_pos_list = []  # type: list[int]
            for i in xrange(len(anno_list)):
                anno_chr, anno_pos = anno_list[i].split("_")[:2]
                anno_pos = int(anno_pos)
                anno_pos_list.append(anno_pos)
                # rebuild org vcf alt
                if anno_pos > vcf_pos:
                    delta = vcf_alt_list[0][:(anno_pos - vcf_pos)]
                    anno_alt_list[i] = (delta + anno_alt_list[i]).strip("-")
                if anno_pos < vcf_pos:
                    should_continue = True
                    break
            if should_continue:
                continue
            for i in xrange(len(anno_alt_list)):
                if anno_alt_list[i] not in vcf_alt_list:
                    candidate_alt_list = filter(lambda x: anno_alt_list[i] in x, vcf_alt_list)
                    if len(candidate_alt_list) != 1:
                        should_continue = True
                        logging.error("more than 1 candidate alt, anno_alt={0} "
                                      "vcf_alt_list={1} candidate_alt_list={2}".format(anno_alt_list[i],
                                                                                       vcf_alt_list,
                                                                                       candidate_alt_list))
                        break
                    else:
                        anno_alt_list[i] == candidate_alt_list[0]
                    # for j in vcf_alt_list:
                    #     if anno_alt_list[i] in j:
                    #         anno_alt_list[i] = j
                    #         break
            if should_continue:
                continue

            # double check
            for i in anno_alt_list:
                if i not in vcf_alt_list:
                    logging.error("anno_alt_list={0} vcf_alt_list={1} "
                                  "id={2} vcf_pos={3} anno_pos={4}".format(anno_alt_list, vcf_alt_list,
                                                                           anno_info.org_id, vcf_pos, anno_pos_list))
                    should_continue = True
                    break
            if should_continue:
                continue
            selector_list = [i in anno_alt_list for i in vcf_alt_list]
            new_vcf_alt_list = list(compress(vcf_alt_list, selector_list))
            vcf_list[4] = ",".join(new_vcf_alt_list)

            # build dict_org2new_alt_num
            org_alt_num_list = list(xrange(1, len(vcf_alt_list) + 1, 1))
            filtered_org_alt_num_list = list(compress(org_alt_num_list, selector_list))
            new_alt_num_list = []
            for i in xrange(len(vcf_alt_list)):
                if not selector_list[i]:
                    new_alt_num_list.append("0")
                else:
                    new_alt_num_list.append(str(filtered_org_alt_num_list.index(org_alt_num_list[i]) + 1))

            dict_org2new_alt_num = dict(zip([str(i) for i in org_alt_num_list], new_alt_num_list))
            for i in xrange(len(vcf_list) - 9):
                gt_list = vcf_list[i + 9].split(":")
                gt_num_list = gt_list[0].split("/")
                for index in xrange(len(gt_num_list)):
                    if gt_num_list[index] in dict_org2new_alt_num:
                        gt_num_list[index] = dict_org2new_alt_num[gt_num_list[index]]

                gt_list[0] = "/".join(gt_num_list)
                vcf_list[i + 9] = ":".join(gt_list)
            fp_out.write("\t".join(vcf_list))


def check_id(vcf_in, id_col):
    logging.basicConfig(filename="check_id.log", level=logging.DEBUG, format=log_format, filemode="w")
    id_dict = {}
    with open(vcf_in, "r") as fp:
        while True:
            vcf_line = fp.readline()
            if not vcf_line:
                break
            if vcf_line.startswith("#"):
                continue
            vcf_list = vcf_line.split("\t")
            vcf_id = vcf_list[int(id_col) - 1]
            if vcf_id in id_dict:
                id_dict[vcf_id] += 1
                print(vcf_id)
            else:
                id_dict[vcf_id] = 1


def split_vcf_line(data_line, result_list):
    """
    Split a row of data of the original vcf into multiple rows of data according to the principle 
    of one alt per row。
    :param result_list: Each row of data of the new vcf
    :type result_list: list[str]
    :type data_line: str
    """
    # print "split_vcf_line begin!"
    data_list = data_line.split("\t")[:7]
    # print len(data_list)
    # print data_list
    if len(data_list) < 5:
        logging.error("This line only has {} element, but at least 5 needed. [{}]".format(data_line))
        return
    alt_list = data_list[4].split(",")
    for alt_index in range(len(alt_list)):
        curr_list = copy.copy(data_list)
        curr_list[4] = alt_list[alt_index]
        result_list.append("\t".join(curr_list))
    return result_list


def split_vcf(vcf_in, output, mode=0):
    logging.basicConfig(filename="split_vcf.log", level=logging.DEBUG, format=log_format, filemode="w")
    if mode == 0:
        ret_list = []
        with open(vcf_in, "r") as fp:
            vcf_data = fp.readlines()
        for vcf_line in vcf_data:
            if vcf_line.startswith("#"):
                ret_list.append(vcf_line.strip())
            split_vcf_line(vcf_line, ret_list)
        with open(output, "w") as fp:
            fp.write("\n".join(ret_list))
    else:
        # for vcf_line in sys.stdin.readlines():
        while True:
            vcf_line = sys.stdin.readline()
            if not vcf_line:
                break
            if vcf_line.startswith("#"):
                sys.stdout.write(vcf_line)
                continue
            vcf_line = vcf_line.strip()
            if len(vcf_line) == 0:
                continue
            ret_list = []
            split_vcf_line(vcf_line, ret_list)
            sys.stdout.write("{}\n".format("\n".join(ret_list)))


def union_intersection(file_list, col_list, output, mode):
    """

    :param mode: "union", "intersection"
    :type output: str
    :type col_list: str
    :type file_list: str
    """

    # logging.basicConfig(filename="union_intersection.log", level=logging.DEBUG, format=log_format, filemode="w")

    def get_key_line(org_line, col_list):
        """

        :type col_list: list[int]
        :type org_line: str
        """
        org_list = org_line.strip().split("\t")
        ret_list = []
        for i in col_list:
            try:
                if i == 1:
                    ret_list.append(org_list[0].strip("chr"))
                else:
                    ret_list.append(org_list[i - 1])
            except:
                print("len(org_list) = {}".format(len(org_list)))
                print("col_list = {}".format(col_list))
        return "\t".join(ret_list)

    file_list = filter(lambda x: len(x.strip()) > 0, file_list.split(","))
    col_list = [int(i) for i in filter(lambda x: len(x.strip()) > 0, col_list.split(","))]
    if mode not in ["union", "intersection"]:
        print("wrong mode: {}".format(mode))
        exit(0)

    # files ---> set lists of key line
    set_list = []
    for each_file in file_list:
        print("loading {}".format(each_file))
        with open(each_file, "r") as fp:
            set_list.append(set(
                [get_key_line(i, col_list) for i in fp.readlines() if len(i.strip()) > 0 and not i.startswith("#")]))

    # union or intersection
    ret_set = set_list[0]
    for index in xrange(len(set_list)):
        if index == 0:
            continue
        if mode == "union":
            ret_set = ret_set | set_list[index]
        elif mode == "intersection":
            ret_set = ret_set & set_list[index]

    # write output
    ret_list = list(ret_set)
    with open(output, "w") as fp:
        fp.write("\n".join(ret_list))
        fp.write("\n")


# def bav_intersection(bystro_vcf, annovar_vcf, vep_vcf, output):
#     with open(bystro_vcf, "r") as fp:
#         bystro_data = [i.strip("\r").strip("\n").strip("\r").split("\t") for i in fp.readlines()]
#         pass


def add_colums(data_file, data_col, lib_file, lib_col, add_cols, output, mode):
    """

    :type add_cols: str
    :type lib_col: str
    :type data_col: str
    """
    logging.basicConfig(filename="add_colums.log", level=logging.DEBUG, format=log_format, filemode="w")

    def get_key(line, col_list):
        """

        :type col_list: list
        :type line: str
        """
        data_list = line.split("\t")
        ret_list = []  # type: list[str]
        for i in col_list:
            ret_list.append(data_list[i - 1])
        if ret_list[0].startswith("chr"):
            ret_list[0] = ret_list[0].strip("chr")
        return "_".join(ret_list)

    assert mode in ["left", "right"]
    add_col_list = [int(i) for i in add_cols.split(",")]
    data_col_list = [int(i) for i in data_col.split(",")]
    lib_col_list = [int(i) for i in lib_col.split(",")]

    # build key2lib_pos_dict
    icounter = 0
    key2lib_pos_dict = {}
    with open(lib_file, "r") as fp_lib:
        while True:
            curr_pos = fp_lib.tell()
            line = fp_lib.readline().strip()
            icounter += 1
            if icounter % 100000 == 0:
                logging.debug("handled {} lib lines".format(icounter))
            if not line:
                break
            if line.startswith("#"):
                continue
            if len(line) == 0:
                continue
            curr_key = get_key(line, lib_col_list)
            if curr_key not in key2lib_pos_dict:
                key2lib_pos_dict[curr_key] = [curr_pos]
            else:
                key2lib_pos_dict[curr_key].append(curr_pos)
    logging.debug("finish building key2lib_pos_dict")

    # handle vcf file line by line
    with open(data_file, "r") as fp_vcf, open(lib_file, "r") as fp_lib, open(output, "w") as fp_out:
        icounter = 0
        while True:
            line = fp_vcf.readline().strip()
            icounter += 1
            if icounter % 10000 == 0:
                logging.debug("handled {} vcf lines".format(icounter))
            if not line:
                break
            if line.startswith("#"):
                # fp_out.write("{}\n".format(line))
                continue
            if len(line) == 0:
                continue
            curr_key = get_key(line, data_col_list)
            if curr_key not in key2lib_pos_dict:
                if mode == "left":
                    fp_out.write("-\t{0}\n".format(line))
                else:
                    fp_out.write("{0}\t-\n".format(line))
                continue
            lib_pos_list = key2lib_pos_dict[curr_key]
            anno_list = []
            for pos in lib_pos_list:
                fp_lib.seek(pos, 0)
                anno_list.append(get_key(fp_lib.readline().strip(), add_col_list))
            anno_list = list(set(anno_list))
            if mode == "left":
                fp_out.write("{0}\t{1}\n".format("|".join(anno_list), line))
            else:
                fp_out.write("{1}\t{0}\n".format("|".join(anno_list), line))
    logging.debug("all done")


def splice_ai_filter(input, score, output, mode=0):
    logging.basicConfig(filename="splice_ai_filter.log", level=logging.DEBUG, format=log_format, filemode="w")

    def handle_element(str):
        ret = tuple(str.split("="))
        if len(ret) > 2:
            logging.error(str)
            exit(0)
        return ret

    def get_max_score(str_in):
        """

        :type str_in: str
        """
        anno_list = str_in.split("|")
        max_score = 0.0
        for anno in anno_list:
            curr_score_dict = dict([handle_element(element) for element in anno.split(";")])
            max_score = max(float(curr_score_dict["DS_AG"]),
                            float(curr_score_dict["DS_AL"]),
                            float(curr_score_dict["DS_DG"]),
                            float(curr_score_dict["DS_DL"]),
                            max_score)
        return max_score

    assert type(input) == str
    assert type(score) == str
    assert type(output) == str

    score_f = float(score)
    if mode == 0:
        with open(input, "r") as fp, open(output, "w") as fp_out:
            while True:
                data_line = fp.readline().strip()
                if not data_line:
                    break
                if data_line.startswith("#") or data_line.startswith("-"):
                    continue
                data_list = data_line.split("\t")
                # score_dict = dict([handle_element(i) for i in data_list[0].split(";")])
                if get_max_score(data_list[0]) >= score_f:
                    fp_out.write("{}\n".format(data_line))
    else:
        while True:
            data_line = sys.stdin.readline().strip()
            if not data_line:
                break
            if data_line.startswith("#") or data_line.startswith("-"):
                continue
            data_list = data_line.split("\t")

            # score_dict = dict([handle_element(i) for i in data_list[0].split(";")])
            if get_max_score(data_list[0]) >= score_f:
                sys.stdout.write("{}\n".format(data_line))


def cut(data_in, num, max_in, min_in, output, mode=0):
    logging.basicConfig(filename="cut.log", level=logging.DEBUG, format=log_format, filemode="w")
    delta = (float(max_in) - float(min_in)) / float(num)
    range_list = [[float(min_in) + i * delta, 0] for i in xrange(int(num))]
    if mode == 0:
        with open(data_in) as fp_in:
            data = [float(i) for i in fp_in.readlines()[1:]]
        data1 = Counter(data).most_common()
        del data
        for element in data1:
            for left_range in range_list:
                if element[0] <= left_range[0] + delta:
                    left_range[1] += element[1]
                    break
        with open(output, "w") as fp_out:
            fp_out.write("range\tnum\n")
            for i in range_list:
                fp_out.write("{0}\t{1}\n".format(i[0], i[1]))
    else:
        while True:
            data_line = sys.stdin.readline()
            if not data_line:
                break
            try:
                element = float(data_line)
            except:
                continue
            if element >= range_list[-1][0] + delta:
                range_list[-1][1] += 1
                continue
            for left_range in range_list:
                if element <= left_range[0] + delta:
                    left_range[1] += 1
                    break
        sys.stdout.write("range\tnum\n")
        for i in range_list:
            sys.stdout.write("{0}\t{1}\n".format(i[0], i[1]))


def add_vcf_head(ver, fai_in, vcf_in, output, org_vcf):
    def load_table_head(org_vcf):
        with open(org_vcf, "r") as fp:
            while True:
                line = fp.readline()
                if not line:
                    break
                if line.startswith("#") and not line.startswith("##"):
                    return line
            return ""

    logging.basicConfig(filename="add_vcf_head.log", level=logging.DEBUG, format=log_format, filemode="w")
    if ver not in ["hg38", "hg19"]:
        print("""version should be hg38 or hg19""")
        exit(0)
    version = "GRCh37/hg19" if ver == "hg19" else "GRCh38/hg38"
    pp = Popen(["awk '{printf(\"##contig=<ID=\"$1\",length=\"$2\">\\n\");}' " + "{0}".format(fai_in)], stdout=PIPE,
               shell=True)
    contigs = pp.stdout.readlines()

    with open(vcf_in, "r") as fp_in, open(output, "w") as fp_out:
        fp_out.write("##fileformat=VCFv4.0\n##assembly={}\n".format(version))
        for contig in contigs:
            fp_out.write(contig)
        fp_out.write(load_table_head(org_vcf))
        while True:
            vcf_line = fp_in.readline()
            if not vcf_line:
                break
            fp_out.write(vcf_line)
            # vcf_list = vcf_line.split("\t")
            # if len(vcf_list) > 8:
            #     vcf_list = vcf_list[:8]
            # if len(vcf_list) < 8:
            #     for i in xrange(8 - len(vcf_list)):
            #         vcf_list.append(".")
            # fp_out.write("{}\n".format("\t".join(vcf_list)))


def select_vcf(vcf_in, num, output, mode=0):
    logging.basicConfig(filename="select_vcf.log", level=logging.DEBUG, format=log_format, filemode="w")
    icounter = 0
    if mode == 0:
        with open(vcf_in, "r") as fp_in, open(output, "w") as fp_out:
            while True:
                vcf_line = fp_in.readline()
                if not vcf_line:
                    break
                if vcf_line.startswith("#"):
                    fp_out.write(vcf_line)
                    continue
                icounter += 1
                if icounter % int(num) == 0:
                    fp_out.write(vcf_line)
                else:
                    continue
    else:
        while True:
            vcf_line = sys.stdin.readline()
            if not vcf_line:
                break
            if vcf_line.startswith("#"):
                sys.stdout.write(vcf_line)
                continue
            icounter += 1
            if icounter % int(num) == 0:
                sys.stdout.write(vcf_line)
            else:
                continue


def component(vcf_in, output_path, need_head):
    logging.basicConfig(filename="component.log", level=logging.DEBUG, format=log_format, filemode="w")

    def is_indel(ref, alt):
        # type: (str, str) -> bool
        alt_list = alt.split(",")
        for one_alt in alt_list:
            if len(one_alt) != len(ref):
                return True
        return False

    with open(vcf_in, "r") as fp_in, open(output_path + "snp1.vcf", "w") as fp_snp1, open(output_path + "snps.vcf",
                                                                                          "w") as fp_snps, open(
        output_path + "indel1.vcf", "w") as fp_indel1, open(output_path + "indels.vcf", "w") as fp_indels:
        while True:
            data_line = fp_in.readline()
            if not data_line:
                break
            if data_line.startswith("#"):
                if need_head:
                    fp_snp1.write(data_line)
                    fp_snps.write(data_line)
                    fp_indel1.write(data_line)
                    fp_indels.write(data_line)
                continue
            data_line_list = data_line.split("\t")
            if is_indel(data_line_list[3], data_line_list[4]):
                if "," in data_line_list[4]:
                    fp_indels.write(data_line)
                    continue
                else:
                    fp_indel1.write(data_line)
                    continue
            else:
                if "," in data_line_list[4]:
                    fp_snps.write(data_line)
                    continue
                else:
                    fp_snp1.write(data_line)
                    continue


def splice_ai_filter_indel(vcf_in, score, output, mode=0):
    logging.basicConfig(filename="splice_ai_filter_indel.log", level=logging.DEBUG, format=log_format,
                        filemode="w")
    score = float(score)
    if mode == 0:
        with open(vcf_in, "r") as fp_in, open(output, "w") as fp_out:
            while True:
                data_line = fp_in.readline()
                if not data_line:
                    break
                data_line = data_line.strip()
                if not data_line or data_line.startswith("#"):
                    continue
                vcf_list = data_line.split("\t")
                anno_list = vcf_list[-1].split("|")
                if len(anno_list) < 6:
                    continue
                try:
                    score_max = max(float(anno_list[2]), float(anno_list[3]), float(anno_list[4]), float(anno_list[5]))
                except:
                    logging.debug("can not get float {}".format(data_line))
                    continue
                if score_max < score:
                    continue
                fp_out.write("{}\n".format(data_line))
    else:
        while True:
            data_line = sys.stdin.readline().strip()
            if not data_line:
                break
            data_line = data_line.strip()
            if data_line.startswith("#") or not data_line:
                continue
            vcf_list = data_line.split("\t")
            anno_list = vcf_list[-1].split("|")
            if len(anno_list) < 6:
                continue
            try:
                score_max = max(float(anno_list[2]), float(anno_list[3]), float(anno_list[4]), float(anno_list[5]))
            except:
                logging.debug("can not get float [{}]".format(data_line))
                continue
            # score_dict = dict([handle_element(i) for i in data_list[0].split(";")])
            if score_max < score:
                continue
            sys.stdout.write("{}\n".format(data_line))


def select_by_cols(vcf_in, vcf_cols, list_file, list_cols, output, not_found):
    """
    Filter vcf according to some columns of list file
    :type list_cols: str
    :type vcf_cols: str
    :param vcf_cols: columns of vcf input. eg:1,2,4,5
    :param vcf_in:
    :param list_file:
    :param output:
    :return:
    """
    logging.basicConfig(filename="select_by_cols.log", level=logging.DEBUG, format=log_format, filemode="w")
    assert type(vcf_in) == str
    assert type(list_file) == str
    assert type(output) == str

    # assert type(not_found) == str

    def get_key_str(data_list, col_list):
        ret = ""
        for col in col_list:
            if not ret:
                ret = data_list[col - 1]
            else:
                ret = ret + "\t" + data_list[col - 1]
        return ret

    vcf_col_list = [int(i) for i in vcf_cols.split(",")]
    list_col_list = [int(i) for i in list_cols.split(",")]
    if len(vcf_col_list) != len(list_col_list):
        logging.warnings("vcf_cols has {0} key columns and list has {1} key columns".format(len(vcf_col_list),
                                                                                            len(list_col_list)))
        exit(0)

    print("building list set")
    with open(list_file, "r") as fp:
        list_set = set([get_key_str(i.strip().strip("chr").split("\t"), list_col_list) for i in
                        fp.readlines() if len(i.strip()) > 0 and not i.startswith("#")])
    icounter = 0
    if not_found is not None:
        with open(vcf_in, "r") as fp_vcf, open(output, "w") as fp_out, open(not_found, "w") as fp_not_found:
            while True:
                vcf_line = fp_vcf.readline()
                if not vcf_line:
                    print("handled {} vcf lines".format(icounter))
                    break
                if vcf_line.startswith("#") or len(vcf_line.strip()) == 0:
                    icounter += 1
                    continue

                if get_key_str(vcf_line.strip().strip("chr").split("\t"), vcf_col_list) in list_set:
                    fp_out.write(vcf_line)
                else:
                    fp_not_found.write(vcf_line)
                icounter += 1
                if icounter > 0 and icounter % 100000 == 0:
                    print("handled {} lines".format(icounter))
        print("all done")
    else:
        with open(vcf_in, "r") as fp_vcf, open(output, "w") as fp_out:
            while True:
                vcf_line = fp_vcf.readline()
                if not vcf_line:
                    print("handled {} vcf lines".format(icounter))
                    break
                if vcf_line.startswith("##") or len(vcf_line.strip()) == 0:
                    icounter += 1
                    continue
                if vcf_line.startswith("#"):
                    fp_out.write(vcf_line)
                    continue
                if get_key_str(vcf_line.strip().strip("chr").split("\t"), vcf_col_list) in list_set:
                    fp_out.write(vcf_line)
                # else:
                #     fp_not_found.write(vcf_line)
                icounter += 1
                if icounter > 0 and icounter % 100000 == 0:
                    print("handled {} lines".format(icounter))
        print("all done")


def splice_ai_merge_snp_indel(snp_in, indel_in, output):
    """
    Merge the filtered snp anno and indel anno into one file
    The comment column and comment format of snp are different from indel and need to be converted
    The annotation of snp is in the first column, format:
    SYMBOL=TUBB8;STRAND=-;TYPE=E;DIST=-53;DS_AG=0.0000;DS_AL=0.0000;DS_DG=0.0000;DS_DL=0.0000;DP_AG=-26;DP_AL=-10;DP_DG=3;DP_DL=35
    The comment of indel is in the eighth column, the format is:
    ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL
    SpliceAI=C|SAMD11|0.00|0.01|0.00|0.00|-126|-131|396|-53
    :param snp_in:
    :param indel_in:
    :param output:
    :return:
    """
    logging.basicConfig(filename="splice_ai_merge_snp_indel.log", level=logging.DEBUG, format=log_format, filemode="w")

    # check whether the file is filtered
    def check_file_ok(file_in):
        print("cheching {}".format(file_in))
        anno_index = -1
        with open(file_in, "r") as fp:
            while True:
                line = fp.readline()
                if not line:
                    break
                if line.startswith("#") or not line.strip():
                    continue

                if "SYMBOL" not in line and "SpliceAI" not in line:
                    return -1
                if anno_index == -1:
                    data_list = line.split("\t")
                    for i in xrange(len(data_list)):
                        if "SYMBOL" in data_list[i] or "SpliceAI" in data_list[i]:
                            anno_index = i
                            break
        return anno_index

    anno_index_snp = check_file_ok(snp_in)
    if anno_index_snp < 0:
        print("Make sure the [{}] file is filtered.")
        exit(0)
    anno_index_indel = check_file_ok(indel_in)
    if anno_index_indel < 0:
        print("Make sure the [{}] file is filtered.")
        exit(0)

    with open(snp_in, "r") as fp_snp, open(indel_in, "r") as fp_indel, open(output, "w") as fp_out:
        # handle snp
        while True:
            line = fp_snp.readline()
            if not line:
                break
            if line.startswith("#") or not line.strip():
                continue
            data_list = line.strip().split("\t")
            ret_str = ""
            anno_format_str = ""
            for i in xrange(len(data_list)):
                if i == anno_index_snp:
                    anno_str = data_list[i]
                    anno_list = anno_str.split("|")
                    for anno in anno_list:
                        try:
                            anno_dict = eval("{'" + anno.replace(";", "','").replace("=", "':'") + "'}")
                        except:
                            print("can not convert to dict [{}]".format(line))
                        if not anno_format_str:
                            anno_format_str = "SYMBOL=" + anno_dict["SYMBOL"] + "|"
                        else:
                            anno_format_str += "SYMBOL=" + anno_dict["SYMBOL"] + "|"
                        anno_format_str += "DS_AG=" + anno_dict["DS_AG"] + "|"
                        anno_format_str += "DS_AL=" + anno_dict["DS_AL"] + "|"
                        anno_format_str += "DS_DG=" + anno_dict["DS_DG"] + "|"
                        anno_format_str += "DS_DL=" + anno_dict["DS_DL"] + "|"
                        anno_format_str += "DP_AG=" + anno_dict["DP_AG"] + "|"
                        anno_format_str += "DP_AL=" + anno_dict["DP_AL"] + "|"
                        anno_format_str += "DP_DG=" + anno_dict["DP_DG"] + "|"
                        anno_format_str += "DP_DL=" + anno_dict["DP_DL"] + ";"
                    anno_format_str = anno_format_str[:-1]

                else:
                    if not ret_str:
                        ret_str = data_list[i]
                    else:
                        ret_str += "\t" + data_list[i]
            ret_str += "\t" + anno_format_str
            fp_out.write("{}\n".format(ret_str))
        # handle indel
        while True:
            ret_str = ""
            anno_format_str = ""
            line = fp_indel.readline()
            if not line:
                break
            if line.startswith("#") or not line.strip():
                continue
            data_list = line.strip().split("\t")
            for i in xrange(len(data_list)):
                if i == anno_index_indel:
                    anno_str = data_list[i]
                    anno_list = anno_str.split("|")
                    anno_format_str = "SYMBOL=" + anno_list[1] + "|"
                    anno_format_str += "DS_AG=" + anno_list[2] + "|"
                    anno_format_str += "DS_AL=" + anno_list[3] + "|"
                    anno_format_str += "DS_DG=" + anno_list[4] + "|"
                    anno_format_str += "DS_DL=" + anno_list[5] + "|"
                    anno_format_str += "DP_AG=" + anno_list[6] + "|"
                    anno_format_str += "DP_AL=" + anno_list[7] + "|"
                    anno_format_str += "DP_DG=" + anno_list[8] + "|"
                    anno_format_str += "DP_DL=" + anno_list[9]
                else:
                    if not ret_str:
                        ret_str = data_list[i]
                    else:
                        ret_str += "\t" + data_list[i]
            ret_str += "\t" + anno_format_str
            fp_out.write("{}\n".format(ret_str))


def build_vcf_index(vcf_in):
    logging.basicConfig(filename="build_vcf_index.log", level=logging.DEBUG, format=log_format, filemode="w")
    with open(vcf_in, "r") as fp, open(vcf_in + "i", "w") as fp_out:
        ret_list = []
        while True:
            curr_pos = fp.tell()
            data_line = fp.readline()
            if not data_line:
                break
            if data_line.startswith("#") or not data_line.strip():
                continue
            if data_line.startswith("chr"):
                data_line = data_line[3:]
            ret_list = []  # type: list[str]
            for new_line in split_vcf_line(data_line, ret_list):
                new_line = line_left_normalization(new_line)
                new_list = new_line.split("\t")
                key = "\t".join([new_list[0], new_list[1], new_list[3], new_list[4]])
                fp_out.write("{0}|{1}\n".format(key, curr_pos))


def db_has_col(cursor, table_name, col_name):
    cmd_str = "select sql from sqlite_master where type = 'table' and name = '{0}'".format(table_name)
    creat_sql = cursor.execute(cmd_str).fetchall()[0][0]
    # creat_sql_list = filter(lambda x: len(x) > 0 and "varchr" not in x, re.split(",| ", creat_sql))
    creat_sql_list = [i.strip("\"") for i in
                      filter(lambda x: len(x) > 0 and "varchar" not in x, re.split(",| |\t|\n", creat_sql))]
    if col_name in creat_sql_list:
        return True
    else:
        return False


def db_add_col(cursor, table_name, col_name, col_type):
    """
    Add column to database table
    @param cursor: cursor
    @param table_name: table name
    @param col_name: column name to add
    @param col_type: column properties
    @return: 0 Success, 1 already has this column, no need to add
    """
    cmd_str = "select sql from sqlite_master where type = 'table' and name = '{0}'".format(table_name)
    creat_sql = cursor.execute(cmd_str).fetchall()[0][0]
    # creat_sql_list = filter(lambda x: len(x) > 0 and "varchr" not in x, re.split(",| ", creat_sql))
    creat_sql_list = [i.strip("\"") for i in
                      filter(lambda x: len(x) > 0 and "varchar" not in x, re.split(",| |\t|\n", creat_sql))]
    if col_name in creat_sql_list:
        return 1
    cmd_str = "ALTER TABLE {0} ADD COLUMN '{1}' {2}".format(table_name, col_name, col_type)
    # print cmd_str
    cursor.execute(cmd_str)
    return 0


def create_db(union01_file, union001_file,
              annovar01, bystro01, dmis01, dsplicing01, spidex01, spliceai01, vep01,
              db_file):
    """
    Create a table with the data of union0.01
    :param union01_file:
    :return:
    """
    logging.basicConfig(filename="create_db.log", level=logging.DEBUG, format=log_format,
                        filemode="w")
    table_name = "variance"

    def db_insert(cursor, table, id, chr, pos, ref, alt, is001, is_annovar,
                  is_bystro, is_dmis, is_dsplicing, is_splidex, is_spliceai, is_vep):
        cursor.execute(
            "insert into {0} (id, chr, pos, ref, alt, is001, annovar, bystro, dmis, dsplicing, spidex, spliceAI, vep) "
            "values ('{1}', '{2}','{3}','{4}','{5}','{6}','{7}','{8}','{9}','{10}','{11}','{12}','{13}')".format(
                table, id, chr, pos, ref, alt, is001,
                is_annovar, is_bystro, is_dmis, is_dsplicing, is_splidex, is_spliceai, is_vep))

    def load_anno_set(anno_set, anno_vcf):
        with open(anno_vcf, "r") as fp:
            while True:
                data_line = fp.readline()
                if not data_line:
                    break
                if data_line.startswith("#") or not data_line.strip():
                    continue
                if data_line.startswith("chr"):
                    data_line = data_line[3:]
                data_list = data_line.strip().split("\t")
                anno_set.add("{0}\t{1}\t{2}\t{3}".format(data_list[0], data_list[1], data_list[3], data_list[4]))

    annovar01_set = set([])
    load_anno_set(annovar01_set, annovar01)
    bystro01_set = set([])
    load_anno_set(bystro01_set, bystro01)
    dmis01_set = set([])
    load_anno_set(dmis01_set, dmis01)
    dsplicing01_set = set([])
    load_anno_set(dsplicing01_set, dsplicing01)
    spidex01_set = set([])
    load_anno_set(spidex01_set, spidex01)
    spliceai01_set = set([])
    load_anno_set(spliceai01_set, spliceai01)
    vep01_set = set([])
    load_anno_set(vep01_set, vep01)
    # print spliceai01_set
    with open(union001_file, "r") as fp001:
        union001_set = set([])
        while True:
            data_line = fp001.readline()
            if not data_line:
                break
            if data_line.startswith("#") or not data_line.strip():
                continue
            if data_line.startswith("chr"):
                data_line = data_line[3:]
            union001_set.add(data_line.strip())
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()
    cursor.execute("DROP TABLE IF EXISTS {0}".format(table_name))
    cursor.execute(
        'create table {0} (id int primary key, chr varchar(20), pos int, ref varchr(40), alt varchr(40), '
        'is001 int, annovar varchar(1), bystro varchar(1), dmis varchar(1), dsplicing varchar(1), '
        'spidex varchar(1), spliceAI varchar(1), vep varchar(1))'.format(
            table_name))
    id = 0
    with open(union01_file, "r") as fp:
        while True:
            data_line = fp.readline()
            if not data_line:
                break
            if data_line.startswith("#") or not data_line.strip():
                continue

            data_line = data_line.strip()
            if data_line.startswith("chr"):
                data_line = data_line[3:]
            id += 1
            data_list = data_line.split("\t")
            # print data_list
            db_insert(cursor, table_name, id, data_list[0], data_list[1],
                      data_list[2], data_list[3], "1" if data_line in union001_set else "0",
                      "1" if data_line in annovar01_set else "0", "1" if data_line in bystro01_set else "0",
                      "1" if data_line in dmis01_set else "0", "1" if data_line in dsplicing01_set else "0",
                      "1" if data_line in spidex01_set else "0", "1" if data_line in spliceai01_set else "0",
                      "1" if data_line in vep01_set else "0")
    logging.debug("done added {} lines to db".format(id))

    cursor.close()
    conn.commit()
    conn.close()


def db_add_sample2table_variance(vcf_file, db_file, table_name):
    logging.basicConfig(filename="db_add_sample2table_variance.log", level=logging.DEBUG, format=log_format,
                        filemode="w")
    # table_name = "variance"
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()

    logging.debug("loading vcf index")
    vcf_index_dict = {}
    icounter = 0
    with open(vcf_file + "i", "r") as fp:
        while True:
            data_line = fp.readline()
            if not data_line:
                break
            if data_line.startswith("#") or not data_line.strip():
                continue
            icounter += 1
            if icounter > 0 and icounter % 1000000 == 0:
                logging.debug("handled {} vcf index".format(icounter))
            data_line = data_line.strip()
            data_list = data_line.split("|")
            vcf_index_dict[data_list[0]] = int(data_list[1])

    # load head list
    logging.debug("loading sample name")
    with open(vcf_file, "r") as fp:
        while True:
            data_line = fp.readline()
            if data_line.startswith("##") or not data_line.strip():
                continue
            if data_line.startswith("#"):
                head_list = data_line[1:].strip().split("\t")
                break
    head_list = ["sample_{}".format(i) for i in filter(lambda x: len(x) > 0, head_list[9:])]

    logging.debug("add sample columns to db")
    db_add_col(cursor, table_name, "vcf_id", "varchr(40)")
    db_add_col(cursor, table_name, "vcf_qual", "varchr(20)")
    db_add_col(cursor, table_name, "vcf_filter", "varchr(20)")
    db_add_col(cursor, table_name, "vcf_info", "varchr(20)")
    db_add_col(cursor, table_name, "vcf_format", "varchr(20)")
    for sample_name in head_list:
        db_add_col(cursor, table_name, sample_name, "varchr(2)")

    logging.debug("begin update sample data")
    icounter = 0
    with open(vcf_file, "r") as fp_vcf:
        sql_list = []
        for element in cursor.execute("SELECT id, chr, pos, ref, alt  FROM {};".format(table_name)):
            if icounter > 0 and icounter % 100 == 0:
                logging.debug("handled {} elements".format(icounter))
            icounter += 1
            id = element[0]
            element_key = "\t".join([str(i) for i in element[1:]])
            ref_alt_key = "{0}\t{1}".format(element[3], element[4])
            fp_vcf.seek(vcf_index_dict[element_key])
            org_vcf_list = fp_vcf.readline().strip().split("\t")
            org_vcf_list_short = org_vcf_list[:5]
            alt_list = org_vcf_list[4].split(",")
            vcf_ref_alt_list = []
            for alt in alt_list:
                org_vcf_list_short[4] = alt
                left_org_vcf_list_short = line_left_normalization("\t".join(org_vcf_list_short)).split("\t")
                vcf_ref_alt_list.append("{0}\t{1}".format(left_org_vcf_list_short[3], left_org_vcf_list_short[4]))
            key2num_dict = dict(zip(vcf_ref_alt_list, [str(i) for i in xrange(1, len(vcf_ref_alt_list) + 1)]))
            element_num = key2num_dict[ref_alt_key]  # type: str

            assert type(element_num) == str
            gene_type_list = [i.split(":")[0] for i in org_vcf_list[9:]]
            sql_list.append("UPDATE {0} SET vcf_id = '{1}' WHERE id = {2}".format(table_name, org_vcf_list[2], id))
            sql_list.append("UPDATE {0} SET vcf_qual = '{1}' WHERE id = {2}".format(table_name, org_vcf_list[5], id))
            sql_list.append("UPDATE {0} SET vcf_filter = '{1}' WHERE id = {2}".format(table_name, org_vcf_list[6], id))
            sql_list.append("UPDATE {0} SET vcf_info = '{1}' WHERE id = {2}".format(table_name, org_vcf_list[7], id))
            sql_list.append("UPDATE {0} SET vcf_format = '{1}' WHERE id = {2}".format(table_name, org_vcf_list[8], id))

            for index in xrange(len(gene_type_list)):
                gene_type_sep_list = gene_type_list[index].split("/")
                if "." in gene_type_list[index]:
                    update = "na"
                elif gene_type_sep_list[0] == element_num and gene_type_sep_list[1] == element_num:
                    update = "2"
                elif gene_type_sep_list[0] == element_num or gene_type_sep_list[1] == element_num:
                    update = "1"
                else:
                    update = "0"
                sql_list.append("UPDATE {0} SET {1} = '{2}' WHERE id = {3}".format(table_name,
                                                                                   head_list[index],
                                                                                   update, id))

        total_num = len(sql_list)
        icounter = 0
        for sql in sql_list:
            cursor.execute(sql)
            icounter += 1
            if icounter % 10000 == 0:
                logging.debug("excuted {0} / {1} sqls".format(icounter, total_num))
    cursor.close()
    conn.commit()
    conn.close()
    logging.debug("all done")


def check_vcf_missing(vcf_in):
    logging.basicConfig(filename="check_vcf_missing.log", level=logging.DEBUG, format=log_format,
                        filemode="w")
    logging.debug("check_vcf_missing begin")
    with open(vcf_in, "r") as fp:
        icounter = 0
        while True:
            vcf_line = fp.readline()
            if not vcf_line:
                break
            if vcf_line.startswith("#") or not vcf_line.strip():
                continue
            vcf_list = vcf_line.strip().split("\t")
            icounter += 1
            for index in xrange(len(vcf_list)):
                if index < 9:
                    continue
                gene_type_list = vcf_list[index].split(":")[0].split("/")
                if (gene_type_list[0] == '.' and gene_type_list[1] != '.') or (
                        gene_type_list[0] != '.' and gene_type_list[1] == '.'):
                    logging.debug(
                        """missing  {0}   {1}    {2}   {3} {4} gene_type_list={5}""".format(vcf_list[0], vcf_list[1],
                                                                                            vcf_list[2], vcf_list[3],
                                                                                            vcf_list[4],
                                                                                            gene_type_list))
                    break
            if icounter % 10000 == 0:
                logging.debug("handled {} vcf lines".format(icounter))
    logging.debug("all done")


def export_vcf_from_db(db_file, table_name, variance_restrict, org_vcf, output, fai_in):
    """

    :param db_file: database file
    :param table_name: variance
    :param org_vcf: The original vcf is used to extract the table header
    :param output: output vcf file
    :return: None
    """
    logging.basicConfig(filename="export_vcf_from_db.log", level=logging.DEBUG, format=log_format,
                        filemode="w")
    logging.debug("db_file=[{0}], table_name=[{1}], org_vcf=[{2}], output=[{3}], fai_in=[{4}], "
                  "variance_restrict=[{5}]".format(db_file, table_name, org_vcf, output, fai_in, variance_restrict))

    # def get_restrict_str(name, org_restrict, value):
    #     assert value in [-1, 0, 1, 2, 3]
    #     if value < 0:
    #         return org_restrict
    #     ret = org_restrict + (
    #         " WHERE ({0}={1}".format(name, value % 2) if not org_restrict else " {0} {1}='{2}'".format(
    #             "AND" if value < 2 else "OR", name, value % 2))
    #     return ret

    # handle condition
    # condition_str = ""
    #
    # condition_str = get_restrict_str("annovar", condition_str, is_annovar)
    # condition_str = get_restrict_str("bystro", condition_str, is_bystro)
    # condition_str = get_restrict_str("dmis", condition_str, is_dmis)
    # condition_str = get_restrict_str("dsplicing", condition_str, is_dsplicing)
    # condition_str = get_restrict_str("spidex", condition_str, is_spidex)
    # condition_str = get_restrict_str("spliceAI", condition_str, is_spliceAI)
    # condition_str = get_restrict_str("vep", condition_str, is_vep)
    #
    # if is_001 >= 0:
    #     condition_str += " WHERE is001={}".format(
    #         is_001 % 2) if not condition_str else ") AND is001='{0}'".format(is_001 % 2)
    # elif condition_str:
    #     condition_str += ")"
    #
    # logging.debug("condition_str=[{}]".format(condition_str))
    # load head from org v
    with open(org_vcf, "r") as fp:
        while True:
            org_vcf_line = fp.readline()
            if not org_vcf_line:
                break
            if org_vcf_line.startswith("##"):
                continue
            if org_vcf_line.startswith("#"):
                head_line = org_vcf_line
                break
    head_list = head_line.strip().split("\t")

    # handle db
    with open(output, "w") as fp_out:

        pp = Popen(["awk '{printf(\"##contig=<ID=\"$1\",length=\"$2\">\\n\");}' " + "{0}".format(fai_in)], stdout=PIPE,
                   shell=True)
        contigs = pp.stdout.readlines()
        fp_out.write("##fileformat=VCFv4.0\n##assembly=GRCh38/hg38\n")
        for contig in contigs:
            fp_out.write(contig)
        fp_out.write(head_line)
        conn = sqlite3.connect(db_file)
        cursor = conn.cursor()
        format_variance_restrict = "" if not variance_restrict else " WHERE {}".format(variance_restrict)
        cmd_str = "SELECT v.chr, v.pos, v.vcf_id, v.ref, v.alt, v.vcf_qual, " \
                  "v.vcf_filter, v.vcf_info, v.vcf_format, {0} FROM " \
                  "{1} AS v{2}".format(", ".join(["v.sample_{}".format(i) for i in head_list[9:]]),
                                       table_name,
                                       format_variance_restrict)
        logging.debug("variance table = {}".format(table_name))
        logging.debug("variance restrict = {}".format(format_variance_restrict))
        # logging.debug("sql=[{}]".format(cmd_str))
        for element in cursor.execute(cmd_str):
            list_element = list(element)
            for index in xrange(len(element)):
                if index < 9:
                    continue
                # logging.debug("list_element[index]=[{0}] type=[{1}]".format(list_element[index],type(list_element[index])))
                if list_element[index] == 0:
                    list_element[index] = "0/0"
                elif list_element[index] == 1:
                    list_element[index] = "0/1"
                elif list_element[index] == 2:
                    list_element[index] = "1/1"
                else:
                    list_element[index] = "./."
            list_element[8] = list_element[8].split(":")[0]  # vcf_format
            tmp_str = "{}\n".format("\t".join([str(i) for i in list_element]))
            if not tmp_str.startswith("chr"):
                tmp_str = "chr{0}".format(tmp_str)
            fp_out.write(tmp_str)
        cursor.close()
        conn.commit()
        conn.close()
    logging.debug("all done")


def db_add_spliceAI_anno(db_file, annotation):
    logging.basicConfig(filename="db_add_spliceAI_anno.log", level=logging.DEBUG, format=log_format,
                        filemode="w")

    def get_score(anno):
        tmp_dict = eval("{\"" + anno.replace("|", "\",\"").replace("=", "\":\"") + "\"}")
        return max([float(i) for i in [tmp_dict["DS_DG"], tmp_dict["DS_AG"], tmp_dict["DS_AL"], tmp_dict["DS_DL"]]])

    anno_dict = {}
    anno_score_dict = {}
    with open(annotation, "r") as fp:
        while True:
            data_line = fp.readline()
            if not data_line:
                break
            if data_line.startswith("#") or not data_line.strip():
                continue
            data_list = data_line.strip().split("\t")
            anno = data_list[7]
            key = "{}\t{}\t{}\t{}".format(data_list[0].strip("chr"), data_list[1], data_list[3], data_list[4])
            logging.debug("key={0} anno={1}".format(key, anno))
            score = max([get_score(i) for i in anno.split(";")])

            anno_dict[key] = anno
            anno_score_dict[key] = score
    table_name = "variance"
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()
    db_add_col(cursor, table_name, "spliceAI_anno", "varchr(500)")
    db_add_col(cursor, table_name, "spliceAI_score", "varchr(20)")
    sql_list = []
    for element in cursor.execute("SELECT id, chr, pos, ref, alt  FROM {} WHERE spliceAI=1".format(table_name)):
        id = element[0]
        element_key = "\t".join([str(i) for i in element[1:]])
        if element_key in anno_dict:
            sql_list.append(
                "UPDATE {0} SET spliceAI_anno = '{1}' WHERE id = {2}".format(table_name, anno_dict[element_key], id))
            sql_list.append(
                "UPDATE {0} SET spliceAI_score = '{1}' WHERE id = {2}".format(table_name, anno_score_dict[element_key],
                                                                              id))
        else:
            logging.debug("{} not in anno_dict")
    for sql in sql_list:
        cursor.execute(sql)
    cursor.close()
    conn.commit()
    conn.close()
    logging.debug("all done")


def db_add_spidex_anno(db_file, spidex_file):
    logging.basicConfig(filename="db_add_spidex_anno.log", level=logging.DEBUG, format=log_format, filemode="w")
    key_dict = {}
    with open(spidex_file, "r") as fp:
        while True:
            data_line = fp.readline()
            if not data_line:
                break
            data_list = data_line.strip().split("\t")
            if len(data_list) < 5 or data_list[0] == "NA":
                continue
            key = "_".join(data_list[:4])
            if key not in key_dict:
                key_dict[key] = [data_list[4]]
            else:
                key_dict[key].append(data_list[4])
                logging.debug("{} repeat".format(key))
    table_name = "variance"
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()

    if not len(cursor.execute("select * from sqlite_master where tbl_name='{}' and sql like '%spidex_score%';".format(
            table_name)).fetchall()):
        db_add_col(cursor, table_name, "spidex_score", "varchr(20)")
    sql_list = []
    for element in cursor.execute("SELECT id, chr, pos, ref, alt  FROM {} WHERE spidex=1".format(table_name)):
        id = element[0]
        element_key = "_".join([str(i) for i in element[1:]])
        sql_list.append(
            "UPDATE {0} SET spidex_score = '{1}' WHERE id = {2}".format(table_name, key_dict[element_key][0], id))

    icounter = 0
    sql_list_len = len(sql_list)
    for sql in sql_list:
        if icounter > 0 and icounter % 10000 == 0:
            logging.debug("executed {0} / {1} sqls".format(icounter, sql_list_len))
        cursor.execute(sql)
        icounter += 1

    cursor.close()
    conn.commit()
    conn.close()
    logging.debug("all done")


def db_add_vcf_pos(db_file, vcf_file):
    logging.basicConfig(filename="db_add_vcf_pos.log", level=logging.DEBUG, format=log_format, filemode="w")
    table_name = "variance"
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()
    logging.debug("loading vcf index")
    vcf_index_dict = {}
    icounter = 0
    with open(vcf_file + "i", "r") as fp:
        while True:
            data_line = fp.readline()
            if not data_line:
                break
            if data_line.startswith("#") or not data_line.strip():
                continue
            icounter += 1
            if icounter > 0 and icounter % 1000000 == 0:
                logging.debug("handled {} vcf index".format(icounter))
            data_line = data_line.strip()
            data_list = data_line.split("|")
            vcf_index_dict[data_list[0]] = int(data_list[1])
    if not len(cursor.execute("select * from sqlite_master where tbl_name='{}' and sql like '%vcf_pos%';".format(
            table_name)).fetchall()):
        db_add_col(cursor, table_name, "vcf_pos", "varchr(20)")
    sql_list = []
    icounter = 0
    with open(vcf_file, "r") as fp_vcf:
        for element in cursor.execute("SELECT id, chr, pos, ref, alt  FROM {};".format(table_name)):
            if icounter > 0 and icounter % 100 == 0:
                logging.debug("handled {} elements".format(icounter))
            icounter += 1
            id = element[0]
            element_key = "\t".join([str(i) for i in element[1:]])
            fp_vcf.seek(vcf_index_dict[element_key])
            org_vcf_list = fp_vcf.readline().strip().split("\t")
            sql_list.append("UPDATE {0} SET vcf_pos = '{1}' WHERE id = {2}".format(table_name, org_vcf_list[1], id))
    for sql in sql_list:
        cursor.execute(sql)

    cursor.close()
    conn.commit()
    conn.close()
    logging.debug("all done")


def db_add_bystro_anno(db_file, tsv_file, table_name):
    def get_max_in_anno(str_in):
        # if "!" in str_in:
        #     logging.debug("str_in= {}".format(str_in))
        anno_list = filter(lambda x: x != "!", re.split('[;|]', str_in))
        if not anno_list:
            return "!"
        else:
            return str(max([float(i) for i in anno_list]))

    logging.basicConfig(filename="db_add_bystro_anno.log", level=logging.DEBUG, format=log_format, filemode="w")
    # table_name = "variance"
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()
    with open(tsv_file, "r") as fp:
        tsv_data = [i.split("\t") for i in fp.readlines()]
        tsv_dict = {}
        for tsv_list in tsv_data[1:]:
            [ref, alt] = bystro_anno_info.build_annovar_ref_alt_key_from_bystro(tsv_list[16], tsv_list[4])
            vcf_key = "{0}_{1}_{2}_{3}".format(tsv_list[0], tsv_list[15], ref, alt)  # chr_vcfpos_annovarref_annovaralt
            if vcf_key not in tsv_dict:
                tsv_dict[vcf_key] = [[get_max_in_anno(tsv_list[14]),
                                      get_max_in_anno(tsv_list[46]),
                                      get_max_in_anno(tsv_list[47]),
                                      get_max_in_anno(tsv_list[48])]]
            else:
                logging.debug("{} more than one".format(vcf_key))
                tsv_dict[vcf_key].append([get_max_in_anno(tsv_list[14]),
                                          get_max_in_anno(tsv_list[46]),
                                          get_max_in_anno(tsv_list[47]),
                                          get_max_in_anno(tsv_list[48])])

    # add columns
    sql = "select * from sqlite_master where tbl_name='{}' and sql like '%bystro_sampleMaf%';".format(table_name)
    if not len(cursor.execute(sql).fetchall()):
        db_add_col(cursor, table_name, "bystro_sampleMaf", "varchr(20)")

    sql = "select * from sqlite_master where tbl_name='{}' and sql like '%bystro_phastCons%';".format(table_name)
    if not len(cursor.execute(sql).fetchall()):
        db_add_col(cursor, table_name, "bystro_phastCons", "varchr(20)")
    sql = "select * from sqlite_master where tbl_name='{}' and sql like '%bystro_phyloP%';".format(table_name)
    if not len(cursor.execute(sql).fetchall()):
        db_add_col(cursor, table_name, "bystro_phyloP", "varchr(20)")
    sql = "select * from sqlite_master where tbl_name='{}' and sql like '%bystro_cadd%';".format(table_name)
    if not len(cursor.execute(sql).fetchall()):
        db_add_col(cursor, table_name, "bystro_cadd", "varchr(20)")

    sql_list = []
    icounter = 0
    for element in cursor.execute("SELECT id, chr, pos, ref, alt  FROM {}".format(table_name)):  # WHERE bystro=1
        if icounter > 0 and icounter % 100 == 0:
            logging.debug("handled {} elements".format(icounter))
        icounter += 1
        id = element[0]
        [new_ref, new_alt, ileft_delete] = left_normalization_one_alt_annovar(element[3], element[4])
        element_key = "chr{0}_{1}_{2}_{3}".format(element[1], element[2], new_ref, new_alt)
        if element_key not in tsv_dict:
            logging.debug("{} not in dict".format(element_key))
            continue
        sql_list.append(
            "UPDATE {0} SET bystro_sampleMaf = '{1}' WHERE id = {2}".format(table_name, tsv_dict[element_key][0][0],
                                                                            id))
        sql_list.append(
            "UPDATE {0} SET bystro_phastCons = '{1}' WHERE id = {2}".format(table_name, tsv_dict[element_key][0][1],
                                                                            id))
        sql_list.append(
            "UPDATE {0} SET bystro_phyloP = '{1}' WHERE id = {2}".format(table_name, tsv_dict[element_key][0][2], id))
        sql_list.append(
            "UPDATE {0} SET bystro_cadd = '{1}' WHERE id = {2}".format(table_name, tsv_dict[element_key][0][3], id))

    icounter = 0
    sql_list_len = len(sql_list)
    for sql in sql_list:
        if icounter > 0 and icounter % 10000 == 0:
            logging.debug("executed {0} / {1} sqls".format(icounter, sql_list_len))
        cursor.execute(sql)
        icounter += 1
    cursor.close()
    conn.commit()
    conn.close()
    logging.debug("all done")


def check_dbNSFP(anno):
    icounter = 0
    inum = 0
    logging.basicConfig(filename="check_dbNSFT.log", level=logging.DEBUG, format=log_format, filemode="w")
    with open(anno, "r") as fp:
        while True:
            data_line = fp.readline()
            icounter += 1
            if icounter > 0 and icounter % 100000 == 0:
                print("handled {} lines".format(icounter))
            if not data_line:
                break
            if data_line.startswith("#") or not data_line.strip():
                continue
            data_list = data_line.strip().split("\t")
            if ";" not in data_list[12]:
                continue
            inum += 1
            target_list = data_list[410:]
            target_str = "\t".join([i.strip(".") for i in target_list]).strip()
            if len(target_str) > 0:
                logging.debug("{0}\t[{1}]".format("_".join(data_list[:4]), "\t".join(data_list[410:])))
    print("{} multiple gene lines".format(inum))


def check_protein_coding_gene(file_in, ):
    logging.basicConfig(filename="check_protein_coding_gene.log", level=logging.DEBUG, format=log_format, filemode="w")
    total_dict = {}
    with open(file_in, "r") as fp:
        while True:
            data_line = fp.readline()
            if not data_line:
                break
            if data_line.startswith("#") or not data_line.strip():
                continue
            data_list = data_line.strip().split("\t")
            # print data_list
            chr = data_list[0]
            start = data_list[3]
            end = data_list[4]
            my_dict = eval("{\"" + data_list[8].replace("; ", ",\"").replace(" ", "\":").strip(";") + "}")
            # print my_dict
            gene_id = my_dict["gene_id"]
            gene_name = my_dict["gene_name"]
            key = "{0}_{1}".format(gene_id, gene_name)
            value = "{0}_{1}_{2}".format(chr, start, end)
            if key not in total_dict:
                total_dict[key] = [value]
            else:
                total_dict[key].append(value)
    icounter = 0
    icounter2 = 0
    for key in total_dict:
        icounter2 += 1
        if len(total_dict[key]) == 1:
            continue
        icounter += 1
        for i in total_dict[key]:
            logging.debug("{0}\t{1}".format(key, i))
        logging.debug("-----------------------")
    print("repeat {0} / {1}".format(icounter, icounter2))


def check_protein_coding_transcript(file_in, log_name="check_protein_coding_transcript.log"):
    def is_overlap(key1_list, key2_list):
        [chr1, start1, end1] = key1_list[:3]
        [chr2, start2, end2] = key2_list[:3]
        if chr1 != chr2:
            return False
        start1 = int(start1)
        end1 = int(end1)
        start2 = int(start2)
        end2 = int(end2)
        if end1 < start2 or start1 > end2:
            # logging.debug("{0}\t{1}\t{2}\t{3} seperate".format(start1,end1,start2,end2))
            return False
        return True

    def append_region(key1_list, key2_list):
        [chr1, start1, end1] = key1_list[:3]
        [chr2, start2, end2] = key2_list[:3]
        start1 = int(start1)
        end1 = int(end1)
        start2 = int(start2)
        end2 = int(end2)
        chr = chr1
        start = min(start1, start2)
        end = max(end1, end2)
        return "_".join([chr, str(start), str(end)])

    logging.basicConfig(filename=log_name, level=logging.DEBUG, format=log_format,
                        filemode="w")
    total_dict = {}
    name_id_dict = {}
    id_name_dict = {}
    with open(file_in, "r") as fp:
        while True:
            data_line = fp.readline()
            if not data_line:
                break
            if data_line.startswith("#") or not data_line.strip():
                continue
            data_list = data_line.strip().split("\t")
            # print data_list
            chr = data_list[0]
            start = data_list[3]
            end = data_list[4]
            my_dict = eval("{\"" + data_list[8].replace("; ", ",\"").replace(" ", "\":").strip(";") + "}")
            # print my_dict
            gene_id = my_dict["gene_id"].split(".")[0]
            gene_name = my_dict["gene_name"]
            value = "{0}_{1}_{2}".format(chr, start, end)

            if gene_id not in total_dict:
                total_dict[gene_id] = [value]
            else:
                appended = False
                for index in xrange(len(total_dict[gene_id])):
                    if is_overlap([chr, start, end], total_dict[gene_id][index].split("_")):
                        total_dict[gene_id][index] = append_region([chr, start, end],
                                                                   total_dict[gene_id][index].split("_"))
                        appended = True
                        break
                if not appended:
                    total_dict[gene_id].append(value)

            if gene_name not in name_id_dict:
                name_id_dict[gene_name] = [gene_id]
            else:
                if gene_id not in name_id_dict[gene_name]:
                    name_id_dict[gene_name].append(gene_id)

            if gene_id not in id_name_dict:
                id_name_dict[gene_id] = [gene_name]
            else:
                if gene_name not in id_name_dict[gene_id]:
                    id_name_dict[gene_id].append(gene_name)
    icounter = 0
    icounter2 = 0
    for key in total_dict:
        icounter2 += 1
        if len(total_dict[key]) == 1:
            continue
        icounter += 1
        for i in total_dict[key]:
            logging.debug("{0}\t{1}".format(key, i))
        logging.debug("-----------------------")
    print("idmultipleintervals{0} / {1}".format(icounter, icounter2))
    icounter = 0
    icounter2 = 0
    for key in name_id_dict:
        icounter2 += 1
        if len(name_id_dict[key]) > 1:
            icounter += 1
            logging.debug("gene_id {0}\t{1}".format(key, name_id_dict[key]))
    print("gene多id {0} / {1}".format(icounter, icounter2))
    icounter = 0
    icounter2 = 0
    for key in id_name_dict:
        icounter2 += 1
        if len(id_name_dict[key]) > 1:
            icounter += 1
            logging.debug("id_gene {0}\t{1}".format(key, id_name_dict[key]))
    print("id多gene {0} / {1}".format(icounter, icounter2))
    return [total_dict, id_name_dict]


def check_dms(dms_in):
    logging.basicConfig(filename="check_dms.log", level=logging.DEBUG, format=log_format, filemode="w")
    logging.debug("begin")
    gene_dict = {}
    # icounter = 0
    with open(dms_in, "r") as fp:
        while True:
            data_line = fp.readline()
            # icounter += 1
            # if icounter % 10000 == 0:
            #     print "handled {} lines".format(icounter)
            if not data_line:
                # print "breaked"
                break
            if data_line.startswith("#") or not data_line.strip():
                continue
            data_list = data_line.strip().split("\t")
            # print data_list

            gene_name = data_list[12]
            tx_name = data_list[1].split(".")[0]
            key = "{0},{1}".format(gene_name, tx_name)

            if key not in gene_dict:
                gene_dict[key] = [data_line]
                # logging.debug(key)
            else:
                gene_dict[key].append(data_line)
    # print 1111
    # exit(0)
    # print gene_dict["NM_001005221"]
    icounter = 0
    icounter2 = 0
    for key in gene_dict:
        icounter2 += 1
        if len(gene_dict[key]) > 1:
            icounter += 1
            for i in gene_dict[key]:
                logging.debug(i)
            logging.debug("--------------")
    print("{0} / {1} 个repeatgene name and tx name's".format(icounter, icounter2))


def split_gene_centric_annotation(file_in, output):
    cmd_str = "sort -n {} | uniq > ii".format(file_in)
    pp = Popen([cmd_str], shell=True)
    pp.wait()
    with open("ii", "r") as fp, open("tmp", "w") as fp_out:

        ret_set = set([])
        icounter = 0
        while True:
            data_line = fp.readline()
            if not data_line:
                break
            if data_line.startswith("#") or data_line.startswith("genename") or data_line.startswith("\n"):
                continue
            ret_list = []
            icounter += 1
            if icounter % 1000 == 0:
                print("handled {} lines".format(icounter))
            data_list = data_line.split("\t")
            # print data_list
            num = len(data_list[1].split(";"))
            for index in xrange(num):
                tmp = []
                for i in data_list:
                    splited_i = i.split(";")
                    if len(splited_i) != num:
                        tmp.append(i)
                    else:
                        tmp.append(splited_i[index])
                ret_str = "\t".join(tmp)
                if ret_str not in ret_set:
                    ret_set.add(ret_str)
                    ret_list.append("\t".join(tmp))
            fp_out.write("".join(ret_list))
    cmd_str = "sort -n tmp | uniq > {}".format(output)
    pp = Popen([cmd_str], shell=True)
    pp.wait()


def db_build_synonymous_snp_table(db_file, table_name, annovar01, bystro01,
                                  vep01, union001_file, union01_file):
    logging.basicConfig(filename="db_build_synonymous_snp_table.log", level=logging.DEBUG, format=log_format,
                        filemode="w")

    # table_name = "synonymous_snp"

    def db_insert(cursor, table, id, chr, pos, ref, alt, is001, is_annovar, is_bystro, is_vep):
        cursor.execute("insert into {0} (id, chr, pos, ref, alt, is001, annovar, bystro, vep) values "
                       "('{1}', '{2}','{3}','{4}','{5}','{6}','{7}','{8}','{9}')".format(table, id, chr, pos,
                                                                                         ref, alt, is001,
                                                                                         is_annovar, is_bystro,
                                                                                         is_vep))

    def load_anno_set(anno_set, anno_key):
        with open(anno_key, "r") as fp:
            while True:
                data_line = fp.readline()
                if not data_line:
                    break
                if data_line.startswith("#") or not data_line.strip():
                    continue
                if data_line.startswith("chr"):
                    data_line = data_line[3:]
                # data_list = data_line.strip().split("\t")
                anno_set.add(data_line.strip())

    # logging.debug("loading vcf index")
    # vcf_index_dict = {}
    # icounter = 0
    # with open(vcf_file + "i", "r") as fp:
    #     while True:
    #         data_line = fp.readline()
    #         if not data_line:
    #             break
    #         if data_line.startswith("#") or not data_line.strip():
    #             continue
    #         icounter += 1
    #         if icounter > 0 and icounter % 1000000 == 0:
    #             logging.debug("handled {} vcf index".format(icounter))
    #         data_line = data_line.strip()
    #         data_list = data_line.split("|")
    #         vcf_index_dict[data_list[0]] = int(data_list[1])

    annovar01_set = set([])
    load_anno_set(annovar01_set, annovar01)
    bystro01_set = set([])
    load_anno_set(bystro01_set, bystro01)
    vep01_set = set([])
    load_anno_set(vep01_set, vep01)
    union001_set = set([])
    load_anno_set(union001_set, union001_file)

    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()
    cursor.execute("DROP TABLE IF EXISTS {0}".format(table_name))
    cursor.execute(
        'create table {0} (id int primary key, chr varchar(20), pos int, '
        'ref varchr(40), alt varchr(40), is001 int, '
        'annovar varchar(1), bystro varchar(1), vep varchar(1))'.format(table_name))
    id = 0
    with open(union01_file, "r") as fp:  # open(vcf_file, "r") as fp_vcf
        while True:
            data_line = fp.readline()
            if not data_line:
                break
            if data_line.startswith("#") or not data_line.strip():
                continue
            data_line = data_line.strip()
            if data_line.startswith("chr"):
                data_line = data_line[3:]
            # fp_vcf.seek(vcf_index_dict[data_line])
            # org_vcf_line = fp_vcf.readline().strip()
            # org_vcf_id = org_vcf_line.split("\t")[2]
            id += 1
            data_list = data_line.split("\t")
            # print data_list
            db_insert(cursor, table_name, id, data_list[0], data_list[1],
                      data_list[2], data_list[3],
                      "1" if data_line in union001_set else "0",
                      "1" if data_line in annovar01_set else "0",
                      "1" if data_line in bystro01_set else "0",
                      "1" if data_line in vep01_set else "0")
            if id % 1000 == 0:
                logging.debug("added {} lines to db".format(id))
    logging.debug("done added {} lines to db".format(id))

    cursor.close()
    conn.commit()
    conn.close()


def db_build_gene_table(protein_coding_gene_transcript, db_file, table_name):
    def db_insert(cursor, table, id, gene_id, gene_name,
                  chr, start_pos, end_pos,
                  chr2=None, start_pos2=-1, end_pos2=-1):
        if chr2:
            cursor.execute(
                "insert into {0} (id, gene_id, gene_name, chr, start_pos, end_pos, chr2, start_pos2, end_pos2) values "
                "('{1}', '{2}','{3}','{4}','{5}','{6}','{7}','{8}','{9}')".format(table, id, gene_id, gene_name,
                                                                                  chr, start_pos, end_pos,
                                                                                  chr2, start_pos2, end_pos2))
        else:
            cursor.execute(
                "insert into {0} (id, gene_id, gene_name, chr, start_pos, end_pos) values "
                "('{1}', '{2}','{3}','{4}','{5}','{6}')".format(table, id, gene_id, gene_name,
                                                                chr, start_pos, end_pos))

    [total_dict, id_name_dict] = check_protein_coding_transcript(protein_coding_gene_transcript,
                                                                 "db_build_gene_table.log")
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()
    cursor.execute("DROP TABLE IF EXISTS {0}".format(table_name))
    cursor.execute('create table {0} (id int primary key, gene_id varchar(40), gene_name varchr(20), '
                   'chr varchr(20), start_pos int, end_pos int, '
                   'chr2 varchr(20), start_pos2 int, end_pos2 int)'.format(table_name))

    id = 0
    for gene_id in total_dict:
        id += 1
        if len(total_dict[gene_id]) == 1:
            [chr, start, end] = total_dict[gene_id][0].split("_")
            chr = chr.strip("chr")
            db_insert(cursor, table_name, id, gene_id, id_name_dict[gene_id][0],
                      chr, start, end)
        elif len(total_dict[gene_id]) == 2:
            [chr, start, end] = total_dict[gene_id][0].split("_")
            [chr2, start2, end2] = total_dict[gene_id][1].split("_")
            chr = chr.strip("chr")
            chr2 = chr2.strip("chr")
            db_insert(cursor, table_name, id, gene_id, id_name_dict[gene_id][0],
                      chr, start, end,
                      chr2, start2, end2)
        else:
            logging.debug("Can not handle it. {0} has more than 2 region. {1}".format(gene_id, total_dict[gene_id]))
        if id % 1000 == 0:
            logging.debug("handled {} gene_ids".format(id))
    logging.debug("all done added {} lines to db".format(id))
    cursor.close()
    conn.commit()
    conn.close()


def db_copy_table(db_file, table_name, table_file):
    logging.basicConfig(filename="db_copy_table.log", level=logging.DEBUG, format=log_format, filemode="w")

    def db_insert(cursor, table, id, head_list, data_list):
        cmd_str = "insert into {0} (id, ".format(table)
        tmp_str = "values ('{}', ".format(id)
        for index in xrange(len(head_list)):
            cmd_str = "{0}{1}, ".format(cmd_str, head_list[index])
            tmp_str = "{0}'{1}', ".format(tmp_str, data_list[index])
        cmd_str = cmd_str.strip(", ")
        tmp_str = tmp_str.strip(", ")
        cmd_str = "{0}){1})".format(cmd_str, tmp_str)
        cursor.execute(cmd_str)

        # cursor.execute(
        #     "insert into {0} (id, gene_id, gene_name, chr, start_pos, end_pos, chr2, start_pos2, end_pos2) values "
        #     "('{1}', '{2}','{3}','{4}','{5}','{6}','{7}','{8}','{9}')".format(table, id, gene_id, gene_name,
        #                                                                       chr, start_pos, end_pos,
        #                                                                           chr2, start_pos2, end_pos2))

    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()
    with open(table_file, "r") as fp_table:
        head_line = fp_table.readline()
        head_list = head_line.strip().split("\t")
        cursor.execute("DROP TABLE IF EXISTS {0}".format(table_name))
        cmd_str = 'create table {0} (id int primary key, '.format(
            table_name)  # gene_id varchar(40), gene_name varchr(20), chr varchr(20), start_pos int, end_pos int, chr2 varchr(20), start_pos2 int, end_pos2 int)
        for col_name in head_list:
            cmd_str = "{0} {1} varchr(20), ".format(cmd_str, col_name)
        cmd_str = cmd_str.strip(", ")
        cmd_str = "{})".format(cmd_str)
        print(cmd_str)
        cursor.execute(cmd_str)
        id = 0
        while True:
            data_line = fp_table.readline()
            if not data_line:
                break
            id += 1
            if id % 100 == 0:
                logging.debug("handled {} lines".format(id))
            data_list = data_line.strip().split("\t")
            db_insert(cursor, table_name, id, head_list, data_list)
    logging.debug("all done {} lines".format(id))
    cursor.close()
    conn.commit()
    conn.close()


spliter_dict = {"tab": "\t", "comma": ","}


def db_add_gene_anno(db_file, table_name, gene_file, gene_id_col, spliter):
    gene_id_col = int(gene_id_col)
    gene_name_col = 3 - gene_id_col
    logging.basicConfig(filename="{}.log".format(sys._getframe().f_code.co_name),
                        level=logging.DEBUG, format=log_format, filemode="w")
    logging.debug("loading gene_file")
    if spliter not in spliter_dict:
        logging.error("spliter should in dict {0}".format(spliter_dict))
        return
    spliter = spliter_dict[spliter]
    gene_dict = {}
    with open(gene_file, "r") as fp:
        head_line = fp.readline().strip()
        head_list = head_line.split(spliter)[2:]
        for i in xrange(len(head_list)):
            if "(" in head_list[i]:
                head_list[i] = head_list[i].replace("(", "_")
            if ")" in head_list[i]:
                head_list[i] = head_list[i].replace(")", "_")
            if "/" in head_list[i]:
                head_list[i] = head_list[i].replace("/", "_")
            if "-" in head_list[i]:
                head_list[i] = head_list[i].replace("-", "_")
            if "." in head_list[i]:
                head_list[i] = head_list[i].replace(".", "_")
            if " " in head_list[i]:
                head_list[i] = head_list[i].replace(" ", "_")
        # print head_list
        icounter = 0
        while True:
            data_line = fp.readline()
            if not data_line:
                break
            icounter += 1
            if icounter % 1000 == 0:
                logging.debug("handled {} lines".format(icounter))
            data_list = data_line.strip().split(spliter)
            gene_id_list = data_list[gene_id_col - 1].split(";")
            for gene_id in gene_id_list:
                if gene_id not in gene_dict:
                    gene_dict[gene_id] = data_list
                else:
                    logging.error("{0} repeat!!!!".format(gene_id))
    logging.debug("load gene file done")
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()
    for col_name in head_list:
        # sql = "select * from sqlite_master where tbl_name='{0}' and sql like '%{1}%';".format(table_name, col_name)
        # if not len(cursor.execute(sql).fetchall()):
        iret = db_add_col(cursor, table_name, col_name, "varchar(1024)")
        if iret != 0:
            logging.info("can not add {0} to table {1}".format(col_name, table_name))

    sql_list = []
    for element in cursor.execute("SELECT id, gene_id, gene_name FROM {0}".format(table_name)):
        id = element[0]
        gene_id = element[1]
        gene_name = element[2]
        if gene_id not in gene_dict:
            # logging.error("{0} not in dict".format(gene_id))
            continue
        data_list = gene_dict[gene_id]
        if gene_name != data_list[gene_name_col - 1]:
            logging.error(
                "different gene name[{0}] [{1}] id=[{2}]".format(gene_name, data_list[gene_name_col - 1], gene_id))
            # continue
        for index in xrange(len(head_list)):
            if "\"" in data_list[index + 2] and "'" in data_list[index + 2]:
                data_list[index + 2] = data_list[index + 2].replace("\"", "'")
            sql_list.append("UPDATE {0} SET {1} = {4}{2}{4} WHERE id = {3}".format(table_name,
                                                                                   head_list[index],
                                                                                   data_list[index + 2],
                                                                                   id,
                                                                                   "\"" if "'" in data_list[
                                                                                       index + 2] else "'"))
    icounter = 0
    sql_list_len = len(sql_list)
    logging.debug("begin to execute sqls {0}".format(sql_list_len))
    for sql in sql_list:
        if icounter > 0 and icounter % 10000 == 0:
            logging.debug("executed {0} / {1} sqls".format(icounter, sql_list_len))
        # print sql
        try:
            cursor.execute(sql)
        except:
            print(sql)
        logging.debug("sql={0}".format(sql))
        icounter += 1
    cursor.close()
    conn.commit()
    conn.close()
    logging.debug("all done")


def db_reduce_col_duplicate(db_file, table_name, col_name):
    logging.basicConfig(filename="{}.log".format(sys._getframe().f_code.co_name),
                        level=logging.DEBUG, format=log_format, filemode="w")
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()
    cursor.execute("SELECT count() FROM {}".format(table_name))
    line_num = cursor.fetchall()[0][0]

    icounter = 0
    for i in xrange(1, line_num, 1):
        cursor.execute("SELECT {0} FROM {1} WHERE id={2}".format(col_name, table_name, i))
        query_ret = cursor.fetchall()
        icounter += 1
        if icounter % 1000 == 0:
            logging.debug("handled {} lines".format(icounter))

        if not query_ret:
            continue
        vcf_id = query_ret[0][0]
        # print vcf_id
        cmd_str = "SELECT id FROM {0} WHERE {1}='{2}'".format(table_name, col_name, vcf_id)
        # print cmd_str
        cursor.execute(cmd_str)
        id_list = cursor.fetchall()
        if len(id_list) == 1:
            continue
        # print id_list
        for j in xrange(len(id_list)):
            cur_vcf_id = "{0}_{1}".format(vcf_id, j + 1)
            cmd_str = "UPDATE {0} SET {1} = '{2}' WHERE id = {3}".format(table_name, col_name, cur_vcf_id,
                                                                         id_list[j][0])
            # print cmd_str
            cursor.execute(cmd_str)
    cursor.close()
    conn.commit()
    conn.close()
    logging.debug("all done")
    # for element in cursor.execute("SELECT id, {1} FROM {0}".format(table_name, col_name)):


def collapse(db_file, variance_table_name, gene_table_name, variance_restrict, gene_restrict, output):
    """
    For each gene, select the snp within its range from the snp table
    variance_restrict eg: v.is001=1 AND v.bystro=1
    gene_restrict eg: g.gene_name IN ('SOCS4', 'RAB22A')
    :type output: str
    :type db_file: str
    :type gene_restrict: str
    :type variance_restrict: str
    :type variance_table_name: str
    """
    logging.basicConfig(filename="{}.log".format(sys._getframe().f_code.co_name),
                        level=logging.DEBUG, format=log_format, filemode="w")
    logging.debug("begin")
    logging.debug("db_file=[{}]".format(db_file))
    logging.debug("variance_table_name=[{}]".format(variance_table_name))
    logging.debug("gene_table_name=[{}]".format(gene_table_name))
    logging.debug("variance_restrict=[{}]".format(variance_restrict))
    logging.debug("gene_restrict=[{}]".format(gene_restrict))
    logging.debug("output=[{}]".format(output))
    if variance_table_name not in ['variance', 'synonymous_snp']:
        logging.debug("illegal table name [{0}]".format(variance_table_name))
        return

    # tempfile = StringIO()
    # source = sqlite3.connect(db_file)
    # for line in source.iterdump():
    #     tempfile.write('%s\n' % line)
    # source.close()
    # tempfile.seek(0)
    #
    # conn = sqlite3.connect(':memory:')
    # conn.cursor().executescript(tempfile.read())
    # conn.commit()
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()
    cursor.execute("SELECT count() FROM {0}".format(variance_table_name))
    tmp_num = cursor.fetchall()[0][0]
    print("all variance number in table {0} is {1}".format(variance_table_name, tmp_num))
    cursor.execute("SELECT count() FROM {0}".format(gene_table_name))
    tmp_num = cursor.fetchall()[0][0]
    print("number of all the gene id in {0} is {1}".format(gene_table_name, tmp_num))

    format_variance_restrict = "" if not variance_restrict else " AND {}".format(variance_restrict)
    format_gene_restrict = "" if not gene_restrict else " AND {}".format(gene_restrict)
    cmd_str = "SELECT g.gene_id, g.gene_name, v.vcf_id, v.id, v.chr, v.pos, v.ref, v.alt FROM {0} AS v " \
              "INNER JOIN {1} AS g " \
              "ON ((v.chr=g.chr AND v.pos>=g.start_pos AND v.pos<=g.end_pos) OR " \
              "(v.chr=g.chr2 AND v.pos>=g.start_pos2 " \
              "AND v.pos<=g.end_pos2)){2}{3}".format(variance_table_name, gene_table_name,
                                                     format_variance_restrict, format_gene_restrict)

    logging.debug("sql={}".format(cmd_str))
    format_variance_restrict = "" if not variance_restrict else " WHERE {}".format(variance_restrict.replace("v.", ""))
    cursor.execute("SELECT count() FROM {0}{1}".format(variance_table_name, format_variance_restrict))
    tmp_num = cursor.fetchall()[0][0]
    print("number of variance under restrict [{0}] is {1}".format(variance_restrict.replace("v.", ""), tmp_num))
    format_gene_restrict = "" if not gene_restrict else " WHERE {}".format(gene_restrict.replace("g.", ""))
    cursor.execute("SELECT count() FROM {0}{1}".format(gene_table_name, format_gene_restrict))
    tmp_num = cursor.fetchall()[0][0]
    print("number of gene under restrict [{0}] is {1}".format(gene_restrict.replace("v.", ""), tmp_num))
    with open(output, "w") as fp:
        fp.write("##db:{0}\n##variance table:\"{1}\"\n##variance_restrict:\"{2}\"\n"
                 "##gene table:\"{3}\"\n##gene restrict:\"{4}\"\n".format(db_file,
                                                                          variance_table_name,
                                                                          variance_restrict.replace("v.", ""),
                                                                          gene_table_name,
                                                                          gene_restrict))
        fp.write("#gene_id\tgene_name\tvcf_id\tvariance_table_id\tchr\tpos\tref\talt\n")
        logging.debug("begin execute")
        cursor.execute(cmd_str)
        sql_ret_list = cursor.fetchall()
        logging.debug("got sql ret")
        gene_variance_dict = {}
        for element in sql_ret_list:
            if element[0] not in gene_variance_dict:
                gene_variance_dict[element[0]] = [element]
            else:
                gene_variance_dict[element[0]].append(element)
        icounter = 0
        total_len = len(sql_ret_list)
        gene_id_set = set([])
        gene_name_set = set([])
        variance_set = set([])
        # for gene_id in set(filter(lambda x: len(gene_variance_dict[x]) > 1, gene_variance_dict)):
        for gene_id in set(gene_variance_dict):
            element_list = gene_variance_dict[gene_id]
            for element in element_list:
                icounter += 1
                if icounter % 1000 == 0:
                    logging.debug("handled {0} / {1}".format(icounter, total_len))
                fp.write("{}\n".format("\t".join([str(i) for i in element])))
                gene_id_set.add(element[0])
                gene_name_set.add(element[1])
                variance_set.add(element[3])
    print("There are {0} gene {1}, {2} {3} and {4} {5} in collapse result" \
          "".format(len(gene_id_set), "id" if len(gene_id_set) < 2 else "ids",
                    len(gene_name_set), "gene name" if len(gene_name_set) < 2 else "gene names",
                    len(variance_set), "variance" if len(variance_set) < 2 else "variances"))
    cursor.close()
    conn.commit()
    conn.close()
    logging.debug("all done")


def build_data_ram_index(line_data, key_col):
    ret = []
    for i in xrange(len(line_data)):
        if i % 100 == 0:
            ret.append([int(line_data[i][key_col - 1]), i])
    return ret


def parse_fai(fai_in):
    """
    fai file ---> chr str:[offset int, last_pos int]
    :param fai_in:
    :return:
    """
    ret_dict = {}
    icounter = 0
    icounter2 = 0
    with open(fai_in, "r") as fp:
        while True:
            fai_line = fp.readline()
            if not fai_line:
                break
            if not fai_line.strip():
                continue
            icounter2 += 1
            fai_list = fai_line.split("\t")
            chr_name = fai_list[0].strip("chr")
            offset = icounter
            icounter += int(fai_list[1])
            ret_dict[chr_name] = [offset, icounter]
            # if icounter2 == 24:
            #     break
    return ret_dict


def chr_pos2absolute_pos(chrom, pos, chr2offset_dict):
    """
    turn chr pos in to absolute position format
    :param chrom: str
    :param pos: int or str
    :param chr2offset_dict:
    :return: absolute position int
    -1 chrom not in dict. may be not mapped or unusual chrom
    """
    if chrom not in chr2offset_dict:
        # logging.debug("not mapped or unusual chrom.")
        return -1
    return int(pos) + chr2offset_dict[chrom][0]


# def snp_in_gene(gene_chr, gene_start_pos, gene_end_pos, gene_chr2, gene_start_pos2, gene_end_pos2, v_chr, v_pos):
#     if gene_chr is not None:
#         gene_chr = str(gene_chr)
#     if gene_chr2 is not None:
#         gene_chr2 = str(gene_chr2)
#     v_chr = str(v_chr)
#     v_pos = int(v_pos)
#     if gene_start_pos is not None:
#         gene_start_pos = int(gene_start_pos)
#     if gene_end_pos is not None:
#         gene_end_pos = int(gene_end_pos)
#     if gene_start_pos2 is not None:
#         gene_start_pos2 = int(gene_start_pos2)
#     if gene_end_pos2 is not None:
#         gene_end_pos2 = int(gene_end_pos2)
#     is_in_region1 = False
#     if gene_chr == v_chr and gene_start_pos <= v_pos <= gene_end_pos:
#         is_in_region1 = True
#     is_in_region2 = False
#     if gene_chr2 == v_chr and gene_start_pos2 <= v_pos <= gene_end_pos2:
#         is_in_region2 = True
#     return is_in_region1 or is_in_region2


def sort_key_last(elem):
    return int(elem[-1])


def sort_key_last_float(elem):
    return float(elem[-1])


def sort_key7(elem):
    return int(elem[6])


def variance_in_region(variance_data, variance_data_index, region, key_col):
    r_min = region[0]
    r_max = region[1]
    ret = []
    if r_min < variance_data_index[0][0]:
        begin_index = 0
    else:
        for i in xrange(len(variance_data_index) - 1):
            try:
                abp, index = variance_data_index[i]
            except:
                logging.debug("can not unpack variance_data_index[i] = {}".format(variance_data_index[i]))
            abp1, index1 = variance_data_index[i + 1]
            if abp <= r_min <= abp1:
                begin_index = index
                break
        else:
            begin_index = variance_data_index[-1][1]
    for j in xrange(begin_index, len(variance_data), 1):
        if variance_data[j][key_col - 1] < r_min:
            continue
        if variance_data[j][key_col - 1] > r_max:
            break
        ret.append(variance_data[j])
    return ret


def collapse_new(db_file, variance_table_name, gene_table_name, variance_restrict, gene_restrict, output, fai_in):
    """
    For each gene, select the snp within its range from the snp table
    variance_restrict eg: v.is001=1 AND v.bystro=1
    gene_restrict eg: g.gene_name IN ('SOCS4', 'RAB22A')
    :type output: str
    :type db_file: str
    :type gene_restrict: str
    :type variance_restrict: str
    :type variance_table_name: str
    """
    chr2offset_dict = parse_fai(fai_in)
    logging.basicConfig(filename="{}.log".format(sys._getframe().f_code.co_name),
                        level=logging.DEBUG, format=log_format, filemode="w")
    logging.debug("begin")
    logging.debug("db_file=[{}]".format(db_file))
    logging.debug("variance_table_name=[{}]".format(variance_table_name))
    logging.debug("gene_table_name=[{}]".format(gene_table_name))
    logging.debug("variance_restrict=[{}]".format(variance_restrict))
    logging.debug("gene_restrict=[{}]".format(gene_restrict))
    logging.debug("output=[{}]".format(output))
    if variance_table_name not in ['variance', 'synonymous_snp']:
        logging.debug("illegal table name [{0}]".format(variance_table_name))
        return
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()
    cursor.execute("SELECT count() FROM {0}".format(variance_table_name))
    tmp_num = cursor.fetchall()[0][0]
    print("all variance number in table {0} is {1}".format(variance_table_name, tmp_num))
    cursor.execute("SELECT count() FROM {0}".format(gene_table_name))
    tmp_num = cursor.fetchall()[0][0]
    print("number of all the gene id in {0} is {1}".format(gene_table_name, tmp_num))

    format_variance_restrict = "" if not variance_restrict else " WHERE {}".format(variance_restrict)
    format_gene_restrict = "" if not gene_restrict else " WHERE {}".format(gene_restrict)
    cmd_str = "SELECT g.gene_id, g.gene_name, g.chr, g.start_pos, g.end_pos, g.chr2, g.start_pos2, g.end_pos2 " \
              "FROM {0} AS g{1}".format(gene_table_name, format_gene_restrict)
    cursor.execute(cmd_str)
    gene_data = cursor.fetchall()
    cmd_str = "SELECT v.vcf_id, v.id, v.chr, v.pos, v.ref, v.alt FROM {0} AS v{1}" \
              "".format(variance_table_name, format_variance_restrict)
    cursor.execute(cmd_str)
    variance_data = cursor.fetchall()
    for i in xrange(len(variance_data)):
        variance_data[i] = list(variance_data[i])
        variance_data[i].append(chr_pos2absolute_pos(str(variance_data[i][2]), variance_data[i][3], chr2offset_dict))
    variance_data.sort(key=sort_key7)
    variance_data_index = build_data_ram_index(variance_data, 7)

    logging.debug("variance_data[0] = {}".format(variance_data[0]))
    logging.debug("variance_data[-1] = {}".format(variance_data[-1]))
    logging.debug("len(variance_data) = {}".format(len(variance_data)))
    # exit(0)
    with open(output, "w") as fp:
        fp.write("##db:{0}\n##variance table:\"{1}\"\n##variance_restrict:\"{2}\"\n"
                 "##gene table:\"{3}\"\n##gene restrict:\"{4}\"\n".format(db_file,
                                                                          variance_table_name,
                                                                          variance_restrict.replace("v.", ""),
                                                                          gene_table_name,
                                                                          gene_restrict))
        fp.write("#gene_id\tgene_name\tvcf_id\tvariance_table_id\tchr\tpos\tref\talt\n")
        icounter = 0
        gene_data_len = len(gene_data)
        tmp_str = ""
        for gene_id, gene_name, gene_chr, gene_start_pos, gene_end_pos, gene_chr2, gene_start_pos2, gene_end_pos2 in gene_data:
            selected_variance = []
            if gene_chr is not None:
                gene_region = [chr_pos2absolute_pos(str(gene_chr), gene_start_pos, chr2offset_dict),
                               chr_pos2absolute_pos(str(gene_chr), gene_end_pos, chr2offset_dict)]
                selected_variance = variance_in_region(variance_data, variance_data_index, gene_region, 7)
            if gene_chr2 is not None:
                gene_region = [chr_pos2absolute_pos(str(gene_chr2), gene_start_pos2, chr2offset_dict),
                               chr_pos2absolute_pos(str(gene_chr2), gene_end_pos2, chr2offset_dict)]
                selected_variance.extend(variance_in_region(variance_data, variance_data_index, gene_region, 7))
            for curr_variance in selected_variance:
                tmp_str = "{0}{1}\t{2}\t{3}\n".format(tmp_str,
                                                      gene_id,
                                                      gene_name,
                                                      "\t".join([str(i) for i in curr_variance[:6]]))

            icounter += 1
            if icounter % 1000 == 0:
                logging.debug("handled {0} / {1} gene datas".format(icounter, gene_data_len))
        fp.write(tmp_str)
    logging.debug("all done")


def build_contingency_table(db_file, collapse_file, phenotype, output, sample_restrict):
    """
           case     control
        ---------------------
        |         |         |
    alt |    A    |    B    |
        |         |         |
        ---------------------
        |         |         |
    ref |    C    |    D    |
        |         |         |
        ---------------------
    :param collapse_file:
    :return:
    """
    logging.basicConfig(filename="{}.log".format(sys._getframe().f_code.co_name),
                        level=logging.DEBUG, format=log_format, filemode="w")
    logging.debug("begin")
    logging.debug("db_file=[{}]".format(db_file))
    logging.debug("collapse_file=[{}]".format(collapse_file))
    logging.debug("phenotype=[{}]".format(phenotype))
    logging.debug("sample_restrict=[{}]".format(sample_restrict))
    logging.debug("output=[{}]".format(output))

    def get_alt_ref(genotype_list):
        """
        从genotype list中计算出alt alelle num和ref alelle num，并返回
        """
        # print genetype_list
        genotype_list = [int(i) for i in filter(lambda x: x.isdigit(), genotype_list)]
        # print genetype_list
        alt = sum(genotype_list)
        ref = sum([2 - i for i in genotype_list])
        # alt1 = len(filter(lambda x: x > 0, genotype_list))
        # ref1 = len(genotype_list) - alt1
        return [alt, ref]

    def format_genotype_list(genotype_list):
        return map(lambda x: int(x) if x.isdigit() else 0, genotype_list)

    if phenotype not in ["heart6", "ps_andor_pa6", "raa6", "iaab6", "pta6",
                         "tof6", "asdall6", "asdalone6", "vsd6", "vsdalone6",
                         "tofall6", "purevsdalone6", "ps_or_pa_and_vsd6", "intracardiac6", "aorticarch6",
                         "heartnoasd6", "tof_or_pta6", "tof_or_pta_or_iaab6", "CTD"]:
        logging.error("illegal phenotype [{}]".format(phenotype))
        return

    # build head_dict and gene_variance_dict
    logging.debug("loading collapse...")
    head_dict = {}
    gene_variance_dict = {}  # gene_id\tgene_name ---> variance_table_id
    ret_head_str = ''
    with open(collapse_file, "r") as fp:
        while True:
            data_line = fp.readline()
            if not data_line:
                break
            if data_line.startswith("##"):
                ret_head_str += data_line
                data_list = re.split("[:<>]", data_line.strip("#").strip())
                head_dict[data_list[0]] = data_list[1]
                continue
            if data_line.startswith("#") or not data_line.strip():
                continue
            data_list = data_line.strip().split("\t")
            key = "{0}\t{1}".format(data_list[0], data_list[1])
            if key not in gene_variance_dict:
                gene_variance_dict[key] = [data_list[3]]
            else:
                gene_variance_dict[key].append(data_list[3])
    variance_table = head_dict["variance table"]

    logging.debug("begin select control list and case list...")
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()
    format_restrict = "" if not sample_restrict else " AND {}".format(sample_restrict)
    cursor.execute("SELECT s.gen_id FROM sampleChdPhenotype AS s WHERE s.{0}='0'{1}".format(phenotype, format_restrict))
    control_id_list = [i[0] for i in cursor.fetchall()]  # m
    print("control number={}".format(len(control_id_list)))
    cursor.execute("SELECT s.gen_id FROM sampleChdPhenotype AS s WHERE s.{0}='1'{1}".format(phenotype, format_restrict))
    case_id_list = [i[0] for i in cursor.fetchall()]  # n
    print("case number={}".format(len(case_id_list)))
    # print len(control_id_list)
    # print len(case_id_list)
    logging.debug("begin handle genes")
    with open(output, "w") as fp:
        icounter = 0
        fp.write("{0}##phenotype:\"{1}\"\n##control number:{2}\n##case number:{3}\n"
                 "##sample restrict:\"{4}\"\n".format(ret_head_str, phenotype,
                                                      len(control_id_list),
                                                      len(case_id_list),
                                                      sample_restrict))
        fp.write("""##                         case     control
##                      ---------------------
##                      |         |         |
##    alt allele number |    A    |    B    |
##                      |         |         |
##                      ---------------------
##                      |         |         |
##    ref allele number |    C    |    D    |
##                      |         |         |
##                      ---------------------
##
##
##                              case     control
##                           ---------------------
##                           |         |         |
##            people has alt |    A1   |    B1   |
##                           |         |         |
##                           ---------------------
##                           |         |         |
##    people doesn't has alt |    C1   |    D1   |
##                           |         |         |
##                           ---------------------
""")
        fp.write("#gene_id\tgene_name\tA\tB\tC\tD\tp_value\todds_ratio\tA1\tB1\tC1\tD1\tp_value1\todds_ratio1\n")
        for key in gene_variance_dict:
            gene_id, gene_name = key.split("\t")
            icounter += 1
            if icounter % 100 == 0:
                logging.debug("handled {0} / {1} genes".format(icounter, len(gene_variance_dict)))
            sql_control_str = "SELECT "
            for control_id in control_id_list:
                sql_control_str = "{0}sample_{1}, ".format(sql_control_str, control_id)
            sql_control_str = sql_control_str.strip(", ")
            sql_control_str = "{0} FROM {1} WHERE id=".format(sql_control_str, variance_table)
            # print "sql_control={}".format(sql_control_str)

            sql_case_str = "SELECT "
            for case_id in case_id_list:
                sql_case_str = "{0}sample_{1}, ".format(sql_case_str, case_id)
            sql_case_str = sql_case_str.strip(", ")
            sql_case_str = "{0} FROM {1} WHERE id=".format(sql_case_str, variance_table)
            # print "sql_case={}".format(sql_case_str)
            AA = 0
            BB = 0
            CC = 0
            DD = 0
            control_people_list = []
            case_people_list = []
            for id in gene_variance_dict[key]:
                # print "id={}".format(id)
                cmd_str = "{0}{1}".format(sql_control_str, id)
                # print cmd_str
                cursor.execute(cmd_str)
                sql_ret = cursor.fetchall()
                # print sql_ret
                control_genetype_list = [str(i) for i in sql_ret[0]]
                alt_num, ref_num = get_alt_ref(control_genetype_list)
                # control_people_list = control_people_list or format_genotype_list(control_genetype_list)  # ?
                if not control_people_list:
                    control_people_list = format_genotype_list(control_genetype_list)
                else:
                    control_people_list = [max(i) for i in
                                           zip(control_people_list, format_genotype_list(control_genetype_list))]
                BB += alt_num
                DD += ref_num
                cmd_str = "{0}{1}".format(sql_case_str, id)
                cursor.execute(cmd_str)
                # print cmd_str
                sql_ret = cursor.fetchall()
                case_genotype_list = [str(i) for i in sql_ret[0]]
                alt_num, ref_num = get_alt_ref(case_genotype_list)
                # case_people_list = case_people_list or format_genotype_list(case_genetype_list)  # ?
                if not case_people_list:
                    case_people_list = format_genotype_list(case_genotype_list)
                else:
                    case_people_list = [max(i) for i in zip(case_people_list, format_genotype_list(case_genotype_list))]
                AA += alt_num
                CC += ref_num
            B = len(filter(lambda x: x > 0, control_people_list))
            D = len(filter(lambda x: x == 0, control_people_list))
            A = len(filter(lambda x: x > 0, case_people_list))
            C = len(filter(lambda x: x == 0, case_people_list))
            if AA + BB < 2:
                continue
            oddsratio, pvalue = stats.fisher_exact([[AA, BB], [CC, DD]])
            oddsratio1, pvalue1 = stats.fisher_exact([[A, B], [C, D]])
            fp.write(
                "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}"
                "\n".format(gene_id, gene_name, AA, BB, CC, DD,
                            pvalue, oddsratio, A, B, C, D,
                            pvalue1, oddsratio1))
            # break
        logging.debug("handled {} genes in total".format(icounter))
    conn.commit()
    conn.close()
    logging.debug("all done")


def build_contingency_table_new(db_file, phenotype,
                                sample_table_name, sample_restrict,
                                gene_table_name, gene_restrict,
                                variance_table, variance_restrict,
                                fai_in, output,
                                job_id, should_log=True):
    if os.path.exists(variance_restrict):
        with open(variance_restrict, "r") as fp:
            vr = fp.readline().strip()
    else:
        vr = variance_restrict
    build_contingency_table_new_inner(db_file, phenotype, sample_table_name, sample_restrict,
                                      gene_table_name, gene_restrict, variance_table, vr,
                                      fai_in, output, job_id, should_log)


def build_contingency_table_new_inner(db_file, phenotype,
                                      sample_table_name, sample_restrict,
                                      gene_table_name, gene_restrict,
                                      variance_table, variance_restrict,
                                      fai_in, output,
                                      job_id, should_log=True):
    """
           case     control
        ---------------------
        |         |         |
    alt |    A    |    B    |
        |         |         |
        ---------------------
        |         |         |
    ref |    C    |    D    |
        |         |         |
        ---------------------
    Carry out fisher exact test in units of genes
    :return:
    """
    if should_log:
        logging.basicConfig(filename="{0}{1}.log".format(sys._getframe().f_code.co_name, job_id),
                            level=logging.DEBUG, format=log_format, filemode="w")
        logging.debug("begin")
        logging.debug("db_file=[{0}]".format(db_file))
        logging.debug("phenotype=[{0}]".format(phenotype))
        logging.debug("sample_table_name=[{0}]".format(sample_table_name))
        logging.debug("sample_restrict=[{0}]".format(sample_restrict))
        logging.debug("output=[{0}]".format(output))

    if phenotype not in ["heart6", "ps_andor_pa6", "raa6", "iaab6", "pta6",
                         "tof6", "asdall6", "asdalone6", "vsd6", "vsdalone6",
                         "tofall6", "purevsdalone6", "ps_or_pa_and_vsd6", "intracardiac6", "aorticarch6",
                         "heartnoasd6", "tof_or_pta6", "tof_or_pta_or_iaab6", "CTD"]:
        if should_log:
            logging.error("illegal phenotype [{}]".format(phenotype))
        return

    chr2offset_dict = parse_fai(fai_in)
    if should_log:
        logging.debug("begin select control list and case list...")
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()
    format_sample_restrict = "" if not sample_restrict else " AND {}".format(sample_restrict)
    cursor.execute("SELECT s.gen_id FROM {2} AS s WHERE s.{0}='0'{1}".format(phenotype,
                                                                             format_sample_restrict,
                                                                             sample_table_name))
    control_id_list = [i[0] for i in cursor.fetchall()]  # m
    print("control number={}".format(len(control_id_list)))
    cursor.execute("SELECT s.gen_id FROM {2} AS s WHERE s.{0}='1'{1}".format(phenotype,
                                                                             format_sample_restrict,
                                                                             sample_table_name))
    case_id_list = [i[0] for i in cursor.fetchall()]  # n
    print("case number={}".format(len(case_id_list)))
    if should_log:
        logging.debug("begin handle genes")
    with open(output, "w") as fp:
        fp.write("##phenotype:\"{0}\"\n##control number:{1}\n##case number:{2}\n"
                 "##sample restrict:\"{3}\"\n".format(phenotype,
                                                      len(control_id_list),
                                                      len(case_id_list),
                                                      sample_restrict))
        fp.write("""##                         case     control
##                      ---------------------
##                      |         |         |
##    alt allele number |    A    |    B    |
##                      |         |         |
##                      ---------------------
##                      |         |         |
##    ref allele number |    C    |    D    |
##                      |         |         |
##                      ---------------------
##
##
##                              case     control
##                           ---------------------
##                           |         |         |
##          subjects have alt|    A1   |    B1   |
##                           |         |         |
##                           ---------------------
##                           |         |         |
##  subjects do not have alt |    C1   |    D1   |
##                           |         |         |
##                           ---------------------
""")
        fp.write("#gene_id\tgene_name\tA\tB\tC\tD\tp_value\todds_ratio\tA1\tB1\tC1\tD1\tp_value1\todds_ratio1\t"
                 "n_variance_in_gene\tvariance_in_gene_name\tn_variance_control\tvariance_in_control_name\t"
                 "n_variance_case\tvariance_in_case_name\n")
        # prepare the variance data
        cmd_str = "SELECT vcf_id, chr, pos, "
        for control_id in control_id_list:
            cmd_str = "{0}sample_{1}, ".format(cmd_str, control_id)
        for case_id in case_id_list:
            cmd_str = "{0}sample_{1}, ".format(cmd_str, case_id)
        cmd_str = cmd_str.strip(", ")
        format_variance_restrict = "" if not variance_restrict else " WHERE {}".format(variance_restrict)
        cmd_str = "{0} FROM {1} as v{2}".format(cmd_str, variance_table, format_variance_restrict)
        if should_log:
            logging.debug("getting the variance data from db under restrict [{}]".format(variance_restrict))
        cursor.execute(cmd_str)
        variance_list = cursor.fetchall()
        if should_log:
            logging.debug("got {} variance data".format(len(variance_list)))
        for i in xrange(len(variance_list)):
            variance_list[i] = list(variance_list[i])
            variance_list[i].append(chr_pos2absolute_pos(str(variance_list[i][1]),
                                                         variance_list[i][2],
                                                         chr2offset_dict))
        variance_list.sort(key=sort_key_last)
        variance_list_index = build_data_ram_index(variance_list, len(variance_list[0]))

        format_gene_restrict = "" if not gene_restrict else " WHERE {}".format(gene_restrict)
        cmd_str = "SELECT g.gene_id, g.gene_name, g.chr, g.start_pos, g.end_pos, g.chr2, g.start_pos2, g.end_pos2 " \
                  "FROM {0} AS g{1}".format(gene_table_name, format_gene_restrict)
        logging.debug("gene sql = [{}]".format(cmd_str))
        cursor.execute(cmd_str)
        gene_data = cursor.fetchall()

        ret_str = ""
        icounter = 0
        gene_data_len = len(gene_data)
        for gene_id, gene_name, gene_chr, gene_start, gene_end, gene_chr2, gene_start2, gene_end2 in gene_data:
            if should_log:
                if icounter % 100 == 0 and icounter > 0:
                    logging.debug("handled {0} / {1} genes".format(icounter, gene_data_len))
            variance_selected_list = []
            if gene_chr is not None:
                region = [chr_pos2absolute_pos(str(gene_chr), gene_start, chr2offset_dict),
                          chr_pos2absolute_pos(str(gene_chr), gene_end, chr2offset_dict)]
                variance_selected_list = variance_in_region(variance_list, variance_list_index, region,
                                                            len(variance_list[0]))
            if gene_chr2 is not None:
                region2 = [chr_pos2absolute_pos(str(gene_chr2), gene_start2, chr2offset_dict),
                           chr_pos2absolute_pos(str(gene_chr2), gene_end2, chr2offset_dict)]
                variance_selected_list.extend(variance_in_region(variance_list, variance_list_index, region2,
                                                                 len(variance_list[0])))
            all_variance_in_gene_num = len(variance_selected_list)
            all_variance_in_gene_name_list = [i[0] for i in variance_selected_list]
            control_data = [i[3:3 + len(control_id_list)] for i in variance_selected_list]
            case_data = [i[3 + len(control_id_list):3 + len(control_id_list) + len(case_id_list)] for i in
                         variance_selected_list]
            tmp_list = []
            control_variance_in_gene_num = 0
            control_variance_in_gene_name_list = []
            for i in xrange(len(control_data)):
                control_line_data = filter(lambda x: type(x) == int, control_data[i])
                tmp_list.extend(control_line_data)
                if sum(control_line_data) > 0:
                    control_variance_in_gene_num += 1
                    control_variance_in_gene_name_list.append(variance_selected_list[i][0])
            B = sum(tmp_list)
            D = 2 * len(tmp_list) - B
            tmp_list = []
            case_variance_in_gene_num = 0
            case_variance_in_gene_name_list = []
            for i in xrange(len(case_data)):
                case_line_data = filter(lambda x: type(x) == int, case_data[i])
                tmp_list.extend(case_line_data)
                if sum(case_line_data) > 0:
                    case_variance_in_gene_num += 1
                    case_variance_in_gene_name_list.append(variance_selected_list[i][0])
            A = sum(tmp_list)
            C = 2 * len(tmp_list) - A
            # if A + B < 2:
            #     icounter += 1
            #     continue
            control_people_data = [sum(filter(lambda x: type(x) == int, people_line)) for people_line in
                                   zip(*control_data)]
            B1 = len(filter(lambda x: x > 0, control_people_data))
            D1 = len(control_people_data) - B1
            case_people_data = [sum(filter(lambda x: type(x) == int, people_line)) for people_line in zip(*case_data)]
            A1 = len(filter(lambda x: x > 0, case_people_data))
            C1 = len(case_people_data) - A1

            if A1 + B1 <= 2:
                icounter += 1
                continue  # the number of people with variance is less than 2. Go to next gene.

            oddsratio, pvalue = stats.fisher_exact([[A, B], [C, D]])
            oddsratio1, pvalue1 = stats.fisher_exact([[A1, B1], [C1, D1]])
            ret_str = "{0}{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t" \
                      "{11}\t{12}\t{13}\t{14}\t{15}\t{16}\t{17}\t{18}\t{19}\t{20}\n" \
                      "".format(ret_str, gene_id, gene_name, A, B, C, D, pvalue, oddsratio, A1, B1, C1, D1,
                                pvalue1, oddsratio1,
                                all_variance_in_gene_num,
                                ";".join(all_variance_in_gene_name_list) if len(
                                    all_variance_in_gene_name_list) > 0 else ".",
                                control_variance_in_gene_num,
                                ";".join(control_variance_in_gene_name_list) if len(
                                    control_variance_in_gene_name_list) > 0 else ".",
                                case_variance_in_gene_num,
                                ";".join(case_variance_in_gene_name_list) if len(
                                    case_variance_in_gene_name_list) > 0 else ".")
            icounter += 1

        fp.write(ret_str)
    if should_log:
        logging.debug("all done")
    # exit(0)


def filter_collapse_and_vcf_with_gene_list(org_vcf, org_collapse, gene_list_file, vcf_out, collapse_out):
    logging.basicConfig(filename="{}.log".format(sys._getframe().f_code.co_name),
                        level=logging.DEBUG, format=log_format, filemode="w")
    logging.debug("begin")

    # load gene set
    gene_set = set([])
    with open(gene_list_file, "r") as fp:
        while True:
            data_line = fp.readline()
            if not data_line:
                break
            if data_line.startswith("#") or not data_line.strip():
                continue
            gene_set.add(data_line.strip().split("\t")[0])

    # load vcf index
    vcf_dict = {}
    with open(org_vcf + "i", "r") as fp:
        data = fp.readlines()
    data = filter(lambda x: len(x.strip()) > 0, data)
    for i in data:
        line_list = i.strip().split("|")
        vcf_dict[line_list[0]] = int(line_list[1])
    # write output collapse and collect vcf_key_list
    vcf_key_list = []
    with open(org_collapse, "r") as fp_co_in, open(collapse_out, "r") as fp_co_out:
        while True:
            data_line = fp_co_in.readline()
            if not data_line:
                break
            if data_line.startswith("#"):
                fp_co_out.write(data_line)
                continue
            if not data_line.strip():
                continue
            data_list = data_line.strip().split("\t")
            if data_list[0] in gene_set or data_list[1] in gene_set:
                fp_co_out.write(data_line)
                vcf_key_list.append("{0}\t{1}\t{2}\t{3}".format(data_list[4], data_list[5], data_list[6], data_list[7]))
    vcf_key_list = list(set(vcf_key_list))

    # write vcf
    with open(vcf_out, "w") as fp_vcf_out, open(org_vcf, "r") as fp_vcf_in:
        while True:
            vcf_line = fp_vcf_in.readline()
            if not vcf_line.startswith("#"):
                break
            fp_vcf_out.write(vcf_line)
        for vcf_key in vcf_key_list:
            if vcf_key not in vcf_dict:
                logging.error("can not find variance {}".format(vcf_key))
                continue
            fp_vcf_in.seek(vcf_dict[vcf_key])
            fp_vcf_out.write(fp_vcf_in.readline())


def get_variance_num_in_sample(db_file, variance_table, variance_restrict,
                               sample_table, sample_restrict, output,
                               need_log=True):
    if need_log:
        logging.basicConfig(filename="{}.log".format(sys._getframe().f_code.co_name),
                            level=logging.DEBUG, format=log_format, filemode="w")
        logging.debug("db_file=[{}]".format(db_file))
        logging.debug("variance_table=[{}]".format(variance_table))
        logging.debug("variance_restrict=[{}]".format(variance_restrict))
        logging.debug("sample_table=[{}]".format(sample_table))
        logging.debug("sample_restrict=[{}]".format(sample_restrict))
        logging.debug("output=[{}]".format(output))
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()
    if need_log:
        logging.debug("selecting sample id...")
    pheno_list = ["heart6", "ps_andor_pa6", "raa6", "iaab6", "pta6",
                  "tof6", "asdall6", "asdalone6", "vsd6", "vsdalone6",
                  "tofall6", "purevsdalone6", "ps_or_pa_and_vsd6", "intracardiac6", "aorticarch6",
                  "heartnoasd6", "tof_or_pta6", "tof_or_pta_or_iaab6"]
    cmd_str = "SELECT gen_id, {0} FROM {1} AS s{2}" \
              "".format(", ".join(pheno_list),
                        sample_table,
                        "" if not sample_restrict else " WHERE {}".format(sample_restrict))
    if need_log:
        logging.debug("sql=[{}]".format(cmd_str))
    cursor.execute(cmd_str)
    sample_data = cursor.fetchall()
    sample_id_list = [int(i[0]) for i in sample_data]

    if need_log:
        logging.debug("selecting genotype...")
    cmd_str = "SELECT {0} FROM {1} AS v{2}" \
              "".format(", ".join(["sample_{}".format(i) for i in sample_id_list]),
                        variance_table,
                        "" if not variance_restrict else " WHERE {}".format(variance_restrict))
    if need_log:
        logging.debug("sql=[{}]".format(cmd_str))
    cursor.execute(cmd_str)
    gene_type_num_list = cursor.fetchall()

    if need_log:
        logging.debug("writing result...")
    with open(output, "w") as fp:
        fp.write("#sample_id\talt_allele_num\t{0}\n".format("\t".join(pheno_list)))
        for i in xrange(len(sample_id_list)):
            fp.write("{0}\t{1}\t{2}\n".format(sample_id_list[i],
                                              sum([int(k) for k in
                                                   filter(lambda x: x.isdigit(),
                                                          [str(j[i]) for j in gene_type_num_list])]),
                                              "\t".join([str(i) for i in sample_data[i][1:]])))
    conn.commit()
    conn.close()
    if need_log:
        logging.debug("all done")
    # cursor


def get_variance_function_number(file_all, file01, file001, col, output):
    logging.basicConfig(filename="{}.log".format(sys._getframe().f_code.co_name),
                        level=logging.DEBUG, format=log_format, filemode="w")
    col = int(col)

    def handle_file(file_name, col, function_num_dict, set_all):
        counter = 0
        with open(file_name, "r") as fp:
            while True:
                line = fp.readline()
                counter += 1
                if counter % 10000 == 0:
                    logging.debug("file {0} handled {1} lines".format(os.path.basename(file_name), counter))
                if not line:
                    break
                name_list = line.split("\t")[col - 1].split(";")
                for name in name_list:
                    set_all.add(name)
                    if name not in function_num_dict:
                        function_num_dict[name] = 1
                    else:
                        function_num_dict[name] += 1

    dict_all = {}
    dict_01 = {}
    dict_001 = {}
    set_all = set([])
    logging.debug("handling {}".format(file_all))
    handle_file(file_all, col, dict_all, set_all)

    logging.debug("handling {}".format(file01))
    handle_file(file01, col, dict_01, set_all)

    logging.debug("handling {}".format(file001))
    handle_file(file001, col, dict_001, set_all)
    with open(output, "w") as fp:
        fp.write("\tall\t01\t001\n")
        for name in set_all:
            fp.write("{0}\t{1}\t{2}\t{3}\n".format(name, dict_all[name] if name in dict_all else "0",
                                                   dict_01[name] if name in dict_01 else "0",
                                                   dict_001[name] if name in dict_001 else "0"))


def test_t(file_in, value_col, group_cols):
    """

    :type value_col: str
    :type group_cols: str
    """
    assert type(value_col) == str
    assert type(group_cols) == str
    assert type(file_in) == str
    with open(file_in, "r") as fp:
        data = fp.readlines()
    data_head_list = data[0].strip().split("\t")

    data = [i.strip().split("\t") for i in data[1:]]
    group_col_list = [int(i) for i in group_cols.strip().split(",")]
    value_col = int(value_col)
    print(data_head_list[value_col - 1])
    value_list = [i[value_col - 1] for i in data]
    for group_col in group_col_list:
        group_name = data_head_list[group_col - 1]
        group = [str(i[group_col - 1]) for i in data]
        rvs1 = [int(value_list[i]) for i in filter(lambda x: group[x] == "0", xrange(len(value_list))) if value_list[i]]
        rvs2 = [int(value_list[i]) for i in filter(lambda x: group[x] == "1", xrange(len(value_list))) if value_list[i]]
        statistic, pvalue = stats.levene(rvs1, rvs2)
        print("{}:".format(group_name))
        print("All input samples are {}from populations with equal variances".format("" if pvalue > 0.05 else "not "))
        print(stats.ttest_ind(rvs1, rvs2, equal_var=(False if pvalue < 0.05 else True)))
        print("""Group 0:
        mean={0}
        median={1}
        extremelybad ptp={2}
        variance var={3}
        standard deviation std={4}""".format(mean(rvs1), median(rvs1), ptp(rvs1), 
        var(rvs1), std(rvs1)))
        print("""Group 1:
        average value mean={0}
        median={1}
        extremely bad ptp={2}
        variance var={3}
        standard deviation std={4}""".format(mean(rvs2), median(rvs2), ptp(rvs2), 
        var(rvs2), std(rvs2)))
        # data = array([rvs1, rvs2])
        # print"""
        # Calculate the covariance of two sets of numbers
        # bias=1 indicates that the result needs to be divided by N, otherwise only 
        the molecular part is calculated
        # The result is a matrix, and the data in row i and column j represents the 
        covariance between group i and group j. The diagonal is the variance
        # covariance matrix cov={0}
        #
        # Calculate the correlation coefficient of two sets of numbers
        # The returned result is a matrix, and the data in the i-th row and the j-th 
        column represents the correlation coefficient between the i-th group and the j-th 
        group. Diagonal is 1
        # correlation coefficient corrcoef={1}""".format(str(cov(data, bias=1)), 
        str(corrcoef(data)))

        print("")


def handle_gt(str_in, target):
    gt_str, gq = str_in.split(":")
    gt_list = gt_str.split("/")
    for i in xrange(len(gt_list)):
        if gt_list[i] == "0" or gt_list[i] == ".":
            continue
        if gt_list[i] == str(target):
            gt_list[i] = "1"
            continue
        else:
            gt_list[i] = "0"
    ret = "/".join(gt_list)
    if ret == "1/0":
        ret = "0/1"
    return "{0}:{1}".format(ret, gq)


def split_vcf_line_with_sample(vcf_line):
    """

    :type vcf_line: str
    """
    ret = []
    assert type(vcf_line) == str
    vcf_list = vcf_line.split("\t")
    alt_list = vcf_list[4].split(",")
    if len(alt_list) == 1:
        return [vcf_line]
    alt_dict = dict(zip(alt_list, xrange(1, len(alt_list) + 1, 1)))
    for alt in alt_dict:
        tmp_list = copy.copy(vcf_list)
        tmp_list[4] = alt
        # print "alt={}".format(alt)
        for index in xrange(9, len(tmp_list), 1):
            # print index
            tmp_list[index] = handle_gt(tmp_list[index], alt_dict[alt])
            # print tmp_list[index]
        ret.append("\t".join(tmp_list))
    return ret


def split_vcf_with_sample():
    while True:
        vcf_line = sys.stdin.readline()
        if not vcf_line:
            break
        if vcf_line.startswith("#"):
            sys.stdout.write(vcf_line)
            continue
        vcf_line = vcf_line.strip()
        if len(vcf_line) == 0:
            continue
        splited_line_list = split_vcf_line_with_sample(vcf_line)
        sys.stdout.write("\n".join(splited_line_list))
        sys.stdout.write("\n")


def cross_genehancer_snp(hancer_file, snp_file, output):
    def snp_is_in_hancer(snp_line, hancer_line):
        snp_chr, snp_pos = snp_line.split("\t")[:2]
        hancer_chr, hancer_start, hancer_end = hancer_line.split("\t")[:3]
        if snp_chr != hancer_chr:
            return False
        isnp_pos = int(snp_pos)
        ihancer_start = int(hancer_start)
        ihancer_end = int(hancer_end)
        if isnp_pos < ihancer_start or isnp_pos > ihancer_end:
            return False
        return True

    logging.basicConfig(filename="{}.log".format(sys._getframe().f_code.co_name),
                        level=logging.DEBUG, format=log_format, filemode="w")
    with open(hancer_file, "r") as fp_han:
        hancer_lines = fp_han.readlines()
    hancer_lines = filter(lambda x: not x.startswith("#"), hancer_lines)
    icounter = 0
    with open(snp_file, "r") as fp_snp, open(output, "w") as fp_out:
        while True:
            snp_line = fp_snp.readline()
            icounter += 1
            # if icounter % 100 == 0:
            logging.debug("handled {} snp lines".format(icounter))
            if not snp_line:
                break
            if snp_line.startswith("#"):
                continue
            if not snp_line.strip():
                continue
            for hancer_line in hancer_lines:
                if snp_is_in_hancer(snp_line, hancer_line):
                    fp_out.write("{0}\t{1}\n".format(hancer_line.strip(), snp_line.strip()))


def cross_genehancer_snp_writer(q2w, output_file):
    icounter = 0
    with open(output_file, "w") as fp_out:
        while True:
            data = q2w.get(block=True)
            if data != "STOP":
                icounter += 1
            if icounter % 1000 == 0:
                logging.debug("handled {} snp results".format(icounter))
            if data == "STOP":
                logging.debug("handled {} snp results".format(icounter))
                break
            if len(data) > 0:
                fp_out.write(data)
    logging.debug("all done")


def get_ip():
    try:
        s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
        s.connect(('8.8.8.8', 80))
        ip = s.getsockname()[0]
    finally:
        s.close()
    return ip


def cross_genehancer_snp_multiple(num, hancer_file, snp_file, output_file):
    logging.basicConfig(filename="{}.log".format(sys._getframe().f_code.co_name),
                        level=logging.DEBUG, format=log_format, filemode="w")
    if float(sys.version[:3]) != 2.7:
        print("CRITICAL: Python version must be 2.7!\n")
        sys.exit(1)
    curr_dir = os.getcwd()
    server_addr = get_ip()
    print("master ip = {}".format(server_addr))
    with open(curr_dir + "/tmp.sh", "w") as fp:
        fp.write("""#!/bin/bash
#$ -j y
#$ -S /bin/bash
#$ -cwd
#$ -N worker
#$ -pe smp 1
#$ -l h_vmem=16G
#$ -m bes
#$ -q all.q
export LC_ALL=C
export MALLOC_ARENA_MAX=4
export PATH=$PATH:$HOME/wyj/.code
wgsa.py cross_genehancer_snp_worker {0} {1} $1""".format(server_addr, snp_file))
        fp.flush()
        fp.close()

    while True:
        if os.access(curr_dir + "/tmp.sh", os.R_OK):
            print("{} read OK".format(curr_dir + "/tmp.sh"))
            break
        time.sleep(1)
    for i in xrange(int(num)):
        Popen(["qsub {0}/tmp.sh {1}".format(curr_dir, i + 1)], shell=True)

    cross_genehancer_snp_master(hancer_file, snp_file, output_file, int(num))
    Popen(["rm {0}/tmp.sh".format(curr_dir)], shell=True)


# def cross_genehancer_snp_worker(q2w, curr_data, snp_file, hancer_lines):
def cross_genehancer_snp_worker(server_addr, snp_file, worker_id):
    def snp_is_in_hancer_(snp_chr, snp_pos, hancer_chr, hancer_start, hancer_end):
        # snp_chr, snp_pos = snp_line.split("\t")[:2]
        if snp_chr != hancer_chr:
            return False
        isnp_pos = int(snp_pos)
        ihancer_start = int(hancer_start)
        ihancer_end = int(hancer_end)
        if isnp_pos < ihancer_start or isnp_pos > ihancer_end:
            return False
        return True

    logging.basicConfig(filename="worker{}.log".format(worker_id),
                        level=logging.DEBUG, format=log_format, filemode="w")
    BaseManager.register("get_task_queue")
    BaseManager.register("get_result_queue")
    BaseManager.register('get_msgc_queue')
    print("Connect to server {}...".format(server_addr))
    logging.debug("worker{0} Connect to server {1}...".format(worker_id, server_addr))
    manager = BaseManager(address=(server_addr, 5000), authkey=b"wgsawgsa")
    # wait for server
    icounter = 0
    while 1:
        try:
            manager.connect()
        except:
            time.sleep(1)
            icounter += 1
            print("waiting for server {}s".format(icounter))
            continue
        break
    print("server {} Connected".format(server_addr))
    # get Queue object
    icounter = 0
    while True:
        try:
            task_queue = manager.get_task_queue()
            result_queue = manager.get_result_queue()
            msgc_queue = manager.get_msgc_queue()
            icounter += 1
        except:
            print("try to get queue object {} times".format(icounter))
            continue
        break
    msgc_queue.put("{}".format(get_ip()), block=False)
    icounter = 0
    while 1:
        try:
            msg_task = task_queue.get(block=True)
            if msg_task == "STOP":
                task_queue.task_done()
                print("I'm done.")
                break
            # handle the task here
            offset_list = msg_task[0]
            hancer_lines = msg_task[1]
            hancer_list = [i.split("\t")[:3] for i in hancer_lines]
            offset_list_len = len(offset_list)
            with open(snp_file, "r") as fp:
                for offset in offset_list:
                    fp.seek(offset)
                    snp_line = fp.readline()
                    snp_chr, snp_pos = snp_line.split("\t")[:2]
                    ret = ""
                    for i in xrange(len(hancer_lines)):
                        if snp_is_in_hancer_(snp_chr, snp_pos, hancer_list[i][0], hancer_list[i][1], hancer_list[i][2]):
                            ret = "{0}{1}\t{2}\n".format(ret, hancer_lines[i].strip(), snp_line.strip())
                    result_queue.put(ret, block=False)
                    icounter += 1
                    if icounter % 1000 == 0:
                        print("handled {} snps".format(icounter))
                        logging.debug("handled {0} / {1} snps".format(icounter, offset_list_len))
            task_queue.task_done()
        except Queue.Empty:
            print("task queue is empty.")
    logging.debug("handled {0} / {1} snps".format(icounter, offset_list_len))
    logging.debug("all done")
    return


if sys.version_info[0] == 2:
    task_queue = Queue.Queue()
    result_queue = Queue.Queue()
    msgc_queue = Queue.Queue()
elif sys.version_info[0] == 3:
    task_queue = queue.Queue()
    result_queue = queue.Queue()
    msgc_queue = queue.Queue()


# msgw_queue = Queue.Queue()

def cross_genehancer_snp_master(hancer_file, snp_file, output_file, num):
    def return_task_queue():
        global task_queue
        return task_queue

    def return_result_queue():
        global result_queue
        return result_queue

    def return_msgc_queue():
        global msgc_queue
        return msgc_queue

    multiprocessing.freeze_support()
    BaseManager.register('get_task_queue', callable=return_task_queue)
    BaseManager.register('get_result_queue', callable=return_result_queue)
    BaseManager.register('get_msgc_queue', callable=return_msgc_queue)
    manager = BaseManager(address=('', 5000), authkey=b"wgsawgsa")
    # start Queue:
    manager.start()
    # get Queue objects:
    task_queue = manager.get_task_queue()
    result_queue = manager.get_result_queue()
    msgc_queue = manager.get_msgc_queue()
    print("starting writer...")
    logging.debug("starting writer...")
    proc_write_result = multiprocessing.Process(target=cross_genehancer_snp_writer, args=(result_queue, output_file))
    proc_write_result.daemon = True
    proc_write_result.start()
    print("loading hancer file...")
    logging.debug("loading hancer file...")
    # Assign task
    with open(hancer_file, "r") as fp_han:
        hancer_lines = fp_han.readlines()
    hancer_lines = filter(lambda x: not x.startswith("#"), hancer_lines)
    print("loading snp index file...")
    logging.debug("loading snp index file...")
    with open(snp_file + "i", "r") as fp:
        snp_index = [int(i.split("|")[1]) for i in fp.readlines()]
    print("assigning task...")
    logging.debug("assigning task...")
    step = len(snp_index) / (num - 1)
    with open(snp_file, "r") as fp:
        for i in xrange(num):
            if i < num - 1:
                curr_index_list = snp_index[i * step:(i * step + step)]
            else:
                curr_index_list = snp_index[i * step:]
            task_queue.put([curr_index_list, hancer_lines])
        # for i in snp_index:
        #    fp.seek(i)
        #    task_queue.put([fp.readline(), hancer_lines])
    print("waiting for tasks...")
    logging.debug("waiting for tasks...")
    # wait for tasks
    task_queue.join()
    result_queue.put("STOP", block=True)
    proc_write_result.join()
    worker_ip_list = []
    while 1:
        if msgc_queue.empty():
            break
        else:
            worker_ip_list.append(msgc_queue.get())
    logging.debug("got {} worker in total".format(len(worker_ip_list)))
    for i in xrange(len(worker_ip_list)):
        task_queue.put("STOP", block=True)
    task_queue.join()
    time.sleep(10)
    manager.shutdown()
    logging.debug("sent {} done message".format(len(worker_ip_list)))
    print('master exit.')


#
# def cross_genehancer_snp_multiprocess(hancer_file, snp_file, output_file, p_num=150, w_num=200):
#     def cross_genehancer_snp_writer(q2w, output_file):
#         icounter = 0
#         with open(output_file, "w") as fp_out:
#             while True:
#                 data = q2w.get(block=True)
#                 if data != "STOP":
#                     icounter += 1
#                 if icounter % 100 == 0:
#                     logging.debug("handled {} snps".format(icounter))
#                 if data == "STOP":
#                     logging.debug("handled {} snps".format(icounter))
#                     break
#                 fp_out.write(data)
#         logging.debug("all done")
#
#     def cross_genehancer_snp_worker(q2w, curr_data, snp_file, hancer_lines):
#         with open(snp_file, "r") as fp:
#             for offset in curr_data:
#                 fp.seek(offset)
#                 snp_line = fp.readline()
#                 for hancer_line in hancer_lines:
#                     if snp_is_in_hancer(snp_line, hancer_line):
#                         q2w.put("{0}\t{1}\n".format(hancer_line.strip(), snp_line.strip()), block=False)
#     q2w = multiprocessing.Manager().Queue()
#     process_pool = multiprocessing.Pool(processes=p_num)
#     proc_write_result = multiprocessing.Process(target=cross_genehancer_snp_writer, args=(q2w, output_file))
#     proc_write_result.daemon = True
#     proc_write_result.start()
#     # setup para
#     with open(hancer_file, "r") as fp_han:
#         hancer_lines = fp_han.readlines()
#     hancer_lines = filter(lambda x: not x.startswith("#"), hancer_lines)
#     with open(snp_file + "i", "r") as fp:
#         snp_index = [int(i.split("|")[1]) for i in fp.readlines()]
#     step = len(snp_index) / w_num - 1
#     for i in xrange(w_num):
#         if i < w_num - 1:
#             curr_data = snp_index[(step * i):(step * i + step)]
#         else:
#             curr_data = snp_index[step * i:]
#         process_pool.apply_async(cross_genehancer_snp_worker, (q2w, curr_data, snp_file, hancer_lines))
#         # setup data
#     process_pool.close()
#     process_pool.join()
#     q2w.put("STOP", block=True)
def check_log(path):
    files = os.listdir(path)  # Get all file names under the folder
    icounter = 0
    icounter2 = 0
    icounter3 = 0
    for file in files:  # traverse folders
        if not os.path.isdir(file) and "worker" in file and ".log" in file:
            log_str = ""
            log_str = Popen(["tail -n 2 {}".format(file)], stdout=PIPE, shell=True).communicate()[0]
            icounter3 += 1
            # print log_str
            if "all done" in log_str:
                num_str_list = re.findall("\[cross_genehancer_snp_worker\]handled (\d+) / ", log_str)
                icounter2 += 1
                if len(num_str_list) > 0:
                    icounter += int(num_str_list[-1])
    print("handled {} snps in total".format(icounter))
    print("{0} / {1} processes quit".format(icounter2, icounter3))


def test_ks(file_in, value_col, group_cols, out_path, p):
    """
    Computes the Kolmogorov-Smirnov statistic on case an control.
    :type out_path: str
    :type value_col: str
    :type group_cols: str
    """
    from scipy.stats import ks_2samp
    assert type(value_col) == str
    assert type(group_cols) == str
    assert type(file_in) == str
    # print "test_ks begin"
    out_path = out_path + ("" if out_path.endswith("/") else "/")
    # print "out_path={}".format(out_path)
    # print "file_in={}".format(file_in)
    # print "value_col={}".format(value_col)
    # print "group_cols={}".format(group_cols)
    # print "p={}".format(p)
    file_name = os.path.basename(file_in)
    p = float(p)
    with open(file_in, "r") as fp:
        data = fp.readlines()
    data_head_list = data[0].strip().split("\t")
    data = [i.strip().split("\t") for i in data[1:]]
    group_col_list = [int(i) for i in group_cols.strip().split(",")]
    value_col = int(value_col)
    value_list = [int(i[value_col - 1]) for i in data]
    for group_col in group_col_list:
        group_name = data_head_list[group_col - 1]
        group = [str(i[group_col - 1]) for i in data]
        rvs1 = [value_list[i] for i in filter(lambda x: group[x] == "0", xrange(len(value_list)))]
        rvs2 = [value_list[i] for i in filter(lambda x: group[x] == "1", xrange(len(value_list)))]
        ret = ks_2samp(rvs1, rvs2)
        if ret.pvalue < p:
            # print "{}:".format(group_name)
            # print "pvalue=[{}]".format(ret.pvalue)
            # print """Group 0:
            # average value mean={0}
            # median={1}
            # extremely bad ptp={2}
            # Variance var={3}
            # standard deviation std={4}""".format(mean(rvs1), median(rvs1), ptp(rvs1), 
            var(rvs1), std(rvs1))
            # print """Group 1:
            # average value mean={0}
            # median={1}
            # extremely bad ptp={2}
            # variance var={3}
            # standard deviation std={4}
            # """.format(mean(rvs2), median(rvs2), ptp(rvs2), var(rvs2), std(rvs2))

            with open(out_path + file_name + ".{0}.ks.{1}.selected".format(group_name, ret.pvalue), "w") as fp:
                fp.write("alt_allele_num\t{}\n".format(group_name))
                fp.write("".join(["{0}\t{1}\n".format(value_list[i], group[i]) for i in xrange(len(value_list))]))


def test_MWW(file_in, value_col, group_cols, out_path, p):
    """
    Compute the Mann-Whitney rank test on samples x and y.
    In statistics, the Mann–Whitney U test (also called the
    Mann–Whitney–Wilcoxon (MWW), Wilcoxon rank-sum test, or
    Wilcoxon–Mann–Whitney test) is a nonparametric test of the
    null hypothesis that it is equally likely that a randomly
    selected value from one population will be less than or
    greater than a randomly selected value from a second population.
    :type value_col: str
    :type group_cols: str
    """
    # print "test_MWW begin"
    from scipy.stats import mannwhitneyu
    assert type(value_col) == str
    assert type(group_cols) == str
    assert type(file_in) == str
    out_path = out_path + ("" if out_path.endswith("/") else "/")
    # print "out_path={}".format(out_path)
    # print "file_in={}".format(file_in)
    # print "value_col={}".format(value_col)
    # print "group_cols={}".format(group_cols)
    # print "p={}".format(p)
    file_name = os.path.basename(file_in)
    p = float(p)
    with open(file_in, "r") as fp:
        data = fp.readlines()
    data_head_list = data[0].strip().split("\t")
    data = [i.strip().split("\t") for i in data[1:]]
    group_col_list = [int(i) for i in group_cols.strip().split(",")]
    value_col = int(value_col)
    value_list = [int(i[value_col - 1]) for i in data]
    for group_col in group_col_list:
        group_name = data_head_list[group_col - 1]
        group = [str(i[group_col - 1]) for i in data]
        rvs1 = [value_list[i] for i in filter(lambda x: group[x] == "0", xrange(len(value_list)))]
        rvs2 = [value_list[i] for i in filter(lambda x: group[x] == "1", xrange(len(value_list)))]
        ret = mannwhitneyu(rvs1, rvs2)
        if ret.pvalue < p:
            # print "{}:".format(group_name)
            # print "pvalue=[{}]".format(ret.pvalue)
            # print """Group 0:
            # average value mean={0}
            # median={1}
            # extremely bad ptp={2}
            # variance var={3}
            # standard deviation std={4}""".format(mean(rvs1), median(rvs1), ptp(rvs1), 
            var(rvs1), std(rvs1))
            # print """Group 1:
            # average value mean={0}
            # median={1}
            # extremely bad ptp={2}
            # variance var={3}
            # standard deviation std={4}
            # """.format(mean(rvs2), median(rvs2), ptp(rvs2), var(rvs2), std(rvs2))
            with open(out_path + file_name + ".{0}.MWW.{1}.selected".format(group_name, ret.pvalue), "w") as fp:
                fp.write("alt_allele_num\t{}\n".format(group_name))
                fp.write("".join(["{0}\t{1}\n".format(value_list[i], group[i]) for i in xrange(len(value_list))]))


def get_geneset_variance_num_in_sample(db_file, variance_table, variance_restrict, sample_table,
                                       sample_restrict, output, gene_id_list):
    logging.basicConfig(filename="{}.log".format(sys._getframe().f_code.co_name),
                        level=logging.DEBUG, format=log_format, filemode="w")
    variance_restrict = variance_restrict.strip()
    logging.debug("db_file=[{}]".format(db_file))
    logging.debug("variance_table=[{}]".format(variance_table))
    logging.debug("variance_restrict=[{}]".format(variance_restrict))
    logging.debug("sample_table=[{}]".format(sample_table))
    logging.debug("sample_restrict=[{}]".format(sample_restrict))
    logging.debug("output=[{}]".format(output))

    gene_id_list = gene_id_list.strip().strip(",").split(",")
    logging.debug("gene_id_list={}".format(gene_id_list))
    assert len(gene_id_list) > 0
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()
    logging.debug("selecting sample id...")
    pheno_list = ["heart6", "ps_andor_pa6", "raa6", "iaab6", "pta6",
                  "tof6", "asdall6", "asdalone6", "vsd6", "vsdalone6",
                  "tofall6", "purevsdalone6", "ps_or_pa_and_vsd6", "intracardiac6", "aorticarch6",
                  "heartnoasd6", "tof_or_pta6", "tof_or_pta_or_iaab6"]
    cmd_str = "SELECT gen_id, {0} FROM {1} AS s{2}" \
              "".format(", ".join(pheno_list),
                        sample_table,
                        "" if not sample_restrict else " WHERE {}".format(sample_restrict))
    logging.debug("sql=[{}]".format(cmd_str))
    cursor.execute(cmd_str)
    sample_data = cursor.fetchall()
    sample_id_list = [int(i[0]) for i in sample_data]

    logging.debug("selecting genotype...")
    cmd_str = ""
    for gene_id in gene_id_list:
        if not variance_restrict:
            variance_restrict_all = " WHERE v.pos BETWEEN (SELECT start_pos FROM gene_table WHERE gene_id=\"{}\") " \
                                    "AND (SELECT end_pos FROM gene_table WHERE gene_id=\"{}\") " \
                                    "AND v.chr=(SELECT chr FROM gene_table WHERE gene_id=\"{}\")" \
                                    "".format(gene_id)
        else:
            variance_restrict_all = " WHERE {0} AND " \
                                    "v.pos BETWEEN (SELECT start_pos FROM gene_table WHERE gene_id=\"{1}\") " \
                                    "AND (SELECT end_pos FROM gene_table WHERE gene_id=\"{1}\") " \
                                    "AND v.chr=(SELECT chr FROM gene_table WHERE gene_id=\"{1}\")" \
                                    "".format(variance_restrict, gene_id)
        if not cmd_str:
            cmd_str = "SELECT {0} FROM {1} AS v{2}".format(", ".join(["sample_{}".format(i) for i in sample_id_list]),
                                                           variance_table,
                                                           variance_restrict_all)
        else:
            cmd_str += "\nUNION\nSELECT {0} FROM {1} AS v{2}" \
                       "".format(", ".join(["sample_{}".format(i) for i in sample_id_list]),
                                 variance_table,
                                 variance_restrict_all)
    logging.debug("sql=[{}]".format(cmd_str))
    cursor.execute(cmd_str)
    gene_type_num_list = cursor.fetchall()

    logging.debug("writing result...")
    with open(output, "w") as fp:
        fp.write("#sample_id\talt_allele_num\t{0}\n".format("\t".join(pheno_list)))
        for i in xrange(len(sample_id_list)):
            fp.write("{0}\t{1}\t{2}\n".format(sample_id_list[i],
                                              sum([int(k) for k in
                                                   filter(lambda x: x.isdigit(),
                                                          [str(j[i]) for j in gene_type_num_list])]),
                                              "\t".join([str(i) for i in sample_data[i][1:]])))
    conn.commit()
    conn.close()
    logging.debug("all done")


def multiple_geneset_variance_num_in_sample(db_file, variance_table, variance_restrict, sample_table, sample_restrict,
                                            output, gene_set_file):
    def one_geneset_variance_num_in_sample(db_file, variance_table, variance_restrict, gene_id_list, sample_id_list):
        # logging.debug("gene_id_list={}".format(gene_id_list))
        conn = sqlite3.connect(db_file)
        cursor = conn.cursor()
        # logging.debug("selecting genotype...")
        variance_id_set = set([])
        cmd_str = ""
        for gene_id in gene_id_list:
            logging.debug("gene_id = {}".format(gene_id))
            if not variance_restrict:
                variance_restrict_all = " WHERE v.pos BETWEEN (SELECT start_pos FROM gene_table WHERE gene_id=\"{0}\" OR " \
                                        "gene_name=\"{0}\") " \
                                        "AND (SELECT end_pos FROM gene_table WHERE gene_id=\"{0}\" OR " \
                                        "gene_name=\"{0}\") " \
                                        "AND v.chr=(SELECT chr FROM gene_table WHERE gene_id=\"{0}\" OR " \
                                        "gene_name=\"{0}\")".format(gene_id)
            else:
                variance_restrict_all = " WHERE {0} AND " \
                                        "v.pos BETWEEN (SELECT start_pos FROM gene_table WHERE gene_id=\"{1}\" OR " \
                                        "gene_name=\"{1}\") " \
                                        "AND (SELECT end_pos FROM gene_table WHERE gene_id=\"{1}\" OR " \
                                        "gene_name=\"{1}\") " \
                                        "AND v.chr=(SELECT chr FROM gene_table WHERE gene_id=\"{1}\" OR " \
                                        "gene_name=\"{1}\")".format(variance_restrict, gene_id)
            # if not cmd_str:
            #     cmd_str = "SELECT {0} FROM {1} AS v{2}".format(
            #         ", ".join(["sample_{}".format(i) for i in sample_id_list]),
            #         variance_table,
            #         variance_restrict_all)
            # else:
            #     cmd_str += "\nUNION\nSELECT {0} FROM {1} AS v{2}" \
            #                "".format(", ".join(["sample_{}".format(i) for i in sample_id_list]),
            #                          variance_table,
            #                          variance_restrict_all)
            cmd_str = "SELECT id from {0} as v{1}".format(variance_table, variance_restrict_all)
            cursor.execute(cmd_str)
            id_data = [i[0] for i in cursor.fetchall()]
            for i in id_data:
                variance_id_set.add(str(i))
        cmd_str = "SELECT {0} from {1} WHERE id in {2}" \
                  "".format(", ".join(["sample_{}".format(i) for i in sample_id_list]),
                            variance_table,
                            "(" + ",".join(variance_id_set) + ")")
        logging.debug("sql=[{}]".format(cmd_str))
        cursor.execute(cmd_str)
        gene_type_num_list = cursor.fetchall()

        # logging.debug("writing result...")
        # with open(output, "w") as fp:
        #     fp.write("#sample_id\talt_allele_num\t{0}\n".format("\t".join(pheno_list)))
        #     for i in xrange(len(sample_id_list)):
        #         fp.write("{0}\t{1}\t{2}\n".format(sample_id_list[i],
        #                                           sum([int(k) for k in
        #                                                filter(lambda x: x.isdigit(),
        #                                                       [str(j[i]) for j in gene_type_num_list])]),
        #                                           "\t".join([str(i) for i in sample_data[i][1:]])))
        conn.commit()
        conn.close()
        return [sum([int(k) for k in filter(lambda x: x.isdigit(), [str(j[i]) for j in gene_type_num_list])]) for i in
                xrange(len(sample_id_list))]

    if os.path.exists(variance_restrict):
        with open(variance_restrict, "r") as fp:
            variance_restrict = fp.readline().strip()

    logging.basicConfig(filename="{}.log".format(sys._getframe().f_code.co_name),
                        level=logging.DEBUG, format=log_format, filemode="w")
    variance_restrict = variance_restrict.strip()
    logging.debug("db_file=[{}]".format(db_file))
    logging.debug("variance_table=[{}]".format(variance_table))
    logging.debug("variance_restrict=[{}]".format(variance_restrict))
    logging.debug("sample_table=[{}]".format(sample_table))
    logging.debug("sample_restrict=[{}]".format(sample_restrict))
    logging.debug("output=[{}]".format(output))

    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()
    logging.debug("selecting sample id...")
    pheno_list = ["heart6", "ps_andor_pa6", "raa6", "iaab6", "pta6",
                  "tof6", "asdall6", "asdalone6", "vsd6", "vsdalone6",
                  "tofall6", "purevsdalone6", "ps_or_pa_and_vsd6", "intracardiac6", "aorticarch6",
                  "heartnoasd6", "tof_or_pta6", "tof_or_pta_or_iaab6", "CTD"]
    cmd_str = "SELECT gen_id, {0} FROM {1} AS s{2}" \
              "".format(", ".join(pheno_list),
                        sample_table,
                        "" if not sample_restrict else " WHERE {}".format(sample_restrict))
    logging.debug("sql=[{}]".format(cmd_str))
    cursor.execute(cmd_str)
    sample_data = cursor.fetchall()
    for i in xrange(len(sample_data)):
        sample_data[i] = list(sample_data[i])
    # print(sample_data[0])
    sample_id_list = [int(i[0]) for i in sample_data]
    conn.commit()
    conn.close()
    head_str = "#sample_id\t" + "\t".join(pheno_list)
    with open(gene_set_file, "r") as fp:
        while True:
            gene_set_line = fp.readline()
            if not gene_set_line:
                break
            if len(gene_set_line.strip()) == 0 or gene_set_line.startswith("#"):
                continue
            gene_set_list = gene_set_line.strip().split("\t")
            gene_set_name = gene_set_list[0]
            print("handling {} ...".format(gene_set_name))
            head_str += "\t" + gene_set_name
            gene_id_list = gene_set_list[1:]
            assert len(gene_id_list) > 0
            sample_variance_in_gene_num = one_geneset_variance_num_in_sample(db_file, variance_table,
                                                                             variance_restrict, gene_id_list,
                                                                             sample_id_list)
            # logging.debug("len(sample_variance_in_gene_num)={}".format())
            assert len(sample_variance_in_gene_num) == len(sample_id_list)
            for i in xrange(len(sample_data)):
                sample_data[i].append(str(sample_variance_in_gene_num[i]))
    head_str += "\n"
    with open(output, "w") as fp:
        fp.write(head_str)
        for i in sample_data:
            fp.write("sample_" + "\t".join([str(k) for k in i]) + "\n")


def test_GLM(file_in, value_col, group_cols, family, out_path, p):
    """
    :type value_col: str
    :type group_cols: str
    """
    import statsmodels.api as sm
    assert type(value_col) == str
    assert type(group_cols) == str
    assert type(file_in) == str
    # print "test_GLM begin"
    out_path = out_path + ("" if out_path.endswith("/") else "/")
    # print "out_path={}".format(out_path)
    # print "file_in={}".format(file_in)
    # print "value_col={}".format(value_col)
    # print "group_cols={}".format(group_cols)
    # print "family={}".format(family)
    # print "p={}".format(p)
    file_name = os.path.basename(file_in)
    p = float(p)
    if family not in ["Binomial", "Gamma", "Gaussian", "InverseGaussian", "NegativeBinomial", "Poisson", "Tweedie"]:
        print("Family should be in [Binomial, Gamma, Gaussian, InverseGaussian, NegativeBinomial, Poisson, Tweedie]." \
              "But it is {} now. \nExit.".format(family))
        exit(0)
    with open(file_in, "r") as fp:
        data = fp.readlines()
    data_head_list = data[0].strip().split("\t")
    data = [i.strip().split("\t") for i in data[1:]]
    group_col_list = [int(i) for i in group_cols.strip().split(",")]
    value_col = int(value_col)
    value_list = [int(i[value_col - 1]) for i in data]
    for group_col in group_col_list:
        group_name = data_head_list[group_col - 1]
        group = [str(i[group_col - 1]) for i in data]
        # rvs1 = [value_list[i] for i in filter(lambda x: group[x] == "0", xrange(len(value_list)))]
        # rvs2 = [value_list[i] for i in filter(lambda x: group[x] == "1", xrange(len(value_list)))]
        exog = [int(i) for i in filter(lambda x: x == "0" or x == "1", group)]
        endog = [value_list[i] for i in filter(lambda x: group[x] == "0" or group[x] == "1", xrange(len(group)))]

        if family == "Binomial":
            my_family = sm.families.Binomial()
        elif family == "Gamma":
            my_family = sm.families.Gamma()
        elif family == "Gaussian":
            my_family = sm.families.Gaussian()
        elif family == "InverseGaussian":
            my_family = sm.families.InverseGaussian()
        elif family == "NegativeBinomial":
            my_family = sm.families.NegativeBinomial()
        elif family == "Poisson":
            my_family = sm.families.Poisson()
        elif family == "Tweedie":
            my_family = sm.families.Tweedie()
        else:
            my_family = sm.families.Gamma()
        model = sm.GLM(array(endog), sm.add_constant(array(exog)), family=my_family)
        result = model.fit()
        # print result.summary()
        if result.pvalues[1] < p:
            # print "{}:".format(group_name)
            # print("pvalue={}".format(result.pvalues[1]))
            # print """Group exog:
            # average value mean={0}
            # median={1}
            # extremely bad ptp={2}
            # variance var={3}
            # standard deviation std={4}""".format(mean(exog), median(exog), ptp(exog), 
            var(exog), std(exog))
            # print """Group endog:
            # average value mean={0}
            # median={1}
            # extremely bad ptp={2}
            # variance var={3}
            # standard deviation std={4}
            # """.format(mean(endog), median(endog), ptp(endog), var(endog), std(endog))

            with open(out_path + file_name + ".{0}.GLM_{1}.{2}.selected".format(group_name, family, result.pvalues[1]),
                      "w") as fp:
                fp.write("alt_allele_num\t{}\n".format(group_name))
                fp.write("".join(["{0}\t{1}\n".format(value_list[i], group[i]) for i in xrange(len(value_list))]))


def check_pvalue_in_table(p_threshold_in, table_file):
    try:
        p_threshold_in = float(p_threshold_in)
        with open(table_file, "r") as fp:
            data_lines = filter(lambda x: (not x.startswith("##")) and len(x) > 3, fp.readlines())
        head_list = data_lines[0].strip().split("\t")
        index_p = head_list.index("p_value")
        index_p1 = head_list.index("p_value1")
        data_list = [i.strip().split("\t") for i in data_lines[1:]]
        p_is_good = False
        p1_is_good = False
        for curr_list in data_list:
            if not p_is_good:
                if float(curr_list[index_p]) < p_threshold_in:
                    p_is_good = True
            if not p1_is_good:
                if float(curr_list[index_p1]) < p_threshold_in:
                    p1_is_good = True
            if p_is_good and p1_is_good:
                break
        table_file_name = os.path.basename(table_file)
        if p_is_good or p1_is_good:
            print(table_file_name + ("\tp" if p_is_good else "\t") + ("\tp1" if p1_is_good else "\t"))
    except:
        pass


def bar_hist_plot(number_file, phenotype):
    import matplotlib.pyplot as plt
    # import pylab
    import scipy
    with open(number_file, "r") as fp:
        data_lines = fp.readlines()
    head_list = data_lines[0].strip().split("\t")
    data_list = [i.strip().split("\t") for i in data_lines[1:]]
    index_value = head_list.index("alt_allele_num")
    if phenotype not in head_list:
        print("please check your phenotype [{}] again.".format(phenotype))
        exit(0)
    index_pheno = head_list.index(phenotype)
    value_list = [int(i[index_value]) for i in data_list]
    group_list = [i[index_pheno] for i in data_list]
    vco_list = [vv for vv, gv in zip(value_list, group_list) if gv == "0"]
    vca_list = [vv for vv, gv in zip(value_list, group_list) if gv == "1"]

    # pylab.font={'family':'fantasy'}
    fig = plt.figure()

    plt.subplot(121)
    x = [0.3, 0.7]
    y = [np.mean(vco_list), np.mean(vca_list)]  # control  case
    # error list
    std_err = [np.std(vco_list), np.std(vca_list)]
    error_params = dict(elinewidth=0.5, ecolor='black', capsize=3, capthick=0.5)  # Set error 
    flag parameters
    # Draw a histogram, set error markers and histogram labels
    plt.bar(x, y,
            color=['#5ecacc', '#eb9390'],
            width=0.25,
            yerr=std_err,
            error_kw=error_params,
            tick_label=['control', 'case'])
    plt.xlim((0, 3))
    # plt.text(0, 1, "mean={0:.2f}\nstd={1:.2f}".format(np.mean(vco_list), np.std(vco_list)), family="Arial")
    # plt.text(0.95, 1, "mean={0:.2f}\nstd={1:.2f}".format(np.mean(vca_list), np.std(vca_list)), family="Arial")

    plt.subplot(122)
    n_co, bins_co, patches_co = plt.hist(vco_list, 50, density=1, facecolor='#5ecacc', alpha=1)
    n_ca, bins_ca, patches_ca = plt.hist(vca_list, 50, density=1, facecolor='#eb9390', alpha=0.6)
    plt.plot(bins_co, scipy.stats.norm.pdf(bins_co, np.mean(vco_list), np.std(vco_list)), 'b')
    plt.plot(bins_ca, scipy.stats.norm.pdf(bins_ca, np.mean(vca_list), np.std(vca_list)), 'r')

    # plt.subplot(223)
    # bplot = plt.boxplot([vco_list, vca_list], patch_artist=True)
    # colors = ['red', 'blue']
    # for patch, color in zip(bplot['boxes'], colors):
    #     patch.set_facecolor(color)  # Fill different colors for different boxplots
    # plt.xticks([1, 2], ['control', 'case'])

    # plt.show()
    file_name = os.path.basename(number_file)
    plt.savefig(os.path.splitext(file_name)[0] + "_{}.pdf".format(phenotype))


def analyze_fisher_test_variance(database, path):
    pp = Popen(["cd {}".format(path)], shell=True)
    pp.wait()

    annotator_list_base = ["annovar = 1", "bystro = 1", "dmis = 1", "dsplicing = 1", "spidex = 1", "spliceAI = 1",
                           "vep = 1"]
    annotator_list = []
    annotator_list.extend(list(itertools.combinations(annotator_list_base, 1)))
    annotator_list.extend(list(itertools.combinations(annotator_list_base, 2)))
    annotator_list.extend(list(itertools.combinations(annotator_list_base, 3)))
    annotator_list.extend(list(itertools.combinations(annotator_list_base, 4)))
    annotator_list.extend(list(itertools.combinations(annotator_list_base, 5)))
    annotator_list.extend(list(itertools.combinations(annotator_list_base, 6)))
    annotator_list.extend(list(itertools.combinations(annotator_list_base, 7)))
    annotator_list = [" AND ({})".format(" OR ".join(i)) for i in annotator_list]
    annotator_list.append("")
    icounter = 1
    icounter2 = 0
    script_str = ''
    for phenotype in ["heart6", "tof6", "aorticarch6", "tof_or_pta_or_iaab6"]:
        for annotator in annotator_list:
            for freq in ["bystro_sampleMaf <= 0.01", "bystro_sampleMaf <= 0.001",
                         "bystro_sampleMaf <= 0.005"]:
                for ph in ["", " AND bystro_phastCons >= 0.2", " AND bystro_phastCons >= 0.3",
                           " AND bystro_phastCons >= 0.4",
                           " AND bystro_phastCons >= 0.5", " AND bystro_phastCons >= 0.6",
                           " AND bystro_phastCons >= 0.7",
                           " AND bystro_phastCons >= 0.8", " AND bystro_phyloP >= -1", " AND bystro_phyloP >= 0",
                           " AND bystro_phyloP >= 1", " AND bystro_phyloP >= 2", " AND bystro_phyloP >= 3",
                           " AND bystro_phyloP >= 4"]:
                    for cadd in ["", " AND bystro_cadd >= 5", " AND bystro_cadd >= 10", " AND bystro_cadd >= 15",
                                 " AND bystro_cadd >= 20", " AND bystro_cadd >= 25", " AND bystro_cadd >= 30"]:
                        variance_restrict = "{0}{1}{2}{3}".format(freq, annotator, ph, cadd)
                        base_name = re.sub("\)$", "", re.sub(" >= ", "", re.sub("bystro_", "", re.sub(" = ", "",
                                                                                                      re.sub(" OR ",
                                                                                                             "_or_",
                                                                                                             re.sub(
                                                                                                                 " AND ",
                                                                                                                 "_and_",
                                                                                                                 re.sub(
                                                                                                                     "\) AND ",
                                                                                                                     "_and_",
                                                                                                                     re.sub(
                                                                                                                         " AND \(",
                                                                                                                         "_and_",
                                                                                                                         re.sub(
                                                                                                                             "bystro_sampleMaf <= 0.",
                                                                                                                             "freq",
                                                                                                                             variance_restrict)))))))))
                        # print variance_restrict
                        # print base_name
                        if script_str == "":
                            script_str = """#!/bin/bash
#$ -j y
#$ -S /bin/bash
#$ -cwd
#$ -N fisher
#$ -pe smp 1
#$ -l h_vmem=8G
#$ -m bes
# -q all.q
export LC_ALL=C
export MALLOC_ARENA_MAX=4
module load python/2.7.15/gcc.4.4.7
export PATH=$PATH:~/wyj/.code/
time=`date`
echo "==START $time =="

fisher.test.sh -d {0} -v variance -p {1} --variance_restrict "{2}" --collapse_out {3}.collapse --vcf_out {3}.vcf --contingency_out {3}.table
time=`date`
echo == $time ==
""".format(database, phenotype, variance_restrict, base_name)
                        else:
                            script_str += "fisher.test.sh -d {0} -v variance -p {1} --variance_restrict \"{2}\" " \
                                          "--collapse_out {3}.collapse --vcf_out {3}.vcf --contingency_out {3}.table" \
                                          "\ntime=`date`\necho == $time ==\n".format(database, phenotype,
                                                                                     variance_restrict, base_name)

                        if icounter == 302:
                            with open("tmp.sh", "w") as fp:
                                fp.write(script_str + "\ndate; echo \"==END==\"")
                            while True:
                                if os.access("tmp.sh", os.R_OK):
                                    break
                                time.sleep(1)
                            pp = Popen(["qsub tmp.sh"], shell=True)
                            pp.wait()
                            script_str = ""
                            icounter = 0
                            icounter2 += 1
                        icounter += 1
    if len(script_str) > 0:
        with open("tmp.sh", "w") as fp:
            fp.write(script_str)
        while True:
            if os.access("tmp.sh", os.R_OK):
                break
            time.sleep(1)
        pp = Popen(["qsub tmp.sh"], shell=True)
        pp.wait()
    pp = Popen(["rm tmp.sh"], shell=True)
    pp.wait()
    print("All the jobs has been submitted. Job number = {}".format(icounter2 + 1))


def whole_genome_variance_test(db_file, variance_table, sample_table, sample_restrict, outputpath, p_threshold):
    logging.basicConfig(filename="whole_genome_variance_test.log", level=logging.INFO, format=log_format, filemode="w")
    pp = Popen(["cd {}".format(outputpath)], shell=True)
    pp.wait()
    logging.info("Begin")
    annotator_list_base = ["annovar = 1", "bystro = 1", "dmis = 1", "dsplicing = 1", "spidex = 1", "spliceAI = 1",
                           "vep = 1"]
    annotator_list = []
    annotator_list.extend(list(itertools.combinations(annotator_list_base, 1)))
    annotator_list.extend(list(itertools.combinations(annotator_list_base, 2)))
    annotator_list.extend(list(itertools.combinations(annotator_list_base, 3)))
    annotator_list.extend(list(itertools.combinations(annotator_list_base, 4)))
    annotator_list.extend(list(itertools.combinations(annotator_list_base, 5)))
    annotator_list.extend(list(itertools.combinations(annotator_list_base, 6)))
    annotator_list.extend(list(itertools.combinations(annotator_list_base, 7)))
    annotator_list = [" AND ({})".format(" OR ".join(i)) for i in annotator_list]
    annotator_list.append("")
    phenotype_list = ["heart6", "tof6", "aorticarch6", "tof_or_pta_or_iaab6"]
    group_str = ""
    icounter = 0
    print(len(annotator_list))
    for annotator in annotator_list:
        for freq in ["bystro_sampleMaf <= 0.01", "bystro_sampleMaf <= 0.001",
                     "bystro_sampleMaf <= 0.005"]:
            for ph in ["", " AND bystro_phastCons >= 0.2", " AND bystro_phastCons >= 0.3",
                       " AND bystro_phastCons >= 0.4",
                       " AND bystro_phastCons >= 0.5", " AND bystro_phastCons >= 0.6",
                       " AND bystro_phastCons >= 0.7",
                       " AND bystro_phastCons >= 0.8", " AND bystro_phyloP >= -1", " AND bystro_phyloP >= 0",
                       " AND bystro_phyloP >= 1", " AND bystro_phyloP >= 2", " AND bystro_phyloP >= 3",
                       " AND bystro_phyloP >= 4"]:
                for cadd in ["", " AND bystro_cadd >= 5", " AND bystro_cadd >= 10", " AND bystro_cadd >= 15",
                             " AND bystro_cadd >= 20", " AND bystro_cadd >= 25", " AND bystro_cadd >= 30"]:
                    variance_restrict = "{0}{1}{2}{3}".format(freq, annotator, ph, cadd)
                    base_name = re.sub("bystro_sampleMaf <= 0.", "freq", variance_restrict)
                    base_name = re.sub(" AND \(", "_and_", base_name)
                    base_name = re.sub("\) AND ", "_and_", base_name)
                    base_name = re.sub(" AND ", "_and_", base_name)
                    base_name = re.sub(" OR ", "_or_", base_name)
                    base_name = re.sub(" = ", "", base_name)
                    base_name = re.sub("bystro_", "", base_name)
                    base_name = re.sub(" >= ", "", base_name)
                    base_name = re.sub("\)$", "", base_name)
                    variance_number_file = base_name + ".varn"
                    logging.info("building varn")
                    get_variance_num_in_sample(db_file, variance_table,
                                               variance_restrict, sample_table,
                                               sample_restrict, variance_number_file,
                                               False)
                    if group_str == "":
                        with open(variance_number_file, "r") as fp:
                            head_list = fp.readline().strip().split("\t")
                        for phenotype in phenotype_list:
                            group_str += "{},".format(str(head_list.index(phenotype) + 1))
                        group_str = group_str.strip(",")
                        print(group_str)
                    logging.info("ks test begin")
                    test_ks(variance_number_file, "2", group_str, outputpath, p_threshold)
                    logging.info("MWW test begin")
                    test_MWW(variance_number_file, "2", group_str, outputpath, p_threshold)
                    logging.info("GLM test begin")
                    test_GLM(variance_number_file, "2", group_str, "Poisson", outputpath, p_threshold)
                    icounter += 1
                    logging.info("handled {} cases".format(icounter))
                    # exit(0)
    logging.info("handled {} cases in total".format(icounter))
    logging.info("all done")


def analyze_fisher_test_variance_new(database, variance_table, sample_restrict, fai_in, path):
    pp = Popen(["cd {}".format(path)], shell=True)
    pp.wait()
    if variance_table == "variance":
        annotator_list_base = ["annovar = 1", "bystro = 1", "dmis = 1", "dsplicing = 1", "spidex = 1", "spliceAI = 1",
                               "vep = 1"]
    elif variance_table == "synonymous_snp":
        annotator_list_base = ["annovar = 1", "bystro = 1", "vep = 1"]
    else:
        print("variance_table could only be variance or synonymous_snp")
        return
    annotator_list = []
    annotator_list.extend(list(itertools.combinations(annotator_list_base, 1)))
    annotator_list.extend(list(itertools.combinations(annotator_list_base, 2)))
    annotator_list.extend(list(itertools.combinations(annotator_list_base, 3)))
    annotator_list.extend(list(itertools.combinations(annotator_list_base, 4)))
    annotator_list.extend(list(itertools.combinations(annotator_list_base, 5)))
    annotator_list.extend(list(itertools.combinations(annotator_list_base, 6)))
    annotator_list.extend(list(itertools.combinations(annotator_list_base, 7)))
    annotator_list = [" AND ({})".format(" OR ".join(i)) for i in annotator_list]
    annotator_list.append("")
    icounter = 1
    icounter2 = 0
    script_str = ''
    for phenotype in ["heart6", "tof6", "aorticarch6", "tof_or_pta_or_iaab6"]:
        for annotator in annotator_list:
            for freq in ["bystro_sampleMaf <= 0.01", "bystro_sampleMaf <= 0.001",
                         "bystro_sampleMaf <= 0.005"]:
                for ph in ["", " AND bystro_phastCons >= 0.2", " AND bystro_phastCons >= 0.3",
                           " AND bystro_phastCons >= 0.4",
                           " AND bystro_phastCons >= 0.5", " AND bystro_phastCons >= 0.6",
                           " AND bystro_phastCons >= 0.7",
                           " AND bystro_phastCons >= 0.8", " AND bystro_phyloP >= -1", " AND bystro_phyloP >= 0",
                           " AND bystro_phyloP >= 1", " AND bystro_phyloP >= 2", " AND bystro_phyloP >= 3",
                           " AND bystro_phyloP >= 4"]:
                    for cadd in ["", " AND bystro_cadd >= 5", " AND bystro_cadd >= 10", " AND bystro_cadd >= 15",
                                 " AND bystro_cadd >= 20", " AND bystro_cadd >= 25", " AND bystro_cadd >= 30"]:
                        variance_restrict = "{0}{1}{2}{3}".format(freq, annotator, ph, cadd)
                        output_name = re.sub("bystro_sampleMaf <= 0.", "freq", variance_restrict)
                        output_name = re.sub(" AND \(", "_and_", output_name)
                        output_name = re.sub("\) AND ", "_and_", output_name)
                        output_name = re.sub(" AND ", "_and_", output_name)
                        output_name = re.sub(" OR ", "_or_", output_name)
                        output_name = re.sub(" = ", "", output_name)
                        output_name = re.sub("bystro_", "", output_name)
                        output_name = re.sub(" >= ", "", output_name)
                        output_name = re.sub("\)$", "", output_name)
                        output_name = "{0}_{1}.table".format(phenotype, output_name)
                        # print output_name
                        if script_str == "":
                            script_str = """#!/bin/bash
#$ -j y
#$ -S /bin/bash
#$ -cwd
#$ -N ct{6}
#$ -pe smp 1
#$ -l h_vmem=8G
#$ -m bes
# -q all.q
export LC_ALL=C
export MALLOC_ARENA_MAX=4
module load python/2.7.15/gcc.4.4.7
module load sqlite3/3.8.11/gcc.4.4.7
time=`date`
echo "==START $time =="
~/miniconda2/bin/python ~/wyj/.code/wgsa.py build_contingency_table_new {0} {1} sampleChdPhenotype '{2}' gene_table '' {7} '{3}' {4} {5} {6}
echo {5} is done
time=`date`
echo == $time ==
""".format(database, phenotype, sample_restrict, variance_restrict, fai_in, output_name, icounter2, variance_table)
                        else:
                            script_str += "~/miniconda2/bin/python ~/wyj/.code/wgsa.py build_contingency_table_new {0} {1} " \
                                          "sampleChdPhenotype '{2}' gene_table '' {6} '{3}' {4} {5} {7}" \
                                          "\necho {5} is done" \
                                          "\ntime=`date`\necho == $time ==\n" \
                                          "".format(database,
                                                    phenotype,
                                                    sample_restrict,
                                                    variance_restrict,
                                                    fai_in,
                                                    output_name,
                                                    variance_table,
                                                    icounter2)

                        if icounter == 100:
                            with open("qsub{}.sh".format(icounter2), "w") as fp:
                                fp.write(script_str + "\ndate; echo \"==END==\"")
                            while True:
                                if os.access("qsub{}.sh".format(icounter2), os.R_OK):
                                    break
                                time.sleep(1)
                            pp = Popen(["qsub qsub{}.sh".format(icounter2)], shell=True)
                            pp.wait()
                            script_str = ""
                            icounter = 0
                            icounter2 += 1
                        icounter += 1
    if len(script_str) > 0:
        with open("qsub{}.sh".format(icounter2), "w") as fp:
            fp.write(script_str + "\ndate; echo \"==END==\"")
        while True:
            if os.access("qsub{}.sh".format(icounter2), os.R_OK):
                break
            time.sleep(1)
        pp = Popen(["qsub qsub{}.sh".format(icounter2)], shell=True)
        pp.wait()

    print(
        "All the jobs has been submitted. Job number = {}\nIf any job has problem, kill it. And qsub the corresponding sh file".format(
            icounter2 + 1))


def rebuild_variance_restrict_from_filename(file_name):
    file_name = file_name.split("\t")[0]

    phenotype = re.findall("aorticarch6|tof6|heart6|tof_or_pta_or_iaab6|ps_andor_pa6|raa6|iaab6|pta6|"
                           "asdall6|asdalone6|vsd6|vsdalone6|tofall6|purevsdalone6|ps_or_pa_and_vsd6|"
                           "intracardiac6|heartnoasd6|tof_or_pta6", file_name)
    if not phenotype:
        return None
    if len(phenotype) != 1:
        return None
    phenotype = phenotype[0]
    file_name = re.sub(phenotype + "_", "", file_name)
    # logging.debug("phenotype = [{0}], data_line = [{1}]".format(phenotype, data_line))

    freq = re.findall("freq(\d+)", file_name)
    if not freq:
        freq = ""
    else:
        file_name = re.sub("freq" + freq[0], "", file_name)
        freq = "bystro_sampleMaf <= 0." + freq[0]
    # logging.debug("freq = [{0}] data_line = [{1}]".format(freq, data_line))

    cadd = re.findall("cadd(\d+).table$", file_name)
    if not cadd:
        cadd = ""
        file_name = re.sub(".table", "", file_name)
    else:
        cadd = " AND bystro_cadd >= " + cadd[0]
        file_name = re.sub("cadd(\d+).table$", "", file_name)
    # logging.debug("cadd = [{0}] data_line = [{1}]".format(cadd, data_line))

    ph = re.findall("phyloP-1|phastCons0.\d|phyloP\d", file_name)
    file_name = re.sub("phyloP-1|phastCons0.\d|phyloP\d", "", file_name)
    if not ph:
        ph = ""
    elif "phastCons0." in ph[0]:
        tmp_list = ph[0].split(".")
        ph = " AND bystro_phastCons >= 0." + tmp_list[1]
    elif "phyloP" in ph[0]:
        ph = " AND bystro_phyloP >= " + re.sub("phyloP", "", ph[0])

    # logging.debug("ph = [{0}], data_line = [{1}]".format(ph, data_line))
    file_name = re.sub("_and_", "", file_name)
    tmp_list = file_name.split("_or_")
    logging.debug("tmp_list={0}".format(tmp_list))
    anno_list = [i[:-1] + " = 1" for i in file_name.split("_or_") if len(i.strip()) > 0]
    # logging.debug("anno_list = {}".format(anno_list))
    annotator = " AND ({})".format(" OR ".join(anno_list)) if len(anno_list) > 0 else ""
    variance_restrict = "{0}{1}{2}{3}".format(freq, annotator, ph, cadd)
    return [phenotype, variance_restrict]


def rerun_selected_combination(db_file, variance_table_name, selected_o_file, path, org_vcf, fai_in, sample_restrict,
                               with_ccrs):
    with_ccrs = int(with_ccrs)
    logging.basicConfig(filename="rerun_selected_combination.log",
                        level=logging.DEBUG, format=log_format, filemode="w")
    logging.debug("begin")
    pp = Popen(["mkdir -p {0} && cd {0}".format(path)], shell=True)
    pp.wait()
    if with_ccrs == 0:
        ccrs_list = [""]
    else:
        ccrs_list = ["",
                     " AND (domain_limbr <= -3.654145493 OR domain_limbr is NULL)",
                     " AND (domain_limbr <= -1.751155961 OR domain_limbr is NULL)",
                     " AND (domain_limbr <= -0.699941137 OR domain_limbr is NULL)",
                     " AND (domain_limbr <= 0.091314981 OR domain_limbr is NULL)",
                     " AND (domain_limbr <= 0.753317939 OR domain_limbr is NULL)",
                     " AND (domain_limbr <= 1.338304674 OR domain_limbr is NULL)",
                     " AND (domain_limbr <= 1.881692228 OR domain_limbr is NULL)",
                     " AND (exone_limbr <= -1.538383601 OR exone_limbr is NULL)",
                     " AND (exone_limbr <= -0.572107415 OR exone_limbr is NULL)",
                     " AND (exone_limbr <= -0.023396951 OR exone_limbr is NULL)",
                     " AND (exone_limbr <= 0.390431738 OR exone_limbr is NULL)",
                     " AND (exone_limbr <= 0.750902136 OR exone_limbr is NULL)",
                     " AND (exone_limbr <= 1.084410106 OR exone_limbr is NULL)",
                     " AND (exone_limbr <= 1.412196324 OR exone_limbr is NULL)",
                     " AND (ccrs >= 95 OR ccrs is NULL)",
                     " AND (ccrs >= 90 OR ccrs is NULL)",
                     " AND (ccrs >= 85 OR ccrs is NULL)",
                     " AND (ccrs >= 80 OR ccrs is NULL)",
                     " AND (ccrs >= 75 OR ccrs is NULL)",
                     " AND (ccrs >= 70 OR ccrs is NULL)",
                     " AND (ccrs >= 65 OR ccrs is NULL)"]
    icounter = 0
    logging.debug("begin main loop")
    with open(selected_o_file, "r") as fp:
        while True:
            data_line = fp.readline()
            if not data_line:
                break
            data_line = data_line.strip().split("\t")[0]
            [phenotype, variance_restrict] = rebuild_variance_restrict_from_filename(data_line)
            for ccrs in ccrs_list:
                variance_restrict_new = variance_restrict + ccrs
                if with_ccrs == 0:
                    output_vcf_name = os.path.join(path, re.sub(".table$", ".vcf", data_line))
                    output_table_name = os.path.join(path, data_line)
                else:
                    output_name = re.sub("bystro_sampleMaf <= 0.", "freq", variance_restrict_new)
                    output_name = re.sub(" AND \(", "_and_", output_name)
                    output_name = re.sub("\) AND ", "_and_", output_name)
                    output_name = re.sub(" AND ", "_and_", output_name)
                    output_name = re.sub(" OR ", "_or_", output_name)
                    output_name = re.sub(" = ", "", output_name)
                    output_name = re.sub("bystro_", "", output_name)
                    output_name = re.sub(" >= ", "", output_name)
                    output_name = re.sub("\)$", "", output_name)
                    output_name = "{0}_{1}".format(phenotype, output_name)
                    output_name = re.sub("_\(", "_", output_name)
                    output_name = re.sub("\)_", "_", output_name)
                    output_name = re.sub(" is ", "", output_name)
                    output_vcf_name = os.path.join(path, output_name + ".vcf")
                    output_table_name = os.path.join(path, output_name + ".table")
                    logging.debug("output_name=[{0}]".format(output_name))
                logging.debug("output = [{0}]".format(output_vcf_name))
                logging.debug("variance_restrict=[{0}]".format(variance_restrict_new))
                logging.debug("phenotype=[{0}]".format(phenotype))
                try:
                    export_vcf_from_db(db_file, variance_table_name, variance_restrict_new, org_vcf, output_vcf_name,
                                       fai_in)
                    build_contingency_table_new(db_file, phenotype, "sampleChdPhenotype", sample_restrict, "gene_table",
                                                "", variance_table_name, variance_restrict_new, fai_in,
                                                output_table_name,
                                                0, False)
                except:
                    logging.debug(traceback.format_exc())
                icounter += 1
                logging.debug("got {0} table".format(icounter))
    print("there should be {0} vcf and {0} table files".format(icounter))
    logging.debug("all done")


def cut_org_exported_vcf_samples(org_vcf, small_vcf):
    def sort_key(elem):
        return elem[0]

    with open(small_vcf, "r") as fp:
        small_head_set = set([i.split("_")[0] for i in fp.readline().strip().split("\t")])
        fp.seek(0, 0)
        small_vcf_data = [[j.split(":")[0] for j in i.strip().split("\t")] for i in fp.readlines()]
        small_vcf_data_T = list(zip(*small_vcf_data))
        small_vcf_data_T_1_9 = small_vcf_data_T[:9]
        small_vcf_data_T_9_ = small_vcf_data_T[9:]
        small_vcf_data_T_9_.sort(key=sort_key)
    with open(org_vcf, "r") as fp:
        org_head_set = set([i.split("_")[0] for i in fp.readline().strip().split("\t")])
        fp.seek(0, 0)
        org_vcf_data = [[j.split(":")[0] for j in i.strip().split("\t")] for i in fp.readlines()]
        for i in xrange(len(org_vcf_data[0])):
            org_vcf_data[0][i] = org_vcf_data[0][i].split("_")[0]
        org_vcf_data_T = list(zip(*org_vcf_data))
        org_vcf_data_T_1_9 = org_vcf_data_T[:9]
        org_vcf_data_T_9_ = org_vcf_data_T[9:]
        org_vcf_data_T_9_.sort(key=sort_key)
    interset_head_set = small_head_set & org_head_set

    filter_small = [i[0].split("_")[0] in interset_head_set for i in small_vcf_data_T_9_]
    small_vcf_data_T_9_ = list(compress(small_vcf_data_T_9_, filter_small))

    filter_org = [i[0].split("_")[0] in interset_head_set for i in org_vcf_data_T_9_]
    org_vcf_data_T_9_ = list(compress(org_vcf_data_T_9_, filter_org))

    small_vcf_data_T_1_9.extend(small_vcf_data_T_9_)
    small_vcf_data = list(zip(*small_vcf_data_T_1_9))

    org_vcf_data_T_1_9.extend(org_vcf_data_T_9_)
    org_vcf_data = list(zip(*org_vcf_data_T_1_9))
    with open(os.path.splitext(org_vcf)[0] + "_cut_sample.vcf", "w") as fp_org, \
            open(os.path.splitext(small_vcf)[0] + "_cut_sample.vcf", "w") as fp_small:
        for i in small_vcf_data:
            tmp_str = "{}\n".format("\t".join(i))
            tmp_str = tmp_str.strip("chr")
            tmp_str = re.sub("NS=1595", "PR", tmp_str)
            tmp_str = re.sub("PASS", ".", tmp_str)
            # tmp_str = re.sub(":\d","",re.sub(":\d.\d+","",tmp_str))
            fp_small.write(tmp_str)
        for i in org_vcf_data:
            tmp_str = "{}\n".format("\t".join(i))
            tmp_str = tmp_str.strip("chr")
            tmp_str = re.sub("NS=1595", "PR", tmp_str)
            tmp_str = re.sub("PASS", ".", tmp_str)
            # tmp_str = re.sub(":\d", "", re.sub(":\d.\d+", "", tmp_str))
            fp_org.write(tmp_str)


def map_vcf_sample_id(vcf_in, id_pair_in, vcf_out):
    with open(id_pair_in, "r") as fp:
        id_pair_list = [i.strip().split("\t") for i in fp.readlines() if len(i.strip()) > 0]
    id_pair_dict = {}
    for SL_id, num_id in id_pair_list:
        id_pair_dict[SL_id] = num_id
    with open(vcf_in, "r") as fp_in, open(vcf_out, "w") as fp_out:
        while True:
            line_data = fp_in.readline()
            if not line_data:
                break
            if line_data.startswith("##"):
                fp_out.write(line_data)
                continue
            if line_data.startswith("#"):
                tmp_str = ""
                line_list = line_data.strip().split("\t")
                print("line_list[0]=[{}]".format(line_list))
                for i in line_list:
                    if i not in id_pair_dict:
                        tmp_str = "{0}\t{1}".format(tmp_str, i)
                    else:
                        tmp_str = "{0}\t{1}".format(tmp_str, id_pair_dict[i])
                tmp_str = tmp_str.strip()
                tmp_str = "{0}{1}".format(tmp_str, "\n")
                fp_out.write(tmp_str)
                continue
            fp_out.write(line_data)


def cut_sample(vcf_file, head_front_list, selected_head_back_list, back_selector):
    def sort_key(elem):
        return elem[0]

    with open(vcf_file, "r") as fp_in, open(os.path.splitext(vcf_file)[0] + "_cut_sample.vcf", "w") as fp_out:
        tmp_list = copy.copy(selected_head_back_list)
        tmp_list.sort()
        head_list = head_front_list + tmp_list
        tmp_str = "{}\n".format("\t".join(head_list))
        fp_out.write(tmp_str)
        fp_in.seek(0, 0)
        fp_in.readline()
        icounter = 0
        ret_str = ""
        while True:
            data_line = fp_in.readline()
            if not data_line:
                # fp_out.write(ret_str)
                # ret_str = ""
                logging.debug("handled {0} data lines in {1}".format(icounter, vcf_file))
                break
            data_list = [i.split(":")[0] for i in data_line.strip().split("\t")]
            data_front_list = data_list[:9]
            data_back_list = data_list[9:]
            selected_back_list = list(compress(data_back_list, back_selector))
            selected_back_list = zip(selected_head_back_list, selected_back_list)
            selected_back_list.sort(key=sort_key)
            selected_back_list = list(zip(*selected_back_list)[1])
            # print data_front_list
            # print selected_back_list
            data_list = data_front_list + selected_back_list
            tmp_str = "{}\n".format("\t".join(data_list))
            tmp_str = tmp_str.strip("chr")
            tmp_str = re.sub("NS=1595", "PR", tmp_str)
            tmp_str = re.sub("PASS", ".", tmp_str)
            fp_out.write(tmp_str)

            icounter += 1
            if icounter % 10000 == 0:
                logging.debug("handled {0} data lines in {1}".format(icounter, vcf_file))
                # fp_out.write(ret_str)
                # ret_str = ""
    logging.debug("done")
    return 0


def cut_org_exported_vcf_samples_new(org_vcf_file, small_vcf_file):
    logging.basicConfig(filename="cut_org_exported_vcf_samples_new.log", level=logging.DEBUG, format=log_format,
                        filemode="w")
    logging.debug("Begin")
    with open(small_vcf_file, "r") as fp:
        small_head_list = [i.split("_")[0] for i in fp.readline().strip().split("\t")]
        small_head_front_list = small_head_list[:9]
        small_head_back_list = small_head_list[9:]
        small_head_set = set(small_head_back_list)
    with open(org_vcf_file, "r") as fp:
        org_head_list = [i.split("_")[0] for i in fp.readline().strip().split("\t")]
        org_head_front_list = org_head_list[:9]
        org_head_back_list = org_head_list[9:]
        org_head_set = set(org_head_back_list)
    interset_head_set = small_head_set & org_head_set
    small_back_selector = [i in interset_head_set for i in small_head_back_list]
    org_back_selector = [i in interset_head_set for i in org_head_back_list]
    small_selected_head_back_list = list(compress(small_head_back_list, small_back_selector))
    org_selected_head_back_list = list(compress(org_head_back_list, org_back_selector))
    logging.debug("handling org vcf...")
    cut_sample(org_vcf_file, org_head_front_list, org_selected_head_back_list, org_back_selector)
    logging.debug("handling small vcf...")
    cut_sample(small_vcf_file, small_head_front_list, small_selected_head_back_list, small_back_selector)


L2PSTEP = 100000


def handle_l2p(file_name):
    icounter = 1
    job_num = 0
    if not os.path.exists(file_name + ".l2p"):
        with open(file_name, "r") as fp, open(file_name + ".l2p", "w") as fp_l2p:
            while True:
                if icounter % L2PSTEP == 0 or icounter == 1:
                    fp_l2p.write("{0}\t{1}\n".format(icounter, fp.tell()))
                    job_num += 1
                data_line = fp.readline()
                icounter += 1
                if icounter % 10000 == 0:
                    print("handling line {} in sam file".format(icounter))

                if not data_line.strip():
                    break
    else:
        print(file_name + ".l2p is OK.")
        with open(file_name + ".l2p", "r") as fp_l2p:
            while len(fp_l2p.readline().strip()) > 0:
                job_num += 1
    return job_num


def sort_key1(elem):
    return elem[0]


def cut_sample_worker(vcf_file, job_str, selected_head_back_list, back_selector, prefix, worker_id):
    logging.debug("worker_id = [{0}] job_str = [{1}]".format(worker_id, job_str.strip()))
    line_num, offset = job_str.split("\t")
    if int(line_num) == 1:
        step = L2PSTEP - 1

    else:
        step = L2PSTEP
    # print("job_str={}".format(job_str))
    with open(vcf_file, "r") as fp_in, open("{0}_{1}.vcf".format(prefix, line_num), "w") as fp_out:
        fp_in.seek(int(offset))
        icounter = 0
        while True:
            data_line = fp_in.readline()
            if not data_line:
                # logging.debug("handled {0} data lines in {1}".format(icounter, vcf_file))
                break
            data_list = [i.split(":")[0] for i in data_line.strip().split("\t")]
            data_front_list = data_list[:9]
            data_back_list = data_list[9:]
            selected_back_list = list(compress(data_back_list, back_selector))
            selected_back_list = zip(selected_head_back_list, selected_back_list)
            selected_back_list.sort(key=sort_key1)
            selected_back_list = list(zip(*selected_back_list)[1])
            data_list = data_front_list + selected_back_list
            tmp_str = "{}\n".format("\t".join(data_list))
            tmp_str = tmp_str.strip("chr")
            tmp_str = re.sub("NS=1595", "PR", tmp_str)
            tmp_str = re.sub("PASS", ".", tmp_str)
            fp_out.write(tmp_str)

            icounter += 1
            if icounter == step:
                # logging.debug("handled {0} data lines in {1}".format(icounter, vcf_file))
                break
    logging.debug("worker{} done".format(worker_id))
    return 0


def cut_org_exported_vcf_samples_multiple(org_vcf_file, small_vcf_file, cpu_num):
    logging.basicConfig(filename="cut_org_exported_vcf_samples_multiple.log", level=logging.DEBUG, format=log_format,
                        filemode="w")
    logging.debug("Begin")
    org_output = os.path.splitext(org_vcf_file)[0] + "_cut_sample_m.vcf"
    small_output = os.path.splitext(small_vcf_file)[0] + "_cut_sample_m.vcf"
    org_job_num = handle_l2p(org_vcf_file)
    small_job_num = handle_l2p(small_vcf_file)
    with open(small_vcf_file, "r") as fp:
        small_head_list = [i.split("_")[0] for i in fp.readline().strip().split("\t")]
        small_head_front_list = small_head_list[:9]
        small_head_back_list = small_head_list[9:]
        small_head_set = set(small_head_back_list)
    with open(org_vcf_file, "r") as fp:
        org_head_list = [i.split("_")[0] for i in fp.readline().strip().split("\t")]
        org_head_front_list = org_head_list[:9]
        org_head_back_list = org_head_list[9:]
        org_head_set = set(org_head_back_list)
    interset_head_set = small_head_set & org_head_set
    small_back_selector = [i in interset_head_set for i in small_head_back_list]
    org_back_selector = [i in interset_head_set for i in org_head_back_list]
    small_selected_head_back_list = list(compress(small_head_back_list, small_back_selector))
    org_selected_head_back_list = list(compress(org_head_back_list, org_back_selector))

    process_pool = multiprocessing.Pool(processes=int(cpu_num))
    icounter = 0
    with open(org_vcf_file + ".l2p", "r") as fp:
        while True:
            line = fp.readline()
            if not line:
                break
            if len(line.strip()) == 0:
                continue
            icounter += 1
            process_pool.apply_async(cut_sample_worker, (org_vcf_file, line,
                                                         org_selected_head_back_list, org_back_selector,
                                                         "cut_multiple_1_", icounter))
    with open(small_vcf_file + ".l2p", "r") as fp:
        while True:
            line = fp.readline()
            if not line:
                break
            if len(line.strip()) == 0:
                continue
            icounter += 1
            process_pool.apply_async(cut_sample_worker, (small_vcf_file, line,
                                                         small_selected_head_back_list, small_back_selector,
                                                         "cut_multiple_2_", icounter))
    process_pool.close()
    process_pool.join()

    logging.debug("handling org vcf...")

    # ph = Popen(["rm -rf {}".format(org_output)], shell=True, stdout=PIPE)
    # ph.wait()
    ph = Popen(["cat cut_multiple_1__* > {}".format(org_output)], shell=True, stdout=PIPE)
    ph.wait()
    ph = Popen(["rm -rf cut_multiple_1_*"], shell=True, stdout=PIPE)
    ph.wait()
    logging.debug("handling small vcf...")

    # ph = Popen(["rm -rf {}".format(small_output)], shell=True, stdout=PIPE)
    # ph.wait()
    ph = Popen(["cat cut_multiple_2__* > {}".format(small_output)], shell=True, stdout=PIPE)
    ph.wait()
    ph = Popen(["rm -rf cut_multiple_2_*"], shell=True, stdout=PIPE)
    ph.wait()


def abp_in_regions(abp, region_list):
    for start, end in region_list:
        if start <= abp <= end:
            return True
    return False


def export_variance_info_with_gene_list_gene_name(db_file, gene_list_file, output_variance, output_sample):
    col_list = ['bystro_gene_name', 'annovar_gene_name', 'vep_gene_id', 'spliceAI_gene_name',
                'dmis_gene_name', 'dsplicing_gene_name']
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()
    cursor.execute('select gene_id, gene_name from gene_table')
    gene_id2name_dict = dict(cursor.fetchall())
    cursor.execute('select gene_name, gene_id from gene_table')
    gene_name2id_dict = dict(cursor.fetchall())

    cursor.execute('select gen_id from sampleChdPhenotype WHERE CTD=1')
    ctd_case_set = set(["sample_{}".format(i[0]) for i in cursor.fetchall()])
    cursor.execute('select gen_id from sampleChdPhenotype WHERE heart6=1')
    heart6_case_set = set(["sample_{}".format(i[0]) for i in cursor.fetchall()])
    cursor.execute('select gen_id from sampleChdPhenotype WHERE heart6=0')
    heart6_control_set = set(["sample_{}".format(i[0]) for i in cursor.fetchall()])
    with open(gene_list_file, 'r') as fp:
        gene_list = [i for i in re.split('[\r\n \t;,|]', fp.read()) if len(i) > 0]
    gene_name_set = set([])
    gene_id_set = set([])
    for gene in gene_list:
        if gene in gene_name2id_dict:
            gene_name_set.add(gene)
            gene_id_set.add(gene_name2id_dict[gene])
        elif gene in gene_id2name_dict:
            gene_name_set.add(gene_id2name_dict[gene])
            gene_id_set.add(gene)
    cursor.execute("PRAGMA table_info([variance])")
    table_info = cursor.fetchall()
    head_str = ""
    for elem in table_info[1:]:
        if not head_str:
            head_str = elem[1]
        elif not elem[1].startswith('sample_'):
            head_str += "\t" + elem[1]
    head_str += "\tctd case num\tctd case samples\theart6 case num\t" \
                "heart6 case samples\theart6 control num\theart6 control samples"
    col2index_dict = dict([[i[1], i[0]] for i in table_info])
    target_variance_list = []
    with open(output_variance, 'w') as fp:
        fp.write("gene_name\tgene_id\t" + head_str + '\n')
        cursor.execute('select * from variance')
        for variance in cursor.fetchall():
            sub_variance_line = [variance[int(i[0])] for i in table_info[1:] if not i[1].startswith('sample_')]
            ctd_num = len([i for i in table_info if (i[1] in ctd_case_set and
                                                     type(variance[int(i[0])]) == int and
                                                     variance[int(i[0])] > 0)])
            ctd_samples = ','.join([i[1].strip('sample_') for i in table_info if (i[1] in ctd_case_set and
                                                                                  type(variance[
                                                                                           int(i[0])]) == int and
                                                                                  variance[int(i[0])] > 0)])
            heart6_case_num = len([i for i in table_info if (i[1] in heart6_case_set and
                                                             type(variance[int(i[0])]) == int and
                                                             variance[int(i[0])] > 0)])
            heart6_case_samples = ','.join(
                [i[1].strip('sample_') for i in table_info if (i[1] in heart6_case_set and
                                                               type(variance[int(i[0])]) == int and
                                                               variance[int(i[0])] > 0)])
            heart6_contrl_num = len([i for i in table_info if (i[1] in heart6_control_set and
                                                               type(variance[int(i[0])]) == int and
                                                               variance[int(i[0])] > 0)])
            heart6_contrl_samples = ','.join(
                [i[1].strip('sample_') for i in table_info if (i[1] in heart6_control_set and
                                                               type(variance[int(i[0])]) == int and
                                                               variance[int(i[0])] > 0)])
            anno_set = set([])  # gene name
            for col_name in col_list:
                anno_str = variance[col2index_dict[col_name]]
                if not anno_str:
                    continue
                # anno_set = anno_set | set([] if not anno_str else anno_str.split(';'))
                for anno in anno_str.split(';'):
                    if anno not in gene_id_set and anno not in gene_name_set:
                        # print('WARNING: unrecognized anno: {}'.format(anno))
                        continue
                    if anno in gene_id_set:
                        anno_set.add(gene_id2name_dict[anno])
                        continue
                    anno_set.add(anno)
                # if anno_set & gene_name_set or anno_set & gene_id_set:
                #     fp.write("\t".join([str(i) for i in variance]) + '\n')
                #     target_variance_list.append(variance)
                #     break
            if anno_set & gene_name_set:
                target_variance_list.append(variance)
                for gene_name in anno_set & gene_name_set:
                    fp.write("{0}\t{1}\t".format(gene_name, gene_name2id_dict[gene_name]) +
                             "\t".join([str(i) for i in sub_variance_line]) + '\t{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n'
                                                                              ''.format(ctd_num, ctd_samples,
                                                                                        heart6_case_num,
                                                                                        heart6_case_samples,
                                                                                        heart6_contrl_num,
                                                                                        heart6_contrl_samples))

    print("len(target_variance_list)={}".format(len(target_variance_list)))
    if len(target_variance_list) == 0:
        print('Do nothing.')
        conn.commit()
        conn.close()
        return
    target_variance_list_T = zip(*target_variance_list)
    tmp = [[i[1].strip('sample_'), sum([k for k in target_variance_list_T[i[0]] if type(k) == int])] for i in table_info
           if i[1].startswith('sample_')]
    print('len(tmp)={}'.format(len(tmp)))
    gen_id_list = [j[0] for j in tmp if j[1] > 0]
    cursor.execute('SELECT * from sampleChdPhenotype where gen_id in ({})'.format(','.join(gen_id_list)))
    sample_data = cursor.fetchall()
    cursor.execute("PRAGMA table_info([sampleChdPhenotype])")
    table_info = cursor.fetchall()
    with open(output_sample, 'w') as fp:
        fp.write('\t'.join([i[1] for i in table_info]) + '\n')  # head
        for sample in sample_data:
            fp.write('\t'.join([str(i) for i in sample]) + '\n')
    conn.commit()
    conn.close()
    print('all done')


def export_vcf_with_gene_list(db_file, variance_table, fai_in, gene_list, output):
    chr2offset_dict = parse_fai(fai_in)
    region_list = []
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()
    print("loading gene_list...")
    with open(gene_list, "r") as fp:
        fp.readline()
        gene_data = [i.strip().split("\t")[0] for i in fp.readlines() if len(i.strip()) > 0]

    for gene_id in gene_data:
        cmd_str = "SELECT chr, start_pos, end_pos, chr2, start_pos2, end_pos2 FROM gene_table WHERE gene_id='{0}'" \
                  "".format(gene_id)
        cursor.execute(cmd_str)
        chr1, start1, end1, chr2, start2, end2 = cursor.fetchall()[0]
        region_list.append([chr_pos2absolute_pos(str(chr1), start1, chr2offset_dict),
                            chr_pos2absolute_pos(str(chr1), end1, chr2offset_dict)])
        if chr2 is not None:
            region_list.append([chr_pos2absolute_pos(str(chr2), start2, chr2offset_dict),
                                chr_pos2absolute_pos(str(chr2), end2, chr2offset_dict)])
    print("got {0} regions".format(len(region_list)))

    with open(output, "w") as fp_out:
        # while True:
        #     org_vcf_line = fp_vcf.readline()
        #     if org_vcf_line.startswith("##"):
        #         fp_out.write(org_vcf_line)
        #         continue
        #     if org_vcf_line.startswith("#"):
        #         vcf_head_list = org_vcf_line.strip().split("\t")
        #         fp_out.write(org_vcf_line)
        #         break
        # cmd_str = "SELECT chr, pos, vcf_id, ref, alt, vcf_qual, vcf_filter, vcf_info, vcf_format"
        # for i in xrange(9, len(vcf_head_list), 1):
        #     cmd_str += ", sample_{0}".format(vcf_head_list[i])
        # cmd_str += " FROM {0}".format(variance_table)
        cmd_str = "PRAGMA table_info([{0}])".format(variance_table)
        cursor.execute(cmd_str)
        table_info = cursor.fetchall()
        head_str = ""
        for elem in table_info:
            if not head_str:
                head_str = elem[1]
            else:
                head_str += "\t" + elem[1]
        head_str += "\n"
        fp_out.write(head_str)
        cmd_str = "SELECT * from {0}".format(variance_table)
        # print("cmd_str={0}".format(cmd_str))
        print("handling variance...")
        icounter = 0
        for element in cursor.execute(cmd_str):
            icounter += 1
            list_element = list(element)
            abp = chr_pos2absolute_pos(str(list_element[1]), list_element[2], chr2offset_dict)
            if not abp_in_regions(abp, region_list):
                print("handled {0} variance".format(icounter))
                continue
            # for index in xrange(len(element)):
            #     if index < 9:
            #         continue
            #     if list_element[index] == 0:
            #         list_element[index] = "0/0"
            #     elif list_element[index] == 1:
            #         list_element[index] = "0/1"
            #     elif list_element[index] == 2:
            #         list_element[index] = "1/1"
            #     else:
            #         list_element[index] = "./."
            # list_element[8] = list_element[8].split(":")[0]  # vcf_format
            # tmp_str = "{}\n".format("\t".join([str(i) for i in list_element]))
            # if not tmp_str.startswith("chr"):
            #     tmp_str = "chr{}".format(tmp_str)
            tmp_str = "\t".join([str(i) for i in list_element])
            fp_out.write(tmp_str + "\n")
            print("handled {0} variance".format(icounter))


def db_add_ccrs_limbr(db_file, ccrs_file, domain_file, exone_file, fai_in):
    def pos_in_region_list(pos, region_list):
        pos = int(pos)
        for i in xrange(len(region_list)):
            start, end = region_list[i]
            if start <= pos <= end:
                return i
        return -1

    def region_in_region_list(region, region_list):
        index_list = []
        for pos in xrange(int(region[0]), int(region[1]) + 1, 1):
            curr_index = pos_in_region_list(pos, region_list)
            if curr_index < 0:
                return []
            index_list.append(curr_index)
        return index_list

    logging.basicConfig(filename="db_add_ccrs_limbr.log", level=logging.DEBUG, format=log_format, filemode="w")
    chr2offset_dict = parse_fai(fai_in)
    ccrs_region_list = []
    ccrs_score_list = []  # type: list[float]
    print("loading ccrs file...")
    logging.debug("loading ccrs file...")
    with open(ccrs_file, "r") as fp:
        while True:
            data_line = fp.readline()
            if not data_line:
                break
            if not data_line.strip():
                continue
            chrom, start, end, score = data_line[3:].strip().split("\t")
            ccrs_region_list.append([chr_pos2absolute_pos(chrom, int(start) + 1, chr2offset_dict),
                                     chr_pos2absolute_pos(chrom, int(end), chr2offset_dict)])
            ccrs_score_list.append(float(score))

    domain_limbr_region_list = []
    domain_limbr_list = []  # type: list[float]
    # domain_limbr_quartile_list = []  # type: list[str]
    print("loading domain file...")
    logging.debug("loading domain file...")
    with open(domain_file, "r") as fp:
        while True:
            data_line = fp.readline()
            if not data_line:
                break
            if not data_line.strip():
                continue
            chrom, start, end, score = data_line[3:].strip().split("\t")
            domain_limbr_region_list.append([chr_pos2absolute_pos(chrom, int(start) + 1, chr2offset_dict),
                                             chr_pos2absolute_pos(chrom, int(end), chr2offset_dict)])
            domain_limbr_list.append(float(score))
            # domain_limbr_quartile_list.append(quartile)

    exone_limbr_region_list = []
    exone_limbr_list = []  # type: list[float]
    # exone_limbr_quartile_list = []  # type: list[str]
    print("loading exone file...")
    logging.debug("loading exone file...")
    with open(exone_file, "r") as fp:
        while True:
            data_line = fp.readline()
            if not data_line:
                break
            if not data_line.strip():
                continue
            chrom, start, end, score = data_line[3:].strip().split("\t")
            exone_limbr_region_list.append([chr_pos2absolute_pos(chrom, int(start) + 1, chr2offset_dict),
                                            chr_pos2absolute_pos(chrom, int(end), chr2offset_dict)])
            exone_limbr_list.append(float(score))
            # exone_limbr_quartile_list.append(quartile)

    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()
    db_add_col(cursor, "variance", "ccrs", "varchr(20)")
    db_add_col(cursor, "variance", "domain_limbr", "varchr(20)")
    db_add_col(cursor, "variance", "exone_limbr", "varchr(20)")
    db_add_col(cursor, "synonymous_snp", "ccrs", "varchr(20)")
    db_add_col(cursor, "synonymous_snp", "domain_limbr", "varchr(20)")
    db_add_col(cursor, "synonymous_snp", "exone_limbr", "varchr(20)")
    sql_list = []
    icounter = 0
    print("begin to handle variance")
    logging.debug("begin to handle variance")
    for id, chrom, pos, ref, alt in cursor.execute("SELECT id, chr, pos, ref, alt  FROM variance"):
        abp = chr_pos2absolute_pos(chrom, pos, chr2offset_dict)
        if len(ref) > len(alt):  # del
            delta = len(ref) - len(alt)
            region = [abp + 1, abp + delta]
        else:  # ins & snv
            region = [abp, abp]

        # ccrs
        index_list = region_in_region_list(region, ccrs_region_list)
        if len(index_list) > 0:
            tmp = [ccrs_score_list[i] for i in index_list]
            if len(Counter(tmp)) == 1:
                sql_list.append("UPDATE variance SET ccrs = '{0}' WHERE id = {1}"
                                "".format(float(sum(tmp)) / len(tmp), id))
            else:
                sql_list.append("UPDATE variance SET ccrs = NULL WHERE id = {0}".format(id))
        else:
            sql_list.append("UPDATE variance SET ccrs = NULL WHERE id = {0}".format(id))

        # domain limbr
        index_list = region_in_region_list(region, domain_limbr_region_list)
        if len(index_list) > 0:
            tmp = [domain_limbr_list[i] for i in index_list]
            sql_list.append("UPDATE variance SET domain_limbr = '{0}' WHERE id = {1}"
                            "".format(float(sum(tmp)) / len(tmp), id))
        else:
            sql_list.append("UPDATE variance SET domain_limbr = NULL WHERE id = {0}".format(id))

        # exone limbr
        index_list = region_in_region_list(region, exone_limbr_region_list)
        if len(index_list) > 0:
            tmp = [exone_limbr_list[i] for i in index_list]
            sql_list.append("UPDATE variance SET exone_limbr = '{0}' WHERE id = {1}"
                            "".format(float(sum(tmp)) / len(tmp), id))
        else:
            sql_list.append("UPDATE variance SET exone_limbr = NULL WHERE id = {0}".format(id))

        icounter += 1
        # if icounter % 1000 == 0:
        print("handled {0} variances".format(icounter))
        logging.debug("handled {0} variances".format(icounter))

    for id, chrom, pos, ref, alt in cursor.execute("SELECT id, chr, pos, ref, alt  FROM synonymous_snp"):
        abp = chr_pos2absolute_pos(chrom, pos, chr2offset_dict)
        if len(ref) > len(alt):
            delta = len(ref) - len(alt)
            region = [abp + 1, abp + delta]
        else:
            region = [abp, abp]

        # ccrs
        index_list = region_in_region_list(region, ccrs_region_list)
        if len(index_list) > 0:
            tmp = [ccrs_score_list[i] for i in index_list]
            if len(Counter(tmp)) == 1:
                sql_list.append("UPDATE synonymous_snp SET ccrs = '{0}' WHERE id = {1}"
                                "".format(float(sum(tmp)) / len(tmp), id))
            else:
                sql_list.append("UPDATE synonymous_snp SET ccrs = NULL WHERE id = {0}".format(id))
        else:
            sql_list.append("UPDATE synonymous_snp SET ccrs = NULL WHERE id = {0}".format(id))

        # domain limbr
        index_list = region_in_region_list(region, domain_limbr_region_list)
        if len(index_list) > 0:
            tmp = [domain_limbr_list[i] for i in index_list]
            sql_list.append("UPDATE synonymous_snp SET domain_limbr = '{0}' WHERE id = {1}"
                            "".format(float(sum(tmp)) / len(tmp), id))
        else:
            sql_list.append("UPDATE synonymous_snp SET domain_limbr = NULL WHERE id = {0}".format(id))

        # exone limbr
        index_list = region_in_region_list(region, exone_limbr_region_list)
        if len(index_list) > 0:
            tmp = [exone_limbr_list[i] for i in index_list]
            sql_list.append("UPDATE synonymous_snp SET exone_limbr = '{0}' WHERE id = {1}"
                            "".format(float(sum(tmp)) / len(tmp), id))
        else:
            sql_list.append("UPDATE synonymous_snp SET exone_limbr = NULL WHERE id = {0}".format(id))
        icounter += 1
        # if icounter % 1000 == 0:
        print("handled {0} variances".format(icounter))
        logging.debug("handled {0} variances".format(icounter))
    print("begin execute sqls")
    logging.debug("begin execute sqls")
    icounter = 0
    for sql in sql_list:
        cursor.execute(sql)
        icounter += 1
        if icounter % 1000 == 0:
            print("executed {0} sqls".format(icounter))
            logging.debug("executed {0} sqls".format(icounter))
    cursor.close()
    print("committing...")
    logging.debug("committing...")
    conn.commit()
    conn.close()
    logging.debug("all done")
    print("all done")


def analyze_fisher_test_variance_2(database, variance_table, sample_restrict, fai_in, path):
    os.chdir(path)
    if variance_table == "variance":
        annotator_list_base = ["annovar = 1", "bystro = 1", "dmis = 1", "dsplicing = 1", "spidex = 1", "spliceAI = 1",
                               "vep = 1"]
    elif variance_table == "synonymous_snp":
        annotator_list_base = ["annovar = 1", "bystro = 1", "vep = 1"]
    else:
        print("variance_table could only be variance or synonymous_snp")
        return
    domain_limbr_list = [" AND domain_limbr <= -3.654145493", " AND domain_limbr <= -1.751155961",
                         " AND domain_limbr <= -0.699941137", " AND domain_limbr <= 0.091314981",
                         " AND domain_limbr <= 0.753317939"]
    exone_limbr_list = [" AND exone_limbr <= -1.538383601", " AND exone_limbr <= -0.572107415",
                        " AND exone_limbr <= -0.023396951", " AND exone_limbr <= 0.390431738",
                        " AND exone_limbr <= 0.750902136"]
    annotator_list = []
    annotator_list.extend(list(itertools.combinations(annotator_list_base, 1)))
    annotator_list.extend(list(itertools.combinations(annotator_list_base, 2)))
    annotator_list.extend(list(itertools.combinations(annotator_list_base, 3)))
    annotator_list.extend(list(itertools.combinations(annotator_list_base, 4)))
    annotator_list.extend(list(itertools.combinations(annotator_list_base, 5)))
    annotator_list.extend(list(itertools.combinations(annotator_list_base, 6)))
    annotator_list.extend(list(itertools.combinations(annotator_list_base, 7)))
    annotator_list = [" AND ({})".format(" OR ".join(i)) for i in annotator_list]
    annotator_list.append("")
    icounter = 1
    icounter2 = 0
    script_str = ''
    for phenotype in ["heart6", "tof6", "aorticarch6", "tof_or_pta_or_iaab6"]:
        for annotator in annotator_list:
            for freq in ["bystro_sampleMaf <= 0.01", "bystro_sampleMaf <= 0.001",
                         "bystro_sampleMaf <= 0.005"]:
                for ph in ["", " AND bystro_phastCons >= 0.2", " AND bystro_phastCons >= 0.3",
                           " AND bystro_phastCons >= 0.4",
                           " AND bystro_phastCons >= 0.5", " AND bystro_phastCons >= 0.6",
                           " AND bystro_phastCons >= 0.7",
                           " AND bystro_phastCons >= 0.8", " AND bystro_phyloP >= -1", " AND bystro_phyloP >= 0",
                           " AND bystro_phyloP >= 1", " AND bystro_phyloP >= 2", " AND bystro_phyloP >= 3",
                           " AND bystro_phyloP >= 4"]:
                    for cadd in ["", " AND bystro_cadd >= 5", " AND bystro_cadd >= 10", " AND bystro_cadd >= 15",
                                 " AND bystro_cadd >= 20", " AND bystro_cadd >= 25", " AND bystro_cadd >= 30"]:
                        for ccrs in [" AND ccrs >= 95", " AND ccrs >= 90", " AND ccrs >= 85",
                                     " AND ccrs >= 80", " AND ccrs >= 75"]:
                            for domain_limbr in domain_limbr_list:
                                for exone_limbr in exone_limbr_list:
                                    variance_restrict = "{0}{1}{2}{3}{4}{5}{6}".format(freq, annotator, ph, cadd,
                                                                                       ccrs, domain_limbr, exone_limbr)
                                    output_name = re.sub("bystro_sampleMaf <= 0.", "freq", variance_restrict)
                                    output_name = re.sub(" AND \(", "_and_", output_name)
                                    output_name = re.sub("\) AND ", "_and_", output_name)
                                    output_name = re.sub(" AND ", "_and_", output_name)
                                    output_name = re.sub(" OR ", "_or_", output_name)
                                    output_name = re.sub(" = ", "", output_name)
                                    output_name = re.sub("bystro_", "", output_name)
                                    output_name = re.sub(" >= ", "", output_name)
                                    output_name = re.sub("\)$", "", output_name)
                                    output_name = "{0}_{1}.table".format(phenotype, output_name)
                                    # print output_name
                                    if script_str == "":
                                        script_str = """#!/bin/bash
#$ -j y
#$ -S /bin/bash
#$ -cwd
#$ -N ct{6}
#$ -pe smp 1
#$ -l h_vmem=8G
#$ -m bes
# -q all.q
export LC_ALL=C
export MALLOC_ARENA_MAX=4
module load python/2.7.15/gcc.4.4.7
module load sqlite3/3.8.11/gcc.4.4.7
time=`date`
echo "==START $time =="
~/miniconda2/bin/python ~/wyj/.code/wgsa.py build_contingency_table_new {0} {1} sampleChdPhenotype '{2}' gene_table '' {7} '{3}' {4} {5} {6}
echo {5} is done
time=`date`
echo == $time ==
""".format(database, phenotype, sample_restrict, variance_restrict, fai_in, output_name, icounter2, variance_table)
                                    else:
                                        script_str += "~/miniconda2/bin/python ~/wyj/.code/wgsa.py build_contingency_table_new {0} {1} " \
                                                      "sampleChdPhenotype '{2}' gene_table '' {6} '{3}' {4} {5} {7}" \
                                                      "\necho {5} is done" \
                                                      "\ntime=`date`\necho == $time ==\n" \
                                                      "".format(database,
                                                                phenotype,
                                                                sample_restrict,
                                                                variance_restrict,
                                                                fai_in,
                                                                output_name,
                                                                variance_table,
                                                                icounter2)

                                    if icounter == 10000:
                                        with open("qsub{}.sh".format(icounter2), "w") as fp:
                                            fp.write(script_str + "\ndate; echo \"==END==\"")
                                        while True:
                                            if os.access("qsub{}.sh".format(icounter2), os.R_OK):
                                                break
                                            time.sleep(1)
                                        pp = Popen(["qsub qsub{}.sh".format(icounter2)], shell=True)
                                        pp.wait()
                                        script_str = ""
                                        icounter = 0
                                        icounter2 += 1
                                    icounter += 1
    if len(script_str) > 0:
        with open("qsub{}.sh".format(icounter2), "w") as fp:
            fp.write(script_str + "\ndate; echo \"==END==\"")
        while True:
            if os.access("qsub{}.sh".format(icounter2), os.R_OK):
                break
            time.sleep(1)
        pp = Popen(["qsub qsub{}.sh".format(icounter2)], shell=True)
        pp.wait()

    print(
        "All the jobs has been submitted. Job number = {}\nIf any job has problem, kill it. And qsub the corresponding sh file".format(
            icounter2 + 1))


domain_limbr_quartile_dict = {10: -3.654145493, 20: -1.751155961, 30: -0.699941137, 40: 0.091314981,
                              50: 0.753317939, 60: 1.338304674, 70: 1.881692228}
exone_limbr_quartile_dict = {10: -1.538383601, 20: -0.572107415, 30: -0.023396951, 40: 0.390431738,
                             50: 0.750902136, 60: 1.084410106, 70: 1.412196324}


def analyze_fisher_test_variance_3(database, variance_table, sample_restrict, fai_in, path):
    os.chdir(path)
    if variance_table == "variance":
        annotator_list_base = ["annovar = 1 OR bystro = 1 OR vep = 1", "dmis = 1",
                               "dsplicing = 1 OR spliceAI = 1 OR spidex = 1"]
    elif variance_table == "synonymous_snp":
        annotator_list_base = ["annovar = 1 OR bystro = 1 OR vep = 1"]
    else:
        print("variance_table could only be variance or synonymous_snp")
        return
    annotator_list = []
    annotator_list.extend(list(itertools.combinations(annotator_list_base, 1)))
    annotator_list.extend(list(itertools.combinations(annotator_list_base, 2)))
    annotator_list.extend(list(itertools.combinations(annotator_list_base, 3)))
    annotator_list = ["({})".format(" OR ".join(i)) for i in annotator_list]

    domain_limbr_list = [" AND (domain_limbr <= -3.654145493 OR domain_limbr is NULL)",
                         " AND (domain_limbr <= -1.751155961 OR domain_limbr is NULL)",
                         " AND (domain_limbr <= -0.699941137 OR domain_limbr is NULL)",
                         " AND (domain_limbr <= 0.091314981 OR domain_limbr is NULL)",
                         " AND (domain_limbr <= 0.753317939 OR domain_limbr is NULL)",
                         " AND (domain_limbr <= 1.338304674 OR domain_limbr is NULL)",
                         " AND (domain_limbr <= 1.881692228 OR domain_limbr is NULL)"]
    ###10%, -3.654145493,20%, -1.751155961,30%, -0.699941137,40%,0.091314981, 50%, 0.753317939,60%  60%,1.338304674, 70% 1.881692228

    exone_limbr_list = [" AND (exone_limbr <= -1.538383601 OR exone_limbr is NULL)",
                        " AND (exone_limbr <= -0.572107415 OR exone_limbr is NULL)",
                        " AND (exone_limbr <= -0.023396951 OR exone_limbr is NULL)",
                        " AND (exone_limbr <= 0.390431738 OR exone_limbr is NULL)",
                        " AND (exone_limbr <= 0.750902136 OR exone_limbr is NULL)",
                        " AND (exone_limbr <= 1.084410106 OR exone_limbr is NULL)",
                        " AND (exone_limbr <= 1.412196324 OR exone_limbr is NULL)",
                        ""]
    all_ccrs_conditions_list = [" AND (ccrs >= 95 OR ccrs is NULL)",
                                " AND (ccrs >= 90 OR ccrs is NULL)",
                                " AND (ccrs >= 85 OR ccrs is NULL)",
                                " AND (ccrs >= 80 OR ccrs is NULL)",
                                " AND (ccrs >= 75 OR ccrs is NULL)",
                                " AND (ccrs >= 70 OR ccrs is NULL)",
                                " AND (ccrs >= 65 OR ccrs is NULL)"]
    ###important info about limbr and ccrs, regions with intolerance are denoted as values below certain cutoffs, which means the smaller the values,
    ###the more constrained the regions are, while CCRS are exactly the opposite. The recommended cutoff for limbr is <=50 percentile, while is >=95 percentile for ccrs
    all_ccrs_conditions_list.extend(domain_limbr_list)
    all_ccrs_conditions_list.extend(exone_limbr_list)

    if variance_table == "variance":
        ph_list = [' AND (bystro_phyloP >= -1 OR bystro_phyloP="!")',
                   ' AND (bystro_phyloP >= 0 OR bystro_phyloP="!")']
    else:
        ph_list = [' AND (bystro_phyloP >= -1 OR bystro_phyloP="!")',
                   ' AND (bystro_phyloP >= 0 OR bystro_phyloP="!")',
                   ""]
    icounter = 1
    icounter2 = 0
    script_str = ''
    for phenotype in ["aorticarch6"]:
        for annotator in annotator_list:
            for freq in [" AND bystro_sampleMaf <= 0.01", ""]:
                for ph in ph_list:
                    for cadd in [""]:
                        for ccrs in all_ccrs_conditions_list:
                            variance_restrict = "{0}{1}{2}{3}{4}".format(annotator, freq, ph, cadd, ccrs)
                            output_name = re.sub("bystro_sampleMaf <= 0.", "freq", variance_restrict)
                            output_name = re.sub(" AND \(", "_and_", output_name)
                            output_name = re.sub("\) AND ", "_and_", output_name)
                            output_name = re.sub(" AND ", "_and_", output_name)
                            output_name = re.sub(" OR ", "_or_", output_name)
                            output_name = re.sub(" = ", "", output_name)
                            output_name = re.sub("bystro_", "", output_name)
                            output_name = re.sub(" >= ", "", output_name)
                            output_name = re.sub("\)$", "", output_name)
                            output_name = "{0}_{1}.table".format(phenotype, output_name)
                            output_name = re.sub("_\(", "_", output_name)
                            output_name = re.sub("\)_", "_", output_name)
                            output_name = re.sub(" is ", "", output_name)
                            logging.debug("output_name=[{0}]".format(output_name))
                            # print output_name
                            if script_str == "":
                                script_str = """#!/bin/bash
#$ -j y
#$ -S /bin/bash
#$ -cwd
#$ -N ct{6}
#$ -pe smp 1
#$ -l h_vmem=8G
#$ -m bes
# -q all.q
export LC_ALL=C
export MALLOC_ARENA_MAX=4
module load python/2.7.15/gcc.4.4.7
module load sqlite3/3.8.11/gcc.4.4.7
time=`date`
echo "==START $time =="
~/miniconda2/bin/python ~/wyj/.code/wgsa.py build_contingency_table_new {0} {1} sampleChdPhenotype '{2}' gene_table '' {7} '{3}' {4} {5} {6}
echo {5} is done
time=`date`
echo == $time ==
""".format(database, phenotype, sample_restrict, variance_restrict, fai_in, output_name, icounter2, variance_table)
                            else:
                                script_str += "~/miniconda2/bin/python ~/wyj/.code/wgsa.py build_contingency_table_new {0} {1} " \
                                              "sampleChdPhenotype '{2}' gene_table '' {6} '{3}' {4} {5} {7}" \
                                              "\necho {5} is done" \
                                              "\ntime=`date`\necho == $time ==\n" \
                                              "".format(database,
                                                        phenotype,
                                                        sample_restrict,
                                                        variance_restrict,
                                                        fai_in,
                                                        output_name,
                                                        variance_table,
                                                        icounter2)

                            if icounter == 100:
                                with open("qsub{}.sh".format(icounter2), "w") as fp:
                                    fp.write(script_str + "\ndate; echo \"==END==\"")
                                while True:
                                    if os.access("qsub{}.sh".format(icounter2), os.R_OK):
                                        break
                                    time.sleep(1)
                                pp = Popen(["qsub qsub{}.sh".format(icounter2)], shell=True)
                                pp.wait()
                                script_str = ""
                                icounter = 0
                                icounter2 += 1
                            icounter += 1
    if len(script_str) > 0:
        with open("qsub{}.sh".format(icounter2), "w") as fp:
            fp.write(script_str + "\ndate; echo \"==END==\"")
        while True:
            if os.access("qsub{}.sh".format(icounter2), os.R_OK):
                break
            time.sleep(1)
        pp = Popen(["qsub qsub{}.sh".format(icounter2)], shell=True)
        pp.wait()

    print(
        "All the jobs has been submitted. Job number = {0}\nIf any job has problem, kill it. And qsub the corresponding sh file".format(
            icounter2 + 1))


def permutation_burden(file_in, out_put, times):
    times = int(times)
    with open(file_in, "r") as fp:
        head_list = fp.readline().strip().split("\t")
        data = [i.strip().split("\t") for i in fp.readlines() if len(i.strip()) > 0]
    data_t = list(zip(*data))
    target = list(data_t[-1])
    num_1 = target.count("1")
    num_2 = target.count("2")
    num_NA = target.count("NA")
    num_all = len(target)
    print("num_1={0} num_2={1} num_NA={2} num_all={3}".format(num_1, num_2, num_NA, num_all))
    # print(target)
    for i in xrange(times):
        shuffled_index = np.random.permutation(num_all)
        tmp = copy.copy(target)
        for j in xrange(len(tmp)):
            tmp[j] = "XXX"
        for j in xrange(0, num_1, 1):
            tmp[shuffled_index[j]] = "1"
        for j in xrange(num_1, num_1 + num_2, 1):
            tmp[shuffled_index[j]] = "2"
        for j in xrange(num_1 + num_2, num_1 + num_2 + num_NA, 1):
            tmp[shuffled_index[j]] = "NA"
        data_t.append(tmp)
        head_list.append("permutation{0}".format(i + 1))
    result = list(zip(*data_t))
    with open(out_put, "w") as fp:
        fp.write("\t".join(head_list) + "\n")
        for i in result:
            fp.write("\t".join(i) + "\n")


def permutation_fisher(db_file, phenotype,
                       sample_table_name, sample_restrict,
                       gene_table_name, gene_restrict,
                       variance_table, variance_restrict,
                       fai_in, path, times=100):
    logging.basicConfig(filename="{}.log".format(sys._getframe().f_code.co_name),
                        level=logging.DEBUG, format=log_format, filemode="w")
    logging.debug("begin")
    logging.debug("db_file=[{0}]".format(db_file))
    logging.debug("phenotype=[{0}]".format(phenotype))
    logging.debug("sample_restrict=[{0}]".format(sample_restrict))
    logging.debug("variance_restrict=[{0}]".format(variance_restrict))
    logging.debug("times=[{0}]".format(times))
    times = int(times)
    if phenotype not in ["heart6", "ps_andor_pa6", "raa6", "iaab6", "pta6",
                         "tof6", "asdall6", "asdalone6", "vsd6", "vsdalone6",
                         "tofall6", "purevsdalone6", "ps_or_pa_and_vsd6", "intracardiac6", "aorticarch6",
                         "heartnoasd6", "tof_or_pta6", "tof_or_pta_or_iaab6", "CTD"]:
        logging.error("illegal phenotype [{}]".format(phenotype))
        return

    chr2offset_dict = parse_fai(fai_in)
    logging.debug("begin select control list and case list...")
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()
    format_sample_restrict = "" if not sample_restrict else " AND {}".format(sample_restrict)
    cursor.execute("SELECT s.gen_id FROM {2} AS s WHERE s.{0}='0'{1}".format(phenotype,
                                                                             format_sample_restrict,
                                                                             sample_table_name))
    # control ids
    control_id_list = [i[0] for i in cursor.fetchall()]
    print("control number={}".format(len(control_id_list)))
    cursor.execute("SELECT s.gen_id FROM {2} AS s WHERE s.{0}='1'{1}".format(phenotype,
                                                                             format_sample_restrict,
                                                                             sample_table_name))
    # case ids
    case_id_list = [i[0] for i in cursor.fetchall()]
    print("case number={}".format(len(case_id_list)))

    # prepare gene data
    format_gene_restrict = "" if not gene_restrict else " WHERE {}".format(gene_restrict)
    cmd_str = "SELECT g.gene_id, g.gene_name, g.chr, g.start_pos, g.end_pos, g.chr2, g.start_pos2, g.end_pos2 " \
              "FROM {0} AS g{1}".format(gene_table_name, format_gene_restrict)
    cursor.execute(cmd_str)
    gene_data = cursor.fetchall()

    # prepare the variance data
    cmd_str = "SELECT id, chr, pos, "
    for control_id in control_id_list:
        cmd_str = "{0}sample_{1}, ".format(cmd_str, control_id)
    for case_id in case_id_list:
        cmd_str = "{0}sample_{1}, ".format(cmd_str, case_id)
    cmd_str = cmd_str.strip(", ")
    format_variance_restrict = "" if not variance_restrict else " WHERE {}".format(variance_restrict)
    cmd_str = "{0} FROM {1}{2}".format(cmd_str, variance_table, format_variance_restrict)
    logging.debug("getting the variance data from db under restrict [{}]".format(variance_restrict))
    cursor.execute(cmd_str)
    org_variance_list = cursor.fetchall()
    logging.debug("len(org_variance_list)={0}".format(len(org_variance_list)))
    # Add absolute position, create sort, index
    for i in xrange(len(org_variance_list)):
        org_variance_list[i] = list(org_variance_list[i])
        org_variance_list[i].append(chr_pos2absolute_pos(str(org_variance_list[i][1]),
                                                         org_variance_list[i][2],
                                                         chr2offset_dict))
    org_variance_list.sort(key=sort_key_last)
    variance_list_index = build_data_ram_index(org_variance_list, len(org_variance_list[0]))
    # permutation
    for times in xrange(1, times + 1, 1):
        print("permutating {}".format(times))
        variance_list = copy.copy(org_variance_list)
        variance_list_t = zip(*variance_list)
        tmp = variance_list_t[3:-1]
        shuffled_index = np.random.permutation(len(case_id_list) + len(control_id_list))
        for i in xrange(3, 3 + len(case_id_list) + len(control_id_list), 1):
            variance_list_t[i] = tmp[shuffled_index[i - 3]]
        variance_list = list(zip(*variance_list_t))
        logging.debug("begin handle genes")
        with open(os.path.join(path, "fisher_permutation{0}".format(times)), "w") as fp:
            fp.write("##phenotype:\"{0}\"\n##control number:{1}\n##case number:{2}\n"
                     "##sample restrict:\"{3}\"\n##sample table:\"{4}\"\n##gene table:\"{5}\"\n"
                     "##gene restrict:\"{6}\"\n##variance table:\"{7}\"\n##variance restreict:\"{8}\"\n"
                     "##data base:\"{9}\""
                     "".format(phenotype,
                               len(control_id_list),
                               len(case_id_list),
                               sample_restrict,
                               sample_table_name,
                               gene_table_name,
                               gene_restrict,
                               variance_table,
                               variance_restrict,
                               db_file))
            fp.write("""
##                              case     control
##                           ---------------------
##                           |         |         |
##          subjects have alt|    A1   |    B1   |
##                           |         |         |
##                           ---------------------
##                           |         |         |
##  subjects do not have alt |    C1   |    D1   |
##                           |         |         |
##                           ---------------------
""")
            fp.write("#gene_id\tgene_name\tA1\tB1\tC1\tD1\tp_value1\todds_ratio1\n")
            ret_str = ""
            icounter = 0
            gene_data_len = len(gene_data)
            for gene_id, gene_name, gene_chr, gene_start, gene_end, gene_chr2, gene_start2, gene_end2 in gene_data:
                if icounter % 100 == 0 and icounter > 0:
                    logging.debug("handled {0} / {1} genes".format(icounter, gene_data_len))
                variance_selected_list = []
                if gene_chr is not None:
                    region = [chr_pos2absolute_pos(str(gene_chr), gene_start, chr2offset_dict),
                              chr_pos2absolute_pos(str(gene_chr), gene_end, chr2offset_dict)]
                    variance_selected_list = variance_in_region(variance_list, variance_list_index, region,
                                                                len(variance_list[0]))
                if gene_chr2 is not None:
                    region2 = [chr_pos2absolute_pos(str(gene_chr2), gene_start2, chr2offset_dict),
                               chr_pos2absolute_pos(str(gene_chr2), gene_end2, chr2offset_dict)]
                    variance_selected_list.extend(variance_in_region(variance_list, variance_list_index, region2,
                                                                     len(variance_list[0])))
                control_data = [i[3:3 + len(control_id_list)] for i in variance_selected_list]
                case_data = [i[3 + len(control_id_list):3 + len(control_id_list) + len(case_id_list)] for i in
                             variance_selected_list]
                control_people_data = [sum(filter(lambda x: type(x) == int, people_line)) for people_line in
                                       zip(*control_data)]
                B1 = len(filter(lambda x: x > 0, control_people_data))
                D1 = len(control_people_data) - B1
                case_people_data = [sum(filter(lambda x: type(x) == int, people_line)) for people_line in
                                    zip(*case_data)]
                A1 = len(filter(lambda x: x > 0, case_people_data))
                C1 = len(case_people_data) - A1
                oddsratio1, pvalue1 = stats.fisher_exact([[A1, B1], [C1, D1]])
                ret_str = "{0}{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n" \
                          "".format(ret_str, gene_id, gene_name, A1, B1, C1, D1,
                                    pvalue1, oddsratio1)
                icounter += 1

            fp.write(ret_str)
    logging.debug("all done")


def qqplot_fisher_permutation(path, fai_in, output):
    def parse_para(file_in):
        tmp_str_list = []
        pvalue_index = 6
        with open(file_in, "r") as fp:
            while True:
                data_line = fp.readline()
                if not data_line:
                    break
                if data_line.startswith("##"):
                    tmp_str_list.append(data_line.strip())
                    continue
                if data_line.startswith("#"):
                    data_list = data_line.split("\t")
                    pvalue_data = filter(lambda x: "p_value" in x, data_list)
                    if len(pvalue_data) == 1:
                        pvalue_index = data_list.index(pvalue_data[0])
                    continue
                break
        phenotype = None
        sample_restrict = None
        sample_table = None
        gene_table = None
        gene_restrict = None
        variance_table = None
        variance_restrict = None
        db_file = None
        for data_line in tmp_str_list:
            print("[{}]".format(data_line))
            data_list = data_line.strip().split(":")
            if len(data_list) != 2:
                continue
            key, value = data_list
            if value.startswith("\""):
                value = value[1:]
            if value.endswith("\""):
                value = value[:-1]
            if "phenotype" in key:
                phenotype = value
            if "sample restrict" in key:
                sample_restrict = value
            if "sample table" in key:
                sample_table = value
            if "gene table" in key:
                gene_table = value
            if "gene restrict" in key:
                gene_restrict = value
            if "variance table" in key:
                variance_table = value
            if "variance restreict" in key:
                variance_restrict = value
            if "data base" in key:
                db_file = value
        if None in [phenotype, sample_table, sample_restrict, gene_table, gene_restrict,
                    variance_table, variance_restrict, db_file, pvalue_index]:
            RuntimeError("illegle input file")
        return [phenotype, sample_table, sample_restrict, gene_table, gene_restrict,
                variance_table, variance_restrict, db_file, pvalue_index]

    def load_pvalue(file_in, pvalue_index):
        tmp_p_list = []
        # logging.debug("file_in={0}".format(file_in))
        # logging.debug("pvalue_index={0}".format(pvalue_index))
        with open(file_in, "r") as fp:
            while True:
                data_line = fp.readline()
                if not data_line:
                    break
                if data_line.startswith("#"):
                    continue
                if not data_line.strip():
                    continue
                data_list = data_line.strip().split("\t")
                # logging.debug(data_list)
                if len(data_list) < pvalue_index + 1:
                    logging.debug(data_list)
                p_value = float(data_list[pvalue_index])
                tmp_p_list.append(p_value)
        tmp_p_list.sort()
        return tmp_p_list

    def get_fisher_p_list(db_file, phenotype,
                          sample_table_name, sample_restrict,
                          gene_table_name, gene_restrict,
                          variance_table, variance_restrict,
                          fai_in):
        chr2offset_dict = parse_fai(fai_in)
        conn = sqlite3.connect(db_file)
        cursor = conn.cursor()
        format_sample_restrict = "" if not sample_restrict else " AND {}".format(sample_restrict)
        cursor.execute("SELECT s.gen_id FROM {2} AS s WHERE s.{0}='0'{1}".format(phenotype,
                                                                                 format_sample_restrict,
                                                                                 sample_table_name))
        # control ids
        control_id_list = [i[0] for i in cursor.fetchall()]
        print("control number={}".format(len(control_id_list)))
        cursor.execute("SELECT s.gen_id FROM {2} AS s WHERE s.{0}='1'{1}".format(phenotype,
                                                                                 format_sample_restrict,
                                                                                 sample_table_name))
        # case ids
        case_id_list = [i[0] for i in cursor.fetchall()]
        print("case number={}".format(len(case_id_list)))

        # prepare gene data
        format_gene_restrict = "" if not gene_restrict else " WHERE {}".format(gene_restrict)
        cmd_str = "SELECT g.gene_id, g.gene_name, g.chr, g.start_pos, g.end_pos, g.chr2, g.start_pos2, g.end_pos2 " \
                  "FROM {0} AS g{1}".format(gene_table_name, format_gene_restrict)
        cursor.execute(cmd_str)
        gene_data = cursor.fetchall()

        # prepare the variance data
        cmd_str = "SELECT id, chr, pos, "
        for control_id in control_id_list:
            cmd_str = "{0}sample_{1}, ".format(cmd_str, control_id)
        for case_id in case_id_list:
            cmd_str = "{0}sample_{1}, ".format(cmd_str, case_id)
        cmd_str = cmd_str.strip(", ")
        format_variance_restrict = "" if not variance_restrict else " WHERE {}".format(variance_restrict)
        cmd_str = "{0} FROM {1}{2}".format(cmd_str, variance_table, format_variance_restrict)
        # logging.debug("getting the variance data from db under restrict [{}]".format(variance_restrict))
        cursor.execute(cmd_str)
        variance_list = cursor.fetchall()
        # logging.debug("len(org_variance_list)={0}".format(len(org_variance_list)))
        # Add absolute position, create sort, index
        for i in xrange(len(variance_list)):
            variance_list[i] = list(variance_list[i])
            variance_list[i].append(chr_pos2absolute_pos(str(variance_list[i][1]),
                                                         variance_list[i][2],
                                                         chr2offset_dict))
        variance_list.sort(key=sort_key_last)
        variance_list_index = build_data_ram_index(variance_list, len(variance_list[0]))
        ret_list = []  # list[list[gene_name, pvalue]]
        for gene_id, gene_name, gene_chr, gene_start, gene_end, gene_chr2, gene_start2, gene_end2 in gene_data:
            variance_selected_list = []
            if gene_chr is not None:
                region = [chr_pos2absolute_pos(str(gene_chr), gene_start, chr2offset_dict),
                          chr_pos2absolute_pos(str(gene_chr), gene_end, chr2offset_dict)]
                variance_selected_list = variance_in_region(variance_list, variance_list_index, region,
                                                            len(variance_list[0]))
            if gene_chr2 is not None:
                region2 = [chr_pos2absolute_pos(str(gene_chr2), gene_start2, chr2offset_dict),
                           chr_pos2absolute_pos(str(gene_chr2), gene_end2, chr2offset_dict)]
                variance_selected_list.extend(variance_in_region(variance_list, variance_list_index, region2,
                                                                 len(variance_list[0])))
            control_data = [i[3:3 + len(control_id_list)] for i in variance_selected_list]
            case_data = [i[3 + len(control_id_list):3 + len(control_id_list) + len(case_id_list)] for i in
                         variance_selected_list]
            control_people_data = [sum(filter(lambda x: type(x) == int, people_line)) for people_line in
                                   zip(*control_data)]
            B1 = len(filter(lambda x: x > 0, control_people_data))
            D1 = len(control_people_data) - B1
            case_people_data = [sum(filter(lambda x: type(x) == int, people_line)) for people_line in
                                zip(*case_data)]
            A1 = len(filter(lambda x: x > 0, case_people_data))
            C1 = len(case_people_data) - A1
            oddsratio1, pvalue1 = stats.fisher_exact([[A1, B1], [C1, D1]])
            ret_list.append([gene_name, pvalue1])
        ret_list.sort(key=sort_key_last_float)
        return ret_list

    print("begin")
    logging.basicConfig(filename="{}.log".format(sys._getframe().f_code.co_name),
                        level=logging.DEBUG, format=log_format, filemode="w")

    main_path = os.path.split(os.path.realpath(__file__))[0]
    # ph = Popen(["mkdir -p {0} && rm -rf {1}".format(path, os.path.join(path, "*"))], shell=True)
    # ph.wait()
    # print("permutating...")
    # permutation_fisher(db_file, phenotype,
    #                    sample_table_name, sample_restrict,
    #                    gene_table_name, gene_restrict,
    #                    variance_table, variance_restrict,
    #                    fai_in, path, permute_times)

    cmd_str = "ls {}".format(os.path.join(path, "*"))
    print("loading file names. cmd={}".format(cmd_str))
    ph = Popen([cmd_str], shell=True, stdout=PIPE)
    # print(1)
    # ph.wait()
    # print(2)
    name_list = [i.strip() for i in ph.stdout.readlines()]
    # print(name_list[0])
    name_list = filter(lambda x: "fisher_permutation" in x, name_list)
    print("name number={}".format(len(name_list)))
    print("parsing para")
    assert len(name_list) > 0
    [phenotype, sample_table, sample_restrict, gene_table, gene_restrict,
     variance_table, variance_restrict, db_file, pvalue_index] = parse_para(name_list[0])

    print("loading p_value")
    # load p_value
    p_matrix = []
    for name in name_list:
        p_matrix.append(load_pvalue(name, pvalue_index))
    p_matrix_t = zip(*p_matrix)

    # list[list[gene_name, pvalue]]
    print("calculating the ob p_value")
    ret_data = get_fisher_p_list(db_file, phenotype, sample_table, sample_restrict, gene_table,
                                 gene_restrict, variance_table, variance_restrict, fai_in)
    # ret_data.sort(key=sort_key_last)
    print(ret_data[:5])
    print("len(ret_data)={}".format(len(ret_data)))
    print("len(p_matrix_t)={}".format(len(p_matrix_t)))
    assert len(ret_data) == len(p_matrix_t)
    for index in xrange(len(ret_data)):
        ret_data[index].extend(p_matrix_t[index])
    with open(os.path.join(main_path, "tmp"), "w") as fp:
        for gene_line in ret_data:
            fp.write("\t".join([str(i) for i in gene_line]) + "\n")

    cmd_str = "Rscript {0} {1} {2} {3}".format(os.path.join(main_path, "qqplot_permutation.R"),
                                               os.path.join(main_path, "tmp"),
                                               output,
                                               len(p_matrix))
    print(cmd_str)
    ph = Popen([cmd_str], shell=True)
    ph.wait()
    # os.remove(os.path.join(main_path, "tmp"))


def modify_gene_name2id(db_file, file_in, gene_table_name, hg19_gene_info, output):
    log_format2 = "%(message)s"
    logging.basicConfig(filename="{}.log".format(sys._getframe().f_code.co_name),
                        level=logging.DEBUG, format=log_format2, filemode="w")
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()
    cursor.execute("SELECT gene_name, gene_id FROM {0} ".format(gene_table_name))

    gene_dict_coding = {}
    id_set = set([])
    gene_set_all = set([])
    gene_set_found = set([])
    for name, id in cursor.fetchall():
        if id is None or name is None:
            continue
        name = name.upper().split(".")[0]
        id = id.upper().split(".")[0]

        id_set.add(id)
        if name not in gene_dict_coding:
            gene_dict_coding[name] = [id]
        else:
            gene_dict_coding[name].append(id)
    old_coding_gene_dict = {}
    old_nocoding_gene_set = set([])
    icounter = 0
    with open(hg19_gene_info, "r") as fp:
        while True:
            data_line = fp.readline()
            icounter += 1
            if not data_line:
                print(data_list)
                break
            data_list = [i.strip("\"") for i in re.split(" |\t", data_line.strip()) if len(i) > 0]
            name = data_list[-1].strip().upper().split(".")[0]
            id = data_list[1].strip().upper().split(".")[0]
            gene_type = data_list[3].strip().upper()
            if gene_type == "PROTEIN_CODING":
                if name not in old_coding_gene_dict:
                    old_coding_gene_dict[name] = [id]
                else:
                    old_coding_gene_dict[name].append(id)
            else:
                old_nocoding_gene_set.add(name)
            if icounter == 1:
                print(data_list)

    # itotal = 0
    # icase1 = 0
    # icase2 = 0
    # icase3 = 0
    # normal_case = 0
    with open(file_in, "r") as fp_in, open(output, "w") as fp_out:
        while True:
            data_line = fp_in.readline()
            if not data_line:
                break
            # data_list = filter(lambda x: len(x) > 0, data_line.strip().split("\t"))
            data_list = [i.strip() for i in data_line.strip().split("\t") if len(i.strip()) > 0]
            out_list = [data_list[0]]
            for ele in data_list[1:]:
                org_ele = ele.split(".")[0]
                ele = ele.upper().split(".")[0]
                # itotal += 1
                if ele in id_set:
                    out_list.append(ele)
                    # icase3 += 1
                    continue

                # ele not in id set
                gene_set_all.add(ele)
                if ele in gene_dict_coding:
                    out_list.extend(gene_dict_coding[ele])
                    gene_set_found.add(ele)
                    continue

                # ele not in gene_dict_coding
                if ele in old_coding_gene_dict:
                    out_list.extend(old_coding_gene_dict[ele])
                    gene_set_found.add(ele)
                    continue

                # ele not in old coding gene dict
                if ele in old_nocoding_gene_set:
                    gene_set_found.add(ele)
                    continue

                # ele not in old nocoding gene set
                out_list.append(ele)
                logging.debug("[{0}] can not find the corresponding id".format(org_ele))

            fp_out.write("\t".join(out_list) + "\n")
    # print("total number of searches = {}".format(itotal))
    # print("gene name total number of searches = {}".format(icase1))
    # print("gene name id not found = {}".format(icase2))
    # print("The number of times the id appears directly = {}".format(icase3))
    # print("gene name An id was found = {}".format(normal_case))
    print("A total of {0} occurred gene name  There are {1} ids found  There are still {2} 
    ids not found"
          "".format(len(gene_set_all), len(gene_set_found),
                    len(gene_set_all) - len(gene_set_found)))
    # print("len(gene_set_all) = {}".format(len(gene_set_all)))
    # print("len(gene_set_found) = {}".format(len(gene_set_found)))


def modify_vcfid(list_in, vcf_in, output):
    print("loading {0} ...".format(list_in))
    sys.stdout.flush()
    key2id_dict = {}
    with open(list_in, "r") as fp:
        while True:
            line_data = fp.readline()
            if not line_data:
                break
            if "-" in line_data:
                continue
            line_list = line_data.strip().split("\t")
            key2id_dict["{0}\t{1}\t{2}\t{3}" \
                        "".format(line_list[0], line_list[1],
                                  line_list[3], line_list[4])] = line_list[5]

    print("handling {0} ...".format(vcf_in))
    sys.stdout.flush()
    icounter = 0
    with open(vcf_in, "r") as fp_in, open(output, "w") as fp_out:
        while True:
            line_data = fp_in.readline()
            icounter += 1
            if not line_data:
                break
            if line_data.startswith("#"):
                fp_out.write(line_data)
                continue
            if len(line_data.strip()) == 0:
                continue
            line_list = line_data.split("\t")
            key = "{0}\t{1}\t{2}\t{3}" \
                  "".format(line_list[0][3:], line_list[1],
                            line_list[3], line_list[4])
            if key not in key2id_dict:
                fp_out.write(line_data)
            else:
                line_list[2] = key2id_dict[key]
                fp_out.write("\t".join(line_list))
            if icounter % 1000:
                print("handled {0} vcf lines".format(icounter))


def add_gene_name(file_in, db_file, gene_table_name, id_col, output):
    logging.basicConfig(filename="{}.log".format(sys._getframe().f_code.co_name),
                        level=logging.DEBUG, format=log_format, filemode="a")
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()
    cursor.execute("SELECT gene_name, gene_id FROM {0} ".format(gene_table_name))
    id2name_dict = {}
    for name, id in cursor.fetchall():
        name = name.upper().split(".")[0]
        id = id.upper().split(".")[0]

        if id not in id2name_dict:
            id2name_dict[id] = [name]
        else:
            id2name_dict[id].append(name)
            logging.warning("id[{0}] has more than 1 name".format(id))
    # print("len(id2name_dict)=[{0}]".format(len(id2name_dict)))
    with open(file_in, "r") as fp, open(output, "w") as fp_out:
        data_head_list = filter(lambda x: len(x) > 0, fp.readline().strip().split(" "))
        data_head_list.append("GENE_NAME")
        fp_out.write("\t".join(data_head_list) + "\n")
        while True:
            data_line = fp.readline()
            if not data_line:
                break
            data_list = filter(lambda x: len(x) > 0, data_line.strip().split(" "))

            id = data_list[int(id_col) - 1].upper().split(".")[0]
            if id not in id2name_dict:
                name = "NA"
                logging.warning("id[{0}] did not find name".format(id))
            else:
                name = id2name_dict[id]
            data_list.extend(name)
            # logging.debug("data_list={0}".format(data_list))
            fp_out.write("\t".join(data_list) + "\n")


def get_id_dup_num(id_str_in):
    if id_str_in.startswith("chr"):
        id_list = id_str_in.split("_")
        if len(id_list) == 2:
            return 0
        return int(id_list[-1])
    else:
        id_list = id_str_in.split("_")
        if len(id_list) == 1:
            return 0
        return int(id_list[-1])


def currID_is_same_with_lastID(curr_id, last_id):
    if last_id.startswith("chr"):
        last_id = "_".join(last_id.split("_")[:2])
        if curr_id == last_id:
            return True
        return False
    spliter_num = last_id.count("_")
    if spliter_num > 0:
        last_id = "_".join(last_id.split("_")[:-1])
    if curr_id == last_id:
        return True
    return False


def modify_duplicateID(vcf_in, output):
    last_id = ""
    with open(vcf_in, "r") as fp_in, open(output, "w") as fp_out:
        while True:
            data_line = fp_in.readline()
            if not data_line:
                break

            # comment ling
            if data_line.startswith("#"):
                fp_out.write(data_line)
                continue
            data_list = data_line.split("\t")

            # the first data line
            if not last_id:
                fp_out.write(data_line)
                last_id = data_list[2]
                continue

            # not the first line. It is a new line
            curr_id = data_list[2]
            if not currID_is_same_with_lastID(curr_id, last_id):
                fp_out.write(data_line)
                last_id = data_list[2]
                continue
            # A line with the same id
            dup_num = get_id_dup_num(last_id)
            data_list[2] = "{0}_{1}".format(data_list[2], dup_num + 1)
            fp_out.write("\t".join(data_list))
            last_id = data_list[2]


def manhattan_direct(file_in, name_col, chr_col, start_col, end_col, p_col, fai_in, output):
    import matplotlib.pyplot as plt
    # color_dict = {"1": (1, 0, 0), "2": (1, 0.2, 0), "3": (1, 0.4, 0), "4": (1, 0.6, 0), "5": (1, 0.8, 0),
    #               "6": (1, 1, 0),
    #               "7": (0.833, 1, 0), "8": (0.66, 1, 0), "9": (0.5, 1, 0), "10": (0.33, 1, 0), "11": (0.166, 1, 0),
    #               "12": (0, 1, 0),
    #               "13": (0, 0.833, 0.166), "14": (0, 0.666, 0.333), "15": (0, 0.5, 0.5), "16": (0, 0.333, 0.666),
    #               "17": (0, 0.166, 0.833),
    #               "18": (0, 0, 1),
    #               "19": (0.167, 0, 1), "20": (0.33, 0, 1), "21": (0.5, 0, 1), "22": (0.67, 0, 1), "X": (0.833, 0.69, 1),
    #               "Y": (1, 0, 1)}
    color_dict = {"1": (0, 0, 0), "2": (0.8, 0.8, 0.8), "3": (0, 0, 0), "4": (0.8, 0.8, 0.8), "5": (0, 0, 0),
                  "6": (0.8, 0.8, 0.8),
                  "7": (0, 0, 0), "8": (0.8, 0.8, 0.8), "9": (0, 0, 0), "10": (0.8, 0.8, 0.8), "11": (0, 0, 0),
                  "12": (0.8, 0.8, 0.8),
                  "13": (0, 0, 0), "14": (0.8, 0.8, 0.8), "15": (0, 0, 0), "16": (0.8, 0.8, 0.8),
                  "17": (0, 0, 0),
                  "18": (0.8, 0.8, 0.8),
                  "19": (0, 0, 0), "20": (0.8, 0.8, 0.8), "21": (0, 0, 0), "22": (0.8, 0.8, 0.8), "X": (0, 0, 0),
                  "Y": (0.8, 0.8, 0.8)}
    fig = plt.figure(figsize=(32, 16))
    ax1 = fig.add_subplot(1, 1, 1)

    chr2offset_dict = parse_fai(fai_in)
    icounter = 0
    font1 = {'family': 'Arial',
             'weight': 'normal',
             'size': 40,
             }
    with open(file_in, "r") as fp:
        while True:
            data_line = fp.readline()
            icounter += 1
            if not data_line:
                break
            if len(data_line.strip()) == 0:
                continue
            data_line = filter(lambda x: len(x) > 0, data_line.strip().split(" "))
            name = data_line[int(name_col) - 1]
            chr = data_line[int(chr_col) - 1]
            start = data_line[int(start_col) - 1]
            end = data_line[int(end_col) - 1]
            p_value = float(data_line[int(p_col) - 1])
            abp_start = chr_pos2absolute_pos(chr, start, chr2offset_dict)
            abp_end = chr_pos2absolute_pos(chr, end, chr2offset_dict)
            ax1.plot([abp_start, abp_end],  # to make it clear set +-5000 as enhancer
                     [-math.log10(p_value), -math.log10(p_value)],
                     color=color_dict[chr], linewidth=5)
            if icounter % 100 == 0:
                print("handled {0} lines".format(icounter))
    ax1.set_ylabel("-log(P)", font1)
    ax1.tick_params(labelsize=35)
    ax1.set_xticks([])
    plt.savefig(output + ".pdf", format="pdf")
    plt.savefig(output + ".png", format="png")
    plt.close(fig)


def get_gene_set_info(db_file, gene_list_file, output, fai_in, phenos, variance_table_name):
    def variance_category(source_list):
        lof_num = sum([int(i) for i in source_list[:3]])
        dmis_num = int(source_list[3])
        dsplicy_num = sum([int(i) for i in source_list[4:7]])
        ret_list = []
        if lof_num > 0:
            ret_list.append("lof")
        if dmis_num > 0:
            ret_list.append("dmis")
        if dsplicy_num > 0:
            ret_list.append("dsplicy")
        return ",".join(ret_list)

    def get_case_control_list(db_cursor, pheno):
        cmd_str = "SELECT gen_id, {0} from sampleChdPhenotype where {0} in ('0','1')".format(pheno)
        db_cursor.execute(cmd_str)
        pheno_data = db_cursor.fetchall()
        case_list = [gen_id for gen_id, pheno_num in pheno_data if pheno_num == 1]
        control_list = [gen_id for gen_id, pheno_num in pheno_data if pheno_num == 0]
        return [case_list, control_list]

    logging.basicConfig(filename="{}.log".format(sys._getframe().f_code.co_name),
                        level=logging.DEBUG, format=log_format, filemode="w")
    gene_table_name = "gene_table"
    variance_restrict = ""
    # variance_table_name = "variance"
    pheno_list = phenos.strip().split(",")

    chr2offset_dict = parse_fai(fai_in)
    with open(gene_list_file, "r") as fp:
        gene_list = [i.strip().split("\t")[0] for i in fp.readlines() if i.strip()]
    gene_list_str = "(\"" + "\",\"".join(gene_list) + "\")"
    gene_restrict = "gene_id in {0} or gene_name in {0}".format(gene_list_str)
    logging.debug("begin")
    logging.debug("db_file=[{}]".format(db_file))
    logging.debug("variance_table_name=[{}]".format(variance_table_name))
    logging.debug("gene_table_name=[{}]".format(gene_table_name))
    logging.debug("variance_restrict=[{}]".format(variance_restrict))
    logging.debug("gene_restrict=[{}]".format(gene_restrict))
    logging.debug("output=[{}]".format(output))
    # if variance_table_name not in ['variance', 'synonymous_snp']:
    #     logging.debug("illegal table name [{0}]".format(variance_table_name))
    #     return
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()

    format_variance_restrict = "" if not variance_restrict else " WHERE {0}".format(variance_restrict)
    format_gene_restrict = "" if not gene_restrict else " WHERE {0}".format(gene_restrict)
    cmd_str = "SELECT g.gene_id, g.gene_name, g.chr, g.start_pos, g.end_pos, g.chr2, g.start_pos2, g.end_pos2 " \
              "FROM {0} AS g{1}".format(gene_table_name, format_gene_restrict)
    cursor.execute(cmd_str)
    gene_data = cursor.fetchall()
    vcf_info_num = 6
    cmd_str = "SELECT v.vcf_id, v.id, v.chr, v.pos, v.ref, v.alt"
    if db_has_col(cursor, variance_table_name, "bystro_cadd"):
        cmd_str += ",v.bystro_cadd"
        vcf_info_num += 1
    if db_has_col(cursor, variance_table_name, "bystro_phyloP"):
        cmd_str += ",v.bystro_phyloP"
        vcf_info_num += 1
    if db_has_col(cursor, variance_table_name, "bystro_phastCons"):
        cmd_str += ",v.bystro_phastCons"
        vcf_info_num += 1
    if db_has_col(cursor, variance_table_name, "annovar"):
        cmd_str += ",v.annovar"
    if db_has_col(cursor, variance_table_name, "bystro"):
        cmd_str += ",v.bystro"
    if db_has_col(cursor, variance_table_name, "vep"):
        cmd_str += ",v.vep"
    if db_has_col(cursor, variance_table_name, "dmis"):
        cmd_str += ",v.dmis"
    if db_has_col(cursor, variance_table_name, "dsplicing"):
        cmd_str += ",v.dsplicing"
    if db_has_col(cursor, variance_table_name, "spidex"):
        cmd_str += ",v.spidex"
    if db_has_col(cursor, variance_table_name, "spliceAI"):
        cmd_str += ",v.spliceAI"
    cmd_str += " FROM {0} AS v{1}".format(variance_table_name, format_variance_restrict)
    cursor.execute(cmd_str)
    variance_data = cursor.fetchall()
    for i in xrange(len(variance_data)):
        variance_data[i] = list(variance_data[i])
        variance_data[i].append(chr_pos2absolute_pos(str(variance_data[i][2]), variance_data[i][3], chr2offset_dict))
    variance_data.sort(key=sort_key_last)
    variance_element_len = len(variance_data[0])
    variance_data_index = build_data_ram_index(variance_data, variance_element_len)

    logging.debug("variance_data[0] = {}".format(variance_data[0]))
    logging.debug("variance_data[-1] = {}".format(variance_data[-1]))
    logging.debug("len(variance_data) = {}".format(len(variance_data)))
    # exit(0)

    head_list = ["#gene_id", "gene_name", "vcf_id", "variance_table_id", "chr", "pos", "ref", "alt"]
    if "bystro_cadd" in cmd_str:
        head_list.append("cadd")
    if "bystro_phyloP" in cmd_str:
        head_list.append("phyloP")
    if "bystro_phastCons" in cmd_str:
        head_list.append("phastCons")
    head_list.append("category")
    pheno_case_control_list = []
    for pheno in pheno_list:
        case_list, control_list = get_case_control_list(cursor, pheno)
        head_list.append("number_of_sample_take_variance_in_{0}_case_group".format(pheno))
        head_list.append("number_of_sample_take_variance_in_{0}_control_group".format(pheno))
        pheno_case_control_list.append([case_list, control_list])

    ret_list = []
    for gene_id, gene_name, gene_chr, gene_start_pos, gene_end_pos, gene_chr2, gene_start_pos2, gene_end_pos2 in gene_data:
        selected_variance = []
        if gene_chr is not None:
            gene_region = [chr_pos2absolute_pos(str(gene_chr), gene_start_pos, chr2offset_dict),
                           chr_pos2absolute_pos(str(gene_chr), gene_end_pos, chr2offset_dict)]
            selected_variance = variance_in_region(variance_data, variance_data_index, gene_region,
                                                   variance_element_len)
        if gene_chr2 is not None:
            gene_region = [chr_pos2absolute_pos(str(gene_chr2), gene_start_pos2, chr2offset_dict),
                           chr_pos2absolute_pos(str(gene_chr2), gene_end_pos2, chr2offset_dict)]
            selected_variance.extend(
                variance_in_region(variance_data, variance_data_index, gene_region, variance_element_len))
        for curr_variance in selected_variance:
            tmp_list = [gene_id, gene_name]
            tmp_list.extend([str(i) for i in curr_variance[:vcf_info_num]])
            if variance_table_name == "variance":
                tmp_list.append(variance_category(curr_variance[9:16]))
            elif variance_table_name == "synonymous_snp":
                tmp_list.append("lof")
            else:
                tmp_list.append("na")
            for i in xrange(len(pheno_list)):
                case_list, control_list = pheno_case_control_list[i]
                cursor.execute("select {0},{1} from {3} where id={2}"
                               "".format(",".join(["sample_{0}".format(x) for x in case_list]),
                                         ",".join(["sample_{0}".format(x) for x in control_list]),
                                         tmp_list[3],
                                         variance_table_name))
                curr_v_p_data = cursor.fetchall()[0]
                v_case_data = filter(lambda x: str(x).isdigit(), curr_v_p_data[:len(case_list)])
                v_control_data = filter(lambda x: str(x).isdigit(), curr_v_p_data[len(case_list):])
                v_case_num = len(filter(lambda x: int(x) > 0, v_case_data))
                v_control_num = len(filter(lambda x: int(x) > 0, v_control_data))
                tmp_list.append("{0}/{1}".format(v_case_num, len(v_case_data)))
                tmp_list.append("{0}/{1}".format(v_control_num, len(v_control_data)))
            ret_list.append("\t".join(tmp_list))

    with open(output, "w") as fp:
        fp.write("##db:{0}\n##variance table:\"{1}\"\n##variance_restrict:\"{2}\"\n"
                 "##gene table:\"{3}\"\n##gene restrict:\"{4}\"\n".format(db_file,
                                                                          variance_table_name,
                                                                          variance_restrict.replace("v.", ""),
                                                                          gene_table_name,
                                                                          gene_restrict))
        fp.write("\t".join(head_list) + "\n")
        fp.write("\n".join(ret_list))
    logging.debug("all done")


def copy_vcf2table(vcf_in, table_name, db_file):
    """
    the vcf should be splited and left normalized first
    @param vcf_in:
    @param table_name:
    @param db_file:
    @return:
    """
    logging.basicConfig(filename="{}.log".format(sys._getframe().f_code.co_name),
                        level=logging.DEBUG, format=log_format, filemode="w")
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()
    cursor.execute("DROP TABLE IF EXISTS {0}".format(table_name))
    cmd_str = 'create table {0} (id int primary key, ' \
              'chr varchr(20), ' \
              'pos varchr(20),' \
              'vcf_id varchr(20),' \
              'ref varchr(20),' \
              'alt varchr(20),' \
              'vcf_qual varchr(20),' \
              'vcf_filter varchr(20),' \
              'vcf_info varchr(20),' \
              'vcf_format varchr(20)'.format(
        table_name)
    sql_list = []
    icounter = 1
    tmp_str = ""
    with open(vcf_in, "r") as fp:
        while True:
            vcf_line = fp.readline()
            if not vcf_line:
                break
            if vcf_line.startswith("##"):
                continue
            if vcf_line.startswith("#"):
                vcf_list = vcf_line.strip().split("\t")
                for sample_num in vcf_list[9:]:
                    sample_num = sample_num.split("_")[0]
                    tmp_str += ",sample_{0}".format(sample_num)
                cmd_str += "{0})".format(tmp_str)
                cursor.execute(cmd_str)
                continue
            left_vcf_list = line_left_normalization(vcf_line).strip().split("\t")
            for i in xrange(len(left_vcf_list)):
                if i < 9:
                    continue
                if "." in left_vcf_list[i]:
                    left_vcf_list[i] = "na"
                    continue
                left_vcf_list[i] = str(sum([int(j) for j in re.split("/|\|", left_vcf_list[i].split(":")[0])]))
            value_str = "','".join(left_vcf_list)
            sql_list.append("insert into {0} (id, chr, pos, vcf_id, ref, alt, vcf_qual, vcf_filter, vcf_info, "
                            "vcf_format{1}) VALUES('{2}','{3}')".format(table_name, tmp_str, icounter, value_str))
            icounter += 1
            if icounter % 100 == 0:
                logging.debug("handled {0} vcf lines".format(icounter))
    icounter = 0
    total_num = len(sql_list)
    for sql in sql_list:
        cursor.execute(sql)
        icounter += 1
        if icounter % 100 == 0:
            logging.debug("excuted {0} / {1} sqls".format(icounter, total_num))
    cursor.close()
    conn.commit()
    conn.close()
    logging.debug("all done")


def db_add_ccds(db_file, ccds_bed, fai_in):
    def abp_is_in_ccds_regions(sorted_ccds_list, target_abp):
        def abp_is_in_region(region, target_abp):
            """
            Judging the absolute position and the inclusion of the area
            @param region: [start abp, end abp]
            @param target_abp:
            @return: 1 target大于region， -1 target小于region， 0 target在region中
            """
            region_start, region_end = region
            if target_abp < region_start:
                return -1
            if target_abp > region_end:
                return 1
            return 0

        for i in xrange(0, len(sorted_ccds_list), 1000):
            if abp_is_in_region(sorted_ccds_list[i], target_abp) > 0:
                continue
            elif abp_is_in_region(sorted_ccds_list[i], target_abp) == 0:
                return True
            else:  # target < sorted_ccds_list[i]
                for j in xrange(i - 1000, i, 1):
                    if abp_is_in_region(sorted_ccds_list[j], target_abp) > 0:
                        continue
                    elif abp_is_in_region(sorted_ccds_list[j], target_abp) == 0:
                        return True
                    else:
                        return False
                else:
                    return False
        else:
            for j in xrange(i, len(sorted_ccds_list), 1):
                if abp_is_in_region(sorted_ccds_list[j], target_abp) > 0:
                    continue
                elif abp_is_in_region(sorted_ccds_list[j], target_abp) == 0:
                    return True
                else:
                    return False
            else:
                return False

    print("add_ccds begin")
    chr2offset_dict = parse_fai(fai_in)
    print("loading ccds data...")
    with open(ccds_bed, "r") as fp:
        # ccds_list = list([
        #     ([chr_pos2absolute_pos(chrom, start, chr2offset_dict), chr_pos2absolute_pos(chrom, end, chr2offset_dict)]
        #      for chrom, start, end in
        #      i.strip().split("\t")) for i in fp.readlines()])
        ccds_data = [i.strip().split("\t") for i in fp.readlines()]
        ccds_list = [
            [chr_pos2absolute_pos(chrom[3:], start, chr2offset_dict),
             chr_pos2absolute_pos(chrom[3:], end, chr2offset_dict)] for
            chrom, start, end in ccds_data]
    # print(ccds_list[0])
    # print(ccds_list[-1])
    # exit(0)
    ccds_list.sort(key=sort_key_last)
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()
    db_add_col(cursor, "variance", "is_ccds", "varchr(1)")
    print("loading variance info from db...")
    cursor.execute("SELECT id, chr, pos FROM variance")
    variance_data = cursor.fetchall()
    sql_list = []
    for id, chrom, pos in variance_data:
        curr_abp = chr_pos2absolute_pos(chrom, pos, chr2offset_dict)
        if not abp_is_in_ccds_regions(ccds_list, curr_abp):
            sql_list.append("UPDATE variance SET is_ccds = '0' WHERE id = {0}".format(id))
        else:
            sql_list.append("UPDATE variance SET is_ccds = '1' WHERE id = {0}".format(id))

    print("loading synonymous info from db...")
    db_add_col(cursor, "synonymous_snp", "is_ccds", "varchr(1)")
    cursor.execute("SELECT id, chr, pos FROM synonymous_snp")
    variance_data = cursor.fetchall()
    for id, chrom, pos in variance_data:
        curr_abp = chr_pos2absolute_pos(chrom, pos, chr2offset_dict)
        if not abp_is_in_ccds_regions(ccds_list, curr_abp):
            sql_list.append("UPDATE synonymous_snp SET is_ccds = '0' WHERE id = {0}".format(id))
        else:
            sql_list.append("UPDATE synonymous_snp SET is_ccds = '1' WHERE id = {0}".format(id))

    print("updating db...")
    sql_num = len(sql_list)
    icounter = 0
    for sql_str in sql_list:
        cursor.execute(sql_str)
        icounter += 1
        if icounter % 1000 == 0:
            print("updating db {0} / {1}\t{2:.2%}".format(icounter, sql_num, icounter / float(sql_num)))
    cursor.close()
    conn.commit()
    conn.close()
    print("all done")


def db_add_gene_number(db_file, phenotype, variance_table, fai_in):
    def sum_list(list_in):
        return sum(filter(lambda x: type(x) == int, list_in))

    print("db_add_gene_number begin")
    print(phenotype)
    print(variance_table)
    case_sample_num_name = "{0}_{1}_case_sample_num".format(variance_table, phenotype)
    case_variance_num_name = "{0}_{1}_case_variance_num".format(variance_table, phenotype)
    case_allele_num_name = "{0}_{1}_case_allele_num".format(variance_table, phenotype)
    control_sample_num_name = "{0}_shared_control_sample_num".format(variance_table)
    control_variance_num_name = "{0}_shared_control_variance_num".format(variance_table)
    control_allele_num_name = "{0}_shared_control_allele_num".format(variance_table)

    chr2offset_dict = parse_fai(fai_in)
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()
    has_control_data = db_has_col(cursor, "gene_table", control_sample_num_name)
    db_add_col(cursor, "gene_table", case_sample_num_name, "varchr(20)")
    db_add_col(cursor, "gene_table", case_variance_num_name, "varchr(20)")
    db_add_col(cursor, "gene_table", case_allele_num_name, "varchr(20)")
    db_add_col(cursor, "gene_table", control_sample_num_name, "varchr(20)")
    db_add_col(cursor, "gene_table", control_variance_num_name, "varchr(20)")
    db_add_col(cursor, "gene_table", control_allele_num_name, "varchr(20)")
    print("loading pheno data...")
    cursor.execute("select gen_id, {0} from sampleChdPhenotype".format(phenotype))
    pheno_data = cursor.fetchall()
    case_list = ["sample_{0}".format(gen_id) for gen_id, group in pheno_data if group == 1]
    control_list = ["sample_{0}".format(gen_id) for gen_id, group in pheno_data if group == 0]
    # print(case_list[0])
    # print(control_list[0])

    print("loading variance data...")
    cmd_str = "SELECT chr, pos, {0}, {1} FROM {2}" \
              "".format(",".join(case_list), ",".join(control_list), variance_table)
    cursor.execute(cmd_str)
    variance_data = cursor.fetchall()
    for i in xrange(len(variance_data)):
        variance_data[i] = list(variance_data[i])
        variance_data[i].append(chr_pos2absolute_pos(str(variance_data[i][0]), variance_data[i][1], chr2offset_dict))
    variance_data.sort(key=sort_key_last)
    variance_element_len = len(variance_data[0])
    variance_data_index = build_data_ram_index(variance_data, variance_element_len)
    print("loading gene data...")
    cmd_str = "SELECT id, gene_id, chr, start_pos, end_pos, chr2, start_pos2, end_pos2 FROM gene_table"
    cursor.execute(cmd_str)
    gene_data = cursor.fetchall()
    sql_list = []
    for id, gene_id, gene_chr, gene_start_pos, gene_end_pos, gene_chr2, gene_start_pos2, gene_end_pos2 in gene_data:
        selected_variance = []
        if gene_chr is not None:
            gene_region = [chr_pos2absolute_pos(str(gene_chr), gene_start_pos, chr2offset_dict),
                           chr_pos2absolute_pos(str(gene_chr), gene_end_pos, chr2offset_dict)]
            selected_variance = variance_in_region(variance_data, variance_data_index, gene_region,
                                                   variance_element_len)
        if gene_chr2 is not None:
            gene_region = [chr_pos2absolute_pos(str(gene_chr2), gene_start_pos2, chr2offset_dict),
                           chr_pos2absolute_pos(str(gene_chr2), gene_end_pos2, chr2offset_dict)]
            selected_variance.extend(
                variance_in_region(variance_data, variance_data_index, gene_region, variance_element_len))
        if len(selected_variance) == 0:
            continue
        case_data = [i[2:2 + len(case_list)] for i in selected_variance]

        case_data_t = zip(*case_data)
        case_sum_list = [sum_list(i) for i in case_data_t]

        # In the phenotype case group, the number of people with variance in a certain gene
        case_sample_num = len(filter(lambda x: x > 0, case_sum_list))

        # In the phenotype case group, the total number of alleles appearing in the 
        variance in a certain gene
        case_allele_num = sum(case_sum_list)

        # In the phenotype case group, the number of variances that appear in a gene
        case_variance_num = len(filter(lambda x: x > 0, [sum_list(i) for i in case_data]))
        sql_list.append("UPDATE {0} SET {1} = {2} WHERE id = {3}"
                        "".format("gene_table", case_sample_num_name, case_sample_num, id))
        sql_list.append("UPDATE {0} SET {1} = {2} WHERE id = {3}"
                        "".format("gene_table", case_allele_num_name, case_allele_num, id))
        sql_list.append("UPDATE {0} SET {1} = {2} WHERE id = {3}"
                        "".format("gene_table", case_variance_num_name, case_variance_num, id))
        if not has_control_data:
            control_data = [i[2 + len(case_list):2 + len(case_list) + len(control_list)] for i in selected_variance]
            control_data_t = zip(*control_data)
            control_sum_list = [sum_list(i) for i in control_data_t]
            control_sample_num = len(filter(lambda x: x > 0, control_sum_list))
            control_allele_num = sum(control_sum_list)
            control_variance_num = len(filter(lambda x: x > 0, [sum_list(i) for i in control_data]))
            sql_list.append("UPDATE {0} SET {1} = {2} WHERE id = {3}"
                            "".format("gene_table", control_sample_num_name, control_sample_num, id))
            sql_list.append("UPDATE {0} SET {1} = {2} WHERE id = {3}"
                            "".format("gene_table", control_allele_num_name, control_allele_num, id))
            sql_list.append("UPDATE {0} SET {1} = {2} WHERE id = {3}"
                            "".format("gene_table", control_variance_num_name, control_variance_num, id))
    print("updating db...")
    icounter = 0
    for sql_str in sql_list:
        cursor.execute(sql_str)
        icounter += 1
        if icounter % 1000 == 0:
            print("handled {0} sqls".format(icounter))
    cursor.close()
    conn.commit()
    conn.close()
    print("all done")


def firth_logistic_regression(db_file, sample_table, phenotype_in, sample_restrict,
                              variance_table, variance_restrict):
    print("begin")
    # chr2offset_dict = parse_fai(fai_in)
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()

    # handle sample
    if sample_restrict.strip():
        cmd_str = "select gen_id, {0},gender,PC1,PC2,PC3,PC4,PC5 from {1} as s where {2}".format(phenotype_in,
                                                                                                 sample_table,
                                                                                                 sample_restrict)
    else:
        cmd_str = "select gen_id, {0},gender,PC1,PC2,PC3,PC4,PC5 from {1}".format(phenotype_in, sample_table)
    print("handling sample data")
    cursor.execute(cmd_str)
    sample_data = cursor.fetchall()
    sample_data_t = zip(*sample_data)

    outcome_list = list([str(i) for i in sample_data_t[1]])
    sample_id_list = list([str(i) for i in sample_data_t[0]])
    gender_list = list([str(i) for i in sample_data_t[2]])
    pc1_list = list([str(i) for i in sample_data_t[3]])
    pc2_list = list([str(i) for i in sample_data_t[4]])
    pc3_list = list([str(i) for i in sample_data_t[5]])
    pc4_list = list([str(i) for i in sample_data_t[6]])
    pc5_list = list([str(i) for i in sample_data_t[7]])

    # handle variance data
    if variance_restrict.strip():
        cmd_str = "select sample_{0} from {1} as v where {2}" \
                  "".format(",sample_".join(sample_id_list), variance_table, variance_restrict)
    else:
        cmd_str = "select sample_{0} from {1}" \
                  "".format(",sample_".join(sample_id_list), variance_table)
    # print cmd_str
    print("handling variance_data")
    cursor.execute(cmd_str)
    variance_data = cursor.fetchall()
    variance_data_t = zip(*variance_data)
    '''
variance_data：
              sample
            1	0	1
variance 	0	0	1
		    0	0	0
		    1	0	0
		    1	0	0

variance_data_t：
             variance
         1	0	0	1	1
sample   0	0	0	0	0
         1	1	0	0	0
    '''
    # print variance_data_t
    pred_list = [str(sum(filter(lambda x: type(x) == int, i))) for i in variance_data_t]

    cmd_r_str = "library(logistf);pred=c({0});outcome=c({1});" \
                "gender=c({2});pc1=c({3});pc2=c({4});pc3=c({5});" \
                "pc4=c({6});pc5=c({7});" \
                "lr2=logistf(outcome ~ pred+gender+pc1+pc2+pc3+pc4+pc5);summary(lr2);" \
                "".format(",".join(pred_list),
                          ",".join(outcome_list),
                          ",".join(gender_list),
                          ",".join(pc1_list),
                          ",".join(pc2_list),
                          ",".join(pc3_list),
                          ",".join(pc4_list),
                          ",".join(pc5_list))
    with open("tmpR.tmpR", "w") as fp:
        fp.write(cmd_r_str)
    # print ("cmd_str=[{0}]".format(cmd_str))
    print("begin to execute R cmd")
    ph = Popen(['R -f tmpR.tmpR > tmpR.result'], shell=True)
    ph.wait()
    print("loading R result")
    with open("tmpR.result", "r") as fp:
        ret = fp.readlines()
    ph = Popen(["rm tmpR*"], shell=True)
    ph.wait()
    print("".join(ret))
    pred_line = filter(lambda k: len(k) > 0, filter(lambda x: x.startswith("pred"), ret)[0].strip().split(" "))
    pred_line2 = filter(lambda k: len(k) > 0, filter(lambda x: x.startswith("pred"), ret)[-1].strip().split(" "))
    # print pred_line
    coef = math.exp(float(pred_line[1]))
    lci = math.exp(float(pred_line[3]))
    uci = math.exp(float(pred_line[4]))
    p_value = float(pred_line2[-1])
    print("p_value = {}".format(p_value))
    print("coef = {}".format(coef))
    print("lower0.95 = {}".format(lci))
    print("upper0.95 = {}".format(uci))
    cursor.close()
    conn.commit()
    conn.close()
    return [coef, lci, uci, p_value]


def firth_logistic_regression_with_gene(db_file, sample_table, phenotype_in, sample_restrict,
                                        variance_table, variance_restrict, gene_table, gene_restrict,
                                        fai_in):
    # print "begin"
    # print "gene_restrict = [{0}]".format(gene_restrict)
    logging.basicConfig(filename="firth_logistic_regression_with_gene.log", level=logging.DEBUG, format=log_format,
                        filemode="w")
    logging.debug("begin")
    logging.debug("db=[{}]".format(db_file))
    logging.debug("sample_table=[{}]".format(sample_table))
    logging.debug("phenotype_in=[{}]".format(phenotype_in))
    logging.debug("sample_restrict=[{}]".format(sample_restrict))
    logging.debug("variance_table=[{}]".format(variance_table))
    logging.debug("variance_restrict=[{}]".format(variance_restrict))
    logging.debug("gene_table=[{}]".format(gene_table))
    logging.debug("gene_restrict=[{}]".format(gene_restrict))
    logging.debug("fai_in=[{}]".format(fai_in))
    chr2offset_dict = parse_fai(fai_in)
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()

    # handle sample
    if sample_restrict.strip():
        cmd_str = "select gen_id, {0},gender,PC1,PC2,PC3,PC4,PC5 from {1} as s where {2}".format(phenotype_in,
                                                                                                 sample_table,
                                                                                                 sample_restrict)
    else:
        cmd_str = "select gen_id, {0},gender,PC1,PC2,PC3,PC4,PC5 from {1}".format(phenotype_in, sample_table)
    # print "handling sample data"
    logging.debug("begin to get sample data. sql = [{}]".format(cmd_str))
    cursor.execute(cmd_str)
    sample_data = cursor.fetchall()
    sample_data_t = zip(*sample_data)

    outcome_list = list([str(i) for i in sample_data_t[1]])
    sample_id_list = list([str(i) for i in sample_data_t[0]])
    gender_list = list([str(i) for i in sample_data_t[2]])
    pc1_list = list([str(i) for i in sample_data_t[3]])
    pc2_list = list([str(i) for i in sample_data_t[4]])
    pc3_list = list([str(i) for i in sample_data_t[5]])
    pc4_list = list([str(i) for i in sample_data_t[6]])
    pc5_list = list([str(i) for i in sample_data_t[7]])

    # handle gene data
    cmd_str = "select chr, start_pos, end_pos, chr2, start_pos2, end_pos2 from {0} as g".format(gene_table)
    if gene_restrict.strip():
        cmd_str += " where {0}".format(gene_restrict.strip())
    logging.debug("handling gene_data. sql=[{0}]".format(cmd_str))
    cursor.execute(cmd_str)
    gene_data = cursor.fetchall()
    region_list = []
    for chr1, start1, end1, chr2, start2, end2 in gene_data:
        region_list.append([chr_pos2absolute_pos(str(chr1), start1, chr2offset_dict),
                            chr_pos2absolute_pos(str(chr1), end1, chr2offset_dict)])
        if chr2 is not None:
            region_list.append([chr_pos2absolute_pos(str(chr2), start2, chr2offset_dict),
                                chr_pos2absolute_pos(str(chr2), end2, chr2offset_dict)])
    logging.debug("len(region_list)={}".format(len(region_list)))
    # handle variance data
    if variance_restrict.strip():
        cmd_str = "select chr, pos, sample_{0} from {1} as v where {2}" \
                  "".format(",sample_".join(sample_id_list), variance_table, variance_restrict)
    else:
        cmd_str = "select chr, pos, sample_{0} from {1}" \
                  "".format(",sample_".join(sample_id_list), variance_table)
    logging.debug("handling variance_data")
    cursor.execute(cmd_str)
    variance_data = cursor.fetchall()
    logging.debug("len(variance_data)={}".format(len(variance_data)))
    # filter variance data with gene data
    logging.debug("filtering variance data with gene data")
    variance_data = filter(lambda x: abp_in_regions(chr_pos2absolute_pos(str(x[0]), x[1], chr2offset_dict),
                                                    region_list), variance_data)
    logging.debug("{} variance data left".format(len(variance_data)))

    variance_data_t = zip(*variance_data)
    variance_data_t = variance_data_t[2:]
    '''
variance_data：
              sample
            1	0	1
variance 	0	0	1
		    0	0	0
		    1	0	0
		    1	0	0

variance_data_t：
             variance
         1	0	0	1	1
sample   0	0	0	0	0
         1	1	0	0	0
    '''
    # print variance_data_t
    pred_list = [str(sum(filter(lambda x: type(x) == int, i))) for i in variance_data_t]

    cmd_r_str = "library(logistf);pred=c({0});outcome=c({1});" \
                "gender=c({2});pc1=c({3});pc2=c({4});pc3=c({5});" \
                "pc4=c({6});pc5=c({7});" \
                "lr2=logistf(outcome ~ pred+gender+pc1+pc2+pc3+pc4+pc5);summary(lr2);" \
                "".format(",".join(pred_list),
                          ",".join(outcome_list),
                          ",".join(gender_list),
                          ",".join(pc1_list),
                          ",".join(pc2_list),
                          ",".join(pc3_list),
                          ",".join(pc4_list),
                          ",".join(pc5_list))
    with open("tmpR.tmpR", "w") as fp:
        fp.write(cmd_r_str)
    logging.debug("begin to execute R cmd")
    ph = Popen(['R -f tmpR.tmpR > tmpR.result'], shell=True)
    ph.wait()
    logging.debug("loading R result")
    with open("tmpR.result", "r") as fp:
        ret = fp.readlines()
    # ph = Popen(["rm tmpR*"], shell=True)
    # ph.wait()
    print("".join(ret))
    pred_line = filter(lambda k: len(k) > 0, filter(lambda x: x.startswith("pred"), ret)[0].strip().split(" "))
    pred_line2 = filter(lambda k: len(k) > 0, filter(lambda x: x.startswith("pred"), ret)[-1].strip().split(" "))
    coef = math.exp(float(pred_line[1]))
    lci = math.exp(float(pred_line[3]))
    uci = math.exp(float(pred_line[4]))
    p_value = float(pred_line2[-1])
    print("p_value = {}".format(p_value))
    print("coef = {}".format(coef))
    print("lower0.95 = {}".format(lci))
    print("upper0.95 = {}".format(uci))
    cursor.close()
    conn.commit()
    conn.close()
    logging.debug("all done")
    return [coef, lci, uci, p_value]


def firth_logistic_regression_genelist(db_file, sample_table, phenotype_in, sample_restrict,
                                       variance_table, variance_restrict, gene_table, gene_restrict,
                                       fai_in, gene_list_file_in, id_cols):
    id_cols_list = [int(i) for i in id_cols.strip().split(",")]
    with open(gene_list_file_in, "r") as fp:
        gene_list_data = [i.strip().split("\t") for i in fp.readlines() if i.strip()]
    gene_list = []
    for col in id_cols_list:
        gene_list.extend(gene_list_data[col - 1])
    # (gene_id in ('ENSG00000180008', 'RAB22A') or gene_name in ('ENSG00000180008', 'RAB22A'))
    gene_restrict2 = "g.gene_id in ({0}) or g.gene_name in ({0})".format("'" + "', '".join(gene_list) + "'")
    if gene_restrict.strip():
        gene_restrict = "(" + gene_restrict.strip() + ") and ({0})".format(gene_restrict2)
    else:
        gene_restrict = gene_restrict2
    return firth_logistic_regression_with_gene(db_file, sample_table, phenotype_in, sample_restrict,
                                               variance_table, variance_restrict, gene_table, gene_restrict, fai_in)


def plot_result(xsize, ysize, file_out, y1, y2, yticks, ytickposes, fontsize, lablesize, marker,
                markersize, linewidth, result_all_list, left, right, top, bottom, hspace):
    """

    @param file_out:
    @param y1:
    @param y2:
    @param yticks:
    @param ytickposes:
    @param fontsize:
    @param lablesize:
    @param marker:
    @param markersize:
    @param linewidth:
    @param result_all_list: [[title, xlable, coef, lci, uci, p_value]]
    @return:
    """
    import matplotlib.pyplot as plt
    font1 = {'family': 'Arial',
             'weight': 'normal',
             'size': int(fontsize),
             }
    ytickposes = [float(i) for i in ytickposes.split(",")]
    yticks = yticks.split(",")
    clist = ["b", "g", "r", "c", "m", "y", "k"]
    x_list = []
    for result_list in result_all_list:
        result_list_t = zip(*(result_list[1:]))
        lci_list = [float(i) for i in result_list_t[1]]
        uci_list = [float(i) for i in result_list_t[2]]
        x_list.append(min(lci_list))
        x_list.append(max(uci_list))
    x1 = min(x_list)
    x2 = max(x_list)
    x2 = x2 + (x2 - x1) * 0.1
    x1 = x1 - (x2 - x1) * 0.1
    del result_list_t
    fig = plt.figure(figsize=(int(xsize), int(ysize)))
    sub_num = len(result_all_list)
    for i in xrange(sub_num):
        ax = plt.subplot(sub_num, 1, i + 1)
        result_list = result_all_list[i]
        title, xlable = result_list[0]
        for j in xrange(len(result_list) - 1):
            coef, lci, uci, p_value = result_list[j + 1]
            coef = float(coef)
            lci = float(lci)
            uci = float(uci)
            p_value = float(p_value)
            plt.scatter(coef, ytickposes[j], s=int(markersize), c=clist[j % len(clist)], marker=marker)
            plt.plot([lci, uci], [ytickposes[j], ytickposes[j]], clist[j % len(clist)] + "-",
                     linewidth=float(linewidth))
        plt.xlim(x1, x2)
        plt.ylim(float(y1), float(y2))
        plt.yticks(ytickposes, yticks)
        plt.xlabel(xlable, font1)
        plt.title(title, font1)
        plt.tick_params(labelsize=int(lablesize))
        labels = ax.get_xticklabels() + ax.get_yticklabels()
        [label.set_fontname('Arial') for label in labels]
    plt.subplots_adjust(left=float(left), bottom=float(bottom), right=float(right), top=float(top),
                        hspace=float(hspace))
    plt.savefig(file_out, format="pdf")


def firth_logistic_regression_plot(script_list_file_in):
    """
    y1	y2	yticks ytickposes fontsize lablesize marker markersize	linewidth	db_file	fai_in left right   top bottom  hspace
    title1	xlable1
    sample_table11	phenotype11	sample_restrict11	variance_table11	variance_restrict11	gene_table11	gene_restrict11
    sample_table12	phenotype12	sample_restrict12	variance_table12	variance_restrict12	gene_table12	gene_restrict12
    @param script_list_file_in:
    @return:
    """
    filename = script_list_file_in + ".pdf"
    script_list = []
    result_all_list = []
    with open(script_list_file_in, "r") as fp:
        [xsize, ysize, y1, y2, yticks, ytickposes, fontsize, lablesize, marker,
         markersize, linewidth, db_file, fai_in, left, bottom, right, top, hspace] = fp.readline().strip("\n").split(
            "\t")
        while True:
            data_line = fp.readline()
            if not data_line:
                break
            if len(data_line.strip()) > 0:
                script_list.append(data_line.strip("\n").split("\t"))
            else:
                print("script_list={0}".format(script_list))
                result_list = [script_list[0]]
                for i in xrange(1, len(script_list), 1):
                    [sample_table, phenotype, sample_restrict, variance_table, variance_restrict,
                     gene_table, gene_restrict] = script_list[i]
                    coef, lci, uci, p_value = firth_logistic_regression_with_gene(db_file, sample_table, phenotype,
                                                                                  sample_restrict, variance_table,
                                                                                  variance_restrict, gene_table,
                                                                                  gene_restrict, fai_in)
                    result_list.append([coef, lci, uci, p_value])
                    # print("result{0}={1}".format(i, [title, xlable, coef, lci, uci, p_value]))
                result_all_list.append(result_list)
                result_list = []
                script_list = []
        if len(script_list) > 0:
            print("script_list={0}".format(script_list))
            result_list = [script_list[0]]
            for i in xrange(1, len(script_list), 1):
                [sample_table, phenotype, sample_restrict, variance_table, variance_restrict,
                 gene_table, gene_restrict] = script_list[i]
                coef, lci, uci, p_value = firth_logistic_regression_with_gene(db_file, sample_table, phenotype,
                                                                              sample_restrict, variance_table,
                                                                              variance_restrict, gene_table,
                                                                              gene_restrict, fai_in)
                result_list.append([coef, lci, uci, p_value])
                # print("result{0}={1}".format(i, [title, xlable, coef, lci, uci, p_value]))
            result_all_list.append(result_list)
            del result_list
            del script_list

    plot_result(xsize, ysize, filename, y1, y2, yticks, ytickposes,
                fontsize, lablesize, marker, markersize, linewidth,
                result_all_list, left, right, top, bottom, hspace)


def build_contingency_table_new_with_gene_list_file(db_file, phenotype, sample_table_name, sample_restrict,
                                                    gene_table_name, gene_list_file, variance_table, variance_restrict,
                                                    fai_in):
    print("db_file=[{0}]\nphenotype=[{1}]\nsample_table_name=[{2}]\nsample_restrict=[{3}]" \
          "".format(db_file, phenotype, sample_table_name, sample_restrict))
    print("gene_table_name=[{0}]\ngene_list_file=[{1}]\nvariance_table=[{2}]\nvariance_restrict=[{3}]" \
          "".format(gene_table_name, gene_list_file, variance_table, variance_restrict))
    print("fai_in=[{}]".format(fai_in))
    with open(gene_list_file, "r") as fp:
        while True:
            gene_data = fp.readline()
            if not gene_data:
                break
            if gene_data.startswith("#"):
                continue
            gene_list = gene_data.strip().split("\t")
            gene_restrict = "g.gene_name IN ('" + "', '".join(gene_list[1:]) + \
                            "') or g.gene_id IN ('" + "', '".join(gene_list[1:]) + "')"
            print(gene_restrict)
            build_contingency_table_new(db_file, phenotype,
                                        sample_table_name, sample_restrict,
                                        gene_table_name, gene_restrict,
                                        variance_table, variance_restrict,
                                        fai_in, gene_list[0], 1, False)


def variance_distribution_plot(db_file, phenotype_in, sample_table, sample_restrict, variance_table,
                               variance_restrict, gene_table, gene_restrict, fai_in, output, h):
    """

    @param db_file:
    @param phenotype:
    @param variance_table:
    @param variance_restrict:
    @param gene_table:
    @param gene_restrict:
    @param fai_in:
    @param output:
    @return:
    """
    import matplotlib.pyplot as plt
    logging.basicConfig(filename="variance_distribution_plot.log", level=logging.DEBUG, format=log_format,
                        filemode="w")
    if os.path.exists(variance_restrict):
        with open(variance_restrict, "r") as fp:
            variance_restrict = fp.readline().strip()
    logging.debug("phenotype_in={}".format(phenotype_in))
    logging.debug("sample_restrict={}".format(sample_restrict))
    logging.debug("variance_table={}".format(variance_table))
    logging.debug("variance_restrict={}".format(variance_restrict))
    logging.debug("gene_restrict={}".format(gene_restrict))
    fontsize = 30
    xsize = 10
    ysize = 10
    font1 = {'family': 'Arial',
             'weight': 'normal',
             'size': int(fontsize),
             }
    fig = plt.figure(figsize=(int(xsize), int(ysize)))
    ax1 = plt.subplot(2, 1, 1)
    ax2 = plt.subplot(2, 1, 2)
    print("gene_restrict=[{0}]".format(gene_restrict))
    chr2offset_dict = parse_fai(fai_in)
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()

    # handle sample with sample restrict
    if sample_restrict.strip():
        cmd_str = "select gen_id, {0} from {1} as s where {2}".format(phenotype_in, sample_table, sample_restrict)
    else:
        cmd_str = "select gen_id, {0} from {1}".format(phenotype_in, sample_table)
    logging.debug("begin to get sample data. sql = [{}]".format(cmd_str))
    cursor.execute(cmd_str)
    sample_data = cursor.fetchall()
    sample_data_t = zip(*sample_data)

    sample_pheno_list = list([str(i) for i in sample_data_t[1]])
    control_num = len([i for i in sample_pheno_list if i == "0"])
    case_num = len([i for i in sample_pheno_list if i == "1"])
    sample_id_list = list([str(i) for i in sample_data_t[0]])

    # get candidate variance
    if variance_restrict.strip():
        cmd_str = "select chr, pos, sample_{0} from {1} as v where {2}" \
                  "".format(",sample_".join(sample_id_list), variance_table, variance_restrict)
    else:
        cmd_str = "select chr, pos, sample_{0} from {1}" \
                  "".format(",sample_".join(sample_id_list), variance_table)
    logging.debug("handling variance_data")
    cursor.execute(cmd_str)
    variance_data = cursor.fetchall()

    # handle gene
    if gene_restrict.startswith("directly"):
        logging.debug("got directly cmd. use all the variance in candidate variance")
    else:
        region_list = []
        if gene_restrict.startswith("restrict:"):
            gene_restrict = gene_restrict[9:].strip()
            if gene_restrict:
                gene_restrict = " where {}".format(gene_restrict)
            else:
                gene_restrict = ""
        elif gene_restrict.startswith("genelist:"):
            gene_list_file = gene_restrict[9:]
            with open(gene_list_file) as fp:
                name_list_str = "('" + "', '".join([i.strip() for i in fp.readlines() if len(i.strip()) > 0]) + "')"
            gene_restrict = " where gene_id in {0} or gene_name in {0}".format(name_list_str)
        else:
            print("illegal gene_restrict. exit.")
            cursor.close()
            conn.commit()
            conn.close()
            exit(0)

        sql_cmd = "select chr, start_pos, end_pos, chr2, start_pos2, end_pos2 from {0} as g{1}" \
                  "".format(gene_table, gene_restrict)
        print("sql_cmd = {}".format(sql_cmd))
        cursor.execute(sql_cmd)
        gene_data = cursor.fetchall()
        cursor.close()
        conn.commit()
        conn.close()
        for chr1, start1, end1, chr2, start2, end2 in gene_data:
            region_list.append([chr_pos2absolute_pos(str(chr1), start1, chr2offset_dict),
                                chr_pos2absolute_pos(str(chr1), end1, chr2offset_dict)])
            if chr2 is not None:
                region_list.append([chr_pos2absolute_pos(str(chr2), start2, chr2offset_dict),
                                    chr_pos2absolute_pos(str(chr2), end2, chr2offset_dict)])
        # filter variance data with gene data
        logging.debug("filtering {} variance data with gene data".format(len(variance_data)))
        variance_data = filter(lambda x: abp_in_regions(chr_pos2absolute_pos(str(x[0]), x[1], chr2offset_dict),
                                                        region_list), variance_data)

    variance_data = filter(lambda x: len(filter(lambda y: type(y) == int and y > 0, x[2:])) <= 20, variance_data)
    v_len = len(variance_data)
    logging.debug("{} variance points in each figure".format(v_len))
    icounter = 0
    max_ax_num = 0
    people_list = []
    allele_list = []
    for one_variance in variance_data:
        control_data = filter(lambda x: type(x) == int, [one_variance[i + 2] for i in xrange(len(sample_pheno_list)) if
                                                         sample_pheno_list[i] == "0"])
        control_people_with_variance_num = len(filter(lambda x: x > 0, control_data)) / float(control_num)
        control_allele_with_variance_num = sum(control_data) / float(control_num)

        case_data = filter(lambda x: type(x) == int,
                           [one_variance[i + 2] for i in xrange(len(sample_pheno_list)) if sample_pheno_list[i] == "1"])
        case_people_with_variance_num = len(filter(lambda x: x > 0, case_data)) / float(case_num)
        case_allele_with_variance_num = sum(case_data) / float(case_num)
        if max_ax_num < max([control_people_with_variance_num, control_allele_with_variance_num,
                             case_people_with_variance_num, case_allele_with_variance_num]):
            max_ax_num = max([control_people_with_variance_num, control_allele_with_variance_num,
                              case_people_with_variance_num, case_allele_with_variance_num])
        people_list.append((case_people_with_variance_num, control_people_with_variance_num))
        allele_list.append((case_allele_with_variance_num, control_allele_with_variance_num))
        icounter += 1
        if icounter % 100 == 0:
            print("{0} / {1} ({2:.2%})".format(icounter, v_len, icounter / float(v_len)))
    # max_ax_num_int = int(math.ceil(max_ax_num))
    people_merge_list = Counter(people_list).most_common()

    people_merge_x = [point[0] for point, times in people_merge_list]  # case
    people_merge_y = [point[1] for point, times in people_merge_list]  # control
    people_times = [times for point, times in people_merge_list]
    people_max_times = float(people_merge_list[0][1])
    print("people_max_times={0} x={1} y={2}".format(people_max_times,
                                                    people_merge_list[0][0][0],
                                                    people_merge_list[0][0][1]))
    print(people_times)
    allele_merge_list = Counter(allele_list).most_common()

    allele_merge_x = [point[0] for point, times in allele_merge_list]  # case
    allele_merge_y = [point[1] for point, times in allele_merge_list]  # control
    allele_times = [times for point, times in allele_merge_list]
    allele_max_times = float(allele_merge_list[0][1])
    print("allele_max_times={0} x={1} y={2}".format(allele_max_times,
                                                    allele_merge_list[0][0][0],
                                                    allele_merge_list[0][0][1]))
    print(allele_times)
    ax1.plot([0, max_ax_num], [0, max_ax_num], color="k", linestyle="--", linewidth=0.5, alpha=0.5)
    ax2.plot([0, max_ax_num], [0, max_ax_num], color="k", linestyle="--", linewidth=0.5, alpha=0.5)
    ax1.plot([0, max_ax_num - float(h)], [float(h), max_ax_num], color="k", linestyle="--", linewidth=0.5, alpha=0.5)
    ax2.plot([0, max_ax_num - float(h)], [float(h), max_ax_num], color="k", linestyle="--", linewidth=0.5, alpha=0.5)
    ax1.plot([float(h), max_ax_num], [0, max_ax_num - float(h)], color="k", linestyle="--", linewidth=0.5, alpha=0.5)
    ax2.plot([float(h), max_ax_num], [0, max_ax_num - float(h)], color="k", linestyle="--", linewidth=0.5, alpha=0.5)
    # sc1 = ax1.scatter(people_merge_x, people_merge_y, s=10, c=[math.log(i, 100) for i in people_times], edgecolors='none')
    # cb1 = plt.colorbar(sc1, ax=ax1)
    # sc2 = ax2.scatter(allele_merge_x, allele_merge_y, s=10, c=[math.log(i, 100) for i in allele_times], edgecolors='none')
    # cb2 = plt.colorbar(sc2, ax=ax2)
    sc1 = ax1.scatter(people_merge_x, people_merge_y, s=[math.log(i, 10) * 20 for i in people_times], c="r",
                      edgecolors='none')
    sc2 = ax2.scatter(allele_merge_x, allele_merge_y, s=[math.log(i, 10) * 20 for i in allele_times], c="r",
                      edgecolors='none')

    # ax1.set_xticks(xrange(0, max_ax_num + 1, 1))
    # ax1.set_yticks(xrange(0, max_ax_num + 1, 1))
    # ax2.set_xticks(xrange(0, max_ax_num + 1, 1))
    # ax2.set_yticks(xrange(0, max_ax_num + 1, 1))

    ax1.set_xlim(-max_ax_num * 0.1, max_ax_num * 1.1)
    ax1.set_ylim(-max_ax_num * 0.1, max_ax_num * 1.1)
    ax2.set_xlim(-max_ax_num * 0.1, max_ax_num * 1.1)
    ax2.set_ylim(-max_ax_num * 0.1, max_ax_num * 1.1)

    ax1.set_title("average sample number")
    ax2.set_title("average allele number")
    ax1.set_xlabel("Averaged nubmer of variants per case")
    ax1.set_ylabel("Averaged nubmer of variants per control")
    ax2.set_xlabel("Averaged nubmer of variants per case")
    ax2.set_ylabel("Averaged nubmer of variants per control")
    plt.show()
    plt.savefig(output, format="pdf")


def select_enriched_variance(db_file, sample_table, sample_restrict, phenotype_in, variance_table_name,
                             variance_restrict, gene_table, gene_restrict, fai_in, bigger_fraction,
                             smaller_fraction, enrich_in_case_out, enrich_in_control_out):
    def gene_overlape_with_variance(gene_chr, gene_start, gene_end, gene_chr2, gene_start2, gene_end2,
                                    v_chr, v_pos, chr2offset_dict):
        # print [gene_chr, gene_start, gene_end, gene_chr2, gene_start2, gene_end2]
        assert type(v_chr) == str
        assert type(gene_chr) == str
        assert type(v_pos) == int
        assert type(gene_start) == int
        assert type(gene_end) == int
        # if v_chr == "16" and v_pos == 284920 and gene_chr == "16" and gene_start == 283152 and gene_end == 287215:
        #     print "got it."
        v_abp = chr_pos2absolute_pos(v_chr, v_pos, chr2offset_dict)
        if gene_chr is not "None":
            ab_start1 = chr_pos2absolute_pos(gene_chr, gene_start, chr2offset_dict)
            ab_end1 = chr_pos2absolute_pos(gene_chr, gene_end, chr2offset_dict)
            if v_abp >= ab_start1 and v_abp <= ab_end1:
                return True
        if gene_chr2 is not "None":
            ab_start2 = chr_pos2absolute_pos(gene_chr2, gene_start2, chr2offset_dict)
            ab_end2 = chr_pos2absolute_pos(gene_chr2, gene_end2, chr2offset_dict)
            if v_abp >= ab_start2 and v_abp <= ab_end2:
                return True
        return False

    def variance_category(source_list):
        if len(source_list) == 7:
            lof_num = sum([int(i) for i in source_list[:3]])
            dmis_num = int(source_list[3])
            dsplicy_num = sum([int(i) for i in source_list[4:7]])
            ret_list = []
            if lof_num > 0:
                ret_list.append("lof")
            if dmis_num > 0:
                ret_list.append("dmis")
            if dsplicy_num > 0:
                ret_list.append("dsplicy")
            return ",".join(ret_list)
        elif len(source_list) == 3:
            return "lof"
        else:
            return "NA"

    if os.path.exists(variance_restrict):
        with open(variance_restrict, "r") as fp:
            variance_restrict = fp.readline().strip()

    bigger_fraction = float(bigger_fraction)
    smaller_fraction = float(smaller_fraction)
    chr2offset_dict = parse_fai(fai_in)
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()
    # handle sample with sample restrict
    if sample_restrict.strip():
        cmd_str = "select gen_id, {0} from {1} as s where {2}".format(phenotype_in, sample_table, sample_restrict)
    else:
        cmd_str = "select gen_id, {0} from {1}".format(phenotype_in, sample_table)
    print("getting sample data. sql = [{}]".format(cmd_str))
    cursor.execute(cmd_str)
    sample_data = cursor.fetchall()
    sample_data_t = zip(*sample_data)

    sample_pheno_list = list([str(i) for i in sample_data_t[1]])
    case_num = float(sample_pheno_list.count("1"))
    control_num = float(sample_pheno_list.count("0"))
    sample_id_list = list([str(i) for i in sample_data_t[0]])

    # get candidate variance
    vcf_info_num = 0
    vcf_category_num = 0

    cmd_str = "SELECT v.vcf_id, v.chr, v.pos, v.ref, v.alt"
    if db_has_col(cursor, variance_table_name, "bystro_cadd"):
        cmd_str += ",v.bystro_cadd"
        vcf_info_num += 1
    if db_has_col(cursor, variance_table_name, "bystro_phyloP"):
        cmd_str += ",v.bystro_phyloP"
        vcf_info_num += 1
    if db_has_col(cursor, variance_table_name, "bystro_phastCons"):
        cmd_str += ",v.bystro_phastCons"
        vcf_info_num += 1
    if db_has_col(cursor, variance_table_name, "ccrs"):
        cmd_str += ",v.ccrs"
        vcf_info_num += 1
    if db_has_col(cursor, variance_table_name, "domain_limbr"):
        cmd_str += ",domain_limbr"
        vcf_info_num += 1
    if db_has_col(cursor, variance_table_name, "exone_limbr"):
        cmd_str += ",exone_limbr"
        vcf_info_num += 1
    if db_has_col(cursor, variance_table_name, "bystro_sampleMaf"):
        cmd_str += ",bystro_sampleMaf"
        vcf_info_num += 1

    if db_has_col(cursor, variance_table_name, "annovar"):
        cmd_str += ",v.annovar"
        vcf_category_num += 1
    if db_has_col(cursor, variance_table_name, "bystro"):
        cmd_str += ",v.bystro"
        vcf_category_num += 1
    if db_has_col(cursor, variance_table_name, "vep"):
        cmd_str += ",v.vep"
        vcf_category_num += 1
    if db_has_col(cursor, variance_table_name, "dmis"):
        cmd_str += ",v.dmis"
        vcf_category_num += 1
    if db_has_col(cursor, variance_table_name, "dsplicing"):
        cmd_str += ",v.dsplicing"
        vcf_category_num += 1
    if db_has_col(cursor, variance_table_name, "spidex"):
        cmd_str += ",v.spidex"
        vcf_category_num += 1
    if db_has_col(cursor, variance_table_name, "spliceAI"):
        cmd_str += ",v.spliceAI"
        vcf_category_num += 1
    cmd_str += ",v.sample_" + ",v.sample_".join(sample_id_list)
    vcf_sample_num = len(sample_id_list)
    if variance_restrict:
        cmd_str += " FROM {0} AS v WHERE {1}".format(variance_table_name, variance_restrict)
    else:
        cmd_str += " FROM {0} AS v".format(variance_table_name)
    print("getting variance data.")
    cursor.execute(cmd_str)
    variance_data = cursor.fetchall()

    # handle gene data
    if len(gene_restrict) == 0:
        sql_cmd = "select gene_id, chr, start_pos, end_pos, chr2, start_pos2, end_pos2, gene_name from {0}" \
                  "".format(gene_table)
    else:
        sql_cmd = "select gene_id, chr, start_pos, end_pos, chr2, start_pos2, end_pos2, gene_name from {0} as g where {1}" \
                  "".format(gene_table, gene_restrict)
    print("gettint gene data. sql_cmd = {}".format(sql_cmd))
    cursor.execute(sql_cmd)
    gene_data = cursor.fetchall()
    cursor.close()
    conn.commit()
    conn.close()
    # print filter(lambda x: x[0] == "ENSG00000185615", gene_data)
    # print filter(lambda x: x[0] == "rs45619835", variance_data)[0][:5]
    # gene_id    vcf_id    chr    pos    ref    alt    category
    # cadd    phylop    phastCons    ccrs    domain_limbr    exon_limbr    Maf
    # case_sample    case_allele    control_sample    control_allele
    with open(enrich_in_case_out, "w") as fp_case, open(enrich_in_control_out, "w") as fp_control:
        icounter_case = 0
        icounter_control = 0
        icounter = 0
        head = "#gene_id\tgene_name\tvcf_id\tchr\tpos\tref\talt\tcategory\tcadd\tphylop\tphastCons\tccrs\tdomain_limbr" \
               "\texon_limbr\tMaf\tcase_sample\tcase_allele\tcontrol_sample\tcontrol_allele\n"
        fp_case.write(head)
        fp_control.write(head)
        variance_len = len(variance_data)
        for one_variance in variance_data:
            icounter += 1
            if icounter % 1000 == 0:
                print("{0} / {1}  ({2:.2%})".format(icounter, variance_len, icounter / float(variance_len)))
            selected_gene = []
            vcf_id, chrom, pos, ref, alt = [str(i) for i in one_variance[:5]]
            vcf_info_list = [str(i) for i in one_variance[5:5 + vcf_info_num]]
            vcf_category_list = one_variance[5 + vcf_info_num:5 + vcf_info_num + vcf_category_num]
            vcf_sample_list = one_variance[
                              5 + vcf_info_num + vcf_category_num:5 + vcf_info_num + vcf_category_num + vcf_sample_num]

            case_data = [vcf_sample_list[i] for i in xrange(vcf_sample_num) if sample_pheno_list[i] == "1"]
            control_data = [vcf_sample_list[i] for i in xrange(vcf_sample_num) if sample_pheno_list[i] == "0"]
            v_in_case = [i for i in case_data if type(i) == int and i > 0]
            v_in_control = [i for i in control_data if type(i) == int and i > 0]
            if len(vcf_info_list) == 7:
                vcf_info_str = "\t" + "\t".join(vcf_info_list)
            else:
                vcf_info_str = "\tNA\tNA\tNA\tNA\tNA\tNA\tNA"
            num_str = "\t{0}/{1}\t{2}/{3}\t{4}/{5}\t{6}/{7}" \
                      "".format(len(v_in_case), case_num,
                                sum(v_in_case), case_num * 2,
                                len(v_in_control), control_num,
                                sum(v_in_control), control_num * 2)

            # handle enrich in case
            if len(v_in_case) / case_num >= bigger_fraction and len(v_in_control) / control_num <= smaller_fraction:
                selected_gene = filter(
                    lambda x: gene_overlape_with_variance(str(x[1]), x[2], x[3], str(x[4]), x[5], x[6], chrom,
                                                          int(pos), chr2offset_dict), gene_data)
                icounter_case += 1
                if len(selected_gene) == 0:
                    fp_case.write("NA\tNA\t" + "\t".join([vcf_id, chrom, pos, ref, alt]) + "\t" +
                                  variance_category(vcf_category_list) + vcf_info_str + num_str + "\n")
                else:
                    tmp_str = "\t" + "\t".join([vcf_id, chrom, pos, ref, alt]) + "\t" + \
                              variance_category(vcf_category_list) + vcf_info_str + num_str + "\n"
                    for one_selected_gene in selected_gene:
                        fp_case.write(one_selected_gene[0] + "\t" + one_selected_gene[7] + "\t" + tmp_str)

            # handle enrich in control
            if len(v_in_control) / control_num >= bigger_fraction and len(v_in_case) / case_num <= smaller_fraction:
                icounter_control += 1
                if len(selected_gene) == 0:
                    selected_gene = filter(
                        lambda x: gene_overlape_with_variance(str(x[1]), x[2], x[3], str(x[4]), x[5], x[6],
                                                              chrom, int(pos), chr2offset_dict),
                        gene_data)
                # if vcf_id == "rs45619835":
                #     print "len(selected_gene)={}".format(len(selected_gene))
                if len(selected_gene) == 0:
                    fp_control.write("NA\tNA\t" + "\t".join([vcf_id, chrom, pos, ref, alt]) + "\t" +
                                     variance_category(vcf_category_list) + vcf_info_str + num_str + "\n")
                else:
                    tmp_str = "\t" + "\t".join([vcf_id, chrom, pos, ref, alt]) + "\t" + \
                              variance_category(vcf_category_list) + vcf_info_str + num_str + "\n"
                    for one_selected_gene in selected_gene:
                        fp_control.write(one_selected_gene[0] + "\t" + one_selected_gene[7] + "\t" + tmp_str)
    print("{} variance enrich in case".format(icounter_case))
    print("{} variance enrich in control".format(icounter_control))


def select_enriched_variance_gene_name(db_file, sample_table, sample_restrict, phenotype_in,
                                       variance_table_name, variance_restrict, gene_table, gene_restrict,
                                       bigger_fraction, smaller_fraction, output, group):
    def variance_category(source_list):

        if len(source_list) == 7:
            # lof_num = sum([int(i) for i in filter(lambda x: x is not None, source_list[:3])])
            # dmis_num = int(source_list[3])
            # dsplicy_num = sum([int(i) for i in filter(lambda x: x is not None, source_list[4:7])])
            # ret_list = []
            # if lof_num > 0:
            #     ret_list.append("lof")
            # if dmis_num > 0:
            #     ret_list.append("dmis")
            # if dsplicy_num > 0:
            #     ret_list.append("dsplicy")
            category_list = ["annovar", "bystro", "vep", "dmis", "dsplicing", "spidex", "spliceAI"]
            return ";".join(list(compress(category_list, [i == '1' for i in source_list])))
        elif len(source_list) == 3:
            category_list = ["annovar", "bystro", "vep"]
            return ";".join(list(compress(category_list, [i == '1' for i in source_list])))
        else:
            return "NA"

    def sub_gene_name(variance_table):
        sub_str_list = []
        if variance_table == "variance":
            # if re.findall("annovar[ =]", variance_restrict):
            sub_str_list.append("annovar_gene_name")
            # if re.findall("bystro[ =]", variance_restrict):
            sub_str_list.append("bystro_gene_name")
            # if re.findall("vep[ =]", variance_restrict):
            sub_str_list.append("vep_gene_id")
            # if re.findall("spliceAI[ =]", variance_restrict):
            sub_str_list.append("spliceAI_gene_name")
            # if re.findall("dmis[ =]", variance_restrict):
            sub_str_list.append("dmis_gene_name")
            # if re.findall("dsplicing[ =]", variance_restrict):
            sub_str_list.append("dsplicing_gene_name")
            # assert len(sub_str_list) > 0

        else:
            # if re.findall("annovar[ =]", variance_restrict):
            sub_str_list.append("annovar_gene_name")
            # if re.findall("bystro[ =]", variance_restrict):
            sub_str_list.append("bystro_gene_name")
            # if re.findall("vep[ =]", variance_restrict):
            sub_str_list.append("vep_gene_id")
            # assert len(sub_str_list) > 0
        # ret_str = ", ".join(sub_str_list)
        return sub_str_list

    gene_name2detail_dict = {"annovar_gene_name": "annovar_function_detail",
                             "bystro_gene_name": "bystro_final_function",
                             "vep_gene_id": "vep_consequence",
                             "spliceAI_gene_name": "spliceAI_anno",
                             "dmis_gene_name": "dmis_HGVSc",
                             "dsplicing_gene_name": "dsplicing_detailed_consequence"}
    gene_name2detail_synonymous_dict = {"annovar_gene_name": "annovar_detail",
                                        "bystro_gene_name": "bystro_final_function",
                                        "vep_gene_id": "vep_consequence"}
    if os.path.exists(variance_restrict):
        with open(variance_restrict, "r") as fp:
            variance_restrict = fp.readline().strip()

    bigger_fraction = float(bigger_fraction)
    smaller_fraction = float(smaller_fraction)
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()
    # handle sample with sample restrict
    if sample_restrict.strip():
        cmd_str = "select gen_id, {0} from {1} as s where {2}".format(phenotype_in, sample_table, sample_restrict)
    else:
        cmd_str = "select gen_id, {0} from {1}".format(phenotype_in, sample_table)
    print("getting sample data. sql = [{}]".format(cmd_str))
    cursor.execute(cmd_str)
    sample_data = cursor.fetchall()
    sample_data_t = zip(*sample_data)

    sample_pheno_list = list([str(i) for i in sample_data_t[1]])
    case_num = float(sample_pheno_list.count("1"))
    control_num = float(sample_pheno_list.count("0"))
    sample_id_list = list([str(i) for i in sample_data_t[0]])

    # get candidate variance
    vcf_info_num = 0
    vcf_category_num = 0

    cmd_str = "SELECT v.vcf_id, v.chr, v.pos, v.ref, v.alt"
    if db_has_col(cursor, variance_table_name, "bystro_cadd"):
        cmd_str += ",v.bystro_cadd"
        vcf_info_num += 1
    if db_has_col(cursor, variance_table_name, "bystro_phyloP"):
        cmd_str += ",v.bystro_phyloP"
        vcf_info_num += 1
    if db_has_col(cursor, variance_table_name, "bystro_phastCons"):
        cmd_str += ",v.bystro_phastCons"
        vcf_info_num += 1
    if db_has_col(cursor, variance_table_name, "ccrs"):
        cmd_str += ",v.ccrs"
        vcf_info_num += 1
    if db_has_col(cursor, variance_table_name, "domain_limbr"):
        cmd_str += ",domain_limbr"
        vcf_info_num += 1
    if db_has_col(cursor, variance_table_name, "exone_limbr"):
        cmd_str += ",exone_limbr"
        vcf_info_num += 1
    if db_has_col(cursor, variance_table_name, "bystro_sampleMaf"):
        cmd_str += ",bystro_sampleMaf"
        vcf_info_num += 1
    if db_has_col(cursor, variance_table_name, "is_ccds"):
        cmd_str += ",is_ccds"
        vcf_info_num += 1

    if db_has_col(cursor, variance_table_name, "annovar"):
        cmd_str += ",v.annovar"
        vcf_category_num += 1
    if db_has_col(cursor, variance_table_name, "bystro"):
        cmd_str += ",v.bystro"
        vcf_category_num += 1
    if db_has_col(cursor, variance_table_name, "vep"):
        cmd_str += ",v.vep"
        vcf_category_num += 1
    if db_has_col(cursor, variance_table_name, "dmis"):
        cmd_str += ",v.dmis"
        vcf_category_num += 1
    if db_has_col(cursor, variance_table_name, "dsplicing"):
        cmd_str += ",v.dsplicing"
        vcf_category_num += 1
    if db_has_col(cursor, variance_table_name, "spidex"):
        cmd_str += ",v.spidex"
        vcf_category_num += 1
    if db_has_col(cursor, variance_table_name, "spliceAI"):
        cmd_str += ",v.spliceAI"
        vcf_category_num += 1

    gene_name_list = sub_gene_name(variance_table_name)
    cmd_str += " ,v." + " ,v.".join(gene_name_list)
    gene_name_num = len(gene_name_list)

    detail_name_list = [
        gene_name2detail_dict[i] if variance_table_name == "variance" else gene_name2detail_synonymous_dict[i]
        for i in gene_name_list]
    cmd_str += " ,v." + " ,v.".join(detail_name_list)
    detail_num = len(detail_name_list)

    cmd_str += ",v.sample_" + ",v.sample_".join(sample_id_list)
    vcf_sample_num = len(sample_id_list)
    if variance_restrict:
        cmd_str += " FROM {0} AS v WHERE {1}".format(variance_table_name, variance_restrict)
    else:
        cmd_str += " FROM {0} AS v".format(variance_table_name)
    print("getting variance data.")
    cursor.execute(cmd_str)
    variance_data = cursor.fetchall()

    with open(output, "w") as fp_out:
        icounter_variance = 0
        icounter = 0
        head = "#case num={0}\n" \
               "#control num={1}\n" \
               "#phenotype={2}\n" \
               "#gene_name\tgene_id\tvcf_id\tchr\tpos\tref\talt\tcategory\tcadd\tphylop\tphastCons\tccrs" \
               "\tdomain_limbr\texon_limbr\tMaf\tis_ccds\tcase_sample\tcase_allele\tcontrol_sample\tcontrol_allele\tdetail\n" \
               "".format(case_num, control_num, phenotype_in)

        fp_out.write(head)
        variance_len = len(variance_data)
        for one_variance in variance_data:
            icounter += 1
            if icounter % 1000 == 0:
                print("{0} / {1}  ({2:.2%})".format(icounter, variance_len, icounter / float(variance_len)))
            selected_gene = []
            vcf_id, chrom, pos, ref, alt = [str(i) for i in one_variance[:5]]
            vcf_info_list = [str(i) for i in one_variance[5:5 + vcf_info_num]]
            offset = 5 + vcf_info_num
            vcf_category_list = one_variance[offset:offset + vcf_category_num]
            offset += vcf_category_num
            gene_name_info_list = one_variance[offset: offset + gene_name_num]
            offset += gene_name_num
            detail_info_list = one_variance[offset: offset + detail_num]
            detail_str = ";".join(filter(lambda x: x is not None, detail_info_list))
            offset += detail_num
            vcf_sample_list = one_variance[offset:offset + vcf_sample_num]

            gene_name_list = []
            for gene_name_info in gene_name_info_list:
                if gene_name_info is not None:
                    gene_name_list.extend(re.split(",|;", gene_name_info))
            gene_name_list = list(set(gene_name_list))

            case_data = [vcf_sample_list[i] for i in xrange(vcf_sample_num) if sample_pheno_list[i] == "1"]
            control_data = [vcf_sample_list[i] for i in xrange(vcf_sample_num) if sample_pheno_list[i] == "0"]
            v_in_case = [i for i in case_data if type(i) == int and i > 0]
            v_in_control = [i for i in control_data if type(i) == int and i > 0]
            if len(vcf_info_list) == 8:
                vcf_info_str = "\t" + "\t".join(vcf_info_list)
            else:
                vcf_info_str = "\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA"
            num_str = "\t{0}\t{1}\t{2}\t{3}".format(len(v_in_case), sum(v_in_case),
                                                    len(v_in_control), sum(v_in_control))

            # handle enrich in case
            if group == "case":
                if len(v_in_case) / case_num >= bigger_fraction and len(v_in_control) / control_num <= smaller_fraction:
                    icounter_variance += 1
                    if len(gene_name_list) == 0:
                        fp_out.write("NA\tNA\t" + "\t".join([vcf_id, chrom, pos, ref, alt]) + "\t" +
                                     variance_category(
                                         vcf_category_list) + vcf_info_str + num_str + "\t" + detail_str + "\n")
                    else:
                        if len(gene_restrict) == 0:
                            cursor.execute("select gene_name, gene_id from {1} where gene_name in ('{0}') "
                                           "or gene_id in ('{0}')"
                                           "".format("', '".join(gene_name_list), gene_table))
                        else:
                            cursor.execute("select gene_name, gene_id from {2} as g where (gene_name in ('{0}') "
                                           "or gene_id in ('{0}')) and {1}"
                                           "".format("', '".join(gene_name_list), gene_restrict, gene_table))
                        selected_gene = cursor.fetchall()
                        tmp_str = "\t".join([vcf_id, chrom, pos, ref, alt]) + "\t" + \
                                  variance_category(
                                      vcf_category_list) + vcf_info_str + num_str + "\t" + detail_str + "\n"
                        for one_selected_gene in selected_gene:
                            fp_out.write(str(one_selected_gene[0]) + "\t" + str(one_selected_gene[1]) + "\t" + tmp_str)

            # handle enrich in control
            elif group == "control":
                if len(v_in_control) / control_num >= bigger_fraction and len(v_in_case) / case_num <= smaller_fraction:
                    icounter_variance += 1
                    if len(gene_name_list) == 0:
                        fp_out.write("NA\tNA\t" + "\t".join([vcf_id, chrom, pos, ref, alt]) + "\t" +
                                     variance_category(
                                         vcf_category_list) + vcf_info_str + num_str + "\t" + detail_str + "\n")
                    else:
                        if len(selected_gene) == 0:
                            if len(gene_restrict) == 0:
                                cursor.execute("select gene_name, gene_id from {1} where gene_name in ('{0}') "
                                               "or gene_id in ('{0}')"
                                               "".format("', '".join(gene_name_list), gene_table))
                            else:
                                cursor.execute("select gene_name, gene_id from {2} as g where (gene_name in ('{0}') "
                                               "or gene_id in ('{0}')) and {1}"
                                               "".format("', '".join(gene_name_list), gene_restrict, gene_table))
                            selected_gene = cursor.fetchall()
                        tmp_str = "\t".join([vcf_id, chrom, pos, ref, alt]) + "\t" + \
                                  variance_category(
                                      vcf_category_list) + vcf_info_str + num_str + "\t" + detail_str + "\n"
                        for one_selected_gene in selected_gene:
                            fp_out.write(str(one_selected_gene[0]) + "\t" + str(one_selected_gene[1]) + "\t" + tmp_str)
            # exit(0)
    print("{0} variance enrich in {1}".format(icounter_variance, group))
    cursor.close()
    conn.commit()
    conn.close()


def gene_distribution_plot(db_file, phenotype_in, sample_table, sample_restrict, variance_table,
                           variance_restrict, gene_table, gene_restrict, fai_in, output, h):
    import matplotlib.pyplot as plt
    logging.basicConfig(filename="gene_distribution_plot.log", level=logging.DEBUG, format=log_format,
                        filemode="w")
    logging.debug("phenotype_in={}".format(phenotype_in))
    logging.debug("sample_restrict={}".format(sample_restrict))
    logging.debug("variance_table={}".format(variance_table))
    logging.debug("variance_restrict={}".format(variance_restrict))
    logging.debug("gene_restrict={}".format(gene_restrict))
    fontsize = 30
    xsize = 10
    ysize = 10
    font1 = {'family': 'Arial',
             'weight': 'normal',
             'size': int(fontsize),
             }
    fig = plt.figure(figsize=(int(xsize), int(ysize)))
    ax1 = plt.subplot(1, 1, 1)
    # ax2 = plt.subplot(2, 1, 2)
    print("gene_restrict=[{0}]".format(gene_restrict))
    chr2offset_dict = parse_fai(fai_in)
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()

    # handle sample with sample restrict
    if sample_restrict.strip():
        cmd_str = "select gen_id, {0} from {1} as s where {2}".format(phenotype_in, sample_table, sample_restrict)
    else:
        cmd_str = "select gen_id, {0} from {1}".format(phenotype_in, sample_table)
    logging.debug("begin to get sample data. sql = [{}]".format(cmd_str))
    cursor.execute(cmd_str)
    sample_data = cursor.fetchall()
    sample_data_t = zip(*sample_data)

    sample_pheno_list = list([str(i) for i in sample_data_t[1]])
    control_num = len([i for i in sample_pheno_list if i == "0"])
    case_num = len([i for i in sample_pheno_list if i == "1"])
    sample_id_list = list([str(i) for i in sample_data_t[0]])

    # get candidate variance
    if variance_restrict.strip():
        cmd_str = "select chr, pos, sample_{0} from {1} as v where {2}" \
                  "".format(",sample_".join(sample_id_list), variance_table, variance_restrict)
    else:
        cmd_str = "select chr, pos, sample_{0} from {1}" \
                  "".format(",sample_".join(sample_id_list), variance_table)
    logging.debug("handling variance_data")
    cursor.execute(cmd_str)
    variance_data = cursor.fetchall()

    # handle gene
    if gene_restrict.startswith("directly"):
        logging.error("You should restrict the gene when plotting the gene plot. Exit.")
        cursor.close()
        conn.commit()
        conn.close()
        exit(0)
    else:
        # region_list = []
        if gene_restrict.startswith("restrict:"):
            gene_restrict = gene_restrict[9:].strip()
            if gene_restrict:
                gene_restrict = " where {}".format(gene_restrict)
            else:
                gene_restrict = ""
        elif gene_restrict.startswith("genelist:"):
            gene_list_file = gene_restrict[9:]
            with open(gene_list_file) as fp:
                name_list_str = "('" + "', '".join([i.strip() for i in fp.readlines() if len(i.strip()) > 0]) + "')"
            gene_restrict = " where gene_id in {0} or gene_name in {0}".format(name_list_str)
        else:
            print("illegal gene_restrict. exit.")
            cursor.close()
            conn.commit()
            conn.close()
            exit(0)

        sql_cmd = "select chr, start_pos, end_pos, chr2, start_pos2, end_pos2 from {0} as g{1}" \
                  "".format(gene_table, gene_restrict)
        print("sql_cmd = {}".format(sql_cmd))
        cursor.execute(sql_cmd)
        gene_data = cursor.fetchall()
        cursor.close()
        conn.commit()
        conn.close()
        # for chr1, start1, end1, chr2, start2, end2 in gene_data:
        #     region_list.append([chr_pos2absolute_pos(str(chr1), start1, chr2offset_dict),
        #                         chr_pos2absolute_pos(str(chr1), end1, chr2offset_dict)])
        #     if chr2 is not None:
        #         region_list.append([chr_pos2absolute_pos(str(chr2), start2, chr2offset_dict),
        #                             chr_pos2absolute_pos(str(chr2), end2, chr2offset_dict)])
        # # filter variance data with gene data
        # logging.debug("filtering {} variance data with gene data".format(len(variance_data)))
        # variance_data = filter(lambda x: abp_in_regions(chr_pos2absolute_pos(str(x[0]), x[1], chr2offset_dict),
        #                                                 region_list), variance_data)

    variance_data = filter(lambda x: len(filter(lambda y: type(y) == int and y > 0, x[2:])) <= 20, variance_data)
    point_list = []
    max_ax_num = 0
    gene_data_len = len(gene_data)
    icounter = 0
    for chr1, start1, end1, chr2, start2, end2 in gene_data:
        icounter += 1
        if icounter % 100 == 0:
            print("handling {0} / {1} gene".format(icounter, gene_data_len))
        region_list = [[chr_pos2absolute_pos(str(chr1), start1, chr2offset_dict),
                        chr_pos2absolute_pos(str(chr1), end1, chr2offset_dict)]]
        if chr2 is not None:
            region_list.append([chr_pos2absolute_pos(str(chr2), start2, chr2offset_dict),
                                chr_pos2absolute_pos(str(chr2), end2, chr2offset_dict)])
        tmp_variance_data = filter(lambda x: abp_in_regions(chr_pos2absolute_pos(str(x[0]), x[1], chr2offset_dict),
                                                            region_list), variance_data)
        if len(tmp_variance_data) == 0:
            continue
        tmp_variance_data_t = zip(*tmp_variance_data)[2:]
        tmp_gene_data = [sum(filter(lambda x: type(x) == int, i)) for i in tmp_variance_data_t]
        control_data = [tmp_gene_data[i] for i in xrange(len(sample_pheno_list)) if sample_pheno_list[i] == "0"]
        control_people_with_variance_in_gene_num = len(filter(lambda x: x > 0, control_data)) / float(control_num)
        case_data = [tmp_gene_data[i] for i in xrange(len(sample_pheno_list)) if sample_pheno_list[i] == "1"]
        case_people_with_variance_in_gene_num = len(filter(lambda x: x > 0, case_data)) / float(case_num)
        point_list.append((case_people_with_variance_in_gene_num, control_people_with_variance_in_gene_num))
        if max_ax_num < max(case_people_with_variance_in_gene_num, control_people_with_variance_in_gene_num):
            max_ax_num = max(case_people_with_variance_in_gene_num, control_people_with_variance_in_gene_num)
    people_merge_list = Counter(point_list).most_common()
    people_merge_x = [point[0] for point, times in people_merge_list]  # case
    people_merge_y = [point[1] for point, times in people_merge_list]  # control
    people_times = [times for point, times in people_merge_list]

    ax1.plot([0, max_ax_num], [0, max_ax_num], color="k", linestyle="--", linewidth=0.5, alpha=0.5)
    ax1.plot([0, max_ax_num - float(h)], [float(h), max_ax_num], color="k", linestyle="--", linewidth=0.5, alpha=0.5)
    ax1.plot([float(h), max_ax_num], [0, max_ax_num - float(h)], color="k", linestyle="--", linewidth=0.5, alpha=0.5)
    sc1 = ax1.scatter(people_merge_x, people_merge_y, s=[math.log(i, 10) * 20 for i in people_times], c="r",
                      edgecolors='none')

    ax1.set_xlim(-max_ax_num * 0.1, max_ax_num * 1.1)
    ax1.set_ylim(-max_ax_num * 0.1, max_ax_num * 1.1)

    ax1.set_title("average sample number")
    ax1.set_xlabel("Averaged nubmer of variants per case")
    ax1.set_ylabel("Averaged nubmer of variants per control")
    plt.show()
    plt.savefig(output, format="pdf")


def gene_set_fisher_test(gene_set_file, db_file, sample_table, sample_restrict, phenotype_in,
                         variance_table, variance_restrict, gene_table, fai_in, output):
    def load_gene_set(gene_set_file):
        gene_set_dict = {}
        with open(gene_set_file, "r") as fp:
            while True:
                data_line = fp.readline()
                if not data_line:
                    break
                if not data_line.strip():
                    continue
                data_list = data_line.strip().split("\t")
                data_list = [i.strip() for i in filter(lambda x: len(x.strip()) > 0, data_list)]
                gene_set_dict[data_list[0]] = data_list[1:]
        return gene_set_dict

    if os.path.exists(variance_restrict):
        with open(variance_restrict, "r") as fp:
            variance_restrict = fp.readline().strip()
    logging.basicConfig(filename="{0}.log".format(sys._getframe().f_code.co_name),
                        level=logging.DEBUG, format=log_format, filemode="w")
    logging.debug("begin")
    logging.debug("gene_set_file=[{0}]".format(gene_set_file))
    logging.debug("db_file=[{0}]".format(db_file))
    logging.debug("phenotype=[{0}]".format(phenotype_in))
    logging.debug("sample_table_name=[{0}]".format(sample_table))
    logging.debug("sample_restrict=[{0}]".format(sample_restrict))
    logging.debug("output=[{0}]".format(output))

    chr2offset_dict = parse_fai(fai_in)
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()

    # handle sample with sample restrict
    if sample_restrict.strip():
        cmd_str = "select gen_id, {0} from {1} as s where {2}".format(phenotype_in, sample_table, sample_restrict)
    else:
        cmd_str = "select gen_id, {0} from {1}".format(phenotype_in, sample_table)
    logging.debug("begin to get sample data. sql = [{}]".format(cmd_str))
    cursor.execute(cmd_str)
    sample_data = cursor.fetchall()
    sample_data_t = zip(*sample_data)

    sample_pheno_list = list([str(i) for i in sample_data_t[1]])
    control_num = len([i for i in sample_pheno_list if i == "0"])
    case_num = len([i for i in sample_pheno_list if i == "1"])
    sample_id_list = list([str(i) for i in sample_data_t[0]])

    # get candidate variance
    if variance_restrict.strip():
        cmd_str = "select vcf_id, chr, pos, sample_{0} from {1} as v where {2}" \
                  "".format(",sample_".join(sample_id_list), variance_table, variance_restrict)
    else:
        cmd_str = "select vcf_id, chr, pos, sample_{0} from {1}" \
                  "".format(",sample_".join(sample_id_list), variance_table)
    logging.debug("handling variance_data")
    cursor.execute(cmd_str)
    variance_data = cursor.fetchall()
    # print(variance_data[0])
    gene_set_dict = load_gene_set(gene_set_file)
    with open(output, "w") as fp:
        fp.write("##phenotype:\"{0}\"\n##control number:{1}\n##case number:{2}\n"
                 "##sample restrict:\"{3}\"\n##test method:fisher exact test\n"
                 "".format(phenotype_in,
                           control_num,
                           case_num,
                           sample_restrict))
        fp.write("""##                         case     control
##                      ---------------------
##                      |         |         |
##    alt allele number |    A    |    B    |
##                      |         |         |
##                      ---------------------
##                      |         |         |
##    ref allele number |    C    |    D    |
##                      |         |         |
##                      ---------------------
##
##
##                              case     control
##                           ---------------------
##                           |         |         |
##          subjects have alt|    A1   |    B1   |
##                           |         |         |
##                           ---------------------
##                           |         |         |
##  subjects do not have alt |    C1   |    D1   |
##                           |         |         |
##                           ---------------------
""")
        fp.write("#gene_set_name\tA\tB\tC\tD\tp_value\todds_ratio\tA1\tB1\tC1\tD1\tp_value1\todds_ratio1\t"
                 "n_variance_in_gene_set\tvariance_in_gene_name\tn_variance_control\tvariance_in_control_name\t"
                 "n_variance_case\tvariance_in_case_name\n")
        for gene_set_name in gene_set_dict:
            logging.debug("handling gene set {}".format(gene_set_name))
            region_list = []
            gene_name_list = gene_set_dict[gene_set_name]
            name_list_str = "('" + "', '".join(gene_name_list) + "')"
            gene_restrict = " where gene_id in {0} or gene_name in {0}".format(name_list_str)
            sql_cmd = "select chr, start_pos, end_pos, chr2, start_pos2, end_pos2 from {0}{1}" \
                      "".format(gene_table, gene_restrict)
            # logging.debug("gene sql = [{}]".format(sql_cmd))
            cursor.execute(sql_cmd)
            gene_data = cursor.fetchall()
            # logging.debug("gene_data={}".format(gene_data))
            for chr1, start1, end1, chr2, start2, end2 in gene_data:
                region_list.append([chr_pos2absolute_pos(str(chr1), start1, chr2offset_dict),
                                    chr_pos2absolute_pos(str(chr1), end1, chr2offset_dict)])
                if chr2 is not None:
                    region_list.append([chr_pos2absolute_pos(str(chr2), start2, chr2offset_dict),
                                        chr_pos2absolute_pos(str(chr2), end2, chr2offset_dict)])

            selected_variance_data = filter(lambda x: abp_in_regions(chr_pos2absolute_pos(str(x[1]),
                                                                                          x[2],
                                                                                          chr2offset_dict),
                                                                     region_list), variance_data)
            if len(selected_variance_data) == 0:
                logging.debug("there is no variance in gene set, {}.".format(gene_set_name))
                logging.debug("region_list={}".format(region_list))
                continue
            all_variance_in_gene_num = len(selected_variance_data)
            all_variance_in_gene_name_list = [i[0] for i in selected_variance_data]

            selected_variance_data_t = zip(*selected_variance_data)[3:]
            assert len(selected_variance_data_t) == len(sample_pheno_list)
            control_data_t = [selected_variance_data_t[i] for i in xrange(len(selected_variance_data_t)) if
                              sample_pheno_list[i] == "0"]
            case_data_t = [selected_variance_data_t[i] for i in xrange(len(selected_variance_data_t)) if
                           sample_pheno_list[i] == "1"]
            control_data = zip(*control_data_t)
            case_data = zip(*case_data_t)

            tmp_list = []
            control_variance_in_gene_num = 0
            control_variance_in_gene_name_list = []
            for i in xrange(len(control_data)):
                control_line_data = filter(lambda x: type(x) == int, control_data[i])
                tmp_list.extend(control_line_data)
                if sum(control_line_data) > 0:
                    control_variance_in_gene_num += 1
                    control_variance_in_gene_name_list.append(selected_variance_data[i][0])
            B = sum(tmp_list)
            D = 2 * len(tmp_list) - B

            tmp_list = []
            case_variance_in_gene_num = 0
            case_variance_in_gene_name_list = []
            for i in xrange(len(case_data)):
                case_line_data = filter(lambda x: type(x) == int, case_data[i])
                tmp_list.extend(case_line_data)
                if sum(case_line_data) > 0:
                    case_variance_in_gene_num += 1
                    case_variance_in_gene_name_list.append(selected_variance_data[i][0])
            A = sum(tmp_list)
            C = 2 * len(tmp_list) - A
            if A + B < 2:
                continue

            control_people_data = [sum(filter(lambda x: type(x) == int, people_line)) for people_line in control_data_t]
            B1 = len(filter(lambda x: x > 0, control_people_data))
            D1 = len(control_people_data) - B1
            case_people_data = [sum(filter(lambda x: type(x) == int, people_line)) for people_line in case_data_t]
            A1 = len(filter(lambda x: x > 0, case_people_data))
            C1 = len(case_people_data) - A1
            oddsratio, pvalue = stats.fisher_exact([[A, B], [C, D]])
            oddsratio1, pvalue1 = stats.fisher_exact([[A1, B1], [C1, D1]])
            fp.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t"
                     "{10}\t{11}\t{12}\t{13}\t{14}\t{15}\t{16}\t{17}\t{18}\n"
                     "".format(gene_set_name, A, B, C, D, pvalue, oddsratio,
                               A1, B1, C1, D1, pvalue1, oddsratio1,
                               all_variance_in_gene_num,
                               ",".join(all_variance_in_gene_name_list) if len(
                                   all_variance_in_gene_name_list) > 0 else ".",
                               control_variance_in_gene_num,
                               ",".join(control_variance_in_gene_name_list) if len(
                                   control_variance_in_gene_name_list) > 0 else ".",
                               case_variance_in_gene_num,
                               ",".join(case_variance_in_gene_name_list) if len(
                                   case_variance_in_gene_name_list) > 0 else "."))
    cursor.close()
    conn.commit()
    conn.close()
    print("all done")
    logging.debug("all done")


def gene_set_fisher_test2(db_file, sample_table, sample_restrict, phenotype_in, variance_table,
                          variance_restrict, gene_table, gene_restrict, gene_set_name, fai_in,
                          output):
    tmp_gene_set_file = "{}.tmp".format(gene_set_name)
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()
    if len(gene_restrict) > 0:
        cmd_str = "select gene_id from {0} as g where {1}".format(gene_table, gene_restrict)
    else:
        cmd_str = "select gene_id from {0}".format(gene_table)
    cursor.execute(cmd_str)
    gene_data = cursor.fetchall()
    gene_id_list = [i[0] for i in gene_data]
    print(gene_id_list)
    with open(tmp_gene_set_file, "w") as fp:
        fp.write("{}\t".format(gene_set_name) + "\t".join(gene_id_list) + "\n")
    gene_set_fisher_test(tmp_gene_set_file, db_file, sample_table, sample_restrict, phenotype_in, variance_table,
                         variance_restrict, gene_table, fai_in, output)


def analyze_fisher_test_variance_May19(database, sample_restrict, fai_in, path):
    table_name_element_dict = {"": "",
                               " AND (annovar = 1)": "annovar",
                               " AND (bystro = 1)": "bystro",
                               " AND (dmis = 1)": "dmis",
                               " AND (dsplicing = 1)": "dsplicing",
                               " AND (spidex = 1)": "spidex",
                               " AND (spliceAI = 1)": "spliceAI",
                               " AND (vep = 1)": "vep",
                               " AND (annovar = 1 or bystro = 1 or vep = 1)": "LOF",
                               " AND (dsplicing = 1 or spidex = 1 or spliceAI = 1)": "dsplicing.all",
                               "tof6": "tof",
                               "CTD": "CTD",
                               "bystro_sampleMaf <= 0.01": "maf.01",
                               "bystro_sampleMaf <= 0.05": "maf.05",
                               "bystro_sampleMaf <= 0.1": "maf.1",
                               # "bystro_sampleMaf <= 0.2": "maf.2",
                               # "bystro_sampleMaf <= 0.3": "maf.3",
                               " AND (bystro_cadd>=10)": "cadd10",
                               " AND (bystro_cadd>=15)": "cadd15",
                               " AND (bystro_cadd>=20)": "cadd20",
                               " AND bystro_phastCons >= 0.4": "Cons.4",
                               " AND bystro_phastCons >= 0.5": "Cons.5",
                               " AND bystro_phastCons >= 0.6": "Cons.6",
                               " AND bystro_phyloP >= -1": "loPn1",
                               " AND bystro_phyloP >= 0": "loP0",
                               " AND bystro_phyloP >= 1": "loP1",
                               " AND (ccrs >= 95)": "ccrs95",
                               " AND (ccrs >= 90)": "ccrs90",
                               " AND (ccrs >= 85)": "ccrs85",
                               " AND (ccrs >= 80)": "ccrs80",
                               " AND (domain_limbr <= {})".format(domain_limbr_quartile_dict[40]): "dlimbr40",
                               " AND (domain_limbr <= {})".format(domain_limbr_quartile_dict[50]): "dlimbr50",
                               " AND (domain_limbr <= {})".format(domain_limbr_quartile_dict[60]): "dlimbr60",
                               " AND (domain_limbr <= {})".format(domain_limbr_quartile_dict[70]): "dlimbr70",
                               " AND (exone_limbr <= {})".format(exone_limbr_quartile_dict[40]): "elimbr40",
                               " AND (exone_limbr <= {})".format(exone_limbr_quartile_dict[50]): "elimbr50",
                               " AND (exone_limbr <= {})".format(exone_limbr_quartile_dict[60]): "elimbr60",
                               " AND (exone_limbr <= {})".format(exone_limbr_quartile_dict[70]): "elimbr70",
                               " AND (is_ccds = 1)": "ccds1"
                               }
    pp = Popen(["cd {}".format(path)], shell=True)
    pp.wait()
    annotator_list = ["",
                      " AND (annovar = 1)",
                      " AND (bystro = 1)",
                      " AND (dmis = 1)",
                      " AND (dsplicing = 1)",
                      " AND (spidex = 1)",
                      " AND (spliceAI = 1)",
                      " AND (vep = 1)",
                      " AND (annovar = 1 or bystro = 1 or vep = 1)",
                      " AND (dsplicing = 1 or spidex = 1 or spliceAI = 1)"]
    icounter = 1
    icounter2 = 0
    script_str = ''
    for phenotype in ["tof6", "CTD"]:
        for annotator in annotator_list:
            for freq in ["bystro_sampleMaf <= 0.05", "bystro_sampleMaf <= 0.1", "bystro_sampleMaf <= 0.01"]:
                for cadd in ["", " AND (bystro_cadd>=10)", " AND (bystro_cadd>=15)", " AND (bystro_cadd>=20)"]:
                    for ph in ["",
                               " AND bystro_phastCons >= 0.4",
                               " AND bystro_phastCons >= 0.5",
                               " AND bystro_phastCons >= 0.6",
                               " AND bystro_phyloP >= -1",
                               " AND bystro_phyloP >= 0",
                               " AND bystro_phyloP >= 1"]:
                        for regional_constraint in ["",
                                                    " AND (ccrs >= 95)",
                                                    " AND (ccrs >= 90)",
                                                    " AND (ccrs >= 85)",
                                                    " AND (ccrs >= 80)",
                                                    " AND (domain_limbr <= {})".format(domain_limbr_quartile_dict[40]),
                                                    " AND (domain_limbr <= {})".format(domain_limbr_quartile_dict[50]),
                                                    " AND (domain_limbr <= {})".format(domain_limbr_quartile_dict[60]),
                                                    " AND (domain_limbr <= {})".format(domain_limbr_quartile_dict[70]),
                                                    " AND (exone_limbr <= {})".format(exone_limbr_quartile_dict[40]),
                                                    " AND (exone_limbr <= {})".format(exone_limbr_quartile_dict[50]),
                                                    " AND (exone_limbr <= {})".format(exone_limbr_quartile_dict[60]),
                                                    " AND (exone_limbr <= {})".format(exone_limbr_quartile_dict[70])]:
                            for is_ccds in ["", " AND (is_ccds = 1)"]:
                                variance_restrict = "{0}{1}{2}{3}{4}{5}".format(freq, annotator, cadd, ph,
                                                                                regional_constraint, is_ccds)
                                name_element_list = [table_name_element_dict[phenotype],
                                                     table_name_element_dict[annotator],
                                                     table_name_element_dict[freq],
                                                     table_name_element_dict[cadd],
                                                     table_name_element_dict[ph],
                                                     table_name_element_dict[regional_constraint],
                                                     table_name_element_dict[is_ccds]]
                                output_name = "{0}.table" \
                                              "".format("_".join(filter(lambda x: len(x) > 0, name_element_list)))
                                # print output_name
                                # print("[{0}]\t[{1}]\t[{2}]\t[{3}]\t[{4}]\t[{5}]\t[{6}]\t[{7}]"
                                #       "".format(output_name, phenotype, annotator, freq, cadd, ph, regional_constraint,
                                #                 is_ccds))
                                output_name = os.path.join(path, output_name)
                                if script_str == "":
                                    script_str = """#!/bin/bash
#$ -j y
#$ -S /bin/bash
#$ -cwd
#$ -N fisher{6}
#$ -pe smp 1
#$ -l h_vmem=8G
#$ -m bes
# -q all.q
export LC_ALL=C
export MALLOC_ARENA_MAX=4
module load python/2.7.15/gcc.4.4.7
module load sqlite3/3.8.11/gcc.4.4.7
time=`date`
echo "==START $time =="
~/miniconda2/bin/python ~/wyj/.code/wgsa.py build_contingency_table_new {0} {1} sampleChdPhenotype '{2}' gene_table '' variance '{3}' {4} {5} {6}
echo {5} is done
time=`date`
echo == $time ==
""".format(database, phenotype, sample_restrict, variance_restrict, fai_in, output_name, icounter2)
                                else:
                                    script_str += "~/miniconda2/bin/python ~/wyj/.code/wgsa.py build_contingency_table_new {0} {1} " \
                                                  "sampleChdPhenotype '{2}' gene_table '' variance '{3}' {4} {5} {6}" \
                                                  "\necho {5} is done" \
                                                  "\ntime=`date`\necho == $time ==\n" \
                                                  "".format(database,
                                                            phenotype,
                                                            sample_restrict,
                                                            variance_restrict,
                                                            fai_in,
                                                            output_name,
                                                            icounter2)

                                if icounter == 30:
                                    shell_file_name = os.path.join(path, "qsub{}.sh".format(icounter2))
                                    with open(shell_file_name, "w") as fp:
                                        fp.write(script_str + "\ndate; echo \"==END==\"")
                                    while True:
                                        if os.access(shell_file_name, os.R_OK):
                                            break
                                        time.sleep(1)
                                    print("qsubing {} batch".format(icounter2 + 1))
                                    pp = Popen(["qsub {}".format(shell_file_name)], shell=True)
                                    pp.wait()
                                    script_str = ""
                                    icounter = 0
                                    icounter2 += 1
                                icounter += 1
    if len(script_str) > 0:
        shell_file_name = os.path.join(path, "qsub{}.sh".format(icounter2))
        with open(shell_file_name, "w") as fp:
            fp.write(script_str + "\ndate; echo \"==END==\"")
        while True:
            if os.access(shell_file_name, os.R_OK):
                break
            time.sleep(1)
        print("qsubing {} batch".format(icounter2 + 1))
        pp = Popen(["qsub {}".format(shell_file_name)], shell=True)
        pp.wait()

    print(
        "All the jobs has been submitted. Job number = {}\nIf any job has problem, kill it. And qsub the corresponding sh file".format(
            icounter2 + 1))


def analyze_geneset_fisher_test(gene_set_file, database, sample_restrict, fai_in, path):
    table_name_element_dict = {"": "",
                               " AND (annovar = 1)": "annovar",
                               " AND (bystro = 1)": "bystro",
                               " AND (dmis = 1)": "dmis",
                               " AND (dsplicing = 1)": "dsplicing",
                               " AND (spidex = 1)": "spidex",
                               " AND (spliceAI = 1)": "spliceAI",
                               " AND (vep = 1)": "vep",
                               " AND (annovar = 1 or bystro = 1 or vep = 1)": "LOF",
                               " AND (dsplicing = 1 or spidex = 1 or spliceAI = 1)": "dsplicing.all",
                               "tof6": "tof",
                               "CTD": "CTD",
                               "bystro_sampleMaf <= 0.01": "maf.01",
                               "bystro_sampleMaf <= 0.05": "maf.05",
                               "bystro_sampleMaf <= 0.1": "maf.1",
                               # "bystro_sampleMaf <= 0.2": "maf.2",
                               # "bystro_sampleMaf <= 0.3": "maf.3",
                               " AND (bystro_cadd>=10)": "cadd10",
                               " AND (bystro_cadd>=15)": "cadd15",
                               " AND (bystro_cadd>=20)": "cadd20",
                               " AND bystro_phastCons >= 0.4": "Cons.4",
                               " AND bystro_phastCons >= 0.5": "Cons.5",
                               " AND bystro_phastCons >= 0.6": "Cons.6",
                               " AND bystro_phyloP >= -1": "loPn1",
                               " AND bystro_phyloP >= 0": "loP0",
                               " AND bystro_phyloP >= 1": "loP1",
                               " AND (ccrs >= 95)": "ccrs95",
                               " AND (ccrs >= 90)": "ccrs90",
                               " AND (ccrs >= 85)": "ccrs85",
                               " AND (ccrs >= 80)": "ccrs80",
                               " AND (domain_limbr <= {})".format(domain_limbr_quartile_dict[40]): "dlimbr40",
                               " AND (domain_limbr <= {})".format(domain_limbr_quartile_dict[50]): "dlimbr50",
                               " AND (domain_limbr <= {})".format(domain_limbr_quartile_dict[60]): "dlimbr60",
                               " AND (domain_limbr <= {})".format(domain_limbr_quartile_dict[70]): "dlimbr70",
                               " AND (exone_limbr <= {})".format(exone_limbr_quartile_dict[40]): "elimbr40",
                               " AND (exone_limbr <= {})".format(exone_limbr_quartile_dict[50]): "elimbr50",
                               " AND (exone_limbr <= {})".format(exone_limbr_quartile_dict[60]): "elimbr60",
                               " AND (exone_limbr <= {})".format(exone_limbr_quartile_dict[70]): "elimbr70",
                               " AND (is_ccds = 1)": "ccds1"
                               }
    pp = Popen(["cd {}".format(path)], shell=True)
    pp.wait()
    annotator_list = ["",
                      " AND (annovar = 1)",
                      " AND (bystro = 1)",
                      " AND (dmis = 1)",
                      " AND (dsplicing = 1)",
                      " AND (spidex = 1)",
                      " AND (spliceAI = 1)",
                      " AND (vep = 1)",
                      " AND (annovar = 1 or bystro = 1 or vep = 1)",
                      " AND (dsplicing = 1 or spidex = 1 or spliceAI = 1)"]
    icounter = 1
    icounter2 = 0
    script_str = ''
    for phenotype in ["tof6", "CTD"]:
        for annotator in annotator_list:
            for freq in ["bystro_sampleMaf <= 0.05", "bystro_sampleMaf <= 0.1", "bystro_sampleMaf <= 0.01"]:
                for cadd in ["", " AND (bystro_cadd>=10)", " AND (bystro_cadd>=15)", " AND (bystro_cadd>=20)"]:
                    for ph in ["",
                               " AND bystro_phastCons >= 0.4",
                               " AND bystro_phastCons >= 0.5",
                               " AND bystro_phastCons >= 0.6",
                               " AND bystro_phyloP >= -1",
                               " AND bystro_phyloP >= 0",
                               " AND bystro_phyloP >= 1"]:
                        for regional_constraint in ["",
                                                    " AND (ccrs >= 95)",
                                                    " AND (ccrs >= 90)",
                                                    " AND (ccrs >= 85)",
                                                    " AND (ccrs >= 80)",
                                                    " AND (domain_limbr <= {})".format(domain_limbr_quartile_dict[40]),
                                                    " AND (domain_limbr <= {})".format(domain_limbr_quartile_dict[50]),
                                                    " AND (domain_limbr <= {})".format(domain_limbr_quartile_dict[60]),
                                                    " AND (domain_limbr <= {})".format(domain_limbr_quartile_dict[70]),
                                                    " AND (exone_limbr <= {})".format(exone_limbr_quartile_dict[40]),
                                                    " AND (exone_limbr <= {})".format(exone_limbr_quartile_dict[50]),
                                                    " AND (exone_limbr <= {})".format(exone_limbr_quartile_dict[60]),
                                                    " AND (exone_limbr <= {})".format(exone_limbr_quartile_dict[70])]:
                            for is_ccds in ["", " AND (is_ccds = 1)"]:
                                variance_restrict = "{0}{1}{2}{3}{4}{5}".format(freq, annotator, cadd, ph,
                                                                                regional_constraint, is_ccds)
                                name_element_list = [table_name_element_dict[phenotype],
                                                     table_name_element_dict[annotator],
                                                     table_name_element_dict[freq],
                                                     table_name_element_dict[cadd],
                                                     table_name_element_dict[ph],
                                                     table_name_element_dict[regional_constraint],
                                                     table_name_element_dict[is_ccds]]
                                output_name = "{0}.table" \
                                              "".format("_".join(filter(lambda x: len(x) > 0, name_element_list)))
                                # print output_name
                                # print("[{0}]\t[{1}]\t[{2}]\t[{3}]\t[{4}]\t[{5}]\t[{6}]\t[{7}]"
                                #       "".format(output_name, phenotype, annotator, freq, cadd, ph, regional_constraint,
                                #                 is_ccds))
                                output_name = os.path.join(path, output_name)
                                if script_str == "":
                                    script_str = """#!/bin/bash
#$ -j y
#$ -S /bin/bash
#$ -cwd
#$ -N gs{7}
#$ -pe smp 1
#$ -l h_vmem=8G
#$ -m bes
# -q all.q
export LC_ALL=C
export MALLOC_ARENA_MAX=4
module load python/2.7.15/gcc.4.4.7
module load sqlite3/3.8.11/gcc.4.4.7
time=`date`
echo "==START $time =="
~/miniconda2/bin/python ~/wyj/.code/wgsa.py gene_set_fisher_test {0} {1} sampleChdPhenotype '{2}' {3} variance '{4}' gene_table {5} {6}
echo {6} is done
time=`date`
echo == $time ==
""".format(gene_set_file, database, sample_restrict, phenotype, variance_restrict, fai_in, output_name, icounter2)
                                    # gene_set_fisher_test(gene_set_file, db_file, sample_table, sample_restrict, phenotype_in,
                                    #                          variance_table, variance_restrict, gene_table, fai_in, output
                                else:
                                    script_str += "~/miniconda2/bin/python ~/wyj/.code/wgsa.py gene_set_fisher_test {0} {1} " \
                                                  "sampleChdPhenotype '{2}' {3} variance '{4}' gene_table {5} {6}" \
                                                  "\necho {5} is done" \
                                                  "\ntime=`date`\necho == $time ==\n" \
                                                  "".format(gene_set_file,
                                                            database,
                                                            sample_restrict,
                                                            phenotype,
                                                            variance_restrict,
                                                            fai_in,
                                                            output_name)

                                if icounter == 30:
                                    shell_file_name = os.path.join(path, "qsub{}.sh".format(icounter2))
                                    with open(shell_file_name, "w") as fp:
                                        fp.write(script_str + "\ndate; echo \"==END==\"")
                                    while True:
                                        if os.access(shell_file_name, os.R_OK):
                                            break
                                        time.sleep(1)
                                    print("qsubing {} batch".format(icounter2))
                                    pp = Popen(["qsub {}".format(shell_file_name)], shell=True)
                                    pp.wait()
                                    script_str = ""
                                    icounter = 0
                                    icounter2 += 1
                                icounter += 1
    if len(script_str) > 0:
        shell_file_name = os.path.join(path, "qsub{}.sh".format(icounter2))
        with open(shell_file_name, "w") as fp:
            fp.write(script_str + "\ndate; echo \"==END==\"")
        while True:
            if os.access(shell_file_name, os.R_OK):
                break
            time.sleep(1)
        print("qsubing {} batch".format(icounter2))
        pp = Popen(["qsub {}".format(shell_file_name)], shell=True)
        pp.wait()

    print(
        "All the jobs has been submitted. Job number = {}\nIf any job has problem, kill it. And qsub the corresponding sh file".format(
            icounter2 + 1))


def fisher_test_synonymous(database, sample_restrict, fai_in, path):
    table_name_element_dict = {"": "",
                               " AND (annovar = 1)": "annovar",
                               " AND (bystro = 1)": "bystro",
                               " AND (dmis = 1)": "dmis",
                               " AND (dsplicing = 1)": "dsplicing",
                               " AND (spidex = 1)": "spidex",
                               " AND (spliceAI = 1)": "spliceAI",
                               " AND (vep = 1)": "vep",
                               " AND (annovar = 1 or bystro = 1 or vep = 1)": "LOF",
                               " AND (dsplicing = 1 or spidex = 1 or spliceAI = 1)": "dsplicing.all",
                               "tof6": "tof",
                               "CTD": "CTD",
                               "bystro_sampleMaf <= 0.01": "maf.01",
                               "bystro_sampleMaf <= 0.05": "maf.05",
                               "bystro_sampleMaf <= 0.1": "maf.1",
                               # "bystro_sampleMaf <= 0.2": "maf.2",
                               # "bystro_sampleMaf <= 0.3": "maf.3",
                               " AND (bystro_cadd>=10)": "cadd10",
                               " AND (bystro_cadd>=15)": "cadd15",
                               " AND (bystro_cadd>=20)": "cadd20",
                               " AND (bystro_cadd<=5)": "caddls5",
                               " AND (bystro_cadd>5 AND bystro_cadd<=10)": "cadd5_10",
                               " AND (bystro_cadd>10 AND bystro_cadd<=15)": "cadd10_15",
                               " AND bystro_phastCons = 0": "Conseq0",
                               " AND bystro_phastCons > 0 AND bystro_phastCons <=0.75": "Con0_.75",
                               " AND bystro_phastCons > 0.75 AND bystro_phastCons<=1": "Con.75_1",
                               " AND bystro_phastCons >= 0.4": "Cons.4",
                               " AND bystro_phastCons >= 0.5": "Cons.5",
                               " AND bystro_phastCons >= 0.6": "Cons.6",
                               " AND bystro_phyloP >= -1": "loPn1",
                               " AND bystro_phyloP >= 0": "loP0",
                               " AND bystro_phyloP >= 1": "loP1",
                               " AND bystro_phyloP <= -2": "loPlsn2",
                               " AND bystro_phyloP > -2 AND bystro_phyloP <= 0": "loPn2_0",
                               " AND bystro_phyloP > 0 AND bystro_phyloP <= 2": "loP0_2",
                               " AND bystro_phyloP >= 2": "loP2",
                               " AND (ccrs >= 95)": "ccrs95",
                               " AND (ccrs >= 90)": "ccrs90",
                               " AND (ccrs >= 85)": "ccrs85",
                               " AND (ccrs >= 80)": "ccrs80",
                               " AND (domain_limbr <= {})".format(domain_limbr_quartile_dict[40]): "dlimbr40",
                               " AND (domain_limbr <= {})".format(domain_limbr_quartile_dict[50]): "dlimbr50",
                               " AND (domain_limbr <= {})".format(domain_limbr_quartile_dict[60]): "dlimbr60",
                               " AND (domain_limbr <= {})".format(domain_limbr_quartile_dict[70]): "dlimbr70",
                               " AND (exone_limbr <= {})".format(exone_limbr_quartile_dict[40]): "elimbr40",
                               " AND (exone_limbr <= {})".format(exone_limbr_quartile_dict[50]): "elimbr50",
                               " AND (exone_limbr <= {})".format(exone_limbr_quartile_dict[60]): "elimbr60",
                               " AND (exone_limbr <= {})".format(exone_limbr_quartile_dict[70]): "elimbr70",
                               " AND (is_ccds = 1)": "ccds1"
                               }
    pp = Popen(["cd {}".format(path)], shell=True)
    pp.wait()
    annotator_list = ["",
                      " AND (annovar = 1)",
                      " AND (bystro = 1)",
                      " AND (vep = 1)"]
    icounter = 1
    icounter2 = 0
    script_str = ''
    for phenotype in ["tof6", "CTD"]:
        for annotator in annotator_list:
            for freq in ["bystro_sampleMaf <= 0.05", "bystro_sampleMaf <= 0.1", "bystro_sampleMaf <= 0.01"]:
                for cadd in ["", " AND (bystro_cadd<=5)", " AND (bystro_cadd>5 AND bystro_cadd<=10)",
                             " AND (bystro_cadd>10 AND bystro_cadd<=15)", " AND (bystro_cadd>=15)"]:
                    for ph in ["",
                               " AND bystro_phastCons = 0",
                               " AND bystro_phastCons > 0 AND bystro_phastCons <=0.75",
                               " AND bystro_phastCons > 0.75 AND bystro_phastCons<=1",
                               " AND bystro_phyloP <= -2",
                               " AND bystro_phyloP > -2 AND bystro_phyloP <= 0",
                               " AND bystro_phyloP > 0 AND bystro_phyloP <= 2",
                               " AND bystro_phyloP >= 2"]:
                        for regional_constraint in ["",
                                                    " AND (ccrs >= 95)",
                                                    " AND (ccrs >= 90)",
                                                    " AND (ccrs >= 85)",
                                                    " AND (ccrs >= 80)",
                                                    " AND (domain_limbr <= {})".format(domain_limbr_quartile_dict[40]),
                                                    " AND (domain_limbr <= {})".format(domain_limbr_quartile_dict[50]),
                                                    " AND (domain_limbr <= {})".format(domain_limbr_quartile_dict[60]),
                                                    " AND (domain_limbr <= {})".format(domain_limbr_quartile_dict[70]),
                                                    " AND (exone_limbr <= {})".format(exone_limbr_quartile_dict[40]),
                                                    " AND (exone_limbr <= {})".format(exone_limbr_quartile_dict[50]),
                                                    " AND (exone_limbr <= {})".format(exone_limbr_quartile_dict[60]),
                                                    " AND (exone_limbr <= {})".format(exone_limbr_quartile_dict[70])]:
                            for is_ccds in ["", " AND (is_ccds = 1)"]:
                                variance_restrict = "{0}{1}{2}{3}{4}{5}".format(freq, annotator, cadd, ph,
                                                                                regional_constraint, is_ccds)
                                name_element_list = [table_name_element_dict[phenotype],
                                                     table_name_element_dict[annotator],
                                                     table_name_element_dict[freq],
                                                     table_name_element_dict[cadd],
                                                     table_name_element_dict[ph],
                                                     table_name_element_dict[regional_constraint],
                                                     table_name_element_dict[is_ccds]]
                                output_name = "{0}.table" \
                                              "".format("_".join(filter(lambda x: len(x) > 0, name_element_list)))
                                # print output_name
                                # print("[{0}]\t[{1}]\t[{2}]\t[{3}]\t[{4}]\t[{5}]\t[{6}]\t[{7}]"
                                #       "".format(output_name, phenotype, annotator, freq, cadd, ph, regional_constraint,
                                #                 is_ccds))
                                output_name = os.path.join(path, output_name)
                                if script_str == "":
                                    script_str = """#!/bin/bash
#$ -j y
#$ -S /bin/bash
#$ -cwd
#$ -N fisher{6}
#$ -pe smp 1
#$ -l h_vmem=8G
#$ -m bes
# -q all.q
export LC_ALL=C
export MALLOC_ARENA_MAX=4
module load python/2.7.15/gcc.4.4.7
module load sqlite3/3.8.11/gcc.4.4.7
time=`date`
echo "==START $time =="
~/miniconda2/bin/python ~/wyj/.code/wgsa.py build_contingency_table_new {0} {1} sampleChdPhenotype '{2}' gene_table '' synonymous_snp '{3}' {4} {5} {6}
echo {5} is done
time=`date`
echo == $time ==
    """.format(database, phenotype, sample_restrict, variance_restrict, fai_in, output_name, icounter2)

                                    # def build_contingency_table_new(db_file, phenotype,
                                    #                                 sample_table_name, sample_restrict,
                                    #                                 gene_table_name, gene_restrict,
                                    #                                 variance_table, variance_restrict,
                                    #                                 fai_in, output,
                                    #                                 job_id, should_log=True)
                                else:
                                    script_str += "~/miniconda2/bin/python ~/wyj/.code/wgsa.py build_contingency_table_new {0} {1} " \
                                                  "sampleChdPhenotype '{2}' gene_table '' synonymous_snp '{3}' {4} {5} {6}" \
                                                  "\necho {5} is done" \
                                                  "\ntime=`date`\necho == $time ==\n" \
                                                  "".format(database,
                                                            phenotype,
                                                            sample_restrict,
                                                            variance_restrict,
                                                            fai_in,
                                                            output_name,
                                                            icounter2)

                                if icounter == 30:
                                    shell_file_name = os.path.join(path, "qsub{}.sh".format(icounter2))
                                    with open(shell_file_name, "w") as fp:
                                        fp.write(script_str + "\ndate; echo \"==END==\"")
                                    while True:
                                        if os.access(shell_file_name, os.R_OK):
                                            break
                                        time.sleep(1)
                                    print("qsubing {} batch".format(icounter2 + 1))
                                    pp = Popen(["qsub {}".format(shell_file_name)], shell=True)
                                    pp.wait()
                                    script_str = ""
                                    icounter = 0
                                    icounter2 += 1
                                icounter += 1
    if len(script_str) > 0:
        shell_file_name = os.path.join(path, "qsub{}.sh".format(icounter2))
        with open(shell_file_name, "w") as fp:
            fp.write(script_str + "\ndate; echo \"==END==\"")
        while True:
            if os.access(shell_file_name, os.R_OK):
                break
            time.sleep(1)
        print("qsubing {} batch".format(icounter2 + 1))
        pp = Popen(["qsub {}".format(shell_file_name)], shell=True)
        pp.wait()

    print(
        "All the jobs has been submitted. Job number = {}\nIf any job has problem, kill it. And qsub the corresponding sh file".format(
            icounter2 + 1))


def prepare_liftover_bed(org_bed_in, value_col_num, bed_out):
    """
         hg19     ------------>        hg38
    chr1  100  103                  chr1  200  205

    It is wrong to convert directly like above.
    The right way to conver should be:

                    prepare bed                   lift over
          hg19     -------------->    hg19     --------------->        hg38
    chr1  100  103               chr1  100  101                  chr1  200  201
                                 chr1  101  102                  chr1  204  205
    @param org_bed_in:
    @param bed_out:
    @return:
    """
    with open(org_bed_in, "r") as fp_in, open(bed_out, "w") as fp_out:
        fp_out.write("#chr\tstart\tend\tvalue\n")
        while True:
            data_line = fp_in.readline()
            if not data_line:
                break
            if data_line.startswith("#"):
                # fp_out.write(data_line)
                continue
            data_list = data_line.strip().split()
            chr = data_list[0] if data_list[0].startswith("chr") else "chr" + data_list[0]
            start = int(data_list[1])
            end = int(data_list[2])
            value = data_list[int(value_col_num) - 1]
            for pos in xrange(start, end, 1):
                fp_out.write("{0}\t{1}\t{2}\t{3}\n".format(chr, pos, pos + 1, value))


def gether_bed_region(input, output):
    with open(input, "r") as fp_in, open(output, "w") as fp_out:
        region_value = ""
        region_start = ""
        region_end = ""
        region_chrom = ""
        while True:
            curr_data_line = fp_in.readline()
            if not curr_data_line:
                fp_out.write("{0}\t{1}\t{2}\t{3}\n".format(region_chrom, region_start, region_end, region_value))
                break
            if curr_data_line.startswith("#"):
                fp_out.write(curr_data_line)
                continue
            curr_chrom, curr_start, curr_end, curr_value = curr_data_line.strip().split("\t")
            if region_chrom == "":
                region_chrom = curr_chrom
                region_start = curr_start
                region_end = curr_end
                region_value = curr_value
                continue
            if curr_start == region_end and curr_value == region_value and curr_chrom == region_chrom:
                region_end = curr_end
                continue
            if curr_value != region_value or curr_start != region_end or curr_chrom != region_chrom:
                fp_out.write("{0}\t{1}\t{2}\t{3}\n".format(region_chrom, region_start, region_end, region_value))
                region_value = curr_value
                region_chrom = curr_chrom
                region_start = curr_start
                region_end = curr_end
                continue


def check_gene_overlap(database, gene_table, fai_in, output):
    def overlap_bp_num(gene_data_line1, gene_data_line2, chr2offset_dict):
        gene_id1, gene_name1, chr11, start11, end11, chr12, start12, end12 = gene_data_line1
        gene_id2, gene_name2, chr21, start21, end21, chr22, start22, end22 = gene_data_line2
        gene_id1 = str(gene_id1)
        gene_id2 = str(gene_id2)
        gene_name1 = str(gene_name1)
        gene_name2 = str(gene_name2)
        if chr11 is not None:
            chr11 = str(chr11)
        if chr12 is not None:
            chr12 = str(chr12)
        if chr21 is not None:
            chr21 = str(chr21)
        if chr22 is not None:
            chr22 = str(chr22)

        region1_list = []
        region2_list = []
        if chr11 is not None:
            region1_list.append([chr_pos2absolute_pos(chr11, start11, chr2offset_dict),
                                 chr_pos2absolute_pos(chr11, end11, chr2offset_dict)])
        if chr12 is not None:
            region1_list.append([chr_pos2absolute_pos(chr12, start12, chr2offset_dict),
                                 chr_pos2absolute_pos(chr12, end12, chr2offset_dict)])
        if chr21 is not None:
            region2_list.append([chr_pos2absolute_pos(chr21, start21, chr2offset_dict),
                                 chr_pos2absolute_pos(chr21, end21, chr2offset_dict)])
        if chr22 is not None:
            region2_list.append([chr_pos2absolute_pos(chr22, start22, chr2offset_dict),
                                 chr_pos2absolute_pos(chr22, end22, chr2offset_dict)])
        if gene_id1 == "ENSG00000204961" and gene_id2 == "ENSG00000250120":
            print(region1_list)
            print(region2_list)
        ret = 0
        for abp11, abp12 in region1_list:
            for abp21, abp22 in region2_list:
                if abp21 <= abp11 <= abp22 and abp12 > abp22:  # 3
                    ret += abp22 - abp11
                if abp21 <= abp12 <= abp22 and abp11 < abp21:  # 1
                    ret += abp12 - abp21
                if abp21 <= abp12 <= abp22 and abp11 >= abp21:  # 2
                    ret += abp12 - abp11
                if abp11 <= abp21 and abp12 >= abp22:
                    ret += abp22 - abp21
        return ret

    chr2offset_dict = parse_fai(fai_in)
    print(chr2offset_dict)
    conn = sqlite3.connect(database)
    cursor = conn.cursor()
    cmd_str = "select gene_id, gene_name, chr, start_pos, end_pos, chr2, start_pos2, end_pos2 from gene_table"
    cursor.execute(cmd_str)
    gene_data = cursor.fetchall()
    with open(output, "w") as fp:
        fp.write("#id1\tname1\tid2\tname2\toverlap_bp_num\n")
        for i in xrange(len(gene_data) - 1):
            # print("handling {} gene".format(i + 1))
            for j in xrange(i + 1, len(gene_data), 1):
                overlap_len = overlap_bp_num(gene_data[i], gene_data[j], chr2offset_dict)
                if overlap_len > 0:
                    fp.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(gene_data[i][0], gene_data[i][1],
                                                                gene_data[j][0], gene_data[j][1],
                                                                overlap_len))


def bystro_build_ref_alt(type, org_alt, org_ref, chrom, vcf_pos, ref_source):
    assert type in ["SNP", "DEL", "INS"]
    if type == "SNP":
        return [org_ref, org_alt]
    elif type == "DEL":
        ref = ref_source.send([chrom, int(vcf_pos), abs(int(org_alt)) + 1])
        alt = ref_source.send([chrom, int(vcf_pos), 1])
        return [ref, alt]
    elif type == "INS":
        assert org_alt.startswith("+")
        org_alt = org_alt[1:]
        ref = ref_source.send([chrom, int(vcf_pos), 1])
        alt = ref + org_alt
        return [ref, alt]


def db_add_bystro_gene_name(db_file, bystro_tsv, ref_in):
    tsv_key2value_dict = {}

    print("loading ref...")
    ref_source = generator_seq_from_ref(ref_in)
    ref_source.next()
    # print(ref_source.send(["chr1", 86023006, 1]))
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()
    print("loading variance in database...")
    cursor.execute("select id, chr, pos, ref, alt from variance")
    db_data = cursor.fetchall()
    db_id2key_dict = dict(
        zip([int(i[0]) for i in db_data], ["{0}\t{1}\t{2}\t{3}".format(i[1], i[2], i[3], i[4]) for i in db_data]))
    print("there are {} variance in database".format(len(db_id2key_dict)))

    with open(bystro_tsv, "r") as fp:
        while True:
            data_line = fp.readline()
            if not data_line:
                break
            if data_line.startswith("chrom\tpos"):
                continue
            data_list = data_line.strip().split("\t")
            # print(len(data_list))
            if len(data_list) != 114:
                continue
            chrom = data_list[0]
            vcf_pos = data_list[15]
            bystro_ref = data_list[16]
            bystro_alt = data_list[4]
            type = data_list[2]
            ref, alt = bystro_build_ref_alt(type, bystro_alt, bystro_ref, chrom, vcf_pos, ref_source)
            # cmd_str = "select id from variance where chr='{0}' and pos={1} and ref='{2}' and alt='{3}'" \
            #           "".format(chrom[3:], vcf_pos, ref,alt)
            # cursor.execute(cmd_str)
            # result = cursor.fetchall()
            # if len(result) != 1:
            #     print(cmd_str)
            #     print(result)
            # exit(0)
            # print("{0}\t{1}".format(chrom, pos))
            # if len(data_list)< 113:
            #     print(data_line)
            final_function = data_list[112]
            gene_name = data_list[113]
            key = "{0}\t{1}\t{2}\t{3}".format(chrom[3:], vcf_pos, ref, alt)
            tsv_key2value_dict[key] = [final_function, gene_name]
    print("there are {} records in tsv file".format(len(tsv_key2value_dict)))
    db_add_col(cursor, "variance", "bystro_final_function", "varchr(256)")
    db_add_col(cursor, "variance", "bystro_gene_name", "varchr(256)")
    sql_list = []
    icounter = 0
    print("collecting sqls... 0%")
    for id in db_id2key_dict:
        db_key = db_id2key_dict[id]
        if db_key not in tsv_key2value_dict:
            sql_list.append("UPDATE variance SET bystro = '0' WHERE id = {0}".format(id))
            sql_list.append("UPDATE variance SET bystro_final_function = NULL WHERE id = {0}".format(id))
            sql_list.append("UPDATE variance SET bystro_gene_name = NULL WHERE id = {0}".format(id))
        else:
            sql_list.append("UPDATE variance SET bystro = '1' WHERE id = {0}".format(id))
            final_function, gene_name = tsv_key2value_dict[db_key]
            sql_list.append(
                "UPDATE variance SET bystro_final_function = '{1}' WHERE id = {0}".format(id, final_function))
            sql_list.append("UPDATE variance SET bystro_gene_name = '{1}' WHERE id = {0}".format(id, gene_name))
        icounter += 1
        if icounter % 10000 == 0:
            print("collecting sqls... {:.2%}".format(float(icounter) / len(db_id2key_dict)))
    print("executing sqls... 0%")
    icounter = 0
    for sql in sql_list:
        cursor.execute(sql)
        icounter += 1
        if icounter % 10000 == 0:
            print("executing sqls... {:.2%}".format(float(icounter) / len(sql_list)))
    cursor.close()
    conn.commit()
    conn.close()
    print("done")


def db_add_annovar_gene_name(db_file, selected_annovar_file, ref_in):
    def build_ref_alt_from_annovar_result(chrom, pos, annovar_ref, annovar_alt, ref_source):
        if annovar_alt == "-":
            delta = ref_source.send([chrom, pos, 1])
            return [delta + annovar_ref, delta]
        elif annovar_ref == "-":
            delta = ref_source.send([chrom, pos, 1])
            return [delta, delta + annovar_alt]
        elif len(annovar_alt) == len(annovar_ref) == 1:
            return [annovar_ref, annovar_alt]
        else:
            RuntimeError("unexpected ref ({0}) and alt ({1})".format(annovar_ref, annovar_alt))

    logging.basicConfig(filename="db_add_annovar_gene_name.log", level=logging.DEBUG, format=log_format, filemode="w")
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()
    print("loading ref...")
    ref_source = generator_seq_from_ref(ref_in)
    ref_source.next()
    print("building id2info dict ...")
    # print(ref_source.send(["chr1", 86023006, 1]))
    sql_list = []
    id2info_dict = {}
    icounter = 0
    with open(selected_annovar_file, "r") as fp:
        while True:
            data_line = fp.readline()
            if not data_line:
                break
            if data_line.startswith("Chr\tStart\t"):
                continue
            if len(data_line.strip()) == 0:
                continue
            data_list = data_line.strip().split("\t")
            chrom = data_list[0]
            pos = data_list[20]
            annovar_ref = data_list[3].upper()
            annovar_alt = data_list[4].upper()
            function = data_list[21]
            detail = data_list[22]
            gene_name = data_list[6]
            ref, alt = build_ref_alt_from_annovar_result(chrom, int(pos), annovar_ref, annovar_alt, ref_source)
            cmd_str = "select id from variance where chr='{0}' and pos='{1}' and ref='{2}' and alt='{3}'" \
                      "".format(chrom[3:], pos, ref, alt)
            cursor.execute(cmd_str)
            sql_result = cursor.fetchall()
            icounter += 1
            if icounter % 100 == 0:
                print("building id2info dict: {0}".format(icounter))
            if len(sql_result) == 0:
                continue
            elif len(sql_result) > 1:
                logging.debug("got {0} result: ".format(len(sql_result)) + cmd_str)
                continue
            else:
                id2info_dict[sql_result[0][0]] = [function, detail, gene_name]
    cursor.execute("select id from variance")
    id_list = [i[0] for i in cursor.fetchall()]

    db_add_col(cursor, "variance", "annovar_gene_name", "varchr(256)")
    db_add_col(cursor, "variance", "annovar_function", "varchr(256)")
    db_add_col(cursor, "variance", "annovar_function_detail", "varchr(256)")
    print("collecting sqls...")
    icounter = 0
    for id in id_list:
        if id not in id2info_dict:
            sql_list.append("UPDATE variance SET annovar = '0' WHERE id = {0}".format(id))
            sql_list.append("UPDATE variance SET annovar_gene_name = NULL WHERE id = {0}".format(id))
            sql_list.append("UPDATE variance SET annovar_function = NULL WHERE id = {0}".format(id))
            sql_list.append("UPDATE variance SET annovar_function_detail = NULL WHERE id = {0}".format(id))
        else:
            annovar_function, annovar_detail, annovar_gene_name = id2info_dict[id]
            sql_list.append("UPDATE variance SET annovar = '1' WHERE id = {0}".format(id))
            sql_list.append("UPDATE variance SET annovar_gene_name = '{1}' WHERE id = {0}"
                            "".format(id, annovar_gene_name))
            sql_list.append("UPDATE variance SET annovar_function = '{1}' WHERE id = {0}"
                            "".format(id, annovar_function))
            sql_list.append("UPDATE variance SET annovar_function_detail = '{1}' WHERE id = {0}"
                            "".format(id, annovar_detail))
        icounter += 1
        if icounter % 1000 == 0:
            print("collecting sqls: {:.2%}".format(icounter / float(len(id_list))))
    print("executing sqls...")
    icounter = 0
    for sql in sql_list:
        cursor.execute(sql)
        icounter += 1
        if icounter % 1000:
            print("executing sqls: {:.2%}".format(float(icounter) / len(sql_list)))
    cursor.close()
    conn.commit()
    conn.close()
    print("done")


def db_add_vep_gene_name(db_file, selected_vep_file):
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()
    id2info_dict = {}
    with open(selected_vep_file, "r") as fp:
        while True:
            data_line = fp.readline()
            if not data_line:
                break
            if data_line.startswith("#chr\t"):
                continue
            if len(data_line.strip()) == 0:
                continue
            chrom, pos, alt, gene_id, feature, consequence = data_line.strip().split("\t")
            alt = alt.upper()
            sql_cmd = "select id from variance where chr='{0}' and pos='{1}' and alt='{2}'".format(chrom, pos, alt)
            cursor.execute(sql_cmd)
            sql_result = cursor.fetchall()
            if len(sql_result) != 1:
                print("Variance (chr={0} pos={1} alt={2}) doesn't exist in variance table. It might be deleted before."
                      "".format(chrom, pos, alt))
                continue
            id2info_dict[sql_result[0][0]] = [gene_id, feature, consequence]
    cursor.execute("select id from variance")
    id_list = [i[0] for i in cursor.fetchall()]

    db_add_col(cursor, "variance", "vep_gene_id", "varchr(256)")
    db_add_col(cursor, "variance", "vep_feature", "varchr(256)")
    db_add_col(cursor, "variance", "vep_consequence", "varchr(256)")

    print("collecting sqls")
    sql_list = []
    for id in id_list:
        if id in id2info_dict:
            gene_id, feature, consequence = id2info_dict[id]
            sql_list.append("UPDATE variance SET vep = '1' WHERE id = {0}".format(id))
            sql_list.append("UPDATE variance SET vep_gene_id = '{1}' WHERE id = {0}".format(id, gene_id))
            sql_list.append("UPDATE variance SET vep_feature = '{1}' WHERE id = {0}".format(id, feature))
            sql_list.append("UPDATE variance SET vep_consequence = '{1}' WHERE id = {0}".format(id, consequence))
        else:
            sql_list.append("UPDATE variance SET vep = '0' WHERE id = {0}".format(id))
            sql_list.append("UPDATE variance SET vep_gene_id = NULL WHERE id = {0}".format(id))
            sql_list.append("UPDATE variance SET vep_feature = NULL WHERE id = {0}".format(id))
            sql_list.append("UPDATE variance SET vep_consequence = NULL WHERE id = {0}".format(id))
    print("executing sql ...")
    icounter = 0
    for sql in sql_list:
        cursor.execute(sql)
        icounter += 1
        if icounter % 1000 == 0:
            print("executing sql ...".format(float(icounter) / len(sql_list)))
    cursor.close()
    conn.commit()
    conn.close()
    print("done")


def db_add_spliceAI_gene_name(db_file, reduced_spliceAI_file):
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()
    id2info_dict = {}
    print("building id2info dict ...")
    with open(reduced_spliceAI_file, "r") as fp:
        while True:
            data_line = fp.readline()
            if not data_line:
                break
            if len(data_line) == 0:
                continue
            chrom, pos, ref, alt, gene_name = data_line.strip().split("\t")
            cursor.execute("select id from variance where chr='{0}' and pos='{1}' and ref='{2}' and alt='{3}'"
                           "".format(chrom, pos, ref, alt))
            sql_result = cursor.fetchall()
            if len(sql_result) != 1:
                print("Variance {} does't exist in data base. It might be deleted before.".format(data_line.strip()))
                continue
            id2info_dict[sql_result[0][0]] = gene_name
    cursor.execute("select id from variance")
    id_list = [i[0] for i in cursor.fetchall()]

    db_add_col(cursor, "variance", "spliceAI_gene_name", "varchr(256)")
    sql_list = []
    print("collecting sql ...")
    for id in id_list:
        if id in id2info_dict:
            sql_list.append("UPDATE variance SET spliceAI = '1' WHERE id = {0}".format(id))
            sql_list.append("UPDATE variance SET spliceAI_gene_name = '{1}' WHERE id = {0}"
                            "".format(id, id2info_dict[id]))
        else:
            sql_list.append("UPDATE variance SET spliceAI = '0' WHERE id = {0}".format(id))
            sql_list.append("UPDATE variance SET spliceAI_gene_name = NULL WHERE id = {0}".format(id))
    print("executing sql ...")
    for sql in sql_list:
        cursor.execute(sql)
    cursor.close()
    conn.commit()
    conn.close()
    print("done")


def db_add_dmis_gene_name(db_file, dmis_file):
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()
    id2info_dict = {}
    print("building id2info dict ...")
    icounter = 0
    with open(dmis_file, "r") as fp:
        while True:
            data_line = fp.readline()
            if not data_line:
                break
            if len(data_line) == 0 or data_line.startswith("#"):
                continue
            chrom, pos, ref, alt, gene_name, Ensembl_geneid, HGVSc, HGVSp = data_line.strip().split("\t")
            cursor.execute("select id from variance where chr='{0}' and pos='{1}' and ref='{2}' and alt='{3}'"
                           "".format(chrom, pos, ref, alt))
            sql_result = cursor.fetchall()
            icounter += 1
            if icounter % 100 == 0:
                print("building id2info dict ... handled {0} lines".format(icounter))
            if len(sql_result) == 0:
                print("Variance {0} does't exist in data base. It might be deleted before.".format(data_line.strip()))
                continue
            elif len(sql_result) > 1:
                print("more than one variance match ({0})".format(data_line.strip()))
                continue
            id2info_dict[sql_result[0][0]] = [gene_name, HGVSc, HGVSp]
    cursor.execute("select id from variance")
    id_list = [i[0] for i in cursor.fetchall()]

    db_add_col(cursor, "variance", "dmis_gene_name", "varchr(256)")
    db_add_col(cursor, "variance", "dmis_HGVSc", "varchr(256)")
    db_add_col(cursor, "variance", "dmis_HGVSp", "varchr(256)")
    sql_list = []
    print("collecting sql ...")
    icounter = 0
    for id in id_list:
        if id in id2info_dict:
            gene_name, HGVSc, HGVSp = id2info_dict[id]
            sql_list.append("UPDATE variance SET dmis = '1' WHERE id = {0}".format(id))
            sql_list.append("UPDATE variance SET dmis_gene_name = '{1}' WHERE id = {0}".format(id, gene_name))
            sql_list.append("UPDATE variance SET dmis_HGVSc = '{1}' WHERE id = {0}".format(id, HGVSc))
            sql_list.append("UPDATE variance SET dmis_HGVSp = '{1}' WHERE id = {0}".format(id, HGVSp))
        else:
            sql_list.append("UPDATE variance SET dmis = '0' WHERE id = {0}".format(id))
            sql_list.append("UPDATE variance SET dmis_gene_name = NULL WHERE id = {0}".format(id))
            sql_list.append("UPDATE variance SET dmis_HGVSc = NULL WHERE id = {0}".format(id))
            sql_list.append("UPDATE variance SET dmis_HGVSp = NULL WHERE id = {0}".format(id))
        icounter += 1
        if icounter % 1000 == 0:
            print("collecting sql ... {:.2%}".format(float(icounter) / len(id_list)))
    print("executing sql ...")
    icounter = 0
    for sql in sql_list:
        cursor.execute(sql)
        icounter += 1
        if icounter % 3000 == 0:
            print("executing sql ... {:.2%}".format(float(icounter) / len(sql_list)))
    cursor.close()
    conn.commit()
    conn.close()
    print("done")


def db_add_dsplicing_geneinfo(db_file, dsplicing_file):
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()
    id2info_dict = {}
    print("building id2info dict ...")
    icounter = 0
    with open(dsplicing_file, "r") as fp:
        while True:
            data_line = fp.readline()
            if not data_line:
                break
            if len(data_line) == 0 or data_line.startswith("chrom\tpos\t"):
                continue
            chrom, pos, ref, alt, region, gene_name, detailed_consequence, ada_score, rf_score = data_line.strip().split(
                "\t")
            cursor.execute("select id from variance where chr='{0}' and pos='{1}' and ref='{2}' and alt='{3}'"
                           "".format(chrom, pos, ref, alt))
            sql_result = cursor.fetchall()
            icounter += 1
            if icounter % 100 == 0:
                print("building id2info dict ... handled {0} lines".format(icounter))
            if len(sql_result) == 0:
                print("Variance {0} does't exist in data base. It might be deleted before.".format(data_line.strip()))
                continue
            elif len(sql_result) > 1:
                print("more than one variance match ({0})".format(data_line.strip()))
                continue
            id2info_dict[sql_result[0][0]] = [region, gene_name, detailed_consequence, ada_score, rf_score]
    cursor.execute("select id from variance")
    id_list = [i[0] for i in cursor.fetchall()]

    db_add_col(cursor, "variance", "dsplicing_gene_name", "varchr(256)")
    db_add_col(cursor, "variance", "dsplicing_region", "varchr(256)")
    db_add_col(cursor, "variance", "dsplicing_detailed_consequence", "varchr(256)")
    db_add_col(cursor, "variance", "dsplicing_ada_score", "varchr(256)")
    db_add_col(cursor, "variance", "dsplicing_rf_score", "varchr(256)")
    sql_list = []
    print("collecting sql ...")
    icounter = 0
    for id in id_list:
        sql_list.append("UPDATE variance SET spidex = NULL WHERE id = {0}".format(id))
        if id in id2info_dict:
            region, gene_name, detailed_consequence, ada_score, rf_score = id2info_dict[id]
            sql_list.append("UPDATE variance SET dsplicing = '1' WHERE id = {0}".format(id))
            sql_list.append("UPDATE variance SET dsplicing_gene_name = '{1}' WHERE id = {0}".format(id, gene_name))
            sql_list.append("UPDATE variance SET dsplicing_region = '{1}' WHERE id = {0}".format(id, region))
            sql_list.append("UPDATE variance SET dsplicing_detailed_consequence = '{1}' WHERE id = {0}"
                            "".format(id, detailed_consequence))
            sql_list.append("UPDATE variance SET dsplicing_ada_score = '{1}' WHERE id = {0}"
                            "".format(id, ada_score))
            sql_list.append("UPDATE variance SET dsplicing_rf_score = '{1}' WHERE id = {0}"
                            "".format(id, rf_score))
        else:
            sql_list.append("UPDATE variance SET dsplicing = '0' WHERE id = {0}".format(id))
            sql_list.append("UPDATE variance SET dsplicing_gene_name = NULL WHERE id = {0}".format(id))
            sql_list.append("UPDATE variance SET dsplicing_region = NULL WHERE id = {0}".format(id))
            sql_list.append("UPDATE variance SET dsplicing_detailed_consequence = NULL WHERE id = {0}".format(id))
            sql_list.append("UPDATE variance SET dsplicing_ada_score = NULL WHERE id = {0}".format(id))
            sql_list.append("UPDATE variance SET dsplicing_rf_score = NULL WHERE id = {0}".format(id))
        icounter += 1
        if icounter % 1000 == 0:
            print("collecting sql ... {:.2%}".format(float(icounter) / len(id_list)))
    print("executing sql ...")
    icounter = 0
    for sql in sql_list:
        cursor.execute(sql)
        icounter += 1
        if icounter % 3000 == 0:
            print("executing sql ... {:.2%}".format(float(icounter) / len(sql_list)))
    cursor.close()
    conn.commit()
    conn.close()
    print("done")


def select_gene_name_in_tsv(tsv_in, tsv_out):
    # 'indel-frameshift', 'indel-nonFrameshift', 'stopGain', 'startLoss', 'stopLoss'
    def check_select(site_type, exonic_function, transcript_ID):
        if transcript_ID.startswith("NR"):
            return False
        if site_type in ["spliceDonor", "spliceAcceptor"]:
            return True
        if site_type == "exonic" and exonic_function in ['indel-frameshift', 'indel-nonFrameshift', 'stopGain',
                                                         'startLoss', 'stopLoss']:
            return True
        return False

    def check_select2(site_type, transcript_ID):
        if transcript_ID.startswith("NR"):
            return False
        if site_type in ["spliceDonor", "spliceAcceptor"]:
            return True
        return False

    def build_key(gene_name, site_type, exonic_allele_function, transcript_ID, mRNA):
        ret = "{0}\t{1}\t{2}\t{3}\t{4}".format(gene_name, site_type, exonic_allele_function, transcript_ID, mRNA)
        return ret

    def build_key2(gene_name, site_type, exonic_allele_function):
        ret = "{0}\t{1}\t{2}".format(gene_name, site_type, exonic_allele_function)
        return ret

    def get_final_function_list(site_type_list, exonic_allele_function_list):
        assert len(site_type_list) == len(exonic_allele_function_list)
        ret_list = []
        for i in xrange(len(site_type_list)):
            if site_type_list[i] == "exonic":
                ret_list.append(exonic_allele_function_list[i])
            else:
                ret_list.append(site_type_list[i])
        return ret_list

    with open(tsv_in, "r") as fp_in, open(tsv_out, "w") as fp_out:
        fp_out.write(fp_in.readline().strip() + "\tfinal_function\tgene_name\n")
        icounter = 2
        while True:
            data_line = fp_in.readline()

            if not data_line:
                break
            if not data_line.strip():
                icounter += 1
                continue
            data_list = data_line.strip().split("\t")

            site_type_list = re.split("[|]", data_list[17])
            exonic_allele_function_list = re.split("[|]", data_list[18])
            gene_name_list = re.split("[|]", data_list[34])
            mRNA_list = re.split("[|]", data_list[27])
            transcript_ID_list = re.split("[|]", data_list[33])

            if len(site_type_list) != len(exonic_allele_function_list) or len(site_type_list) != len(
                    gene_name_list) or len(site_type_list) != len(mRNA_list) or len(site_type_list) != len(
                transcript_ID_list):
                print(icounter)
                print(site_type_list)
                print(exonic_allele_function_list)
                print(gene_name_list)
                print(mRNA_list)
                print(transcript_ID_list)
            assert len(site_type_list) == len(exonic_allele_function_list) == len(gene_name_list) == len(
                mRNA_list) == len(transcript_ID_list)
            ret_set = set([])
            ret_set2 = set([])
            for index in xrange(len(site_type_list)):
                if exonic_allele_function_list[index] != "!":
                    sub_site_type_list = site_type_list[index].split(";")
                    sub_exonic_allele_function_list = exonic_allele_function_list[index].split(";")
                    sub_gene_name_list = gene_name_list[index].split(";")
                    sub_mRNA_list = mRNA_list[index].split(";")
                    sub_transcript_ID_list = transcript_ID_list[index].split(";")

                    if len(sub_site_type_list) != len(sub_exonic_allele_function_list) or len(
                            sub_site_type_list) != len(sub_gene_name_list) or len(sub_site_type_list) != len(
                        sub_mRNA_list) or len(sub_site_type_list) != len(sub_transcript_ID_list):
                        # print exonic_allele_function_list
                        print(sub_site_type_list)
                        print(sub_exonic_allele_function_list)
                        print(sub_gene_name_list)
                        print(sub_mRNA_list)
                        print(sub_transcript_ID_list)
                    assert len(sub_site_type_list) == len(sub_exonic_allele_function_list) == len(
                        sub_gene_name_list) == len(sub_mRNA_list) == len(sub_transcript_ID_list)
                    selector = [check_select(sub_site_type_list[i],
                                             sub_exonic_allele_function_list[i],
                                             sub_transcript_ID_list[i]) for i in
                                xrange(len(sub_gene_name_list))]
                    # print(sub_site_type_list)
                    # print(sub_exonic_allele_function_list)
                    # print(sub_transcript_ID_list)
                    # print(selector)
                    sub_site_type_list = list(compress(sub_site_type_list, selector))
                    sub_exonic_allele_function_list = list(compress(sub_exonic_allele_function_list, selector))
                    sub_gene_name_list = list(compress(sub_gene_name_list, selector))
                    sub_mRNA_list = list(compress(sub_mRNA_list, selector))
                    sub_transcript_ID_list = list(compress(sub_transcript_ID_list, selector))

                    for i in xrange(len(sub_gene_name_list)):
                        ret_set.add(build_key(sub_gene_name_list[i],
                                              sub_site_type_list[i],
                                              sub_exonic_allele_function_list[i],
                                              sub_transcript_ID_list[i],
                                              sub_mRNA_list[i]))
                        ret_set2.add(build_key2(sub_gene_name_list[i],
                                                sub_site_type_list[i],
                                                sub_exonic_allele_function_list[i]))
                else:
                    sub_site_type_list = site_type_list[index].split(";")
                    sub_gene_name_list = gene_name_list[index].split(";")
                    sub_mRNA_list = mRNA_list[index].split(";")
                    sub_transcript_ID_list = transcript_ID_list[index].split(";")
                    if len(sub_site_type_list) != len(sub_gene_name_list) or len(sub_site_type_list) != len(
                            sub_mRNA_list) or len(sub_site_type_list) != len(sub_transcript_ID_list):
                        print(sub_site_type_list)
                        print(sub_gene_name_list)
                        print(sub_mRNA_list)
                        print(sub_transcript_ID_list)
                    assert len(sub_site_type_list) == len(sub_gene_name_list) == len(sub_mRNA_list) == len(
                        sub_transcript_ID_list)
                    selector = [check_select2(sub_site_type_list[i],
                                              sub_transcript_ID_list[i]) for i in xrange(len(sub_gene_name_list))]
                    # print(sub_site_type_list)
                    # print(sub_transcript_ID_list)
                    # print(selector)
                    sub_site_type_list = list(compress(sub_site_type_list, selector))
                    sub_gene_name_list = list(compress(sub_gene_name_list, selector))
                    sub_mRNA_list = list(compress(sub_mRNA_list, selector))
                    sub_transcript_ID_list = list(compress(sub_transcript_ID_list, selector))
                    # if len(sub_gene_name_list) > 0:
                    #     ret_gene_name_list.append(";".join(sub_gene_name_list))
                    #     ret_site_type_list.append(";".join(sub_site_type_list))
                    #     ret_exonic_allele_function_list.append("!")
                    for i in xrange(len(sub_gene_name_list)):
                        ret_set.add(build_key(sub_gene_name_list[i],
                                              sub_site_type_list[i],
                                              "!",
                                              sub_transcript_ID_list[i],
                                              sub_mRNA_list[i]))
                        ret_set2.add(build_key2(sub_gene_name_list[i],
                                                sub_site_type_list[i],
                                                "!"))
            data_list[34] = ";".join([i.split("\t")[0] for i in ret_set2])  # name2: gene name
            data_list[17] = ";".join([i.split("\t")[1] for i in ret_set2])  # siteType
            data_list[18] = ";".join([i.split("\t")[2] for i in ret_set2])  # exonicAlleleFunction
            data_list[33] = ";".join([i.split("\t")[3] for i in ret_set])  # name: transcript ID
            data_list[27] = ";".join([i.split("\t")[4] for i in ret_set])  # mRNA: the transcript ID starting with NM_

            final_function_list = get_final_function_list([i.split("\t")[1] for i in ret_set2],
                                                          [i.split("\t")[2] for i in ret_set2])
            name2_list = [i.split("\t")[0] for i in ret_set2]
            if len(set(final_function_list)) == 1 or len(set(name2_list)) == 1:
                final_function_list = list(set(final_function_list))
                name2_list = list(set(name2_list))
            data_list.append(";".join(final_function_list))
            data_list.append(";".join(name2_list))
            fp_out.write("\t".join(data_list) + "\n")
            icounter += 1


def select_gene_name_in_synonymous_tsv(tsv_in, tsv_out):
    def check_select(site_type, exonic_function, transcript_ID):
        if transcript_ID.startswith("NR"):
            return False
        if site_type == "exonic" and exonic_function in ['synonymous']:
            return True
        return False

    def build_key(gene_name, site_type, exonic_allele_function, transcript_ID, mRNA):
        ret = "{0}\t{1}\t{2}\t{3}\t{4}".format(gene_name, site_type, exonic_allele_function, transcript_ID, mRNA)
        return ret

    def build_key2(gene_name, site_type, exonic_allele_function):
        ret = "{0}\t{1}\t{2}".format(gene_name, site_type, exonic_allele_function)
        return ret

    def get_final_function_list(site_type_list, exonic_allele_function_list):
        assert len(site_type_list) == len(exonic_allele_function_list)
        ret_list = []
        for i in xrange(len(site_type_list)):
            if site_type_list[i] == "exonic":
                ret_list.append(exonic_allele_function_list[i])
            else:
                ret_list.append(site_type_list[i])
        return ret_list

    with open(tsv_in, "r") as fp_in, open(tsv_out, "w") as fp_out:
        fp_out.write(fp_in.readline().strip() + "\tfinal_function\tgene_name\n")
        icounter = 2
        while True:
            data_line = fp_in.readline()

            if not data_line:
                break
            if not data_line.strip():
                icounter += 1
                continue
            data_list = data_line.strip().split("\t")

            site_type_list = re.split("[|]", data_list[17])
            exonic_allele_function_list = re.split("[|]", data_list[18])
            gene_name_list = re.split("[|]", data_list[34])
            mRNA_list = re.split("[|]", data_list[27])
            transcript_ID_list = re.split("[|]", data_list[33])

            if len(site_type_list) != len(exonic_allele_function_list) or len(site_type_list) != len(
                    gene_name_list) or len(site_type_list) != len(mRNA_list) or len(site_type_list) != len(
                transcript_ID_list):
                print(icounter)
                print(site_type_list)
                print(exonic_allele_function_list)
                print(gene_name_list)
                print(mRNA_list)
                print(transcript_ID_list)
            assert len(site_type_list) == len(exonic_allele_function_list) == len(gene_name_list) == len(
                mRNA_list) == len(transcript_ID_list)
            ret_set = set([])
            ret_set2 = set([])
            for index in xrange(len(site_type_list)):
                sub_site_type_list = site_type_list[index].split(";")
                sub_exonic_allele_function_list = exonic_allele_function_list[index].split(";")
                sub_gene_name_list = gene_name_list[index].split(";")
                sub_mRNA_list = mRNA_list[index].split(";")
                sub_transcript_ID_list = transcript_ID_list[index].split(";")

                if len(sub_site_type_list) != len(sub_exonic_allele_function_list) or len(
                        sub_site_type_list) != len(sub_gene_name_list) or len(sub_site_type_list) != len(
                    sub_mRNA_list) or len(sub_site_type_list) != len(sub_transcript_ID_list):
                    print(sub_site_type_list)
                    print(sub_exonic_allele_function_list)
                    print(sub_gene_name_list)
                    print(sub_mRNA_list)
                    print(sub_transcript_ID_list)
                assert len(sub_site_type_list) == len(sub_exonic_allele_function_list) == len(
                    sub_gene_name_list) == len(sub_mRNA_list) == len(sub_transcript_ID_list)
                selector = [check_select(sub_site_type_list[i],
                                         sub_exonic_allele_function_list[i],
                                         sub_transcript_ID_list[i]) for i in
                            xrange(len(sub_gene_name_list))]
                sub_site_type_list = list(compress(sub_site_type_list, selector))
                sub_exonic_allele_function_list = list(compress(sub_exonic_allele_function_list, selector))
                sub_gene_name_list = list(compress(sub_gene_name_list, selector))
                sub_mRNA_list = list(compress(sub_mRNA_list, selector))
                sub_transcript_ID_list = list(compress(sub_transcript_ID_list, selector))

                for i in xrange(len(sub_gene_name_list)):
                    ret_set.add(build_key(sub_gene_name_list[i],
                                          sub_site_type_list[i],
                                          sub_exonic_allele_function_list[i],
                                          sub_transcript_ID_list[i],
                                          sub_mRNA_list[i]))
                    ret_set2.add(build_key2(sub_gene_name_list[i],
                                            sub_site_type_list[i],
                                            sub_exonic_allele_function_list[i]))

            data_list[34] = ";".join([i.split("\t")[0] for i in ret_set2])  # name2: gene name
            data_list[17] = ";".join([i.split("\t")[1] for i in ret_set2])  # siteType
            data_list[18] = ";".join([i.split("\t")[2] for i in ret_set2])  # exonicAlleleFunction
            data_list[33] = ";".join([i.split("\t")[3] for i in ret_set])  # name: transcript ID
            data_list[27] = ";".join([i.split("\t")[4] for i in ret_set])  # mRNA: the transcript ID starting with NM_

            final_function_list = get_final_function_list([i.split("\t")[1] for i in ret_set2],
                                                          [i.split("\t")[2] for i in ret_set2])
            name2_list = [i.split("\t")[0] for i in ret_set2]
            if len(set(final_function_list)) == 1 or len(set(name2_list)) == 1:
                final_function_list = list(set(final_function_list))
                name2_list = list(set(name2_list))
            data_list.append(";".join(final_function_list))
            data_list.append(";".join(name2_list))
            fp_out.write("\t".join(data_list) + "\n")
            icounter += 1


def select_gene_name_in_annovar_result(file_in, file_out):
    def get_vcf_pos(start, org_ref, org_alt):
        if org_alt == "-":
            return int(start) - 1
        elif org_ref == "-":
            return int(start)
        elif len(org_ref) == len(org_alt):
            return int(start)
        else:
            RuntimeError("意外的ref和alt。ref={0} alt={1}".format(org_ref, org_alt))

    def get_function(func, exonic_func):
        if exonic_func == ".":
            return func
        return re.sub("exonic", "{0}".format(exonic_func), func)

    def get_detail(gene_detail, aachange):
        if gene_detail == ".":
            return aachange
        if aachange == ".":
            return gene_detail
        return "GeneDetail({0}) AAChange({1})".format(gene_detail, aachange)

    with open(file_in, "r") as fp_in, open(file_out, "w") as fp_out:
        fp_out.write(fp_in.readline().strip() + "\tpos\tfunction\tdetail\n")
        while True:
            data_line = fp_in.readline()
            if not data_line:
                break
            if not data_line.strip():
                continue
            data_list = data_line.strip().split("\t")
            data_list.append(str(get_vcf_pos(data_list[1], data_list[3], data_list[4])))
            data_list.append(get_function(data_list[5], data_list[8]))
            data_list.append(get_detail(data_list[7], data_list[9]))
            fp_out.write("\t".join(data_list) + "\n")


def reduce_dsplicing_result(file_in, output):
    with open(file_in, "r") as fp_in, open(output, "w") as fp_out:
        fp_out.write(
            "chrom\tpos\tref\talt\tRefSeq_region_dsplicing\tRefSeq_gene_dsplicing\tRefSeq_detailed.consequence_dsplicing\tada_score\trf_score\n")
        while True:
            data_line = fp_in.readline()
            if not data_line:
                break
            if data_line.startswith("chr\t"):
                continue
            data_list = data_line.strip().split("\t")
            ref = data_list[2]
            alt = data_list[3]
            chrom = data_list[4]
            pos = data_list[5]
            region = data_list[8]
            gene = data_list[9]
            detail = data_list[11]
            ada_score = data_list[16]
            rf_score = data_list[17]

            if float(ada_score) < 0.6:
                continue

            tmp_detail = ";".join(re.findall("\((\S*?)\)", gene))
            gene = re.sub("\((\S*?)\)", "", gene)
            if detail == "." and len(tmp_detail) > 0:
                detail = tmp_detail
            elif detail != "." and len(tmp_detail) > 0:
                detail = "{0}{1}".format(detail, tmp_detail)
            else:
                pass

            fp_out.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n".format(chrom, pos, ref, alt, region, gene,
                                                                                detail, ada_score, rf_score))


def generator_seq_from_ref(fasta_file_name):
    chr_name = "1"
    start = 0
    length = 0
    with open(fasta_file_name, "r") as fp:
        fasta_data = fp.read()
    fasta_list = fasta_data.split(">")
    del fasta_data
    fasta_list = [i.split("\n") for i in fasta_list]
    fasta_dict = dict(
        zip([i[0].split(" ")[0] for i in fasta_list], ["".join(i[1:]) for i in fasta_list]))  # chromesome: sequence
    del fasta_list
    while True:
        assert type(chr_name) == str
        assert type(start) == int
        assert type(length) == int
        if chr_name in fasta_dict:
            [chr_name, start, length] = (yield fasta_dict[chr_name][start - 1:start - 1 + length])
        else:
            [chr_name, start, length] = yield None


def reduce_vep_result(file_in, output, ref_in):
    def get_chr_pos_alt(location, allele, ref_source):
        if allele == "-":
            # print(location)
            chrom = re.findall("^chr(\S+):", location)[0]
            # pos = int(re.findall(":(\S+)-", location)[0]) - 1
            pos = int(re.split("[:-]", location)[1]) - 1
            return [chrom,
                    pos,
                    ref_source.send(["chr" + chrom, pos, 1])]
        elif re.match("^\S+:\S+-\S+$", location):
            chrom = re.findall("^chr(\S+):", location)[0]
            pos = int(re.findall(":(\S+)-", location)[0])
            return [chrom, pos, ref_source.send(["chr" + chrom, pos, 1]) + allele]
        else:
            # print(location)
            return [re.findall("^chr(\S+):", location)[0],
                    int(re.findall(":(\S+)$", location)[0]),
                    allele]

    print("loading ref...")
    ref_source = generator_seq_from_ref(ref_in)
    ref_source.next()
    # print(ref_source.send(["chr1", 86023006, 1]))
    gene_id_dict = {}
    feature_dict = {}
    consequence_dict = {}
    with open(file_in, "r") as fp_in, open(output, "w") as fp_out:
        print("loading data...")
        while True:
            data_line = fp_in.readline()
            if not data_line:
                break
            if data_line.startswith("#"):
                continue
            data_list = data_line.strip().split("\t")
            gene = data_list[3]
            feature = data_list[4]
            consequence = data_list[6]
            location = data_list[1]
            allele = data_list[2]
            chrom, pos, alt = get_chr_pos_alt(location, allele, ref_source)

            key = "{0}\t{1}\t{2}".format(chrom, pos, alt)
            if key not in gene_id_dict:
                gene_id_dict[key] = set([gene])
            else:
                gene_id_dict[key].add(gene)

            if key not in feature_dict:
                feature_dict[key] = set([feature])
            else:
                feature_dict[key].add(feature)
            if key not in consequence_dict:
                consequence_dict[key] = set([consequence])
            else:
                consequence_dict[key].add(consequence)

        print("writing result...")
        fp_out.write("#chr\tpos\talt\tGene\tFeature\tConsequence\n")
        for key in gene_id_dict:
            gene = ";".join(gene_id_dict[key])
            feature = ";".join(feature_dict[key])
            consequence = ";".join(consequence_dict[key])
            chrom, pos, alt = key.split("\t")
            fp_out.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(chrom, pos, alt, gene, feature, consequence))


def reduce_splice_ai(file_in, output):
    with open(file_in, "r") as fp_in, open(output, "w") as fp_out:
        while True:
            data_line = fp_in.readline()
            if not data_line:
                break
            data_list = data_line.strip().split("\t")
            # print(data_list)
            chrom = data_list[0][3:] if data_list[0].startswith("chr") else data_list[0]
            pos = data_list[1]
            ref = data_list[3]
            alt = data_list[4]
            # print(data_list[8])
            symbol = re.findall("SYMBOL=(\S+?)\|", data_list[7])[0]
            fp_out.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(chrom, pos, ref, alt, symbol))


def db_add_synonymous_annovar_gene_name(db_file, file_in):
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()
    id2info_dict = {}
    print("building id2info_dict ...")
    with open(file_in, "r") as fp:
        icounter = 0
        while True:
            data_line = fp.readline()
            if not data_line:
                break
            data_list = data_line.split("\t")
            tmp_list = filter(lambda x: len(x) > 0, data_list[2].split(","))
            gene_name_list = list(set([re.findall("(\S+?):\S+", i)[0] for i in tmp_list]))
            gene_name_str = ";".join(gene_name_list)
            detail = data_list[2]
            chrom = data_list[3][3:]
            pos = data_list[4]
            ref = data_list[6]
            alt = data_list[7]
            sql_cmd = "select id from synonymous_snp where chr='{0}' and pos='{1}' and ref='{2}' and alt='{3}'" \
                      "".format(chrom, pos, ref, alt)
            cursor.execute(sql_cmd)
            sql_result = cursor.fetchall()
            if len(sql_result) == 0:
                print("variance {0} {1} {2} {3} doesn't exist in data base.".format(chrom, pos, ref, alt))
                icounter += 1
                continue
            elif len(sql_result) > 1:
                print("variance {0} {1} {2} {3} has more than 1 record in data base.".format(chrom, pos, ref, alt))
                icounter += 1
                continue
            id2info_dict[sql_result[0][0]] = [gene_name_str, detail]
            icounter += 1
            if icounter % 1000 == 0:
                print("building id2info_dict ... handled {} lines".format(icounter))
    cursor.execute("select id from synonymous_snp")
    id_list = [i[0] for i in cursor.fetchall()]
    db_add_col(cursor, "synonymous_snp", "annovar_gene_name", "varchr(256)")
    db_add_col(cursor, "synonymous_snp", "annovar_detail", "varchr(256)")

    sql_list = []
    icounter = 0
    for id in id_list:
        if id in id2info_dict:
            gene_name_str, detail = id2info_dict[id]
            sql_list.append("UPDATE synonymous_snp SET annovar = '1' WHERE id = {0}".format(id))
            sql_list.append(
                "UPDATE synonymous_snp SET annovar_gene_name = '{1}' WHERE id = {0}".format(id, gene_name_str))
            sql_list.append("UPDATE synonymous_snp SET annovar_detail = '{1}' WHERE id = {0}".format(id, detail))
        else:
            sql_list.append("UPDATE synonymous_snp SET annovar = '0' WHERE id = {0}".format(id))
            sql_list.append("UPDATE synonymous_snp SET annovar_gene_name = NULL WHERE id = {0}".format(id))
            sql_list.append("UPDATE synonymous_snp SET annovar_detail = NULL WHERE id = {0}".format(id))
        icounter += 1
        if icounter % 1000 == 0:
            print("collecting sql ... {:.2%}".format(float(icounter) / len(id_list)))
    print("executing sql ...")
    icounter = 0
    for sql in sql_list:
        cursor.execute(sql)
        icounter += 1
        if icounter % 3000 == 0:
            print("executing sql ... {:.2%}".format(float(icounter) / len(sql_list)))
    cursor.close()
    conn.commit()
    conn.close()
    print("done")


def db_add_synonymous_vep_gene_name(db_file, selected_vep_file):
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()
    id2info_dict = {}
    with open(selected_vep_file, "r") as fp:
        while True:
            data_line = fp.readline()
            if not data_line:
                break
            if data_line.startswith("#chr\t"):
                continue
            if len(data_line.strip()) == 0:
                continue
            chrom, pos, alt, gene_id, feature, consequence = data_line.strip().split("\t")
            alt = alt.upper()
            sql_cmd = "select id from synonymous_snp where chr='{0}' and pos='{1}' and alt='{2}'".format(chrom, pos,
                                                                                                         alt)
            cursor.execute(sql_cmd)
            sql_result = cursor.fetchall()
            if len(sql_result) == 0:
                print("Variance (chr={0} pos={1} alt={2}) doesn't exist in variance table. It might be deleted before."
                      "".format(chrom, pos, alt))
                continue
            if len(sql_result) > 1:
                print("Variance (chr={0} pos={1} alt={2}) has more than 1 record in data base.".format(chrom, pos, alt))
                continue
            id2info_dict[sql_result[0][0]] = [gene_id, feature, consequence]
    cursor.execute("select id from synonymous_snp")
    id_list = [i[0] for i in cursor.fetchall()]

    db_add_col(cursor, "synonymous_snp", "vep_gene_id", "varchr(256)")
    db_add_col(cursor, "synonymous_snp", "vep_feature", "varchr(256)")
    db_add_col(cursor, "synonymous_snp", "vep_consequence", "varchr(256)")

    print("collecting sqls")
    sql_list = []
    for id in id_list:
        if id in id2info_dict:
            gene_id, feature, consequence = id2info_dict[id]
            sql_list.append("UPDATE synonymous_snp SET vep = '1' WHERE id = {0}".format(id))
            sql_list.append("UPDATE synonymous_snp SET vep_gene_id = '{1}' WHERE id = {0}".format(id, gene_id))
            sql_list.append("UPDATE synonymous_snp SET vep_feature = '{1}' WHERE id = {0}".format(id, feature))
            sql_list.append("UPDATE synonymous_snp SET vep_consequence = '{1}' WHERE id = {0}".format(id, consequence))
        else:
            sql_list.append("UPDATE synonymous_snp SET vep = '0' WHERE id = {0}".format(id))
            sql_list.append("UPDATE synonymous_snp SET vep_gene_id = NULL WHERE id = {0}".format(id))
            sql_list.append("UPDATE synonymous_snp SET vep_feature = NULL WHERE id = {0}".format(id))
            sql_list.append("UPDATE synonymous_snp SET vep_consequence = NULL WHERE id = {0}".format(id))
    print("executing sql ...")
    icounter = 0
    for sql in sql_list:
        cursor.execute(sql)
        icounter += 1
        if icounter % 1000 == 0:
            print("executing sql ...{:.2%}".format(float(icounter) / len(sql_list)))
    cursor.close()
    conn.commit()
    conn.close()
    print("done")


def db_add_synonymous_bystro_gene_name(db_file, bystro_tsv, ref_in):
    tsv_key2value_dict = {}

    print("loading ref...")
    ref_source = generator_seq_from_ref(ref_in)
    ref_source.next()
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()
    print("loading variance in database...")
    cursor.execute("select id, chr, pos, ref, alt from synonymous_snp")
    db_data = cursor.fetchall()
    db_id2key_dict = dict(
        zip([int(i[0]) for i in db_data], ["{0}\t{1}\t{2}\t{3}".format(i[1], i[2], i[3], i[4]) for i in db_data]))
    print("there are {} variance in database".format(len(db_id2key_dict)))

    with open(bystro_tsv, "r") as fp:
        while True:
            data_line = fp.readline()
            if not data_line:
                break
            if data_line.startswith("chrom\tpos"):
                continue
            data_list = data_line.strip().split("\t")
            # print(len(data_list))
            if len(data_list) != 114:
                continue
            chrom = data_list[0]
            vcf_pos = data_list[15]
            bystro_ref = data_list[16]
            bystro_alt = data_list[4]
            type = data_list[2]
            ref, alt = bystro_build_ref_alt(type, bystro_alt, bystro_ref, chrom, vcf_pos, ref_source)

            final_function = data_list[112]
            gene_name = data_list[113]
            key = "{0}\t{1}\t{2}\t{3}".format(chrom[3:], vcf_pos, ref, alt)
            tsv_key2value_dict[key] = [final_function, gene_name]
    print("there are {} records in tsv file".format(len(tsv_key2value_dict)))
    db_add_col(cursor, "synonymous_snp", "bystro_final_function", "varchr(256)")
    db_add_col(cursor, "synonymous_snp", "bystro_gene_name", "varchr(256)")
    sql_list = []
    icounter = 0
    print("collecting sqls... 0%")
    for id in db_id2key_dict:
        db_key = db_id2key_dict[id]
        if db_key not in tsv_key2value_dict:
            sql_list.append("UPDATE synonymous_snp SET bystro = '0' WHERE id = {0}".format(id))
            sql_list.append("UPDATE synonymous_snp SET bystro_final_function = NULL WHERE id = {0}".format(id))
            sql_list.append("UPDATE synonymous_snp SET bystro_gene_name = NULL WHERE id = {0}".format(id))
        else:
            sql_list.append("UPDATE synonymous_snp SET bystro = '1' WHERE id = {0}".format(id))
            final_function, gene_name = tsv_key2value_dict[db_key]
            sql_list.append(
                "UPDATE synonymous_snp SET bystro_final_function = '{1}' WHERE id = {0}".format(id, final_function))
            sql_list.append("UPDATE synonymous_snp SET bystro_gene_name = '{1}' WHERE id = {0}".format(id, gene_name))
        icounter += 1
        if icounter % 10000 == 0:
            print("collecting sqls... {:.2%}".format(float(icounter) / len(db_id2key_dict)))
    print("executing sqls... 0%")
    icounter = 0
    for sql in sql_list:
        cursor.execute(sql)
        icounter += 1
        if icounter % 10000 == 0:
            print("executing sqls... {:.2%}".format(float(icounter) / len(sql_list)))
    cursor.close()
    conn.commit()
    conn.close()
    print("done")


def check_gene_name(db_file):
    print("checking variance table")
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()
    cursor.execute(
        "select chr, pos,  ref, alt, bystro_gene_name, annovar_gene_name, vep_gene_id, spliceAI_gene_name, dmis_gene_name, dsplicing_gene_name from variance")
    variance_data = cursor.fetchall()
    for record in variance_data:
        chrom, pos, ref, alt, bystro_gene_name, annovar_gene_name, vep_gene_id, spliceAI_gene_name, dmis_gene_name, dsplicing_gene_name = record
        # print("chrom={0} pos={1} ref={2} alt={3} bystro_gene_name={4} annovar_gene_name={5}".format(chrom, pos, ref, alt, bystro_gene_name, annovar_gene_name))
        if bystro_gene_name is not None:
            bystro_gene_name_list = bystro_gene_name.split(";")
            for name in bystro_gene_name_list:
                cursor.execute("select id from gene_table where gene_name='{0}'".format(name))
                sql_result = cursor.fetchall()
                if len(sql_result) == 0:
                    print("variance: {0}\t{1}\t{2}\t{3}\tbystro gene name ({4}) doesn't have record in gene table"
                          "".format(chrom, pos, ref, alt, name))
                elif len(sql_result) > 1:
                    print("variance: {0}\t{1}\t{2}\t{3}\tbystro gene name ({4}) has more than one records in gene table"
                          "".format(chrom, pos, ref, alt, name))
        if annovar_gene_name is not None:
            # print("annovar_gene_name = {}".format(annovar_gene_name))
            annovar_gene_name_list = annovar_gene_name.split(";")
            for name in annovar_gene_name_list:
                cursor.execute("select id from gene_table where gene_name='{0}'".format(name))
                sql_result = cursor.fetchall()
                if len(sql_result) == 0:
                    print("variance: {0}\t{1}\t{2}\t{3}\tannovar gene name ({4}) doesn't have record in gene table"
                          "".format(chrom, pos, ref, alt, name))
                elif len(sql_result) > 1:
                    print(
                        "variance: {0}\t{1}\t{2}\t{3}\tannovar gene name ({4}) has more than one records in gene table"
                        "".format(chrom, pos, ref, alt, name))
        if vep_gene_id is not None and vep_gene_id != "None":
            vep_gene_id_list = vep_gene_id.split(";")
            for id in vep_gene_id_list:
                cursor.execute("select id from gene_table where gene_id='{0}'".format(id))
                sql_result = cursor.fetchall()
                if len(sql_result) == 0:
                    print("variance: {0}\t{1}\t{2}\t{3}\tvep gene id ({4}) doesn't have record in gene table"
                          "".format(chrom, pos, ref, alt, id))
                elif len(sql_result) > 1:
                    print("variance: {0}\t{1}\t{2}\t{3}\tvep gene id ({4}) has more than one records in gene table"
                          "".format(chrom, pos, ref, alt, id))
        if spliceAI_gene_name is not None:
            spliceAI_gene_name_list = spliceAI_gene_name.split(";")
            for name in spliceAI_gene_name_list:
                cursor.execute("select id from gene_table where gene_name='{0}'".format(name))
                sql_result = cursor.fetchall()
                if len(sql_result) == 0:
                    print("variance: {0}\t{1}\t{2}\t{3}\tspliceAI gene name ({4}) doesn't have record in gene table"
                          "".format(chrom, pos, ref, alt, name))
                elif len(sql_result) > 1:
                    print(
                        "variance: {0}\t{1}\t{2}\t{3}\tspliceAI gene name ({4}) has more than one records in gene table"
                        "".format(chrom, pos, ref, alt, name))
        if dmis_gene_name is not None:
            dmis_gene_name_list = dmis_gene_name.split(";")
            for name in dmis_gene_name_list:
                cursor.execute("select id from gene_table where gene_name='{0}'".format(name))
                sql_result = cursor.fetchall()
                if len(sql_result) == 0:
                    print("variance: {0}\t{1}\t{2}\t{3}\tdmis gene name ({4}) doesn't have record in gene table"
                          "".format(chrom, pos, ref, alt, name))
                elif len(sql_result) > 1:
                    print("variance: {0}\t{1}\t{2}\t{3}\tdmis gene name ({4}) has more than one records in gene table"
                          "".format(chrom, pos, ref, alt, name))
        if dsplicing_gene_name is not None:
            dsplicing_gene_name_list = dsplicing_gene_name.split(";")
            for name in dsplicing_gene_name_list:
                cursor.execute("select id from gene_table where gene_name='{0}'".format(name))
                sql_result = cursor.fetchall()
                if len(sql_result) == 0:
                    print("variance: {0}\t{1}\t{2}\t{3}\tdsplicing gene name ({4}) doesn't have record in gene table"
                          "".format(chrom, pos, ref, alt, name))
                elif len(sql_result) > 1:
                    print(
                        "variance: {0}\t{1}\t{2}\t{3}\tdsplicing gene name ({4}) has more than one records in gene table"
                        "".format(chrom, pos, ref, alt, name))
        # exit(0)

    print("checking synonymous_snp table")
    cursor.execute("select chr, pos,  ref, alt, annovar_gene_name, vep_gene_id, bystro_gene_name from synonymous_snp")
    synonymous_data = cursor.fetchall()
    for record in synonymous_data:
        chrom, pos, ref, alt, annovar_gene_name, vep_gene_id, bystro_gene_name = record
        if bystro_gene_name is not None:
            bystro_gene_name_list = bystro_gene_name.split(";")
            for name in bystro_gene_name_list:
                cursor.execute("select id from gene_table where gene_name='{0}'".format(name))
                sql_result = cursor.fetchall()
                if len(sql_result) == 0:
                    print("synonymous: {0}\t{1}\t{2}\t{3}\tbystro gene name ({4}) doesn't have record in gene table"
                          "".format(chrom, pos, ref, alt, name))
                elif len(sql_result) > 1:
                    print(
                        "synonymous: {0}\t{1}\t{2}\t{3}\tbystro gene name ({4}) has more than one records in gene table"
                        "".format(chrom, pos, ref, alt, name))
        if annovar_gene_name is not None:
            annovar_gene_name_list = re.split("[;,]", annovar_gene_name)
            for name in annovar_gene_name_list:
                cursor.execute("select id from gene_table where gene_name='{0}'".format(name))
                sql_result = cursor.fetchall()
                if len(sql_result) == 0:
                    print("synonymous: {0}\t{1}\t{2}\t{3}\tannovar gene name ({4}) doesn't have record in gene table"
                          "".format(chrom, pos, ref, alt, name))
                elif len(sql_result) > 1:
                    print(
                        "synonymous: {0}\t{1}\t{2}\t{3}\tannovar gene name ({4}) has more than one records in gene table"
                        "".format(chrom, pos, ref, alt, name))
        if vep_gene_id is not None:
            vep_gene_id_list = vep_gene_id.split(";")
            for id in vep_gene_id_list:
                cursor.execute("select id from gene_table where gene_id='{0}'".format(id))
                sql_result = cursor.fetchall()
                if len(sql_result) == 0:
                    print("synonymous: {0}\t{1}\t{2}\t{3}\tvep gene id ({4}) doesn't have record in gene table"
                          "".format(chrom, pos, ref, alt, id))
                elif len(sql_result) > 1:
                    print("synonymous: {0}\t{1}\t{2}\t{3}\tvep gene id ({4}) has more than one records in gene table"
                          "".format(chrom, pos, ref, alt, id))


def modify_gene_name(db_file, search_list, variance_not_in_gene_table):
    """
    Process each row in variance_not_in_gene_table.
            1 If the gene name has multiple records in the gene table, no processing will be done
            2 If the gene name is in the specified reserved set, do not process
            3 If the gene name is not in the search_list. If there is no comment, you need to 
              set the corresponding annotator to 0, and delete the comment information
            4 If the gene name is in the search_list, check it again with the id of the search_list
              If found, update the corresponding gene name in the variance table and synonymous
              table to keep the gene name in the gene_table
              If not found, print out
    @param db_file:
    @param search_list: enter
    @param variance_not_in_gene_table: enter
    @return:
    """

    def sql_del_gene_info(table_name, annotator, chrom, pos, ref, alt):
        if table_name == "variance":
            if annotator == "annovar":
                return "UPDATE variance SET annovar = '0', annovar_gene_name=NULL, " \
                       "annovar_function=NULL, annovar_function_detail=NULL " \
                       "WHERE chr='{0}' and pos='{1}' and ref='{2}' and alt='{3}'" \
                       "".format(chrom, pos, ref, alt)
            elif annotator == "bystro":
                return "UPDATE variance SET bystro = '0', bystro_final_function=NULL, " \
                       "bystro_gene_name=NULL WHERE chr='{0}' and pos='{1}' and ref='{2}' and alt='{3}'" \
                       "".format(chrom, pos, ref, alt)
            elif annotator == "vep":
                return "UPDATE variance SET vep = '0', vep_gene_id=NULL, " \
                       "vep_feature=NULL, vep_consequence=NULL " \
                       "WHERE chr='{0}' and pos='{1}' and ref='{2}' and alt='{3}'" \
                       "".format(chrom, pos, ref, alt)
            elif annotator == "spliceAI":
                return "UPDATE variance SET spliceAI = '0', spliceAI_gene_name=NULL " \
                       "WHERE chr='{0}' and pos='{1}' and ref='{2}' and alt='{3}'" \
                       "".format(chrom, pos, ref, alt)
            elif annotator == "dmis":
                return "UPDATE variance SET dmis = '0', dmis_gene_name=NULL, dmis_HGVSc=NULL, dmis_HGVSp=NULL  " \
                       "WHERE chr='{0}' and pos='{1}' and ref='{2}' and alt='{3}'" \
                       "".format(chrom, pos, ref, alt)
            elif annotator == "dsplicing":
                return "UPDATE variance SET dsplicing = '0', dsplicing_gene_name=NULL, " \
                       "dsplicing_region=NULL, dsplicing_detailed_consequence=NULL, " \
                       "dsplicing_ada_score=NULL, dsplicing_rf_score=NULL " \
                       "WHERE chr='{0}' and pos='{1}' and ref='{2}' and alt='{3}'" \
                       "".format(chrom, pos, ref, alt)
            else:
                RuntimeError("unexpect annotator=[{}]".format(annotator))
        elif table_name == "synonymous_snp":
            if annotator == "annovar":
                return "UPDATE synonymous_snp SET annovar = '0',annovar_gene_name=NULL, " \
                       "annovar_detail=NULL " \
                       "WHERE chr='{0}' and pos='{1}' and ref='{2}' and alt='{3}'" \
                       "".format(chrom, pos, ref, alt)
            elif annotator == "vep":
                return "UPDATE synonymous_snp SET vep = '0',vep_gene_id=NULL, " \
                       "vep_feature=NULL, vep_consequence=NULL " \
                       "WHERE chr='{0}' and pos='{1}' and ref='{2}' and alt='{3}'" \
                       "".format(chrom, pos, ref, alt)
            elif annotator == "bystro":
                return "UPDATE synonymous_snp SET bystro = '0', bystro_final_function=NULL, " \
                       "bystro_gene_name=NULL " \
                       "WHERE chr='{0}' and pos='{1}' and ref='{2}' and alt='{3}'" \
                       "".format(chrom, pos, ref, alt)
            else:
                RuntimeError("unexpect annotator=[{}]".format(annotator))

    def update_gene_name(table_name, annotator, chrom, pos, ref, alt, new_gene_name):
        if table_name == "variance":
            if annotator == "annovar":
                return "UPDATE variance SET annovar_gene_name='{4}' " \
                       "WHERE chr='{0}' and pos='{1}' and ref='{2}' and alt='{3}'" \
                       "".format(chrom, pos, ref, alt, new_gene_name)
            elif annotator == "bystro":
                return "UPDATE variance SET bystro_gene_name='{4}' " \
                       "WHERE chr='{0}' and pos='{1}' and ref='{2}' and alt='{3}'" \
                       "".format(chrom, pos, ref, alt, new_gene_name)
            elif annotator == "vep":
                RuntimeError("should not appear vep annotator here")
            elif annotator == "spliceAI":
                return "UPDATE variance SET spliceAI_gene_name='{4}' " \
                       "WHERE chr='{0}' and pos='{1}' and ref='{2}' and alt='{3}'" \
                       "".format(chrom, pos, ref, alt, new_gene_name)
            elif annotator == "dmis":
                return "UPDATE variance SET spliceAI_gene_name='{4}' " \
                       "WHERE chr='{0}' and pos='{1}' and ref='{2}' and alt='{3}'" \
                       "".format(chrom, pos, ref, alt, new_gene_name)
            elif annotator == "dsplicing":
                return "UPDATE variance SET dsplicing_gene_name='{4}' " \
                       "WHERE chr='{0}' and pos='{1}' and ref='{2}' and alt='{3}'" \
                       "".format(chrom, pos, ref, alt, new_gene_name)
        elif table_name == "synonymous_snp":
            if annotator == "annovar":
                return "UPDATE synonymous_snp SET annovar_gene_name='{4}' " \
                       "WHERE chr='{0}' and pos='{1}' and ref='{2}' and alt='{3}'" \
                       "".format(chrom, pos, ref, alt, new_gene_name)
            elif annotator == "vep":
                RuntimeError("vep annotator should not apear here")
            elif annotator == "bystro":
                return "UPDATE synonymous_snp SET bystro_gene_name='{4}' " \
                       "WHERE chr='{0}' and pos='{1}' and ref='{2}' and alt='{3}'" \
                       "".format(chrom, pos, ref, alt, new_gene_name)

    def insert2gene_table(gene_name, chrom, start_pos, end_pos, cursor):
        cursor.execute("select id from gene_table where gene_name='{}'".format(gene_name))
        sql_result = cursor.fetchall()
        if len(sql_result) == 0:
            cursor.execute("INSERT INTO gene_table (gene_name, chr, start_pos, end_pos) "
                           "VALUES ('{0}', '{1}', {2}, {3})".format(gene_name, chrom, start_pos, end_pos))

    map_table_name_dict = {"variance": "variance", "synonymous": "synonymous_snp"}
    # print("loading search list ...")
    search_dict = {}
    with open(search_list, "r") as fp:
        while True:
            data_line = fp.readline()
            if not data_line:
                break
            gene_name, gene_id = data_line.strip().split("\t")
            search_dict[gene_name] = gene_id

    keep_set = {"PYDC5", "CNPY3-GNMT", "MMP24-AS1-EDEM2", "PRH1-TAS2R14", "SETDB2-PHF11",
                "TBC1D7-LOC100130357", "TNFAIP8L2-SCNM1", "GPR75-ASB3", "PRR25", "OR9G9",
                "OR8U8", "C8orf44-SGK3", "CRIPAK"}
    del_set = {"AHSA2", "OR51J1", "ZNF788", "KIAA0125", "C22orf24", "C7orf69", "C20orf78", "C11orf44",
               "C10orf111", "C20orf197", "C15orf56", "C8orf31", "C9orf62", "C10orf91", "C9orf139",
               "C17orf77", "C15orf32", "C17orf82", "C1orf204", "C22orf34", "C9orf163", "C17orf102",
               "C5orf56", "MDS2", "C6orf99", "C10orf126", "C9orf170", "C1orf195", "C8orf44", "C7orf65",
               "C1orf229", "COL4A2-AS2", "C12orf77", "C17orf112", "C3orf79", "C6orf183", "C8orf87",
               "C8orf49", "THEG5", "AP002884.2", "AC135048.1", "TMEM105", "CTD-2574D22.6", "AC005779.2",
               "KIR3DX1", "KM-PA-2", "SERF2-C15ORF63", "PRSS46"}
    map_gene_name_dict = {"HMHA1": "ARHGAP45"}
    insert2gene_table_dict = {"PYDC5": ["1", 158999971, 159000312],
                              "TBC1D7-LOC100130357": ["6", 13266542, 13328583],
                              "CNPY3-GNMT": ["6", 42929480, 42963880],
                              "TNFAIP8L2-SCNM1": ["1", 151156649, 151170296],
                              "MMP24-AS1-EDEM2": ["20", 35115364, 35278122],
                              "GPR75-ASB3": ["2", 53670316, 53860033],
                              "PRH1-TAS2R14": ["12", 10937410, 11171611],
                              "PRR25": ["16", 805443, 813861],
                              "OR9G9": ["11", 56700388, 56701305],
                              "SETDB2-PHF11": ["13", 49444274, 49528976],
                              "OR8U8": ["11", 56375624, 56376499],
                              "C8orf44-SGK3": ["8", 66667596, 66862022],
                              "CRIPAK": ["4", 1391552, 1395989]}
    sql_list = []
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()
    for gene_name in insert2gene_table_dict:
        chrom, start_pos, end_pos = insert2gene_table_dict[gene_name]
        insert2gene_table(gene_name, chrom, start_pos, end_pos, cursor)
    with open(variance_not_in_gene_table, "r") as fp:
        while True:
            data_line = fp.readline()
            if not data_line:
                break
            if data_line.startswith("checking"):
                continue
            data_list = re.findall(
                "^(variance|synonymous): (\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+) gene name \((\S+)\) doesn't have record in gene table\n",
                data_line)
            if len(data_list) == 0:
                continue
            table_name, chrom, pos, ref, alt, annotator, gene_name = data_list[0]
            table_name = map_table_name_dict[table_name]
            if gene_name in keep_set:
                continue

            if gene_name not in search_dict:
                sql_list.append(sql_del_gene_info(table_name, annotator, chrom, pos, ref, alt))
            else:
                if gene_name in del_set:
                    sql_list.append(sql_del_gene_info(table_name, annotator, chrom, pos, ref, alt))
                    continue
                if gene_name in map_gene_name_dict:
                    sql_list.append(update_gene_name(table_name, annotator, chrom, pos, ref, alt,
                                                     map_gene_name_dict[gene_name]))
                    continue
                gene_id = search_dict[gene_name]
                cursor.execute("select gene_name from gene_table where gene_id='{0}'".format(gene_id))
                sql_ret = cursor.fetchall()
                if len(sql_ret) == 0:
                    print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\tgene id ({7}) doesn't have record in gene table"
                          "".format(table_name, chrom, pos, ref, alt, annotator, gene_name, gene_id))
                elif len(sql_ret) > 1:
                    print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\tgene id ({7}) has {8} records in gene table"
                          "".format(table_name, chrom, pos, ref, alt, annotator, gene_name, gene_id, len(sql_ret)))
                else:
                    sql_list.append(update_gene_name(table_name, annotator, chrom, pos, ref, alt, sql_ret[0][0]))
    for sql in sql_list:
        cursor.execute(sql)

    cursor.close()
    conn.commit()
    conn.close()


def build_contingency_table_gene_name_synonymous(db_file, phenotype,
                                                 sample_table_name, sample_restrict,
                                                 gene_table_name, gene_restrict,
                                                 variance_restrict, synonymous_restrict,
                                                 output, job_id, should_log=True):
    build_contingency_table_gene_name_inner("synonymous_snp", db_file, phenotype,
                                            sample_table_name, sample_restrict,
                                            gene_table_name, gene_restrict,
                                            variance_restrict, synonymous_restrict,
                                            output, job_id, should_log)


def build_contingency_table_gene_name(db_file, phenotype,
                                      sample_table_name, sample_restrict,
                                      gene_table_name, gene_restrict,
                                      variance_restrict, synonymous_restrict,
                                      output, job_id, should_log=True):
    """
           case     control
        ---------------------
        |         |         |
    alt |    A    |    B    |
        |         |         |
        ---------------------
        |         |         |
    ref |    C    |    D    |
        |         |         |
        ---------------------
    以基因为单位进行fisher exact test
    :return:
    """
    build_contingency_table_gene_name_inner("PE_GATK_variant", db_file, phenotype,
                                            sample_table_name, sample_restrict,
                                            gene_table_name, gene_restrict,
                                            variance_restrict, synonymous_restrict,
                                            output, job_id, should_log)


def build_contingency_table_gene_name_inner(variance_table, db_file, phenotype,
                                            sample_table_name, sample_restrict,
                                            gene_table_name, gene_restrict,
                                            variance_restrict, synonymous_restrict,
                                            output, job_id, should_log=True):
    """
           case     control
        ---------------------
        |         |         |
    alt |    A    |    B    |
        |         |         |
        ---------------------
        |         |         |
    ref |    C    |    D    |
        |         |         |
        ---------------------
    Carry out fisher exact test in units of genes
    :return:
    """

    def sub_gene_name(variance_table, variance_restrict):
        sub_str_list = []
        if variance_table == "variance" or variance_table == 'PE_GATK_variant':
            if re.findall("annovar[ =]", variance_restrict):
                sub_str_list.append("annovar_gene_name")
            if re.findall("bystro[ =]", variance_restrict):
                sub_str_list.append("bystro_gene_name")
            if re.findall("vep[ =]", variance_restrict):
                sub_str_list.append("vep_gene_id")
            if re.findall("spliceAI[ =]", variance_restrict):
                sub_str_list.append("spliceAI_gene_name")
            if re.findall("dmis[ =]", variance_restrict):
                sub_str_list.append("dmis_gene_name")
            if re.findall("dsplicing[ =]", variance_restrict):
                sub_str_list.append("dsplicing_gene_name")
            assert len(sub_str_list) > 0

        else:
            if re.findall("annovar[ =]", variance_restrict):
                sub_str_list.append("annovar_gene_name")
            if re.findall("bystro[ =]", variance_restrict):
                sub_str_list.append("bystro_gene_name")
            if re.findall("vep[ =]", variance_restrict):
                sub_str_list.append("vep_gene_id")
            assert len(sub_str_list) > 0
        ret_str = ", ".join(sub_str_list)
        # print(ret_str)
        return ret_str

    def build_gene_name_list(sql_data):
        ret_list = []
        for i in sql_data:
            ret_list.extend(re.split("[,;]", i))
        return ret_list

    synonymous_table = "synonymous_snp"
    if should_log:
        logging.basicConfig(filename="{0}{1}.log".format(sys._getframe().f_code.co_name, job_id),
                            level=logging.DEBUG, format=log_format, filemode="w")
        logging.debug("begin")
        logging.debug("db_file=[{0}]".format(db_file))
        logging.debug("phenotype=[{0}]".format(phenotype))
        logging.debug("sample_table_name=[{0}]".format(sample_table_name))
        logging.debug("sample_restrict=[{0}]".format(sample_restrict))
        logging.debug("variance_restrict=[{0}]".format(variance_restrict))
        logging.debug("synonymous_restrict=[{0}]".format(synonymous_restrict))
        logging.debug("output=[{0}]".format(output))

    if phenotype not in ["heart6", "ps_andor_pa6", "raa6", "iaab6", "pta6",
                         "tof6", "asdall6", "asdalone6", "vsd6", "vsdalone6",
                         "tofall6", "purevsdalone6", "ps_or_pa_and_vsd6", "intracardiac6", "aorticarch6",
                         "heartnoasd6", "tof_or_pta6", "tof_or_pta_or_iaab6", "CTD"]:
        if should_log:
            logging.error("illegal phenotype [{}]".format(phenotype))
        return

    if should_log:
        logging.debug("begin select control list and case list...")
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()

    gene_table_info = cursor.execute("PRAGMA table_info([gene_table])")
    gene_set_col_name_list = [i[1] for i in gene_table_info if i[1].startswith('gene_set_')]
    format_sample_restrict = "" if not sample_restrict else " AND {}".format(sample_restrict)
    cursor.execute("SELECT s.gen_id FROM {2} AS s WHERE s.{0}='0'{1}".format(phenotype,
                                                                             format_sample_restrict,
                                                                             sample_table_name))
    control_id_list = [i[0] for i in cursor.fetchall()]  # m
    print("control number={}".format(len(control_id_list)))
    cursor.execute("SELECT s.gen_id FROM {2} AS s WHERE s.{0}='1'{1}".format(phenotype,
                                                                             format_sample_restrict,
                                                                             sample_table_name))
    case_id_list = [i[0] for i in cursor.fetchall()]  # n
    print("case number={}".format(len(case_id_list)))
    if should_log:
        logging.debug("begin handle genes")
    with open(output, "w") as fp:
        fp.write("##phenotype:\"{0}\"\n##control number:{1}\n##case number:{2}\n"
                 "##sample restrict:\"{3}\"\n##version: 2 (select variance with gene name instead of gene region)\n"
                 "".format(phenotype, len(control_id_list), len(case_id_list), sample_restrict))
        fp.write("""##                         case     control
##                      ---------------------
##                      |         |         |
##    alt allele number |    A    |    B    |
##                      |         |         |
##                      ---------------------
##                      |         |         |
##    ref allele number |    C    |    D    |
##                      |         |         |
##                      ---------------------
##
##
##                              case     control
##                           ---------------------
##                           |         |         |
##          subjects have alt|    A1   |    B1   |
##                           |         |         |
##                           ---------------------
##                           |         |         |
##  subjects do not have alt |    C1   |    D1   |
##                           |         |         |
##                           ---------------------
""")
        fp.write(
            "#gene_id\tgene_name\tA\tB\tC\tD\tp_value\todds_ratio\tsynonymous_case_alt(C)\tsynonymous_control_alt(D)\todds_ratio_synonymous\tA1\tB1\tC1\tD1\tp_value1\todds_ratio1\t"
            "n_variance_in_gene\tvariance_in_gene_name\tn_variance_control\tvariance_in_control_name\t"
            "n_variance_case\tvariance_in_case_name\tmouse_mean_WT_OFT_TPM_2\tmouse_mean_WT_LV_TPM_3\tmouse_mean_WT_PA_TPM_3\t"
            "Interactions_IntAct\tInteractions_BioGRID\tInteractions_ConsensusPathDB\tHIPred\tGHIS\tgnomAD_pLI\t"
            "GDI_Phred\tGene_damage_prediction__all_disease_causing_genes\tEssential_gene\tEssential_gene_CRISPR\t"
            "Gene_indispensability_pred\tgevir_percentile\tloeuf_percentile\tvirlof_percentile\t"
            "HHE_E14.5\t{}\n".format("\t".join(gene_set_col_name_list)))
        # gene data
        format_gene_restrict = "" if not gene_restrict else " WHERE {}".format(gene_restrict)
        cmd_str = "SELECT g.gene_id, g.gene_name, mouse_mean_WT_OFT_TPM_2, mouse_mean_WT_LV_TPM_3, mouse_mean_WT_PA_TPM_3, " \
                  "Interactions_IntAct_, Interactions_BioGRID_, Interactions_ConsensusPathDB_, " \
                  "HIPred, GHIS, gnomAD_pLI, GDI_Phred, Gene_damage_prediction__all_disease_causing_genes_ , " \
                  "Essential_gene, Essential_gene_CRISPR , Gene_indispensability_pred , gevir_percentile, " \
                  "loeuf_percentile, virlof_percentile, " \
                  "\"HHE_E14.5\", {2} FROM {0} AS g{1}".format(gene_table_name, format_gene_restrict,
                                                               ", ".join(gene_set_col_name_list))
        logging.debug("gene sql = [{}]".format(cmd_str))
        cursor.execute(cmd_str)
        gene_data = cursor.fetchall()
        print("gene number = {0}".format(len(gene_data)))
        # prepare the variance data
        print("build id data ...")
        cmd_str = "select id, {0} from {1} where {2}".format(sub_gene_name(variance_table, variance_restrict),
                                                             variance_table,
                                                             variance_restrict)
        cursor.execute(cmd_str)
        variance_id_data = [[i[0], build_gene_name_list(filter(lambda x: x is not None, i[1:]))] for i in
                            cursor.fetchall()]

        cmd_str = "select id, {0} from {1} where {2}".format(sub_gene_name(synonymous_table, synonymous_restrict),
                                                             synonymous_table,
                                                             synonymous_restrict)
        cursor.execute(cmd_str)
        synonymous_id_data = [[i[0], build_gene_name_list(filter(lambda x: x is not None, i[1:]))] for i in
                              cursor.fetchall()]
        cmd_str = "SELECT vcf_id, chr, pos, "
        for control_id in control_id_list:
            cmd_str = "{0}sample_{1}, ".format(cmd_str, control_id)
        for case_id in case_id_list:
            cmd_str = "{0}sample_{1}, ".format(cmd_str, case_id)
        cmd_str = cmd_str.strip(", ")
        assert variance_restrict != ""
        # format_variance_restrict = " WHERE {0}".format(variance_restrict)
        cmd_synonymous_str = "{0} FROM {1} as v".format(cmd_str, synonymous_table)
        cmd_str = "{0} FROM {1} as v".format(cmd_str, variance_table)

        ret_str = ""
        icounter = 0
        gene_data_len = len(gene_data)
        for gene_data_list in gene_data:
            gene_id = gene_data_list[0]
            gene_name = gene_data_list[1]
            mouse_mean_WT_OFT_TPM_2 = gene_data_list[2]
            mouse_mean_WT_LV_TPM_3 = gene_data_list[3]
            mouse_mean_WT_PA_TPM_3 = gene_data_list[4]
            Interactions_IntAct_ = gene_data_list[5]
            Interactions_BioGRID_ = gene_data_list[6]
            Interactions_ConsensusPathDB_ = gene_data_list[7]
            HIPred = gene_data_list[8]
            GHIS = gene_data_list[9]
            gnomAD_pLI = gene_data_list[10]
            GDI_Phred = gene_data_list[11]
            Gene_damage_prediction__all_disease_causing_genes_ = gene_data_list[12]
            Essential_gene = gene_data_list[13]
            Essential_gene_CRISPR = gene_data_list[14]
            Gene_indispensability_pred = gene_data_list[15]
            gevir_percentile = gene_data_list[16]
            loeuf_percentile = gene_data_list[17]
            virlof_percentile = gene_data_list[18]
            hhe = gene_data_list[19]
            if should_log:
                if icounter % 100 == 0 and icounter > 0:
                    logging.debug("handled {0} / {1} genes".format(icounter, gene_data_len))
            variance_id_list = [str(i[0]) for i in variance_id_data if gene_id is not None and (
                    gene_name in i[1] or gene_id in i[1]) or gene_id is None and gene_name in i[1]]
            synonymous_id_list = [str(i[0]) for i in synonymous_id_data if gene_id is not None and (
                    gene_name in i[1] or gene_id in i[1]) or gene_id is None and gene_name in i[1]]
            if len(variance_id_list) == 0:
                icounter += 1
                continue
            tmp_cmd_str = cmd_str + " where v.id in (" + ", ".join(variance_id_list) + ")"
            # print(tmp_cmd_str)
            cursor.execute(tmp_cmd_str)

            variance_selected_list = [list(i) for i in cursor.fetchall()]

            tmp_cmd_str = cmd_synonymous_str + " where v.id in (" + ", ".join(synonymous_id_list) + ")"
            cursor.execute(tmp_cmd_str)
            synonymous_selected_list = [list(i) for i in cursor.fetchall()]

            all_variance_in_gene_num = len(variance_selected_list)
            all_variance_in_gene_name_list = [i[0] for i in variance_selected_list]  # vcf_id
            control_data = [i[3:3 + len(control_id_list)] for i in variance_selected_list]
            case_data = [i[3 + len(control_id_list):3 + len(control_id_list) + len(case_id_list)] for i in
                         variance_selected_list]

            synonymous_control_data = [i[3:3 + len(control_id_list)] for i in synonymous_selected_list]
            synonumous_case_data = [i[3 + len(control_id_list):3 + len(control_id_list) + len(case_id_list)] for i in
                                    synonymous_selected_list]
            tmp_list = []
            control_variance_in_gene_num = 0
            control_variance_in_gene_name_list = []
            for i in xrange(len(control_data)):
                control_line_data = filter(lambda x: type(x) == int, control_data[i])
                tmp_list.extend(control_line_data)
                if sum(control_line_data) > 0:
                    control_variance_in_gene_num += 1
                    control_variance_in_gene_name_list.append(variance_selected_list[i][0])
            B = sum(tmp_list)
            D = 2 * len(tmp_list) - B

            tmp_list = []
            for i in xrange(len(synonymous_control_data)):
                synonymous_control_line_data = filter(lambda x: type(x) == int, synonymous_control_data[i])
                tmp_list.extend(synonymous_control_line_data)
            DD = sum(tmp_list)

            tmp_list = []
            case_variance_in_gene_num = 0
            case_variance_in_gene_name_list = []
            for i in xrange(len(case_data)):
                case_line_data = filter(lambda x: type(x) == int, case_data[i])
                tmp_list.extend(case_line_data)
                if sum(case_line_data) > 0:
                    case_variance_in_gene_num += 1
                    case_variance_in_gene_name_list.append(variance_selected_list[i][0])
            A = sum(tmp_list)
            C = 2 * len(tmp_list) - A

            tmp_list = []
            for i in xrange(len(synonumous_case_data)):
                synonymous_case_line_data = filter(lambda x: type(x) == int, synonumous_case_data[i])
                tmp_list.extend(synonymous_case_line_data)
            CC = sum(tmp_list)

            control_people_data = [sum(filter(lambda x: type(x) == int, people_line)) for people_line in
                                   zip(*control_data)]
            B1 = len(filter(lambda x: x > 0, control_people_data))
            D1 = len(control_people_data) - B1
            case_people_data = [sum(filter(lambda x: type(x) == int, people_line)) for people_line in zip(*case_data)]
            A1 = len(filter(lambda x: x > 0, case_people_data))
            C1 = len(case_people_data) - A1

            # if A1 + B1 <= 2:
            #     icounter += 1
            #     continue  # the number of people with variance is less than 2. Go to next gene.

            oddsratio, pvalue = stats.fisher_exact([[A, B], [C, D]])
            oddsratio1, pvalue1 = stats.fisher_exact([[A1, B1], [C1, D1]])
            or_new = "inf" if B == 0 or CC == 0 else float(A * DD) / (B * CC)
            ret_str = "{0}{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t" \
                      "{11}\t{12}\t{13}\t{14}\t{15}\t{16}\t{17}\t{18}\t{19}\t{20}\t" \
                      "{21}\t{22}\t{23}\t{24}\t{25}\t{26}\t{27}\t{28}\t{29}\t{30}\t" \
                      "{31}\t{32}\t{33}\t{34}\t{35}\t{36}\t{37}\t{38}\t{39}\t{40}\t" \
                      "{41}\t{42}\n" \
                      "".format(ret_str, gene_id, gene_name, A, B, C, D, pvalue, oddsratio, CC, DD, or_new,
                                A1, B1, C1, D1,
                                pvalue1, oddsratio1,
                                all_variance_in_gene_num,
                                ";".join(all_variance_in_gene_name_list) if len(
                                    all_variance_in_gene_name_list) > 0 else ".",
                                control_variance_in_gene_num,
                                ";".join(control_variance_in_gene_name_list) if len(
                                    control_variance_in_gene_name_list) > 0 else ".",
                                case_variance_in_gene_num,
                                ";".join(case_variance_in_gene_name_list) if len(
                                    case_variance_in_gene_name_list) > 0 else ".",
                                mouse_mean_WT_OFT_TPM_2, mouse_mean_WT_LV_TPM_3, mouse_mean_WT_PA_TPM_3,
                                Interactions_IntAct_, Interactions_BioGRID_,
                                Interactions_ConsensusPathDB_, HIPred,
                                GHIS, gnomAD_pLI,
                                GDI_Phred, Gene_damage_prediction__all_disease_causing_genes_,
                                Essential_gene, Essential_gene_CRISPR,
                                Gene_indispensability_pred, gevir_percentile,
                                loeuf_percentile, virlof_percentile,
                                hhe, "\t".join([str(i) for i in gene_data_list[20:]]))
            icounter += 1
            print("handled {0} / {1} genes".format(icounter, gene_data_len))

        print("handled {0} / {1} genes in total".format(icounter, gene_data_len))
        fp.write(ret_str)
    if should_log:
        logging.debug("all done")


if sys.version_info[0] == 2:
    import urllib2
    import ssl

    ssl._create_default_https_context = ssl._create_unverified_context


def get_ensembl_mouse_gene_name(file_in, output):
    def grep_gene_name(gene_id):
        url = "https://useast.ensembl.org/Mus_musculus/Gene/Summary?db=core;g={0}".format(gene_id)
        try:
            response1 = urllib2.urlopen(url)
        except:
            return "please handle {0} manualy".format(gene_id)

        context = response1.read()
        resp_code = response1.getcode()
        if resp_code != 200:
            return "error {0} please handle {1} manualy".format(resp_code, gene_id)
        if "This identifier is not in the current EnsEMBL database" in context:
            return "This identifier is not in the current EnsEMBL database"
        return re.findall(
            "<h1 class=\"summary-heading\">Gene: (\S+) <span class=\"summary-subhead\">{0}</span></h1>".format(gene_id),
            context)[0]

    icounter = 0
    with open(file_in, "r") as fp_in, open(output, "w") as fp_out:
        while True:
            data_line = fp_in.readline()
            if not data_line:
                break
            if data_line.startswith("#"):
                continue
            data_list = data_line.strip().split("\t")
            gene_id = data_list[0]
            icounter += 1
            print("handling {0} gene_id= {1}".format(icounter, gene_id))
            new_gene_name = grep_gene_name(gene_id)
            data_list.append(new_gene_name)
            fp_out.write("\t".join(data_list) + "\n")


def db_update_mouse_data3(db_file, gene_table, mouse_data, mouse2human_file):
    """
    Find the homologous human gene name corresponding to the mouse gene name through ncbi, and
    then add the mouse data to the database
    @param db_file:
    @param gene_table: gene_table
    @param mouse_data: mouse data
    @param mouse2human_file: The gene mapping relationship between mice and humans 
    detected by ncbi
    @return:
    """

    def load_mouse2human_file(mouse2human_file):
        mouse2human_dict = {}
        with open(mouse2human_file, "r") as fp:
            while True:
                data_line = fp.readline()
                if not data_line:
                    break
                data_list = data_line.strip().split("\t")
                if len(data_list) >= 2:
                    mouse2human_dict[data_list[0]] = [data_list[1]]
                if len(data_list) == 3:
                    mouse2human_dict[data_list[0]].extend(data_list[2].split(";"))
        return mouse2human_dict

    mouse2human_dict = load_mouse2human_file(mouse2human_file)
    with open(mouse_data, "r") as fp:
        mouse_data_list = [i.strip().split("\t") for i in fp.readlines() if not i.startswith("#")]
    human_gene2mouse_data_dict = {}
    for mouse_data_line in mouse_data_list:
        mouse_gene_name = mouse_data_line[6]
        pa = mouse_data_line[7]
        oft = mouse_data_line[8]
        lv = mouse_data_line[9]
        if mouse_gene_name in mouse2human_dict:
            human_gene_list = mouse2human_dict[mouse_gene_name]
            for human_gene in human_gene_list:
                human_gene2mouse_data_dict[human_gene] = [pa, oft, lv, mouse_gene_name]

    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()
    cursor.execute("select id, gene_name from {0}".format(gene_table))
    gene_data = [[id, name] for id, name in cursor.fetchall()]

    db_add_col(cursor, gene_table, "mouse_mean_WT_OFT_TPM_2", "FLOAT")
    db_add_col(cursor, gene_table, "mouse_mean_WT_LV_TPM_3", "FLOAT")
    db_add_col(cursor, gene_table, "mouse_mean_WT_PA_TPM_3", "FLOAT")
    sql_list = []
    for id, human_gene_name in gene_data:
        if human_gene_name not in human_gene2mouse_data_dict:
            sql_list.append("UPDATE {0} SET MGI_mouse_gene=NULL, mouse_mean_WT_LV_TPM_3=NULL, "
                            "mouse_mean_WT_OFT_TPM_2=NULL, mouse_mean_WT_PA_TPM_3=NULL WHERE id='{1}'"
                            "".format(gene_table, id))
            print("{0} doesn't have homology gene name in mouse genom".format(human_gene_name))
            continue
        pa, oft, lv, mouse_gene_name = human_gene2mouse_data_dict[human_gene_name]
        sql_list.append("UPDATE {0} SET MGI_mouse_gene='{1}', mouse_mean_WT_LV_TPM_3='{2}', "
                        "mouse_mean_WT_OFT_TPM_2='{3}', mouse_mean_WT_PA_TPM_3='{4}' WHERE id='{5}'"
                        "".format(gene_table, mouse_gene_name, lv, oft, pa, id))
    for sql in sql_list:
        cursor.execute(sql)

    cursor.close()
    conn.commit()
    conn.close()


def db_update_mouse_data2(db_file, gene_table, mouse_data, target_mouse_gene_file):
    """
    the gene name in db is human gene. human gene name ----> target mouse gene name ----> target mouse data
    @param db_file:
    @param gene_table:
    @param mouse_data: Go to ensamble to find the current mouse gene name
    @param target_mouse_gene_file: Use biomart find the homology gene name in mouse genom.
    @return:
    """

    def load_target_mouse_gene(target_mouse_gene_file):
        with open(target_mouse_gene_file, "r") as fp:
            data = [currnt_line.strip().split("\t") for currnt_line in fp.readlines()]
        human2target_dict = {}
        for human_name, mouse_name in data:
            if human_name not in human2target_dict:
                human2target_dict[human_name] = [mouse_name]
            else:
                human2target_dict[human_name].append(mouse_name)
        return human2target_dict

    def load_mouse_data(mouse_data):
        name2value_dict = {}
        with open(mouse_data, "r") as fp:
            while True:
                mouse_line = fp.readline()
                if not mouse_line:
                    break
                if mouse_line.startswith("#"):
                    continue
                mouse_list = mouse_line.strip().split("\t")
                org_name, pa, oft, lv, new_name = mouse_list[6:11]
                if org_name == new_name:
                    name2value_dict[org_name] = [float(oft), float(lv), float(pa)]
                else:
                    name2value_dict[org_name] = [float(oft), float(lv), float(pa)]
                    name2value_dict[new_name] = [float(oft), float(lv), float(pa)]
        return name2value_dict

    human2target_dict = load_target_mouse_gene(target_mouse_gene_file)
    mouse_name2value_dict = load_mouse_data(mouse_data)
    sql_list = []
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()
    cursor.execute("select id, gene_name from {0}".format(gene_table))
    gene_data = [[id, name] for id, name in cursor.fetchall()]

    db_add_col(cursor, gene_table, "mouse_mean_WT_OFT_TPM_2", "FLOAT")
    db_add_col(cursor, gene_table, "mouse_mean_WT_LV_TPM_3", "FLOAT")
    db_add_col(cursor, gene_table, "mouse_mean_WT_PA_TPM_3", "FLOAT")

    for id, human_gene_name in gene_data:
        if human_gene_name not in human2target_dict:
            sql_list.append("UPDATE {0} SET MGI_mouse_gene=NULL, mouse_mean_WT_LV_TPM_3=NULL, "
                            "mouse_mean_WT_OFT_TPM_2=NULL, mouse_mean_WT_PA_TPM_3=NULL WHERE id='{1}'"
                            "".format(gene_table, id))
            print("{0} doesn't have homology gene name in mouse genom".format(human_gene_name))
            continue
        if len(human2target_dict[human_gene_name]) > 1:
            sql_list.append("UPDATE {0} SET MGI_mouse_gene=NULL, mouse_mean_WT_LV_TPM_3=NULL, "
                            "mouse_mean_WT_OFT_TPM_2=NULL, mouse_mean_WT_PA_TPM_3=NULL WHERE id='{1}'"
                            "".format(gene_table, id))
            print("{0} has more than 1 ({1}) homology gene names in mouse genom".format(human_gene_name,
                                                                                        len(human2target_dict[
                                                                                                human_gene_name])))
            continue
        target_mouse_gene_name = human2target_dict[human_gene_name][0]
        if target_mouse_gene_name not in mouse_name2value_dict:
            sql_list.append("UPDATE {0} SET MGI_mouse_gene='{1}', mouse_mean_WT_LV_TPM_3=NULL, "
                            "mouse_mean_WT_OFT_TPM_2=NULL, mouse_mean_WT_PA_TPM_3=NULL WHERE id='{2}'"
                            "".format(gene_table, target_mouse_gene_name, id))
            print("{0}({1}) doesn't have mouse data".format(human_gene_name, target_mouse_gene_name))
        else:
            oft, lv, pa = mouse_name2value_dict[target_mouse_gene_name]
            sql_list.append("UPDATE {0} SET MGI_mouse_gene='{1}', mouse_mean_WT_LV_TPM_3='{2}', "
                            "mouse_mean_WT_OFT_TPM_2='{3}', mouse_mean_WT_PA_TPM_3='{4}' WHERE id='{5}'"
                            "".format(gene_table, target_mouse_gene_name, lv, oft, pa, id))
    for sql in sql_list:
        cursor.execute(sql)

    cursor.close()
    conn.commit()
    conn.close()


def db_update_mouse_data(db_file, gene_table, mouse_data):
    print("load mouse data ...")
    name2value_dict = {}
    with open(mouse_data, "r") as fp:
        while True:
            mouse_line = fp.readline()
            if not mouse_line:
                break
            if mouse_line.startswith("#"):
                continue
            mouse_list = mouse_line.strip().split("\t")
            org_name, oft, lv, new_name, pa = mouse_list[6:11]
            if org_name == new_name:
                name2value_dict[org_name] = [float(oft), float(lv), float(pa)]
            else:
                name2value_dict[org_name] = [float(oft), float(lv), float(pa)]
                name2value_dict[new_name] = [float(oft), float(lv), float(pa)]
    sql_list = []
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()
    cursor.execute("select id, MGI_mouse_gene from {0}".format(gene_table))
    gene_data = [[id, name] for id, name in cursor.fetchall() if name is not None]
    db_add_col(cursor, gene_table, "mouse_mean_WT_OFT_TPM_2", "FLOAT")
    db_add_col(cursor, gene_table, "mouse_mean_WT_LV_TPM_3", "FLOAT")
    db_add_col(cursor, gene_table, "mouse_mean_WT_PA_TPM_3", "FLOAT")
    print("collecting sqls...")
    for id, name in gene_data:
        if name in name2value_dict:
            oft, lv, pa = name2value_dict[name]
            sql_list.append(
                "UPDATE {0} SET mouse_mean_WT_OFT_TPM_2={1}, mouse_mean_WT_LV_TPM_3={2}, mouse_mean_WT_PA_TPM_3={3} WHERE id='{4}'"
                "".format(gene_table, oft, lv, pa, id))
        else:
            print("name({0}) which is from db is not in mouse data".format(name))

    print("executing sql...")
    for sql in sql_list:
        cursor.execute(sql)

    cursor.close()
    conn.commit()
    conn.close()


def gene_set2gene_set_with_name(db_file, gene_set_file):
    id2name_dict = {}

    def id2key(cursor, id, id2name_dict):
        if id in id2name_dict:
            name = id2name_dict[id]
            ret = "{0}_{1}".format(name, id)
            print(ret)
            return ret
        cmd_str = "select gene_name from gene_table where gene_id='{}'".format(id)
        cursor.execute(cmd_str)
        sql_ret = cursor.fetchall()
        if len(sql_ret) == 0:
            ret = "None_{0}".format(id)
            print(ret)
            id2name_dict[id] = "None"
            return ret
        name = sql_ret[0][0]
        id2name_dict[id] = name
        ret = "{0}_{1}".format(name, id)
        print(ret)
        return ret

    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()
    with open(gene_set_file, "r") as fp, open(gene_set_file + ".withname", "w") as fp_out:
        while True:
            gene_set_line = fp.readline()
            if not gene_set_line:
                break
            if len(gene_set_line.strip()) == 0:
                continue
            gene_set_list = gene_set_line.strip().split("\t")
            name_id_list = [id2key(cursor, i, id2name_dict) for i in filter(lambda x: len(x) > 0, gene_set_list[1:])]
            fp_out.write("{0}\t{1}\n".format(gene_set_list[0], "\t".join(name_id_list)))
    cursor.close()
    conn.commit()
    conn.close()


def test_gene_binom_test(gene_set_file, fet_result, mode):
    # http://yulab-smu.top/clusterProfiler-book/chapter2.html#over-representation-analysis
    from scipy.stats import binom_test
    # logging.basicConfig(filename="test_gene_binom_test.log", level=logging.DEBUG, format=log_format, filemode="a")
    # logging.debug("gene_set_file={}".format(gene_set_file))
    # logging.debug("fet_result={}".format(fet_result))
    # logging.debug("mode={}".format(mode))
    # logging.debug("output={}".format(output))
    N = 19995.0
    # print("loading gene set data...")
    # conn = sqlite3.connect(db_file)
    # cursor = conn.cursor()
    gene_set_dict = {}
    with open(gene_set_file, "r") as fp:
        while True:
            gene_set_line = fp.readline()
            if not gene_set_line:
                break
            if len(gene_set_line.strip()) == 0:
                continue
            gene_set_list = gene_set_line.strip().split("\t")
            gene_set_dict[gene_set_list[0]] = filter(lambda x: len(x) > 0, gene_set_list[1:])
    # cursor.close()
    # conn.commit()
    # conn.close()
    # print(gene_set_dict)
    # print("loading fet result data...")
    gene_list = []
    NN = 0
    gene_N_set = set([])
    with open(fet_result, "r") as fp:
        while True:
            data_line = fp.readline()
            if not data_line:
                break
            if data_line.startswith("#"):
                continue
            data_list = data_line.strip().split("\t")
            gene_id = data_list[0]
            gene_name = data_list[1]
            gene_id = "{0}_{1}".format(gene_name, gene_id)
            B1 = int(data_list[12])
            A1 = int(data_list[11])
            if A1 + B1 > 0:
                NN += 1
                gene_N_set.add(gene_id)

            p_value1 = float(data_list[15])
            or1 = float(data_list[16])
            if mode == "A3.2.B1.0":
                if (B1 <= 0 and A1 >= 2) or (B1 <= 1 and A1 >= 3):
                    gene_list.append(gene_id)
            if mode == "A4.3.B1.0":
                if (B1 <= 0 and A1 >= 3) or (B1 <= 1 and A1 >= 4):
                    gene_list.append(gene_id)
            if mode == "A3.B0":
                if B1 <= 0 and A1 >= 3:
                    gene_list.append(gene_id)
            if mode == "A4.B0":
                if B1 <= 0 and A1 >= 4:
                    gene_list.append(gene_id)
            if mode == "A3.B1":
                if B1 <= 1 and A1 >= 3:
                    gene_list.append(gene_id)
            if mode == "A4.B1":
                if B1 <= 1 and A1 >= 4:
                    gene_list.append(gene_id)
            if mode == "A2.B0":
                if B1 <= 0 and A1 >= 2:
                    gene_list.append(gene_id)
            elif mode == "p1_case":
                if p_value1 <= 0.05 and or1 > 1:
                    gene_list.append(gene_id)
            if mode == "A1.0.B3.2":
                if (A1 <= 0 and B1 >= 2) or (A1 <= 1 and B1 >= 3):
                    gene_list.append(gene_id)
            if mode == "A1.0.B4.3":
                if (A1 <= 0 and B1 >= 3) or (A1 <= 1 and B1 >= 4):
                    gene_list.append(gene_id)
            if mode == "A0.B3":
                if A1 <= 0 and B1 >= 3:
                    gene_list.append(gene_id)
            if mode == "A1.B4":
                if A1 <= 1 and B1 >= 4:
                    gene_list.append(gene_id)
            if mode == "A1.B3":
                if A1 <= 1 and B1 >= 3:
                    gene_list.append(gene_id)
            if mode == "A0.B2":
                if A1 <= 0 and B1 >= 2:
                    gene_list.append(gene_id)
            if mode == "A0.B4":
                if A1 <= 0 and B1 >= 4:
                    gene_list.append(gene_id)
            elif mode == "p1_control":
                if p_value1 <= 0.05 and or1 < 1:
                    gene_list.append(gene_id)
            else:
                RuntimeError("mode should be B1 or p1")
    n = len(gene_list)
    for gene_set_name in gene_set_dict:
        gene_set_list = gene_set_dict[gene_set_name]
        gene_set_set = set(gene_set_list)
        gene_set_and_N_set = gene_set_set & gene_N_set
        MM = len(gene_set_and_N_set)
        m_list = list(set(gene_list) & gene_set_set)
        mm = len(m_list)
        # mmm = len(set(gene_list) & gene_set_and_N_set)
        M = len(gene_set_list)
        m = len(filter(lambda x: x in gene_set_list, gene_list))
        p_value = binom_test(m, n, float(M) / N, alternative='greater')
        p1_value = binom_test(mm, n, float(MM) / NN, alternative='greater')
        ret_str = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}" \
                  "".format(os.path.basename(fet_result), gene_set_name, p_value, m, n, M, N,
                            p1_value, mm, n, MM, NN, ";".join(m_list), ";".join(gene_list))
        print(ret_str)


def fisher_gene_name_binom_test(database, sample_restrict, gene_set, path, script):
    table_name_element_dict = {"": "",
                               " AND (annovar = 1)": "annovar",
                               " AND (bystro = 1)": "bystro",
                               " AND (dmis = 1)": "dmis",
                               " AND (dsplicing = 1)": "dsplicing",
                               " AND (bystro=1 or annovar=1 or vep=1 or dmis=1 or dsplicing=1 or spliceAI=1)": "all.6",
                               " AND (spliceAI = 1)": "spliceAI",
                               " AND (vep = 1)": "vep",
                               " AND (annovar = 1 or bystro = 1 or vep = 1)": "LOF",
                               " AND (dsplicing = 1 or spliceAI = 1)": "dsplicing.all",
                               " AND (annovar = 1 and bystro = 1)": "annovar.bystro",
                               "tof6": "tof",
                               "CTD": "CTD",
                               "heart6": "heart6",
                               "bystro_sampleMaf <= 0.01": "maf.01",
                               "bystro_sampleMaf <= 0.05": "maf.05",
                               "bystro_sampleMaf <= 0.1": "maf.1",
                               "bystro_sampleMaf <= 1": "",
                               # "bystro_sampleMaf <= 0.2": "maf.2",
                               # "bystro_sampleMaf <= 0.3": "maf.3",
                               " AND (bystro_cadd>=10)": "cadd10",
                               " AND (bystro_cadd>=15)": "cadd15",
                               " AND (bystro_cadd>=20)": "cadd20",
                               " AND bystro_phastCons >= 0.4": "Cons.4",
                               " AND bystro_phastCons >= 0.5": "Cons.5",
                               " AND bystro_phastCons >= 0.6": "Cons.6",
                               " AND bystro_phyloP >= -1": "loPn1",
                               " AND bystro_phyloP >= 0": "loP0",
                               " AND bystro_phyloP >= 1": "loP1",
                               " AND (ccrs >= 95)": "ccrs95",
                               " AND (ccrs >= 90)": "ccrs90",
                               " AND (ccrs >= 85)": "ccrs85",
                               " AND (ccrs >= 80)": "ccrs80",
                               " AND (domain_limbr <= {})".format(domain_limbr_quartile_dict[40]): "dlimbr40",
                               " AND (domain_limbr <= {})".format(domain_limbr_quartile_dict[50]): "dlimbr50",
                               " AND (domain_limbr <= {})".format(domain_limbr_quartile_dict[60]): "dlimbr60",
                               " AND (domain_limbr <= {})".format(domain_limbr_quartile_dict[70]): "dlimbr70",
                               " AND (exone_limbr <= {})".format(exone_limbr_quartile_dict[40]): "elimbr40",
                               " AND (exone_limbr <= {})".format(exone_limbr_quartile_dict[50]): "elimbr50",
                               " AND (exone_limbr <= {})".format(exone_limbr_quartile_dict[60]): "elimbr60",
                               " AND (exone_limbr <= {})".format(exone_limbr_quartile_dict[70]): "elimbr70",
                               " AND (is_ccds = 1)": "ccds1"}

    pp = Popen(["mkdir -p {}".format(path)], shell=True)
    pp.wait()
    annotator_list = [" AND (bystro=1 or annovar=1 or vep=1 or dmis=1 or dsplicing=1 or spliceAI=1)",
                      " AND (annovar = 1)",
                      " AND (bystro = 1)",
                      " AND (dmis = 1)",
                      " AND (vep = 1)",
                      " AND (annovar = 1 or bystro = 1 or vep = 1)",
                      " AND (dsplicing = 1 or spliceAI = 1)",
                      " AND (annovar = 1 and bystro = 1)"]
    synonymous_annotator_list = ["bystro=1 or vep=1 or annovar=1",
                                 "annovar=1",
                                 "bystro=1",
                                 "annovar=1",
                                 "vep=1",
                                 "bystro=1 or vep=1 or annovar=1",
                                 "annovar=1",
                                 "annovar=1 and bystro=1"]
    annotator2synonymous_annotator_dict = dict(zip(annotator_list, synonymous_annotator_list))
    icounter = 1
    icounter2 = 0
    script_str = ''
    for phenotype in ["tof6", "CTD", "heart6"]:
        for annotator in annotator_list:
            for freq in ["bystro_sampleMaf <= 1"]:
                for cadd in ["", " AND (bystro_cadd>=10)", " AND (bystro_cadd>=15)", " AND (bystro_cadd>=20)"]:
                    for ph in ["",
                               " AND bystro_phastCons >= 0.4",
                               " AND bystro_phastCons >= 0.5",
                               " AND bystro_phastCons >= 0.6",
                               " AND bystro_phyloP >= -1",
                               " AND bystro_phyloP >= 0",
                               " AND bystro_phyloP >= 1"]:
                        for regional_constraint in ["",
                                                    " AND (ccrs >= 95)",
                                                    " AND (ccrs >= 90)",
                                                    " AND (ccrs >= 85)",
                                                    " AND (ccrs >= 80)",
                                                    " AND (domain_limbr <= {})".format(domain_limbr_quartile_dict[40]),
                                                    " AND (domain_limbr <= {})".format(domain_limbr_quartile_dict[50]),
                                                    " AND (domain_limbr <= {})".format(domain_limbr_quartile_dict[60]),
                                                    " AND (domain_limbr <= {})".format(domain_limbr_quartile_dict[70]),
                                                    " AND (exone_limbr <= {})".format(exone_limbr_quartile_dict[40]),
                                                    " AND (exone_limbr <= {})".format(exone_limbr_quartile_dict[50]),
                                                    " AND (exone_limbr <= {})".format(exone_limbr_quartile_dict[60]),
                                                    " AND (exone_limbr <= {})".format(exone_limbr_quartile_dict[70])]:
                            for is_ccds in ["", " AND (is_ccds = 1)"]:
                                variance_restrict = "{0}{1}{2}{3}{4}{5}".format(freq, annotator, cadd, ph,
                                                                                regional_constraint, is_ccds)
                                name_element_list = [table_name_element_dict[phenotype],
                                                     table_name_element_dict[annotator],
                                                     table_name_element_dict[freq],
                                                     table_name_element_dict[cadd],
                                                     table_name_element_dict[ph],
                                                     table_name_element_dict[regional_constraint],
                                                     table_name_element_dict[is_ccds]]
                                output_name = "{0}.fet" \
                                              "".format("_".join(filter(lambda x: len(x) > 0, name_element_list)))
                                synonymous_annotator = annotator2synonymous_annotator_dict[annotator]
                                # output_name = os.path.join(path, output_name)

                                if script_str == "":
                                    script_str = """#!/bin/bash
#$ -j y
#$ -S /bin/bash
#$ -cwd
#$ -N fg{6}
#$ -pe smp 1
#$ -l h_vmem=8G
#$ -m bes
# -q all.q
export LC_ALL=C
export MALLOC_ARENA_MAX=4
module load python/2.7.15/gcc.4.4.7
module load sqlite3/3.8.11/gcc.4.4.7
time=`date`
echo "==START $time =="
{9} {0} {1} sampleChdPhenotype "{2}" gene_table "" "{3}" "{4}" {5} {6} {7} {8}
echo {5} is done
time=`date`
echo == $time ==
""".format(database, phenotype, sample_restrict, variance_restrict, synonymous_annotator, output_name, icounter2,
           gene_set, path, script)
                                else:
                                    script_str += "{9} {0} {1} sampleChdPhenotype \"{2}\" gene_table \"\" " \
                                                  "\"{3}\" \"{4}\" {5} {6} {7} {8}" \
                                                  "\necho {5} is done" \
                                                  "\ntime=`date`\necho == $time ==\n" \
                                                  "".format(database, phenotype,
                                                            sample_restrict, variance_restrict,
                                                            synonymous_annotator, output_name,
                                                            icounter2, gene_set,
                                                            path, script)

                                if icounter == 53:
                                    shell_file_name = os.path.join(path, "qsub{}.sh".format(icounter2))
                                    with open(shell_file_name, "w") as fp:
                                        fp.write(script_str + "\ndate; echo \"==END==\"")
                                    while True:
                                        if os.access(shell_file_name, os.R_OK):
                                            break
                                        time.sleep(1)
                                    print("qsubing {} batch".format(icounter2 + 1))
                                    pp = Popen(["qsub {}".format(shell_file_name)], shell=True)
                                    pp.wait()
                                    script_str = ""
                                    icounter = 0
                                    icounter2 += 1
                                icounter += 1
    if len(script_str) > 0:
        shell_file_name = os.path.join(path, "qsub{}.sh".format(icounter2))
        with open(shell_file_name, "w") as fp:
            fp.write(script_str + "\ndate; echo \"==END==\"")
        while True:
            if os.access(shell_file_name, os.R_OK):
                break
            time.sleep(1)
        print("qsubing {} batch".format(icounter2 + 1))
        pp = Popen(["qsub {}".format(shell_file_name)], shell=True)
        pp.wait()

    print(
        "All the jobs has been submitted. Job number = {}\nIf any job has problem, kill it. And qsub the corresponding sh file".format(
            icounter2 + 1))

def fisher_gene_name_binom_test_slurm(database, sample_restrict, gene_set, path, script):
    table_name_element_dict = {"": "",
                               " AND (annovar = 1)": "annovar",
                               " AND (bystro = 1)": "bystro",
                               " AND (dmis = 1)": "dmis",
                               " AND (dsplicing = 1)": "dsplicing",
                               " AND (bystro=1 or annovar=1 or vep=1 or dmis=1 or dsplicing=1 or spliceAI=1)": "all.6",
                               " AND (spliceAI = 1)": "spliceAI",
                               " AND (vep = 1)": "vep",
                               " AND (annovar = 1 or bystro = 1 or vep = 1)": "LOF",
                               " AND (dsplicing = 1 or spliceAI = 1)": "dsplicing.all",
                               " AND (annovar = 1 and bystro = 1)": "annovar.bystro",
                               "tof6": "tof",
                               "CTD": "CTD",
                               "heart6": "heart6",
                               "bystro_sampleMaf <= 0.01": "maf.01",
                               "bystro_sampleMaf <= 0.05": "maf.05",
                               "bystro_sampleMaf <= 0.1": "maf.1",
                               "bystro_sampleMaf <= 1": "",
                               # "bystro_sampleMaf <= 0.2": "maf.2",
                               # "bystro_sampleMaf <= 0.3": "maf.3",
                               " AND (bystro_cadd>=10)": "cadd10",
                               " AND (bystro_cadd>=15)": "cadd15",
                               " AND (bystro_cadd>=20)": "cadd20",
                               " AND bystro_phastCons >= 0.4": "Cons.4",
                               " AND bystro_phastCons >= 0.5": "Cons.5",
                               " AND bystro_phastCons >= 0.6": "Cons.6",
                               " AND bystro_phyloP >= -1": "loPn1",
                               " AND bystro_phyloP >= 0": "loP0",
                               " AND bystro_phyloP >= 1": "loP1",
                               " AND (ccrs >= 95)": "ccrs95",
                               " AND (ccrs >= 90)": "ccrs90",
                               " AND (ccrs >= 85)": "ccrs85",
                               " AND (ccrs >= 80)": "ccrs80",
                               " AND (domain_limbr <= {})".format(domain_limbr_quartile_dict[40]): "dlimbr40",
                               " AND (domain_limbr <= {})".format(domain_limbr_quartile_dict[50]): "dlimbr50",
                               " AND (domain_limbr <= {})".format(domain_limbr_quartile_dict[60]): "dlimbr60",
                               " AND (domain_limbr <= {})".format(domain_limbr_quartile_dict[70]): "dlimbr70",
                               " AND (exone_limbr <= {})".format(exone_limbr_quartile_dict[40]): "elimbr40",
                               " AND (exone_limbr <= {})".format(exone_limbr_quartile_dict[50]): "elimbr50",
                               " AND (exone_limbr <= {})".format(exone_limbr_quartile_dict[60]): "elimbr60",
                               " AND (exone_limbr <= {})".format(exone_limbr_quartile_dict[70]): "elimbr70",
                               " AND (is_ccds = 1)": "ccds1"}

    pp = Popen(["mkdir -p {}".format(path)], shell=True)
    pp.wait()
    annotator_list = [" AND (bystro=1 or annovar=1 or vep=1 or dmis=1 or dsplicing=1 or spliceAI=1)",
                      " AND (annovar = 1)",
                      " AND (bystro = 1)",
                      " AND (dmis = 1)",
                      " AND (vep = 1)",
                      " AND (annovar = 1 or bystro = 1 or vep = 1)",
                      " AND (dsplicing = 1 or spliceAI = 1)",
                      " AND (annovar = 1 and bystro = 1)"]
    synonymous_annotator_list = ["bystro=1 or vep=1 or annovar=1",
                                 "annovar=1",
                                 "bystro=1",
                                 "annovar=1",
                                 "vep=1",
                                 "bystro=1 or vep=1 or annovar=1",
                                 "annovar=1",
                                 "annovar=1 and bystro=1"]
    annotator2synonymous_annotator_dict = dict(zip(annotator_list, synonymous_annotator_list))
    icounter = 1
    icounter2 = 0
    script_str = ''
    for phenotype in ["tof6", "CTD", "heart6"]:
        for annotator in annotator_list:
            for freq in ["bystro_sampleMaf <= 1"]:
                for cadd in ["", " AND (bystro_cadd>=10)", " AND (bystro_cadd>=15)", " AND (bystro_cadd>=20)"]:
                    for ph in ["",
                               " AND bystro_phastCons >= 0.4",
                               " AND bystro_phastCons >= 0.5",
                               " AND bystro_phastCons >= 0.6",
                               " AND bystro_phyloP >= -1",
                               " AND bystro_phyloP >= 0",
                               " AND bystro_phyloP >= 1"]:
                        for regional_constraint in ["",
                                                    " AND (ccrs >= 95)",
                                                    " AND (ccrs >= 90)",
                                                    " AND (ccrs >= 85)",
                                                    " AND (ccrs >= 80)",
                                                    " AND (domain_limbr <= {})".format(domain_limbr_quartile_dict[40]),
                                                    " AND (domain_limbr <= {})".format(domain_limbr_quartile_dict[50]),
                                                    " AND (domain_limbr <= {})".format(domain_limbr_quartile_dict[60]),
                                                    " AND (domain_limbr <= {})".format(domain_limbr_quartile_dict[70]),
                                                    " AND (exone_limbr <= {})".format(exone_limbr_quartile_dict[40]),
                                                    " AND (exone_limbr <= {})".format(exone_limbr_quartile_dict[50]),
                                                    " AND (exone_limbr <= {})".format(exone_limbr_quartile_dict[60]),
                                                    " AND (exone_limbr <= {})".format(exone_limbr_quartile_dict[70])]:
                            for is_ccds in ["", " AND (is_ccds = 1)"]:
                                variance_restrict = "{0}{1}{2}{3}{4}{5}".format(freq, annotator, cadd, ph,
                                                                                regional_constraint, is_ccds)
                                name_element_list = [table_name_element_dict[phenotype],
                                                     table_name_element_dict[annotator],
                                                     table_name_element_dict[freq],
                                                     table_name_element_dict[cadd],
                                                     table_name_element_dict[ph],
                                                     table_name_element_dict[regional_constraint],
                                                     table_name_element_dict[is_ccds]]
                                output_name = "{0}.fet" \
                                              "".format("_".join(filter(lambda x: len(x) > 0, name_element_list)))
                                synonymous_annotator = annotator2synonymous_annotator_dict[annotator]
                                # output_name = os.path.join(path, output_name)

                                if script_str == "":
                                    script_str = """#!/bin/bash
#SBATCH -p unlimited
#SBATCH -J fetbinom
#SBATCH -n 1
#SBATCH --mem=30GB
#SBATCH -o fetbinom_%A.log
#SBATCH -t UNLIMITED

time=`date`
echo "==START $time =="
{9} {0} {1} sampleChdPhenotype "{2}" gene_table "" "{3}" "{4}" {5} {6} {7} {8}
echo {5} is done
time=`date`
echo == $time ==
""".format(database, phenotype, sample_restrict, variance_restrict, synonymous_annotator, output_name, icounter2,
           gene_set, path, script)
                                else:
                                    script_str += "{9} {0} {1} sampleChdPhenotype \"{2}\" gene_table \"\" " \
                                                  "\"{3}\" \"{4}\" {5} {6} {7} {8}" \
                                                  "\necho {5} is done" \
                                                  "\ntime=`date`\necho == $time ==\n" \
                                                  "".format(database, phenotype,
                                                            sample_restrict, variance_restrict,
                                                            synonymous_annotator, output_name,
                                                            icounter2, gene_set,
                                                            path, script)

                                if icounter == 180:
                                    shell_file_name = os.path.join(path, "sbatch{}.sh".format(icounter2))
                                    with open(shell_file_name, "w") as fp:
                                        fp.write(script_str + "\ndate; echo \"==END==\"")
                                    while True:
                                        if os.access(shell_file_name, os.R_OK):
                                            break
                                        time.sleep(1)
                                    print("submitting {} batch".format(icounter2 + 1))
                                    pp = Popen(["sbatch {}".format(shell_file_name)], shell=True)
                                    pp.wait()
                                    script_str = ""
                                    icounter = 0
                                    icounter2 += 1
                                icounter += 1
    if len(script_str) > 0:
        shell_file_name = os.path.join(path, "sbatch{}.sh".format(icounter2))
        with open(shell_file_name, "w") as fp:
            fp.write(script_str + "\ndate; echo \"==END==\"")
        while True:
            if os.access(shell_file_name, os.R_OK):
                break
            time.sleep(1)
        print("submitting {} batch".format(icounter2 + 1))
        pp = Popen(["sbatch {}".format(shell_file_name)], shell=True)
        pp.wait()

    print(
        "All the jobs has been submitted. Job number = {}\nIf any job has problem, kill it. And qsub the corresponding sh file".format(
            icounter2 + 1))

def fisher_gene_name_binom_test_synonymous(database, sample_restrict, gene_set, path):
    table_name_element_dict = {"": "",
                               " AND (annovar = 1)": "annovar",
                               " AND (bystro = 1)": "bystro",
                               " AND (dmis = 1)": "dmis",
                               " AND (dsplicing = 1)": "dsplicing",
                               " AND (bystro=1 or annovar=1 or vep=1 or dmis=1 or dsplicing=1 or spliceAI=1)": "all.6",
                               " AND (spliceAI = 1)": "spliceAI",
                               " AND (vep = 1)": "vep",
                               " AND (annovar = 1 or bystro = 1 or vep = 1)": "LOF",
                               " AND (dsplicing = 1 or spliceAI = 1)": "dsplicing.all",
                               " AND (annovar = 1 and bystro = 1)": "annovar.bystro",
                               "tof6": "tof",
                               "CTD": "CTD",
                               "heart6": "heart6",
                               "bystro_sampleMaf <= 0.01": "maf.01",
                               "bystro_sampleMaf <= 0.05": "maf.05",
                               "bystro_sampleMaf <= 0.1": "maf.1",
                               "bystro_sampleMaf <= 1": "",
                               # "bystro_sampleMaf <= 0.2": "maf.2",
                               # "bystro_sampleMaf <= 0.3": "maf.3",
                               " AND (bystro_cadd>=10)": "cadd10",
                               " AND (bystro_cadd>=15)": "cadd15",
                               " AND (bystro_cadd>=20)": "cadd20",
                               " AND bystro_phastCons >= 0.4": "Cons.4",
                               " AND bystro_phastCons >= 0.5": "Cons.5",
                               " AND bystro_phastCons >= 0.6": "Cons.6",
                               " AND bystro_phyloP >= -1": "loPn1",
                               " AND bystro_phyloP >= 0": "loP0",
                               " AND bystro_phyloP >= 1": "loP1",
                               " AND (ccrs >= 95)": "ccrs95",
                               " AND (ccrs >= 90)": "ccrs90",
                               " AND (ccrs >= 85)": "ccrs85",
                               " AND (ccrs >= 80)": "ccrs80",
                               " AND (domain_limbr <= {})".format(domain_limbr_quartile_dict[40]): "dlimbr40",
                               " AND (domain_limbr <= {})".format(domain_limbr_quartile_dict[50]): "dlimbr50",
                               " AND (domain_limbr <= {})".format(domain_limbr_quartile_dict[60]): "dlimbr60",
                               " AND (domain_limbr <= {})".format(domain_limbr_quartile_dict[70]): "dlimbr70",
                               " AND (exone_limbr <= {})".format(exone_limbr_quartile_dict[40]): "elimbr40",
                               " AND (exone_limbr <= {})".format(exone_limbr_quartile_dict[50]): "elimbr50",
                               " AND (exone_limbr <= {})".format(exone_limbr_quartile_dict[60]): "elimbr60",
                               " AND (exone_limbr <= {})".format(exone_limbr_quartile_dict[70]): "elimbr70",
                               " AND (is_ccds = 1)": "ccds1"}

    pp = Popen(["mkdir -p {}".format(path)], shell=True)
    pp.wait()
    annotator_list = [" AND (annovar = 1)",
                      " AND (bystro = 1)",
                      " AND (vep = 1)",
                      " AND (annovar = 1 or bystro = 1 or vep = 1)",
                      " AND (annovar = 1 and bystro = 1)"]
    synonymous_annotator_list = ["annovar=1",
                                 "bystro=1",
                                 "vep=1",
                                 "bystro=1 or vep=1 or annovar=1",
                                 "annovar=1 and bystro=1"]
    annotator2synonymous_annotator_dict = dict(zip(annotator_list, synonymous_annotator_list))
    icounter = 1
    icounter2 = 0
    script_str = ''
    for phenotype in ["tof6", "CTD", "heart6"]:
        for annotator in annotator_list:
            for freq in ["bystro_sampleMaf <= 1"]:
                for cadd in ["", " AND (bystro_cadd>=10)", " AND (bystro_cadd>=15)", " AND (bystro_cadd>=20)"]:
                    for ph in ["",
                               " AND bystro_phastCons >= 0.4",
                               " AND bystro_phastCons >= 0.5",
                               " AND bystro_phastCons >= 0.6",
                               " AND bystro_phyloP >= -1",
                               " AND bystro_phyloP >= 0",
                               " AND bystro_phyloP >= 1"]:
                        for regional_constraint in ["",
                                                    " AND (ccrs >= 95)",
                                                    " AND (ccrs >= 90)",
                                                    " AND (ccrs >= 85)",
                                                    " AND (ccrs >= 80)",
                                                    " AND (domain_limbr <= {})".format(domain_limbr_quartile_dict[40]),
                                                    " AND (domain_limbr <= {})".format(domain_limbr_quartile_dict[50]),
                                                    " AND (domain_limbr <= {})".format(domain_limbr_quartile_dict[60]),
                                                    " AND (domain_limbr <= {})".format(domain_limbr_quartile_dict[70]),
                                                    " AND (exone_limbr <= {})".format(exone_limbr_quartile_dict[40]),
                                                    " AND (exone_limbr <= {})".format(exone_limbr_quartile_dict[50]),
                                                    " AND (exone_limbr <= {})".format(exone_limbr_quartile_dict[60]),
                                                    " AND (exone_limbr <= {})".format(exone_limbr_quartile_dict[70])]:
                            for is_ccds in ["", " AND (is_ccds = 1)"]:
                                variance_restrict = "{0}{1}{2}{3}{4}{5}".format(freq, annotator, cadd, ph,
                                                                                regional_constraint, is_ccds)
                                name_element_list = [table_name_element_dict[phenotype],
                                                     table_name_element_dict[annotator],
                                                     table_name_element_dict[freq],
                                                     table_name_element_dict[cadd],
                                                     table_name_element_dict[ph],
                                                     table_name_element_dict[regional_constraint],
                                                     table_name_element_dict[is_ccds]]
                                output_name = "{0}.fet" \
                                              "".format("_".join(filter(lambda x: len(x) > 0, name_element_list)))
                                synonymous_annotator = annotator2synonymous_annotator_dict[annotator]
                                # output_name = os.path.join(path, output_name)

                                if script_str == "":
                                    script_str = """#!/bin/bash
#$ -j y
#$ -S /bin/bash
#$ -cwd
#$ -N fg{6}
#$ -pe smp 1
#$ -l h_vmem=8G
#$ -m bes
# -q all.q
export LC_ALL=C
export MALLOC_ARENA_MAX=4
module load python/2.7.15/gcc.4.4.7
module load sqlite3/3.8.11/gcc.4.4.7
time=`date`
echo "==START $time =="
fet_binom_synonymous.sh {0} {1} sampleChdPhenotype "{2}" gene_table "" "{3}" "{4}" {5} {6} {7} {8}
echo {5} is done
time=`date`
echo == $time ==
""".format(database, phenotype, sample_restrict, variance_restrict, synonymous_annotator, output_name, icounter2,
           gene_set, path)
                                else:
                                    script_str += "fet_binom_synonymous.sh {0} {1} sampleChdPhenotype \"{2}\" gene_table \"\" " \
                                                  "\"{3}\" \"{4}\" {5} {6} {7} {8}" \
                                                  "\necho {5} is done" \
                                                  "\ntime=`date`\necho == $time ==\n" \
                                                  "".format(database, phenotype,
                                                            sample_restrict, variance_restrict,
                                                            synonymous_annotator, output_name,
                                                            icounter2, gene_set, path)

                                if icounter == 53:
                                    shell_file_name = os.path.join(path, "qsub{}.sh".format(icounter2))
                                    with open(shell_file_name, "w") as fp:
                                        fp.write(script_str + "\ndate; echo \"==END==\"")
                                    while True:
                                        if os.access(shell_file_name, os.R_OK):
                                            break
                                        time.sleep(1)
                                    print("qsubing {} batch".format(icounter2 + 1))
                                    pp = Popen(["qsub {}".format(shell_file_name)], shell=True)
                                    pp.wait()
                                    script_str = ""
                                    icounter = 0
                                    icounter2 += 1
                                icounter += 1
    if len(script_str) > 0:
        shell_file_name = os.path.join(path, "qsub{}.sh".format(icounter2))
        with open(shell_file_name, "w") as fp:
            fp.write(script_str + "\ndate; echo \"==END==\"")
        while True:
            if os.access(shell_file_name, os.R_OK):
                break
            time.sleep(1)
        print("qsubing {} batch".format(icounter2 + 1))
        pp = Popen(["qsub {}".format(shell_file_name)], shell=True)
        pp.wait()

    print(
        "All the jobs has been submitted. Job number = {}\nIf any job has problem, kill it. And qsub the corresponding sh file".format(
            icounter2 + 1))


def binom_pileup(path, output):
    def parse_filename(filename):
        # type: (str) -> list[str]
        phenotype = ""
        filename = filename[:-6]
        if filename.startswith("CTD_"):
            phenotype = "CTD"
            filename = filename[4:]
        elif filename.startswith("heart6_"):
            phenotype = "heart6"
            filename = filename[7:]
        elif filename.startswith("tof_"):
            phenotype = "tof"
            filename = filename[4:]
        else:
            return None
        # print(filename)
        binom_c = ""
        fet_c = ""
        if filename.endswith("A3.2.B1.0"):
            binom_c = "A3.2.B1.0"
            fet_c = filename[:-10]
        elif filename.endswith("A3.B0"):
            binom_c = "A3.B0"
            fet_c = filename[:-6]
        elif filename.endswith("A3.B1"):
            binom_c = "A3.B1"
            fet_c = filename[:-6]
        elif filename.endswith("A2.B0"):
            binom_c = "A2.B0"
            fet_c = filename[:-6]
        elif filename.endswith("A1.0.B3.2"):
            binom_c = "A1.0.B3.2"
            fet_c = filename[:-10]
        elif filename.endswith("A1.B3"):
            binom_c = "A1.B3"
            fet_c = filename[:-6]
        elif filename.endswith("A0.B3"):
            binom_c = "A0.B3"
            fet_c = filename[:-6]
        elif filename.endswith("A0.B2"):
            binom_c = "A0.B2"
            fet_c = filename[:-6]
        else:
            return None
        return [phenotype, binom_c, fet_c]

    g = os.walk(path)
    with open(output, "w") as fp_out:
        fp_out.write("#https://yulab-smu.github.io/clusterProfiler-book/chapter2.html\n"
                     "#gene_set_name\tfet_c\tphenotype\tbinom_c\tp1\tm\tn\tM1\tN1\tp2\tm\tn\tM2\tN2\tm_name\tn_name\n")
        for path, dir_list, file_list in g:
            for file_name in file_list:
                if not file_name.endswith(".binom"):
                    continue
                ret = parse_filename(file_name)
                if ret is None:
                    continue
                phenotype, binom_c, fet_c = ret
                with open(os.path.join(path, file_name), "r") as fp_in:
                    while True:
                        data_line = fp_in.readline()
                        if not data_line:
                            break
                        if len(data_line.strip()) == 0:
                            continue
                        data_list = data_line.strip("\n").split("\t")
                        if len(data_list) < 3:
                            continue
                        gene_set_name = data_list[1]
                        value_list = data_list[2:12]
                        tmp_str = "{0}\t{1}\t{2}\t{3}\t{4}".format(gene_set_name, fet_c, phenotype,
                                                                   binom_c, "\t".join(value_list))
                        if len(data_list) != 14:
                            print(data_list)
                        assert len(data_list) == 14
                        tmp_str += "\t{0}\t{1}\n".format(data_list[12], data_list[13])
                        fp_out.write(tmp_str)
                        # print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}"
                        #       "".format(os.path.basename(fet_result), gene_set_name, p_value, m, n, M, N,
                        #                 p1_value, mm, n, MM, NN, ";".join(m_list), ";".join(gene_list)))


def file_line_num(file_name, comment_chr):
    # cmd_str = 'grep -v "^{0}" {1} | wc -l'.format(comment_chr, file_name)
    if not os.path.exists(file_name):
        return 0
    return int(os.popen('grep -v "^{0}" {1} | wc -l'.format(comment_chr, file_name)).read().split()[0])


def flat_pileup(pileup_in, out_put, heart6_binom_c1, ctd_binom_c, tof_binom_c, heart6_binom_c2):
    logging.basicConfig(filename="flat_pileup.log", level=logging.DEBUG, format=log_format, filemode="w")
    logging.debug("pileup_in={}".format(pileup_in))
    logging.debug("ctd_binom_c={}".format(ctd_binom_c))
    logging.debug("tof_binom_c={}".format(tof_binom_c))
    logging.debug("heart6_binom_c2={}".format(heart6_binom_c2))
    logging.debug("heart6_binom_c1={}".format(heart6_binom_c1))
    pre_key = ""
    heart6_value1 = ""
    ctd_value = ""
    heart6_value2 = ""
    tof_value = ""
    logging.debug("calculating total line num of pileup_in ...")
    total_line = float(file_line_num(pileup_in, "+"))
    with open(pileup_in, "r") as fp_in, open(out_put, "w") as fp_out:
        ret_str = "##heart6_binom_c1={4}\n" \
                  "##ctd_binom_c={0}\n" \
                  "##tof_binom_c={1}\n" \
                  "##heart6_binom_c2={2}\n" \
                  "##binom_pileup={3}\n" \
                  "#fet_c\tgene_set_name\t" \
                  "heart6_p1_1\theart6_m_1\theart6_n_1\theart6_M1_1\theart6_N1_1\theart6_p2_1\t" \
                  "heart6_m_1\theart6_n_1\theart6_M2_1\theart6_N2_1\theart6_m_name_1\theart6_n_name_1\t" \
                  "ctd_p1\tctd_m\tctd_n\tctd_M1\t" \
                  "ctd_N1\tctd_p2\tctd_m\tctd_n\tctd_M2\tctd_N2\tctd_m_name\tctd_n_name\t" \
                  "tof_p1\ttof_m\ttof_n\ttof_M1\ttof_N1\ttof_p2\ttof_m\ttof_n\t" \
                  "tof_M2\ttof_N2\ttof_m_name\ttof_n_name\t" \
                  "heart6_p1_2\theart6_m_2\theart6_n_2\theart6_M1_2\theart6_N1_2\theart6_p2_2\t" \
                  "heart6_m_2\theart6_n_2\theart6_M2_2\theart6_N2_2\theart6_m_name_2\theart6_n_name_2\n" \
                  "".format(ctd_binom_c, tof_binom_c, heart6_binom_c2, pileup_in, heart6_binom_c1)
        fp_out.write(ret_str)
        # print(ret_str)
        logging.debug("main loop begin ...")
        icounter = 0
        while True:
            if icounter % 100000 == 0:
                logging.debug("handled {0:.2%}\t\t{1} / {2} ".format(icounter / total_line, icounter, int(total_line)))
            data_line = fp_in.readline()
            icounter += 1
            # logging.debug("data line=[{}]".format(data_line))
            if not data_line:
                assert len(heart6_value1) > 0
                assert len(ctd_value) > 0
                assert len(heart6_value2) > 0
                assert len(tof_value) > 0
                ret_str = "{0}\t{1}\t{2}\t{3}\t{4}\n" \
                          "".format(pre_key, heart6_value1,
                                    ctd_value, tof_value, heart6_value2)
                # logging.debug("write {0}".format(ret_str))
                fp_out.write(ret_str)
                # print(ret_str)
                break
            if data_line.startswith("#"):
                # logging.debug("continue")
                continue
            data_list = data_line.strip("\n").split("\t")
            key = "{0}\t{1}".format(data_list[1], data_list[0])
            phenotype = data_list[2]
            binom_c = data_list[3]
            values = data_list[4:]
            # logging.debug("phenotype={0} binom_c={1}".format(phenotype, binom_c))
            if key != pre_key:
                if len(pre_key) > 0:
                    if len(ctd_value) == 0:
                        print("key = {}".format(pre_key))
                        print("ctd_binom_c={}".format(ctd_binom_c))
                    if len(heart6_value2) == 0:
                        print("key = {}".format(pre_key))
                        print("heart6_binom_c2={}".format(heart6_binom_c2))
                    if len(tof_value) == 0:
                        print("key = {}".format(pre_key))
                        print("tof_binom_c={}".format(tof_binom_c))
                    if len(heart6_value1) == 0:
                        print('key = {}'.format(pre_key))
                        print('heart6_binom_c1={}'.format(heart6_binom_c1))
                    assert len(heart6_value1) > 0
                    assert len(ctd_value) > 0
                    assert len(heart6_value2) > 0
                    assert len(tof_value) > 0
                    ret_str = "{0}\t{1}\t{2}\t{3}\t{4}\n" \
                              "".format(pre_key, heart6_value1, ctd_value, tof_value, heart6_value2)
                    # logging.debug("write {}".format(ret_str))
                    fp_out.write(ret_str)
                    # print(ret_str)
                heart6_value1 = ''
                ctd_value = ""
                tof_value = ""
                heart6_value2 = ""
                pre_key = key
            if phenotype == "CTD" and binom_c == ctd_binom_c:
                ctd_value = "\t".join(values)
            if phenotype == "tof" and binom_c == tof_binom_c:
                tof_value = "\t".join(values)
            if phenotype == "heart6" and binom_c == heart6_binom_c2:
                heart6_value2 = "\t".join(values)
            if phenotype == 'heart6' and binom_c == heart6_binom_c1:
                heart6_value1 = '\t'.join(values)
    logging.debug("all done")


def geneid2genename(db_file, file_in, cols, output):
    # type: (str, str, str, str) -> None
    col_list = [int(i) for i in cols.split(",")]
    print(col_list)
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()
    with open(file_in, "r") as fp_in, open(output, "w") as fp_out:
        while True:
            data_line = fp_in.readline()
            if not data_line:
                break
            if data_line.startswith("#"):
                fp_out.write(data_line)
                continue
            data_list = data_line.strip("\n").split("\t")
            for col in col_list:
                gene_id_list = data_list[col - 1].split()


def db_update_gene_table(db_file, file_in, key_col):
    logging.basicConfig(filename="db_update_gene_table.log", level=logging.DEBUG, format=log_format, filemode="w")
    logging.debug("Begin")
    logging.debug("db_file={}".format(db_file))
    logging.debug("file_in={}".format(file_in))
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()
    data_list_all = []
    with open(file_in, "r") as fp:
        while True:
            data_line = fp.readline()
            if not data_line:
                break
            if len(data_line.strip()) == 0:
                continue
            data_list = data_line.strip("\n").strip("\r").strip("\n").split("\t")
            data_list_all.append(data_list)
    data_list_all_T = zip(*data_list_all)
    cursor.execute("select {0} from gene_table".format(key_col))
    db_id_list = [i[0] for i in cursor.fetchall()]
    sql_list = []
    logging.debug("collecting sqls")
    print("collecting sqls")
    for one_set_data in data_list_all_T:
        title = one_set_data[0]
        gene_id_set = set(filter(lambda x: len(x) > 0, one_set_data[1:]))
        logging.debug("add col title={}".format(title))
        db_add_col(cursor, "gene_table", title, "int")

        for db_id in db_id_list:
            if db_id in gene_id_set:
                cmd_str = "UPDATE gene_table SET '{0}'=1 WHERE {1}='{2}'".format(title, key_col, db_id)
            else:
                cmd_str = "UPDATE gene_table SET '{0}'=0 WHERE {1}='{2}'".format(title, key_col, db_id)
            sql_list.append(cmd_str)
    logging.debug("execute sqls")
    icounter = 0
    for sql in sql_list:
        logging.debug("sql = {}".format(sql))
        cursor.execute(sql)
        icounter += 1
        if icounter % 100 == 0:
            print("execute sqls\t{0:.2%}\t{1} / {2}".format(float(icounter) / len(sql_list), icounter, len(sql_list)))
    cursor.close()
    conn.commit()
    conn.close()


def modify_mouse_data(db_file):
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()
    cursor.execute("select id, mouse_mean_WT_LV_TPM_3 from gene_table")
    db_data = cursor.fetchall()
    sql_list = []
    for id, value in db_data:
        if value is not None:
            new_value = value / 3
            sql_list.append("UPDATE gene_table SET mouse_mean_WT_LV_TPM_3={0} WHERE id='{1}'".format(new_value, id))
    for sql in sql_list:
        cursor.execute(sql)
    cursor.close()
    conn.commit()
    conn.close()


def db_add_gene_info(db_file, file_in):
    with open(file_in, "r") as fp:
        head_list = fp.readline().strip().split("\t")
        data_dict = {}
        while True:
            data_line = fp.readline()
            if not data_line:
                break
            if len(data_line.strip()) == 0:
                continue
            data_list = data_line.strip("\n").strip("\r").strip("\n").split("\t")
            data_dict[data_list[0]] = data_list

    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()
    for i in head_list[1:]:
        try:
            db_add_col(cursor, "gene_table", i, "varchr(256)")
        except:
            print(i)
            print(traceback.format_exc())
            exit(0)

    cursor.execute("select gene_name from gene_table")
    db_gene_name_list = [i[0] for i in cursor.fetchall()]
    sql_list = []
    for gene_name in db_gene_name_list:
        if gene_name in data_dict:
            target_data_list = data_dict[gene_name]
            cmd_str = "UPDATE gene_table SET "
            for i in range(len(head_list) - 1):
                cmd_str += "'{0}'='{1}', ".format(head_list[i + 1], target_data_list[i + 1])
            cmd_str = cmd_str[:-2] + " where gene_name='{0}'".format(gene_name)
            sql_list.append(cmd_str)
    sql_list_len = len(sql_list)
    icounter = 0
    for sql in sql_list:
        # print(sql)
        try:
            cursor.execute(sql)
        except:
            print(sql)
            print(traceback.format_exc())
            exit(0)
        icounter += 1
        if icounter % 100 == 0:
            print("{0} / {1}\t{2:.2%}".format(icounter, sql_list_len, float(icounter) / sql_list_len))
    cursor.close()
    conn.commit()
    conn.close()


def selected_flat2variance_info(phenotype, selected_flat, fet_path, db_file, output):
    def load_fet(fet_file):
        """
        build gene id 2 variance name dict
        @param fet_file:
        @return: gene id --> [variance_in_gene_name list, variance_in_control_name list, variance_in_case_name list]
        """
        ret_dict = {}
        with open(fet_file, "r") as fp:
            fet_data = [i.strip().split("\t") for i in filter(lambda x: not x.startswith("#"), fp.readlines()) if
                        len(i.strip()) > 0]
        for fet_line in fet_data:
            ret_dict[fet_line[0]] = [filter(lambda x: x != ".", fet_line[18].split(";")),  # variance_in_gene_name
                                     filter(lambda x: x != ".", fet_line[20].split(";")),  # variance_in_control_name
                                     filter(lambda x: x != ".", fet_line[22].split(";"))]  # variance_in_case_name
        return ret_dict

    logging.basicConfig(filename="selected_flat2variance_info.log", level=logging.DEBUG, format=log_format,
                        filemode="w")
    phenotype2gene_name_col_dict = {"CTD": 26, "heart6": 14, "tof": 38, "control": 50}
    assert phenotype in phenotype2gene_name_col_dict
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()
    # get field names
    cursor.execute("PRAGMA table_info([PE_GATK_variant])")
    field_data = [i[1] for i in cursor.fetchall() if
                  (not i[1].startswith("sample") and i[1] != "id" and i[1] != "vcf_id")]
    field_str = ",".join(field_data)
    head_str = "##selected_flat={0}\n##phenotype={1}\n#gene_name\tgene_id\t" \
               "vcf_id\tin_case\tin_control\t".format(selected_flat, phenotype) + \
               "\t".join(field_data) + "\n"
    job_size = file_line_num(selected_flat, "#")
    icounter1 = 0
    with open(selected_flat, "r") as fp_flat, open(output, "w") as fp_out:
        fp_out.write(head_str)
        while True:
            flat_line = fp_flat.readline()
            if not flat_line:
                break
            if len(flat_line.strip()) == 0 or flat_line.startswith("#"):
                continue
            icounter1 += 1

            flat_list = flat_line.strip().split("\t")
            fet_file_name = os.path.join(fet_path, "{0}_{1}.fet".format(phenotype if phenotype != 'control' else 'CTD',
                                                                        flat_list[0]))
            fet_dict = load_fet(fet_file_name)
            gene_list = [i.split("_") for i in flat_list[phenotype2gene_name_col_dict[phenotype] - 1].split(";")]
            for gene_name, gene_id in gene_list:
                # print("handling {}".format(gene_name))
                assert gene_id in fet_dict
                variance_in_control_set = set(fet_dict[gene_id][1])
                variance_in_case_set = set(fet_dict[gene_id][2])
                variance_num = len(fet_dict[gene_id][0])
                icounter = 1
                for variance_id in fet_dict[gene_id][0]:
                    print("handling {0} / {1}\t{2}\t{3} / {4}".format(icounter1, job_size, gene_name, icounter,
                                                                      variance_num))
                    icounter += 1
                    ret_str = "{0}\t{1}".format(gene_name, gene_id)
                    ret_str += "\t{}".format(variance_id)
                    ret_str += "\t1" if variance_id in variance_in_case_set else "\t0"
                    ret_str += "\t1" if variance_id in variance_in_control_set else "\t0"
                    cmd_str = "select {0} from PE_GATK_variant where vcf_id='{1}'".format(field_str, variance_id)
                    logging.debug("cmd_str={}".format(cmd_str))
                    cursor.execute(cmd_str)
                    sql_result = [str(i) for i in cursor.fetchall()[0]]
                    ret_str += "\t" + "\t".join(sql_result)
                    fp_out.write(ret_str + "\n")
    cursor.close()
    conn.commit()
    conn.close()
    print("all done")


def selected_flat2sample_info(phenotype, selected_flat, fet_path, db_file, output):
    def load_fet(fet_file):
        """
        build gene id 2 variance name dict
        @param fet_file:
        @return: gene id --> [variance_in_gene_name list, variance_in_control_name list, variance_in_case_name list]
        """
        ret_dict = {}
        with open(fet_file, "r") as fp:
            fet_data = [i.strip().split("\t") for i in filter(lambda x: not x.startswith("#"), fp.readlines()) if
                        len(i.strip()) > 0]
        for fet_line in fet_data:
            ret_dict[fet_line[0]] = [filter(lambda x: x != ".", fet_line[18].split(";")),  # variance_in_gene_name
                                     filter(lambda x: x != ".", fet_line[20].split(";")),  # variance_in_control_name
                                     filter(lambda x: x != ".", fet_line[22].split(";"))]  # variance_in_case_name
        return ret_dict

    logging.basicConfig(filename="selected_flat2sample_info.log", level=logging.DEBUG, format=log_format, filemode="w")
    phenotype2gene_name_col_dict = {"CTD": 26, "heart6": 14, "tof": 38, "control": 50}
    assert phenotype in phenotype2gene_name_col_dict
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()
    # get field names
    cursor.execute("PRAGMA table_info([sampleChdPhenotype])")
    field_data = [i[1] for i in cursor.fetchall() if i[1] != "id"]
    field_str = ",".join(field_data)
    head_str = "##selected_flat={0}\n##phenotype={1}\n#gene_name\tgene_id\t" \
               "vcf_id\tin_case\tin_control\tgenotype\t".format(selected_flat, phenotype) + \
               "\t".join(field_data) + "\n"
    job_size = file_line_num(selected_flat, "#")
    icounter1 = 0
    print("loading gen_id from db ...")
    cursor.execute("select gen_id from sampleChdPhenotype")
    sample_id_list = [str(i[0]) for i in cursor.fetchall()]
    print("loading sample info from db ...")
    sample_id2sample_info_dict = {}
    cursor.execute("select {} from sampleChdPhenotype".format(field_str))
    sql_result = cursor.fetchall()
    for sample_line in sql_result:
        sample_id2sample_info_dict[str(sample_line[0])] = "\t".join([str(i) for i in sample_line])

    print("loading genotype info from db ...")
    genotype_dict = {}
    cursor.execute(
        "select vcf_id,{} from PE_GATK_variant".format(",".join(["sample_{}".format(i) for i in sample_id_list])))
    sql_result = cursor.fetchall()
    for genotype_line in sql_result:
        genotype_dict[genotype_line[0]] = dict(zip(sample_id_list, [str(i) for i in genotype_line[1:]]))
    # print(sample_id_list[:2])
    with open(selected_flat, "r") as fp_flat, open(output, "w") as fp_out:
        fp_out.write(head_str)
        while True:
            flat_line = fp_flat.readline()
            if not flat_line:
                break
            if len(flat_line.strip()) == 0 or flat_line.startswith("#"):
                continue
            icounter1 += 1

            flat_list = flat_line.strip().split("\t")
            fet_file_name = os.path.join(fet_path, "{0}_{1}.fet".format(phenotype if phenotype != 'control' else 'CTD',
                                                                        flat_list[0]))
            fet_dict = load_fet(fet_file_name)
            gene_list = [i.split("_") for i in flat_list[phenotype2gene_name_col_dict[phenotype] - 1].split(";")]
            for gene_name, gene_id in gene_list:
                # print("handling {}".format(gene_name))
                assert gene_id in fet_dict
                variance_in_control_set = set(fet_dict[gene_id][1])
                variance_in_case_set = set(fet_dict[gene_id][2])
                variance_num = len(fet_dict[gene_id][0])
                icounter = 1
                for variance_id in fet_dict[gene_id][0]:
                    print("handling flat line: {0} / {1}\tcurrent gene: {2}\tcurrent variance: {3} / {4}"
                          "".format(icounter1, job_size, gene_name, icounter, variance_num))
                    icounter += 1
                    ret_str = "{0}\t{1}".format(gene_name, gene_id)
                    ret_str += "\t{}".format(variance_id)

                    ret_str += "\t1" if variance_id in variance_in_case_set else "\t0"
                    ret_str += "\t1" if variance_id in variance_in_control_set else "\t0"
                    for sample_id in sample_id_list:
                        # cmd_str = "select {0} from variance where vcf_id='{1}'" \
                        #           "".format("sample_{}".format(sample_id), variance_id)
                        # cursor.execute(cmd_str)
                        # sql_result = cursor.fetchall()
                        # genotype = str(sql_result[0][0])
                        genotype = genotype_dict[variance_id][sample_id]
                        if genotype.isdigit() and int(genotype) > 0:
                            # cursor.execute("select {0} from sampleChdPhenotype where gen_id={1}"
                            #                "".format(field_str, sample_id))
                            # fp_out.write(ret_str + "\t" + "\t".join([str(i) for i in cursor.fetchall()[0]]) + "\n")
                            fp_out.write(
                                ret_str + "\t" + genotype + "\t" + sample_id2sample_info_dict[sample_id] + "\n")
    cursor.close()
    conn.commit()
    conn.close()
    print("all done")


def build_contingency_table_gene_name_with_gene_list_file(db_file, phenotype, sample_restrict, gene_set_file,
                                                          variance_restrict, output):
    """
    Analyze in units of gene list, use annotation to determine whether variance belongs to gene
    @param db_file:
    @param phenotype:
    @param sample_restrict:
    @param gene_set_file: It needs to be the version with gene name, the format is 
    set name    name1_id1    name2_id2
    @param variance_restrict:
    @param output:
    @return:
    """
    sample_table_name = 'sampleChdPhenotype'

    def sub_gene_name(variance_table, variance_restrict):
        sub_str_list = []
        if variance_table == "variance":
            if re.findall("annovar[ =]", variance_restrict):
                sub_str_list.append("annovar_gene_name")
            if re.findall("bystro[ =]", variance_restrict):
                sub_str_list.append("bystro_gene_name")
            if re.findall("vep[ =]", variance_restrict):
                sub_str_list.append("vep_gene_id")
            if re.findall("spliceAI[ =]", variance_restrict):
                sub_str_list.append("spliceAI_gene_name")
            if re.findall("dmis[ =]", variance_restrict):
                sub_str_list.append("dmis_gene_name")
            if re.findall("dsplicing[ =]", variance_restrict):
                sub_str_list.append("dsplicing_gene_name")
            assert len(sub_str_list) > 0

        else:
            if re.findall("annovar[ =]", variance_restrict):
                sub_str_list.append("annovar_gene_name")
            if re.findall("bystro[ =]", variance_restrict):
                sub_str_list.append("bystro_gene_name")
            if re.findall("vep[ =]", variance_restrict):
                sub_str_list.append("vep_gene_id")
            assert len(sub_str_list) > 0
        ret_str = ", ".join(sub_str_list)
        # print(ret_str)
        return ret_str

    def build_gene_name_list(sql_data):
        ret_list = []
        for i in sql_data:
            ret_list.extend(re.split("[,;]", i))
        return ret_list

    print("db_file=[{0}]\nphenotype=[{1}]\nsample_restrict=[{2}]".format(db_file, phenotype, sample_restrict))
    print("gene_list_file=[{0}]\nvariance_restrict=[{1}]".format(gene_set_file, variance_restrict))
    with open(gene_set_file, "r") as fp:
        gene_data = [i.strip().split('\t') for i in
                     filter(lambda x: len(x.strip()) > 0 and not x.startswith('#'), fp.readlines())]

    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()
    cursor.execute("select id, {0} from variance as v where {1}"
                   "".format(sub_gene_name('variance', variance_restrict), variance_restrict))
    variance_id_data = [[i[0], build_gene_name_list(filter(lambda x: x is not None, i[1:]))] for i in
                        cursor.fetchall()]  # [id, gene_name list(lable)]

    format_sample_restrict = "" if not sample_restrict else " AND {0}".format(sample_restrict)
    cursor.execute("SELECT s.gen_id FROM {2} AS s WHERE s.{0}='0'{1}".format(phenotype,
                                                                             format_sample_restrict,
                                                                             sample_table_name))
    control_id_list = [i[0] for i in cursor.fetchall()]  # m
    print("control number={}".format(len(control_id_list)))
    cursor.execute("SELECT s.gen_id FROM {2} AS s WHERE s.{0}='1'{1}".format(phenotype,
                                                                             format_sample_restrict,
                                                                             sample_table_name))
    case_id_list = [i[0] for i in cursor.fetchall()]  # n
    print("case number={}".format(len(case_id_list)))

    # write output head
    with open(output, "w") as fp:
        fp.write("##phenotype:\"{0}\"\n##control number:{1}\n##case number:{2}\n"
                 "##sample restrict:\"{3}\"\n##version: 3 (gene list based analysis, "
                 "select variance with gene name instead of gene region)\n"
                 "".format(phenotype, len(control_id_list), len(case_id_list), sample_restrict))
        fp.write("""##                         case     control
##                      ---------------------
##                      |         |         |
##    alt allele number |    A    |    B    |
##                      |         |         |
##                      ---------------------
##                      |         |         |
##    ref allele number |    C    |    D    |
##                      |         |         |
##                      ---------------------
##
##
##                              case     control
##                           ---------------------
##                           |         |         |
##          subjects have alt|    A1   |    B1   |
##                           |         |         |
##                           ---------------------
##                           |         |         |
##  subjects do not have alt |    C1   |    D1   |
##                           |         |         |
##                           ---------------------
""")
        fp.write(
            "#gene_set_name\tA\tB\tC\tD\tp_value\todds_ratio\tA1\tB1\tC1\tD1\tp_value1\todds_ratio1\t"
            "n_variance_in_gene_set\tvariance_in_gene_set_name\tn_variance_control\tvariance_in_control_name\t"
            "n_variance_case\tvariance_in_case_name\n")
        for gene_data_line in gene_data:
            gene_set_name = gene_data_line[0]
            variance_id2data_dict = {}
            gene_set = set(re.split('[_\t]', "\t".join(gene_data_line[1:])))  # gene name and gene id
            selected_variance_id = [str(id) for id, lable_list in variance_id_data if set(lable_list) & gene_set]

            n_variance_in_gene_set = len(selected_variance_id)
            cmd_str = "SELECT id, vcf_id, "
            for control_id in control_id_list:
                cmd_str = "{0}sample_{1}, ".format(cmd_str, control_id)
            for case_id in case_id_list:
                cmd_str = "{0}sample_{1}, ".format(cmd_str, case_id)
            cmd_str = cmd_str.strip(", ")
            cmd_str = "{0} FROM variance as v".format(cmd_str)
            cursor.execute(cmd_str + " where v.id in (" + ", ".join(selected_variance_id) + ")")
            sql_result = cursor.fetchall()
            if len(sql_result) == 0:
                continue
            # print("gene_set_name=[{0}] len(sql_result)={1}".format(gene_set_name, len(sql_result)))
            plain_data = zip(*sql_result)
            variance_in_gene_set_name = ",".join(zip(*sql_result)[1])
            control_plain_data = plain_data[2:2 + len(control_id_list)]
            case_plain_data = plain_data[-len(case_id_list):]
            A1 = len([i for i in case_plain_data if sum([int(j) for j in i if type(j) == int]) > 0])
            C1 = len(case_plain_data) - A1
            B1 = len([i for i in control_plain_data if sum([int(j) for j in i if type(j) == int]) > 0])
            D1 = len(control_plain_data) - B1
            oddsratio1, pvalue1 = stats.fisher_exact([[A1, B1], [C1, D1]])

            variance_name_in_control_list = [i[1] for i in zip(*plain_data[:2 + len(control_id_list)]) if
                                             sum([int(j) for j in i[2:] if type(j) == int]) > 0]
            n_variance_in_control = len(variance_name_in_control_list)
            case_data = plain_data[:2]  # id vcf_id
            case_data.extend(case_plain_data)
            variance_name_in_case_list = [i[1] for i in zip(*case_data) if
                                          sum([int(j) for j in i[2:] if type(j) == int]) > 0]
            n_variance_in_case = len(variance_name_in_case_list)

            # calculate ABCD
            for i in sql_result:
                variance_id2data_dict[int(i[0])] = i[2:]
            ret_data = []
            for target_gene in gene_data_line[1:]:
                selected_variance_id_list = [int(id) for id, lable_list in variance_id_data if
                                             set(lable_list) & set(target_gene.split('_'))]
                for selected_variance_id in selected_variance_id_list:
                    ret_data.append(variance_id2data_dict[selected_variance_id])
            control_data = zip(*ret_data)[:len(control_id_list)]
            case_data = zip(*ret_data)[-len(case_id_list):]
            B = sum([sum([int(j) for j in i if type(j) == int]) for i in control_data])
            D = sum([len(i) for i in control_data]) * 2 - B
            A = sum([sum([int(j) for j in i if type(j) == int]) for i in case_data])
            C = sum([len(i) for i in case_data]) * 2 - A
            oddsratio, pvalue = stats.fisher_exact([[A, B], [C, D]])
            fp.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}"
                     "\t{10}\t{11}\t{12}\t{13}\t{14}\t{15}\t{16}\t{17}\t{18}\n"
                     "".format(gene_set_name, A, B, C, D, pvalue, oddsratio, A1, B1, C1, D1, pvalue1, oddsratio1,
                               n_variance_in_gene_set, variance_in_gene_set_name,
                               n_variance_in_control, ",".join(variance_name_in_control_list),
                               n_variance_in_case, ','.join(variance_name_in_case_list)))

    cursor.close()
    conn.commit()
    conn.close()


def analyze_geneset_fisher_test_v3(database, sample_restrict, gene_set_file, path):
    table_name_element_dict = {"": "",
                               " AND (annovar = 1)": "annovar",
                               " AND (bystro = 1)": "bystro",
                               " AND (dmis = 1)": "dmis",
                               " AND (dsplicing = 1)": "dsplicing",
                               " AND (bystro=1 or annovar=1 or vep=1 or dmis=1 or dsplicing=1 or spliceAI=1)": "all.6",
                               " AND (spliceAI = 1)": "spliceAI",
                               " AND (vep = 1)": "vep",
                               " AND (annovar = 1 or bystro = 1 or vep = 1)": "LOF",
                               " AND (dsplicing = 1 or spliceAI = 1)": "dsplicing.all",
                               " AND (annovar = 1 and bystro = 1)": "annovar.bystro",
                               "tof6": "tof",
                               "CTD": "CTD",
                               "heart6": "heart6",
                               "bystro_sampleMaf <= 0.01": "maf.01",
                               "bystro_sampleMaf <= 0.05": "maf.05",
                               "bystro_sampleMaf <= 0.1": "maf.1",
                               "bystro_sampleMaf <= 1": "",
                               # "bystro_sampleMaf <= 0.2": "maf.2",
                               # "bystro_sampleMaf <= 0.3": "maf.3",
                               " AND (bystro_cadd>=10)": "cadd10",
                               " AND (bystro_cadd>=15)": "cadd15",
                               " AND (bystro_cadd>=20)": "cadd20",
                               " AND bystro_phastCons >= 0.4": "Cons.4",
                               " AND bystro_phastCons >= 0.5": "Cons.5",
                               " AND bystro_phastCons >= 0.6": "Cons.6",
                               " AND bystro_phyloP >= -1": "loPn1",
                               " AND bystro_phyloP >= 0": "loP0",
                               " AND bystro_phyloP >= 1": "loP1",
                               " AND (ccrs >= 95)": "ccrs95",
                               " AND (ccrs >= 90)": "ccrs90",
                               " AND (ccrs >= 85)": "ccrs85",
                               " AND (ccrs >= 80)": "ccrs80",
                               " AND (domain_limbr <= {})".format(domain_limbr_quartile_dict[40]): "dlimbr40",
                               " AND (domain_limbr <= {})".format(domain_limbr_quartile_dict[50]): "dlimbr50",
                               " AND (domain_limbr <= {})".format(domain_limbr_quartile_dict[60]): "dlimbr60",
                               " AND (domain_limbr <= {})".format(domain_limbr_quartile_dict[70]): "dlimbr70",
                               " AND (exone_limbr <= {})".format(exone_limbr_quartile_dict[40]): "elimbr40",
                               " AND (exone_limbr <= {})".format(exone_limbr_quartile_dict[50]): "elimbr50",
                               " AND (exone_limbr <= {})".format(exone_limbr_quartile_dict[60]): "elimbr60",
                               " AND (exone_limbr <= {})".format(exone_limbr_quartile_dict[70]): "elimbr70",
                               " AND (is_ccds = 1)": "ccds1"}

    pp = Popen(["mkdir -p {}".format(path)], shell=True)
    pp.wait()
    annotator_list = [" AND (bystro=1 or annovar=1 or vep=1 or dmis=1 or dsplicing=1 or spliceAI=1)",
                      " AND (annovar = 1)",
                      " AND (bystro = 1)",
                      " AND (dmis = 1)",
                      " AND (vep = 1)",
                      " AND (annovar = 1 or bystro = 1 or vep = 1)",
                      " AND (dsplicing = 1 or spliceAI = 1)",
                      " AND (annovar = 1 and bystro = 1)"]
    icounter = 1
    icounter2 = 0
    script_str = ''
    for phenotype in ["tof6", "CTD", "heart6"]:
        for annotator in annotator_list:
            for freq in ["bystro_sampleMaf <= 1"]:
                for cadd in ["", " AND (bystro_cadd>=10)", " AND (bystro_cadd>=15)", " AND (bystro_cadd>=20)"]:
                    for ph in ["",
                               " AND bystro_phastCons >= 0.4",
                               " AND bystro_phastCons >= 0.5",
                               " AND bystro_phastCons >= 0.6",
                               " AND bystro_phyloP >= -1",
                               " AND bystro_phyloP >= 0",
                               " AND bystro_phyloP >= 1"]:
                        for regional_constraint in ["",
                                                    " AND (ccrs >= 95)",
                                                    " AND (ccrs >= 90)",
                                                    " AND (ccrs >= 85)",
                                                    " AND (ccrs >= 80)",
                                                    " AND (domain_limbr <= {})".format(domain_limbr_quartile_dict[40]),
                                                    " AND (domain_limbr <= {})".format(domain_limbr_quartile_dict[50]),
                                                    " AND (domain_limbr <= {})".format(domain_limbr_quartile_dict[60]),
                                                    " AND (domain_limbr <= {})".format(domain_limbr_quartile_dict[70]),
                                                    " AND (exone_limbr <= {})".format(exone_limbr_quartile_dict[40]),
                                                    " AND (exone_limbr <= {})".format(exone_limbr_quartile_dict[50]),
                                                    " AND (exone_limbr <= {})".format(exone_limbr_quartile_dict[60]),
                                                    " AND (exone_limbr <= {})".format(exone_limbr_quartile_dict[70])]:
                            for is_ccds in ["", " AND (is_ccds = 1)"]:
                                variance_restrict = "{0}{1}{2}{3}{4}{5}".format(freq, annotator, cadd, ph,
                                                                                regional_constraint, is_ccds)
                                name_element_list = [table_name_element_dict[phenotype],
                                                     table_name_element_dict[annotator],
                                                     table_name_element_dict[freq],
                                                     table_name_element_dict[cadd],
                                                     table_name_element_dict[ph],
                                                     table_name_element_dict[regional_constraint],
                                                     table_name_element_dict[is_ccds]]
                                output_name = "{0}.geneset.fet" \
                                              "".format("_".join(filter(lambda x: len(x) > 0, name_element_list)))
                                output_name = os.path.join(path, output_name)
                                if script_str == "":
                                    script_str = """#!/bin/bash
#$ -j y
#$ -S /bin/bash
#$ -cwd
#$ -N fg{6}
#$ -pe smp 1
#$ -l h_vmem=8G
#$ -m bes
# -q all.q
export LC_ALL=C
export MALLOC_ARENA_MAX=4
module load python/2.7.15/gcc.4.4.7
module load sqlite3/3.8.11/gcc.4.4.7
time=`date`
echo "==START $time =="
wgsa.py build_contingency_table_gene_name_with_gene_list_file {0} {1} \"{2}\" {3} \"{4}\" {5}
echo {5} is done
time=`date`
echo == $time ==
""".format(database, phenotype, sample_restrict, gene_set_file, variance_restrict, output_name, icounter2)
                                # build_contingency_table_gene_name_with_gene_list_file(db_file, phenotype, sample_restrict, gene_set_file,
                                # variance_restrict, output)
                                else:
                                    script_str += "wgsa.py build_contingency_table_gene_name_with_gene_list_file " \
                                                  "{0} {1} \"{2}\" {3} \"{4}\" {5}" \
                                                  "\necho {5} is done" \
                                                  "\ntime=`date`\necho == $time ==\n" \
                                                  "".format(database, phenotype,
                                                            sample_restrict, gene_set_file,
                                                            variance_restrict, output_name)

                                if icounter == 53:
                                    shell_file_name = os.path.join(path, "qsub{}.sh".format(icounter2))
                                    with open(shell_file_name, "w") as fp:
                                        fp.write(script_str + "\ndate; echo \"==END==\"")
                                    while True:
                                        if os.access(shell_file_name, os.R_OK):
                                            break
                                        time.sleep(1)
                                    print("qsubing {} batch".format(icounter2 + 1))
                                    pp = Popen(["qsub {}".format(shell_file_name)], shell=True)
                                    pp.wait()
                                    script_str = ""
                                    icounter = 0
                                    icounter2 += 1
                                icounter += 1
    if len(script_str) > 0:
        shell_file_name = os.path.join(path, "qsub{}.sh".format(icounter2))
        with open(shell_file_name, "w") as fp:
            fp.write(script_str + "\ndate; echo \"==END==\"")
        while True:
            if os.access(shell_file_name, os.R_OK):
                break
            time.sleep(1)
        print("qsubing {} batch".format(icounter2 + 1))
        pp = Popen(["qsub {}".format(shell_file_name)], shell=True)
        pp.wait()

    print(
        "All the jobs has been submitted. Job number = {}\nIf any job has problem, kill it. And qsub the corresponding sh file".format(
            icounter2 + 1))


def db_add_human_gene_set2gene_table(db_file, human_gene_set_file):
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()
    sql_list = []
    with open(human_gene_set_file, "r") as fp:
        while 1:
            gene_set_line = fp.readline().strip()
            if not gene_set_line:
                break
            gene_set_name = gene_set_line.split("\t")[0]
            gene_set = set([i.split('_')[0] for i in gene_set_line.split("\t")[1:]])
            col_name = "gene_set_" + gene_set_name
            db_add_col(cursor, "gene_table", col_name, "varch(256)")
            cursor.execute("select id, gene_id, gene_name from gene_table")
            for id, gene_id, gene_name in cursor.fetchall():
                if (gene_id in gene_set) or (gene_name in gene_set):
                    sql_list.append("update gene_table set {0}=1 where id={1}".format(col_name, id))
                else:
                    sql_list.append("update gene_table set {0}=0 where id={1}".format(col_name, id))
    print("executing sql...")
    total_len = len(sql_list)
    icounter = 0
    print("executing {0} / {1}".format(icounter, total_len))
    for sql in sql_list:
        cursor.execute(sql)
        icounter += 1
        if icounter % 1000 == 0:
            print("executing {0} / {1}".format(icounter, total_len))
    cursor.close()
    conn.commit()
    conn.close()
    print("all done")


def load_mapping_data(mouse2human):
    mouse2human_dict = {}
    with open(mouse2human, "r") as fp:
        while True:
            data_line = fp.readline()
            if not data_line:
                break
            data_list = data_line.strip().split("\t")
            if len(data_list) == 1:
                continue
            if len(data_list) == 2:
                mouse2human_dict[data_list[0]] = [data_list[1]]
                continue
            if len(data_list) == 3:
                mouse_name = data_list[0]
                official_human_name = data_list[1]
                aka_list = data_list[2].split(";")
                mouse2human_dict[mouse_name] = [official_human_name]
                mouse2human_dict[mouse_name].extend(aka_list)
    return mouse2human_dict


def db_add_mouse_gene_set2gene_table(db_file, mouse_set_file, mouse2human):
    mouse2human_dict = load_mapping_data(mouse2human)
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()
    sql_list = []
    with open(mouse_set_file, "r") as fp:
        while 1:
            gene_set_line = fp.readline().strip()
            if not gene_set_line:
                break
            gene_set_name = gene_set_line.split("\t")[0]
            # gene_set = set(gene_set_line.split("\t")[1:])
            gene_set = []
            for mouse_name in gene_set_line.split("\t")[1:]:
                if mouse_name in mouse2human_dict:
                    gene_set.extend(mouse2human_dict[mouse_name])
            gene_set = set(gene_set)
            col_name = "gene_set_" + gene_set_name
            db_add_col(cursor, "gene_table", col_name, "varch(256)")
            cursor.execute("select id, gene_id, gene_name from gene_table")
            for id, human_gene_id, human_gene_name in cursor.fetchall():
                if (human_gene_id in gene_set) or (human_gene_name in gene_set):
                    sql_list.append("update gene_table set {0}=1 where id={1}".format(col_name, id))
                else:
                    sql_list.append("update gene_table set {0}=0 where id={1}".format(col_name, id))
    print("executing sql...")
    total_len = len(sql_list)
    icounter = 0
    print("executing {0} / {1}".format(icounter, total_len))
    for sql in sql_list:
        cursor.execute(sql)
        icounter += 1
        if icounter % 1000 == 0:
            print("executing {0} / {1}".format(icounter, total_len))
    cursor.close()
    conn.commit()
    conn.close()
    print("all done")


def transform_mouse2human_gene_set(mouse_gene_set_file, mapping_file, output):
    """
    把老鼠gene set转换成人类的gene set。因为有曾用名，因此转换完了每一个gene set会多出很多gene。
    @param mouse_gene_set_file:
    @param mapping_file:
    @return:
    """
    mouse2human_dict = load_mapping_data(mapping_file)
    with open(mouse_gene_set_file, 'r') as fp:
        mouse_gene_set = [i.strip().split('\t') for i in fp.readlines()]
    with open(output, "w") as fp_out:
        for mouse_gene_line in mouse_gene_set:
            gene_set_name = mouse_gene_line[0]
            human_gene_set_list = []
            for mouse_gene in mouse_gene_line[1:]:
                if mouse_gene in mouse2human_dict:
                    human_gene_set_list.extend(mouse2human_dict[mouse_gene])
            fp_out.write("{0}\t{1}\n".format(gene_set_name, "\t".join(human_gene_set_list)))


def gene_set2gene_set_with_id(db_file, gene_set_file):
    """
    把只有gene name的gene set，转换成 geneName_geneId的gene set
    @param db_file:
    @param gene_set_file:
    @return:
    """
    name2id_dict = {}

    def name2key(cursor, name, name2id_dict):
        # print('name={}'.format(name))
        if name in name2id_dict:
            id = name2id_dict[name]
            ret = "{0}_{1}".format(name, id)
            # print(ret)
            return ret
        cmd_str = """select gene_id from gene_table where gene_name=\"{0}\"""".format(name)
        # print(cmd_str)
        cursor.execute(cmd_str)
        sql_ret = cursor.fetchall()
        if len(sql_ret) == 0:
            ret = "{0}_None".format(name)
            # print(ret)
            name2id_dict[name] = "None"
            return ret
        gene_id = sql_ret[0][0]
        name2id_dict[name] = gene_id
        ret = "{0}_{1}".format(name, gene_id)
        # print(ret)
        sys.stdout.flush()
        return ret

    print("db_file={}".format(db_file))
    print("gene set file={}".format(gene_set_file))
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()
    with open(gene_set_file, "r") as fp, open(gene_set_file + ".withid", "w") as fp_out:
        while True:
            gene_set_line = fp.readline()
            if not gene_set_line:
                break
            if len(gene_set_line.strip()) == 0:
                continue

            gene_set_list = gene_set_line.strip().split("\t")
            print("handling {}".format(gene_set_list[0]))
            name_id_list = [name2key(cursor, i, name2id_dict) for i in filter(lambda x: len(x) > 0, gene_set_list[1:])]
            fp_out.write("{0}\t{1}\n".format(gene_set_list[0], "\t".join(name_id_list)))
    cursor.close()
    conn.commit()
    conn.close()


def fet_pileup(selected_list_file, fet_path, threshold, mode):
    logging.basicConfig(filename="fet_pileup.log", level=logging.DEBUG, format=log_format, filemode="w")
    assert mode in ['p', 'p1']
    threshold = float(threshold)
    with open(selected_list_file, "r") as fp:
        selected_list = [i.strip("\n").split("\t") for i in fp.readlines() if
                         len(i.strip()) > 0 and not i.startswith("#")]
    if mode == 'p':
        selected_file_name_list = [os.path.join(fet_path, i[0]) for i in filter(lambda x: x[1] == 'p', selected_list)]
    else:
        selected_file_name_list = [os.path.join(fet_path, i[0]) for i in filter(lambda x: x[2] == 'p1', selected_list)]
    icounter = 1
    for selected_file_name in selected_file_name_list:
        with open(selected_file_name, "r") as fp:
            logging.debug("handling {0} ... {1} / {2}".format(selected_file_name,
                                                              icounter,
                                                              len(selected_file_name_list)))
            icounter += 1
            index_p = 0
            index_p1 = 0
            while True:
                data_line = fp.readline()
                if not data_line:
                    break
                if data_line.startswith("##"):
                    continue
                if data_line.startswith("#"):
                    data_list = data_line.strip().split("\t")
                    index_p = data_list.index("p_value")
                    index_p1 = data_list.index("p_value1")
                    continue
                data_list = data_line.strip().split("\t")

                target_p = float(data_list[index_p]) if mode == 'p' else float(data_list[index_p1])
                if target_p <= threshold:
                    print("{0}\t{1}".format(selected_file_name, data_line))
    logging.debug("all done")


def variance_restrict2variance_info_sample_info(db_file, variance_restrict, variance_info_out, sample_info_out):
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()
    cursor.execute('select gen_id from sampleChdPhenotype')
    gen_id_list = [i[0] for i in cursor.fetchall()]
    sample_list = ["sample_{}".format(i) for i in gen_id_list]
    # variance info
    cursor.execute('select * from variance where {}'.format(variance_restrict))
    with open(variance_info_out, "w") as fp:
        for i in cursor.fetchall():
            fp.write("\t".join([str(j) for j in i]) + '\n')
    cursor.execute('select {0} from variance where {1}'.format(",".join(sample_list), variance_restrict))
    sql_ret = cursor.fetchall()

    gen_id_selector = [sum([int(j) for j in i if j != 'na']) > 0 for i in zip(*sql_ret)]
    sql = 'select * from sampleChdPhenotype where gen_id in ({})'.format(
        ",".join([str(i) for i in list(compress(gen_id_list, gen_id_selector))]))
    cursor.execute(sql)
    with open(sample_info_out, "w") as fp:
        for i in cursor.fetchall():
            fp.write("\t".join([str(j) for j in i]) + '\n')
    conn.commit()
    conn.close()


def variance_restrict2variance_info(db_file, variance_restrict, variance_info_out, variance_table='variance'):
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()
    cursor.execute('select gene_id,gene_name from gene_table')
    sql_ret = cursor.fetchall()
    gene_id2name_dict = dict(sql_ret)
    cursor.execute('select gene_name, gene_id from gene_table')
    sql_ret = cursor.fetchall()
    gene_name2id_dict = dict(sql_ret)

    cursor.execute("PRAGMA table_info([{}])".format(variance_table))
    field_data = [i[1] for i in cursor.fetchall() if
                  (not i[1].startswith("sample") and i[1] != "id" and i[1] != "vcf_id")]
    field_str = ",".join(field_data)

    # variance info
    cursor.execute('select id, vcf_id, bystro_gene_name, annovar_gene_name, '
                   'vep_gene_id, spliceAI_gene_name, dmis_gene_name, '
                   'dsplicing_gene_name from {0} where {1}'.format(variance_table, variance_restrict))
    sql_ret = cursor.fetchall()

    with open(variance_info_out, "w") as fp:
        fp.write("#method=wgsa.py variance_restrict2variance_info\n#database={0}\n#variance_table={2}\n"
                 "#variance_restrict={1}\n#gene_name\tgene_id\tvcf_id\t{3}\n"
                 "".format(db_file, variance_restrict, variance_table, "\t".join(field_data)))
        for sql_line in sql_ret:
            table_id = sql_line[0]
            vcf_id = sql_line[1]
            key_set = set([])
            for name_id_str in sql_line[2:]:
                name_id_list = name_id_str.split(";") if name_id_str is not None else []
                for name_id in name_id_list:
                    if name_id in gene_name2id_dict:
                        key_set.add("{0}\t{1}".format(name_id, gene_name2id_dict[name_id]))
                    elif name_id in gene_id2name_dict:
                        key_set.add("{0}\t{1}".format(gene_id2name_dict[name_id], name_id))
                    else:
                        print("can't convert {}".format(name_id))
            cursor.execute('select {0} from {1} where id = {2}'.format(field_str, variance_table, table_id))
            variance_info = "\t".join([str(i) for i in cursor.fetchone()])
            for key in key_set:
                fp.write("{0}\t{1}\t{2}\n".format(key, vcf_id, variance_info))

    conn.commit()
    conn.close()


def variance_restrict2sample_info(db_file, variance_restrict, sample_info_out, variance_table):
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()
    cursor.execute('select gene_id,gene_name from gene_table')
    sql_ret = cursor.fetchall()
    gene_id2name_dict = dict(sql_ret)
    cursor.execute('select gene_name, gene_id from gene_table')
    sql_ret = cursor.fetchall()
    gene_name2id_dict = dict(sql_ret)

    cursor.execute('select gen_id from sampleChdPhenotype')
    gen_id_list = [i[0] for i in cursor.fetchall()]
    sample_list = ["sample_{}".format(i) for i in gen_id_list]

    cursor.execute("PRAGMA table_info([sampleChdPhenotype])")
    field_data = [i[1] for i in cursor.fetchall() if i[1] != "id"]
    field_str = ",".join(field_data)

    # variance info
    cursor.execute('select id, vcf_id, bystro_gene_name, annovar_gene_name, '
                   'vep_gene_id, spliceAI_gene_name, dmis_gene_name, '
                   'dsplicing_gene_name from {0} where {1}'.format(variance_table, variance_restrict))
    sql_ret = cursor.fetchall()

    with open(sample_info_out, "w") as fp:
        fp.write("#method=wgsa.py variance_restrict2sample_info\n"
                 "#database={0}\n#variance_table={1}\n"
                 "#variance_restrict={2}\n#gene_name\tgene_id\tvcf_id\t{3}\n"
                 "".format(db_file, variance_table, variance_restrict, "\t".join(field_data)))
        for sql_line in sql_ret:
            table_id = sql_line[0]
            vcf_id = sql_line[1]
            key_set = set([])
            for name_id_str in sql_line[2:]:
                name_id_list = name_id_str.split(";") if name_id_str is not None else []
                for name_id in name_id_list:
                    if name_id in gene_name2id_dict:
                        key_set.add("{0}\t{1}".format(name_id, gene_name2id_dict[name_id]))
                    elif name_id in gene_id2name_dict:
                        key_set.add("{0}\t{1}".format(gene_id2name_dict[name_id], name_id))
                    else:
                        print("can't convert {}".format(name_id))
            cursor.execute('select {0} from {1} where id = {2}'.format(','.join(sample_list),
                                                                       variance_table,
                                                                       table_id))
            genotype_list = cursor.fetchone()
            for key in key_set:
                for index in range(len(genotype_list)):
                    if genotype_list[index] != 'na' and int(genotype_list[index]) > 0:
                        cursor.execute('select {0} from sampleChdPhenotype where gen_id={1}'
                                       ''.format(field_str, gen_id_list[index]))
                        sample_info_list = [str(i) for i in cursor.fetchone()]
                        fp.write("{0}\t{1}\t{2}\n".format(key, vcf_id, '\t'.join(sample_info_list)))

    conn.commit()
    conn.close()


def parse_fet_name(fet_name):
    """
    parse the file name of fet into phenotype and variance restrict
    @param fet_name:
    @return:
    """
    # print("input={}".format(fet_name))
    assert fet_name.endswith('.fet')
    fet_name = fet_name[:-4]
    mouse_cut = None
    mouse_list = re.findall('mouse(\d+)$', fet_name)
    assert len(mouse_list) < 2
    if len(mouse_list) == 1:
        mouse_cut = int(mouse_list[0])
    # print('aa={}'.format(fet_name))
    fet_name = re.sub('\.mouse\d+', "", fet_name)
    # print('aaa={}'.format(fet_name))
    table_name_element_dict = {"": "",
                               " AND (annovar = 1)": "annovar",
                               " AND (bystro = 1)": "bystro",
                               " AND (dmis = 1)": "dmis",
                               " AND (dsplicing = 1)": "dsplicing",
                               " AND (bystro=1 or annovar=1 or vep=1 or dmis=1 or dsplicing=1 or spliceAI=1)": "all.6",
                               " AND (spliceAI = 1)": "spliceAI",
                               " AND (vep = 1)": "vep",
                               " AND (annovar = 1 or bystro = 1 or vep = 1)": "LOF",
                               " AND (dsplicing = 1 or spliceAI = 1)": "dsplicing.all",
                               " AND (annovar = 1 and bystro = 1)": "annovar.bystro",
                               "tof6": "tof",
                               "CTD": "CTD",
                               "heart6": "heart6",
                               "bystro_sampleMaf <= 0.01": "maf.01",
                               "bystro_sampleMaf <= 0.05": "maf.05",
                               "bystro_sampleMaf <= 0.1": "maf.1",
                               "bystro_sampleMaf <= 1": "",
                               # "bystro_sampleMaf <= 0.2": "maf.2",
                               # "bystro_sampleMaf <= 0.3": "maf.3",
                               " AND (bystro_cadd>=10)": "cadd10",
                               " AND (bystro_cadd>=15)": "cadd15",
                               " AND (bystro_cadd>=20)": "cadd20",
                               " AND bystro_phastCons >= 0.4": "Cons.4",
                               " AND bystro_phastCons >= 0.5": "Cons.5",
                               " AND bystro_phastCons >= 0.6": "Cons.6",
                               " AND bystro_phyloP >= -1": "loPn1",
                               " AND bystro_phyloP >= 0": "loP0",
                               " AND bystro_phyloP >= 1": "loP1",
                               " AND (ccrs >= 95)": "ccrs95",
                               " AND (ccrs >= 90)": "ccrs90",
                               " AND (ccrs >= 85)": "ccrs85",
                               " AND (ccrs >= 80)": "ccrs80",
                               " AND (domain_limbr <= {})".format(domain_limbr_quartile_dict[40]): "dlimbr40",
                               " AND (domain_limbr <= {})".format(domain_limbr_quartile_dict[50]): "dlimbr50",
                               " AND (domain_limbr <= {})".format(domain_limbr_quartile_dict[60]): "dlimbr60",
                               " AND (domain_limbr <= {})".format(domain_limbr_quartile_dict[70]): "dlimbr70",
                               " AND (exone_limbr <= {})".format(exone_limbr_quartile_dict[40]): "elimbr40",
                               " AND (exone_limbr <= {})".format(exone_limbr_quartile_dict[50]): "elimbr50",
                               " AND (exone_limbr <= {})".format(exone_limbr_quartile_dict[60]): "elimbr60",
                               " AND (exone_limbr <= {})".format(exone_limbr_quartile_dict[70]): "elimbr70",
                               " AND (is_ccds = 1)": "ccds1"}
    name_ele2condition_dict = {}
    for con in table_name_element_dict:
        ele = table_name_element_dict[con]
        name_ele2condition_dict[ele] = con
    contex_str = os.path.basename(fet_name)
    contex_list = contex_str.split("_")
    # print("contex_list=[{}]".format(contex_list))
    assert len(contex_list) > 1
    phenotype = contex_list[0]
    variance_restrict = "bystro_sampleMaf <= 1"
    for i in range(1, len(contex_list)):
        variance_restrict += name_ele2condition_dict[contex_list[i]]
    return [phenotype, variance_restrict, mouse_cut]


def gene_binom_selected(A1, B1, p_value1, or1, binomc):
    assert type(A1) == int
    assert type(B1) == int
    if binomc == "A3.2.B1.0":
        return (B1 <= 0 and A1 >= 2) or (B1 <= 1 and A1 >= 3)
    if binomc == "A4.3.B1.0":
        return (B1 <= 0 and A1 >= 3) or (B1 <= 1 and A1 >= 4)
    if binomc == "A3.B0":
        return B1 <= 0 and A1 >= 3
    if binomc == "A4.B0":
        return B1 <= 0 and A1 >= 4
    if binomc == "A3.B1":
        return B1 <= 1 and A1 >= 3
    if binomc == "A4.B1":
        return B1 <= 1 and A1 >= 4
    if binomc == "A2.B0":
        return B1 <= 0 and A1 >= 2
    if binomc == "A1.B0":
        return B1 <= 0 and A1 >= 1
    if binomc == "p1_case":
        return p_value1 <= 0.05 and or1 > 1
    if binomc == "A1.0.B3.2":
        return (A1 <= 0 and B1 >= 2) or (A1 <= 1 and B1 >= 3)
    if binomc == "A1.0.B4.3":
        return (A1 <= 0 and B1 >= 3) or (A1 <= 1 and B1 >= 4)
    if binomc == "A0.B3":
        return A1 <= 0 and B1 >= 3
    if binomc == "A1.B4":
        return A1 <= 1 and B1 >= 4
    if binomc == "A1.B3":
        return A1 <= 1 and B1 >= 3
    if binomc == "A0.B2":
        return A1 <= 0 and B1 >= 2
    if binomc == "A0.B1":
        return A1 <= 0 and B1 >= 1
    if binomc == "A0.B4":
        return A1 <= 0 and B1 >= 4
    if binomc == "p1_control":
        return p_value1 <= 0.05 and or1 < 1
    if binomc == 'none':
        return True
    else:
        RuntimeError("wrong binom condition {}".format(binomc))


def fet_binomc2sample_variance_info(db_file, fet_file, binom_c, variance_table, sample_info_out, variance_info_out):
    phenotype, variance_restrict, mouse_cut = parse_fet_name(fet_file)
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()
    cursor.execute('select gen_id from sampleChdPhenotype')
    gen_id_list = [i[0] for i in cursor.fetchall()]
    sample_list = ["sample_{}".format(i) for i in gen_id_list]

    cursor.execute("PRAGMA table_info([{}])".format(variance_table))
    variance_field_data = [i[1] for i in cursor.fetchall() if
                           (not i[1].startswith("sample") and i[1] != "id" and i[1] != "vcf_id")]
    variance_field_str = ",".join(variance_field_data)

    cursor.execute("PRAGMA table_info([sampleChdPhenotype])")
    sample_field_data = [i[1] for i in cursor.fetchall() if i[1] != "id"]
    sample_field_str = ",".join(sample_field_data)

    with open(sample_info_out, "w") as fp_sample_out, open(fet_file, 'r') as fp_in, open(variance_info_out,
                                                                                         'w') as fp_variance_out:
        fp_sample_out.write("#method=wgsa.py fet_binomc2sample_variance_info\n"
                            "#database={0}\n#variance_table={1}\n"
                            "#variance_restrict={2}\n#gene_name\tgene_id\tvcf_id\t{3}\n"
                            "".format(db_file, variance_table, variance_restrict, "\t".join(sample_field_data)))
        fp_variance_out.write("#method=wgsa.py fet_binomc2sample_variance_info\n"
                              "#database={0}\n#variance_table={1}\n"
                              "#variance_restrict={2}\n#gene_name\tgene_id\tvcf_id\t{3}\n"
                              "".format(db_file, variance_table, variance_restrict, "\t".join(variance_field_data)))
        counter = 1
        while True:
            fet_line = fp_in.readline()
            if not fet_line:
                break
            if fet_line.startswith("#"):
                continue
            fet_list = fet_line.strip().split('\t')
            gene_id = fet_list[0]
            gene_name = fet_list[1]
            A1 = int(fet_list[11])
            B1 = int(fet_list[12])
            p_value1 = float(fet_list[15])
            or1 = float(fet_list[16])
            if not gene_binom_selected(A1, B1, p_value1, or1, binom_c):
                continue
            print("gene {0} selected {1}".format(gene_name, counter))
            counter += 1
            vcf_id_list = fet_list[18].split(';')
            for vcf_id in vcf_id_list:
                cursor.execute(
                    'select {0} from {1} where vcf_id = "{2}"'.format(variance_field_str, variance_table, vcf_id))

                variance_info_data = [str(i) for i in cursor.fetchone()]
                fp_variance_out.write("{0}\t{1}\t{2}\t{3}\n".format(gene_name,
                                                                    gene_id,
                                                                    vcf_id,
                                                                    '\t'.join(variance_info_data)))
                cursor.execute(
                    'select {0} from {1} where vcf_id = "{2}"'.format(','.join(sample_list), variance_table, vcf_id))

                genotype_list = cursor.fetchone()
                for index in range(len(genotype_list)):
                    if genotype_list[index] != 'na' and int(genotype_list[index]) > 0:
                        cursor.execute('select {0} from sampleChdPhenotype where gen_id={1}'
                                       ''.format(sample_field_str, gen_id_list[index]))
                        sample_info_list = [str(i) for i in cursor.fetchone()]
                        fp_sample_out.write(
                            "{0}\t{1}\t{2}\t{3}\n".format(gene_name, gene_id, vcf_id, '\t'.join(sample_info_list)))
    conn.commit()
    conn.close()


def fet_genelist2sample_variance_info(db_file, fet_file, gene_list, variance_table, sample_info_out, variance_info_out):
    def gene_selected(gene_id, gene_name, gene_list):
        if gene_id in gene_list:
            return True
        if gene_name in gene_list:
            return True
        return False

    print('db={}'.format(db_file))
    print('fet={}'.format(fet_file))
    print('gene_list={}'.format(gene_list))
    print('variance_table={}'.format(variance_table))

    gene_list = re.split('[;,]', gene_list)
    phenotype, variance_restrict, mouse_cut = parse_fet_name(fet_file)
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()
    cursor.execute('select gen_id from sampleChdPhenotype')
    gen_id_list = [i[0] for i in cursor.fetchall()]
    sample_list = ["sample_{}".format(i) for i in gen_id_list]

    cursor.execute("PRAGMA table_info([{}])".format(variance_table))
    variance_field_data = [i[1] for i in cursor.fetchall() if
                           (not i[1].startswith("sample") and i[1] != "id" and i[1] != "vcf_id")]
    variance_field_str = ",".join(variance_field_data)

    cursor.execute("PRAGMA table_info([sampleChdPhenotype])")
    sample_field_data = [i[1] for i in cursor.fetchall() if i[1] != "id"]
    sample_field_str = ",".join(sample_field_data)

    with open(sample_info_out, "w") as fp_sample_out, open(fet_file, 'r') as fp_in, open(variance_info_out,
                                                                                         'w') as fp_variance_out:
        fp_sample_out.write("##method=wgsa.py fet_genelist2sample_variance_info\n"
                            "##database={0}\n##variance_table={1}\n"
                            "##variance_restrict={2}\n##gene_list={3}\n#gene_name\tgene_id\tvcf_id\t{4}\n"
                            "".format(db_file, variance_table, variance_restrict, ','.join(gene_list),
                                      "\t".join(sample_field_data)))
        fp_variance_out.write("##method=wgsa.py fet_genelist2sample_variance_info\n"
                              "##database={0}\n##variance_table={1}\n"
                              "##variance_restrict={2}\n##gene_list={3}\n#gene_name\tgene_id\tvcf_id\t{4}\n"
                              "".format(db_file, variance_table, variance_restrict, ','.join(gene_list),
                                        "\t".join(variance_field_data)))
        counter = 1
        while True:
            fet_line = fp_in.readline()
            if not fet_line:
                break
            if fet_line.startswith("#"):
                continue
            fet_list = fet_line.strip().split('\t')
            gene_id = fet_list[0]
            gene_name = fet_list[1]

            if not gene_selected(gene_id, gene_name, gene_list):
                continue
            print("gene {0} selected {1}".format(gene_name, counter))
            counter += 1
            vcf_id_list = fet_list[18].split(';')
            for vcf_id in vcf_id_list:
                cursor.execute(
                    'select {0} from {1} where vcf_id = "{2}"'.format(variance_field_str, variance_table, vcf_id))
                variance_info_data = [str(i) for i in cursor.fetchone()]
                fp_variance_out.write("{0}\t{1}\t{2}\t{3}\n".format(gene_name,
                                                                    gene_id,
                                                                    vcf_id,
                                                                    '\t'.join(variance_info_data)))
                cursor.execute(
                    'select {0} from {1} where vcf_id = "{2}"'.format(','.join(sample_list), variance_table, vcf_id))
                genotype_list = cursor.fetchone()
                for index in range(len(genotype_list)):
                    if genotype_list[index] != 'na' and int(genotype_list[index]) > 0:
                        cursor.execute('select {0} from sampleChdPhenotype where gen_id={1}'
                                       ''.format(sample_field_str, gen_id_list[index]))
                        sample_info_list = [str(i) for i in cursor.fetchone()]
                        fp_sample_out.write(
                            "{0}\t{1}\t{2}\t{3}\n".format(gene_name, gene_id, vcf_id, '\t'.join(sample_info_list)))
    conn.commit()
    conn.close()


def grep_human_exon_count_from_ncbi(gene_symbol):
    gene_name_head = 'Name/Gene ID'
    gene_desciption_head = 'Description'
    gene_aliases_head = 'Aliases'
    print("handling [{}]".format(gene_symbol))
    url = "https://www.ncbi.nlm.nih.gov/gene/?term={0}".format(gene_symbol)
    driver = webdriver.Chrome("/usr/local/bin/chromedriver")
    driver.get(url)
    time.sleep(1)
    soup = BeautifulSoup(driver.page_source, "html.parser")

    page_no_list = soup.find_all(id="pageno")
    if len(page_no_list) > 0:
        last_page_num = int(page_no_list[0]['last'])
    else:
        last_page_num = 1
    human_uri = ''

    # looking for the mouse_uri in the search list
    for page_index in range(last_page_num):
        soup = BeautifulSoup(driver.page_source, "html.parser")
        tables = soup.find_all('table', id="ui-ncbigrid-7")
        if len(tables) != 1:
            break
        table = tables[0]
        head_list = [i.getText() for i in table.findAll('th')]
        gene_name_col_idx = head_list.index(gene_name_head)
        gene_description_col_idx = head_list.index(gene_desciption_head)
        gene_aliases_col_idx = head_list.index(gene_aliases_head)
        tr_list = table.findAll('tr')
        for tr in tr_list:
            # print(tr)
            td_list = tr.findAll('td')
            if len(td_list) == 0:
                continue
            gene_name_a = td_list[gene_name_col_idx].find('a')
            gene_name_str = gene_name_a.getText()
            href_str = gene_name_a.get('href')
            description_str = td_list[gene_description_col_idx].getText()
            aliases_str = td_list[gene_aliases_col_idx].getText()
            if ('Homo sapiens' not in description_str) or (
                    gene_symbol not in aliases_str and gene_symbol != str(gene_name_str)):
                continue
            human_uri = str(href_str)
            break

        if human_uri:
            break
        next_link = soup.find('a', class_="active page_link next")
        # print("next_link={}".format(next_link))
        if next_link:
            driver.find_element_by_partial_link_text("Next").click()
        time.sleep(1)
    # if did not find the mouse uri in the list, parse the summary
    if not human_uri:
        info_icon_list = soup.find('ul', id='msgportlet').find_all('li', class_='info icon')
        if info_icon_list:
            message = info_icon_list[0].span.getText()
            print('message=[{}]'.format(message))
            if message == 'No items found.':
                return 'na'
            else:
                print('unknow message')
        else:
            print('no message')

    else:
        print("got the human uri ({})in the list".format(human_uri))
        url = "https://www.ncbi.nlm.nih.gov{}".format(human_uri)
        driver.get(url)
        soup = BeautifulSoup(driver.page_source, "html.parser")
        exon_dl_list = soup.find_all('dl', class_='dl-chr-info exon-count')
        if exon_dl_list:
            return exon_dl_list[0].dd.getText()
        else:
            print("{} no exon dl".format(gene_symbol))
            exit(0)


def add_exon_count2file(raw_file, gene_symbol_col, excon_count_col, output):
    gene_symbol2exon_count_dict = {}
    with open(raw_file, "r") as fp, open(output, "w") as fp_out:
        fp_out.write(fp.readline())
        while True:
            data = fp.readline()
            if not data:
                break
            if len(data.strip()) == 0:
                continue
            data_list = data.strip("\n").split("\t")
            gene_symbol_list = [i.strip() for i in data_list[int(gene_symbol_col) - 1].strip().split(",") if
                                len(i.strip()) > 0]
            exon_count_list = []
            for symbol in gene_symbol_list:
                if symbol in gene_symbol2exon_count_dict:
                    count = gene_symbol2exon_count_dict[symbol]
                else:
                    count = grep_human_exon_count_from_ncbi(symbol)
                    gene_symbol2exon_count_dict[symbol] = count
                exon_count_list.append(count)
            data_list[int(excon_count_col) - 1] = ','.join(exon_count_list)
            fp_out.write('\t'.join(data_list) + "\n")


def add_gene_size2file(raw_file, gene_symbol_col, gene_size_col, db_file, output):
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()
    with open(raw_file, 'r') as fp, open(output, 'w') as fp_out:
        fp_out.write(fp.readline())
        while True:
            data = fp.readline()
            if not data:
                break
            data_list = data.strip("\n").split("\t")
            symbol_list = [i.strip() for i in data_list[int(gene_symbol_col) - 1].strip().split(",") if
                           len(i.strip()) > 0]
            size_list = []
            for symbol in symbol_list:
                sql_str = 'select start_pos, end_pos, start_pos2, end_pos2 from gene_table where gene_name=\'{}\''.format(
                    symbol)
                cursor.execute(sql_str)
                sql_ret = cursor.fetchall()
                if sql_ret:
                    start, end, start2, end2 = sql_ret[0]
                    assert start2 is None
                    assert end2 is None
                    size_list.append(str(abs(int(start) - int(end))))
                else:
                    size_list.append('na')
                    print('not found {}'.format(symbol))
            data_list[int(gene_size_col) - 1] = ','.join(size_list)
            fp_out.write("\t".join(data_list) + '\n')
    conn.commit()
    conn.close()


def build_contingency_table_gene_name_cnv(variance_table, db_file, phenotype,
                                          sample_table_name, sample_restrict,
                                          gene_table_name, gene_restrict,
                                          variance_restrict, synonymous_restrict,
                                          output, job_id, hg38_cnv_file, should_log=True):
    """
           case     control
        ---------------------
        |         |         |
    alt |    A    |    B    |
        |         |         |
        ---------------------
        |         |         |
    ref |    C    |    D    |
        |         |         |
        ---------------------
    Carry out fisher exact test in units of genes
    hg38_cnv_file: chrom start end sample_id heart6 CTD
    :return:
    """

    def sub_gene_name(variance_table, variance_restrict):
        sub_str_list = []
        if variance_table == "variance":
            if re.findall("annovar[ =]", variance_restrict):
                sub_str_list.append("annovar_gene_name")
            if re.findall("bystro[ =]", variance_restrict):
                sub_str_list.append("bystro_gene_name")
            if re.findall("vep[ =]", variance_restrict):
                sub_str_list.append("vep_gene_id")
            if re.findall("spliceAI[ =]", variance_restrict):
                sub_str_list.append("spliceAI_gene_name")
            if re.findall("dmis[ =]", variance_restrict):
                sub_str_list.append("dmis_gene_name")
            if re.findall("dsplicing[ =]", variance_restrict):
                sub_str_list.append("dsplicing_gene_name")
            assert len(sub_str_list) > 0

        else:
            if re.findall("annovar[ =]", variance_restrict):
                sub_str_list.append("annovar_gene_name")
            if re.findall("bystro[ =]", variance_restrict):
                sub_str_list.append("bystro_gene_name")
            if re.findall("vep[ =]", variance_restrict):
                sub_str_list.append("vep_gene_id")
            assert len(sub_str_list) > 0
        ret_str = ", ".join(sub_str_list)
        # print(ret_str)
        return ret_str

    def build_gene_name_list(sql_data):
        ret_list = []
        for i in sql_data:
            ret_list.extend(re.split("[,;]", i))
        return ret_list

    def is_cnv_gene_overlap(gene_line, cnv_line):
        def region_overlap(chr1, start1, end1, chr2, start2, end2):
            if chr1 != chr2:
                return False
            if end1 < start2 or start1 > end2:
                return False
            return True

        chrom, start_pos, end_pos, chrom2, start_pos2, end_pos2 = gene_line[2:]
        cnv_chrom, cnv_start, cnv_end = cnv_line[:3]
        start_pos = int(start_pos)
        end_pos = int(end_pos)
        if chrom2 is not None:
            start_pos2 = int(start_pos2)
            end_pos2 = int(end_pos2)
        cnv_chrom = cnv_chrom.strip('chr')
        cnv_start = int(cnv_start)
        cnv_end = int(cnv_end)
        if region_overlap(chrom, start_pos, end_pos, cnv_chrom, cnv_start, cnv_end):
            return True
        if chrom2:
            if region_overlap(chrom2, start_pos2, end_pos2, cnv_chrom, cnv_start, cnv_end):
                return True
        return False

    def get_gene_region_dict(cursor):
        ret_dict = {}
        sql_str = 'SELECT gene_id, gene_name, chr, start_pos, end_pos, chr2, start_pos2, end_pos2 FROM gene_table'
        cursor.execute(sql_str)
        gene_data = cursor.fetchall()
        for i in gene_data:
            ret_dict["{0}_{1}".format(i[0], i[1])] = i
        return ret_dict

    # variance_table, db_file, phenotype,
    # sample_table_name, sample_restrict,
    # gene_table_name, gene_restrict,
    # variance_restrict, synonymous_restrict,
    # output, job_id, hg38_cnv_file
    print("variance_table={}".format(variance_table))
    print("db_file={}".format(db_file))
    print("phenotype={}".format(phenotype))
    print('sample_table_name={}'.format(sample_table_name))
    print("sample_restrict={}".format(sample_restrict))
    print("gene_table_name={}".format(gene_table_name))
    print("gene_restrict={}".format(gene_restrict))
    print("variance_resctict={}".format(variance_restrict))
    print("synonymous_restrict={}".format(synonymous_restrict))
    print("output={}".format(output))
    print("hg38_cnv_file={}".format(hg38_cnv_file))
    with open(hg38_cnv_file, 'r') as fp:
        cnv_data = [i.strip().split('\t') for i in fp.readlines()]
    cnv_case = [i for i in cnv_data if i[5] == '1']
    cnv_control = [i for i in cnv_data if i[5] == '0']

    gene_set_col_name_list = ["gene_set_Cluster3_CMs_STable1_Figure1_top500marker",
                              "gene_set_Cluster10_ST_STable1_Figure1_top500marker",
                              "gene_set_Cluster15_MLPs_STable1_Figure1_439marker",
                              "gene_set_Mesp1_KO_vs_WT_STable4_cluster7_pSHF_top500_up",
                              "gene_set_Mesp1_KO_vs_WT_STable4_top500_down_cluster4_MLPs",
                              "gene_set_Mesp1_KO_vs_WT_STable4_top500_down_cluster5_aSHF_SOM",
                              "gene_set_Mesp1_KO_vs_WT_STable4_top500_up_cluster4_MLPs",
                              "gene_set_Mesp1_KO_vs_WT_STable4_top500_up_cluster5_aSHF_SOM",
                              "gene_set_Mesp1_KO_vs_WT_STable4_Cluster7_pSHF_top500_down",
                              "gene_set_Mesp1_KO_vs_WT_STable4_Cluster8_CMs_top500_down",
                              "gene_set_Mesp1_KO_vs_WT_STable4_Cluster8_CMs_top500_up",
                              "gene_set_Tbx1_KO_vs_Het_STable5_Cluster5_MLPs_top500_down",
                              "gene_set_Tbx1_KO_vs_Het_STable5_Cluster5_MLPs_top500_up",
                              "gene_set_Tbx1Cre_KO_vs_Het_STable5_Cluster4_aSHF_SOM_top500_down",
                              "gene_set_Tbx1Cre_KO_vs_Het_STable5_Cluster4_aSHF_SOM_top500_up",
                              "gene_set_Tbx1Cre_KO_vs_Het_STable5_Cluster13_CMs_top500_up",
                              "gene_set_Tbx1Cre_KO_vs_WT_STable5_Cluster12_pSHF_top500_down",
                              "gene_set_Tbx1Cre_KO_vs_WT_STable5_Cluster12_pSHF_top500_up",
                              "gene_set_Tbx1Cre_KO_vs_WT_STable5_Cluster13_CMs_top500_down",
                              "gene_set_Tbx1Cre_WT_E9point5_MLPs_cluster7_marker",
                              "gene_set_Tbx1Cre_WT_E9point5_pSHF_cluster9_marker",
                              "gene_set_SysCilia_genes_302",
                              "gene_set_Cilia_genes_669",
                              "gene_set_FoxJ1_genes_116",
                              "gene_set_CHD_genes_no_cilia_genes_402",
                              "gene_set_Chromatin_163",
                              "gene_set_Notch_associated_130",
                              "gene_set_TGF_B_431",
                              "gene_set_Cytoskeletal_excluding_cilia_791",
                              "gene_set_Ser_Thr_kinases_47",
                              "gene_set_Brain_expressed_400",
                              "gene_set_Liver_expressed_400",
                              "gene_set_Mature_heart_left_vent_expressed_400",
                              "gene_set_Autism_Spectrum_Disorder_genes_86",
                              "gene_set_Fibroblast_growth_factor_signaling_genes_87",
                              "gene_set_Hedgehog_signaling_genes_149",
                              "gene_set_Platelet_derived_growth_factor_signaling_genes_116",
                              "gene_set_WNT_signaling_genes_297",
                              "gene_set_Housekeeping_genes_923",
                              "gene_set_Manually_curated_set_of_known_CHD_genes_by_tier_and_inheritance_mode_185",
                              "gene_set_Curated_known_Human_Mouse_CHD_genes_253",
                              "gene_set_folate_metabolism_pathway_34",
                              "gene_set_vitamin_b12_metabolism_21",
                              "gene_set_List_of_non_canonical_PCP_genes_54_",
                              "gene_set_TBX1_downregulated_Jun_chip_data_288",
                              "gene_set_abnormal_cardiac_outflow_tract_development_MP_0006126",
                              "gene_set_abnormal_heart_ventricle_outflow_tract_morphology_MP_0010224_84",
                              "gene_set_abnormal_heart_right_ventricle_outflow_tract_morphology_MP_0010428_48",
                              "gene_set_abnormal_heart_left_ventricle_outflow_tract_morphology_MP_0010429_62",
                              "gene_set_Deyou_Module_4_579",
                              "gene_set_Core_CellEG_957",
                              "gene_set_Dang_et_al_haploinsufficiency_262",
                              "gene_set_MGI_heterozygous_phenotype_313",
                              "gene_set_GDI__all_disease_causing_genes_807",
                              "gene_set_Apoptosis_136",
                              "gene_set_Axon_guidance_191",
                              "gene_set_Cell_cycle_124",
                              "gene_set_ECM_receptor_interaction_88",
                              "gene_set_ErbB_signaling_pathway_85",
                              "gene_set_Focal_adhesion_201",
                              "gene_set_Hedgehog_signaling_pathway_51",
                              "gene_set_HIF_1_signaling_pathway_109",
                              "gene_set_Hippo_signaling_pathway_157",
                              "gene_set_Insulin_signaling_pathway_139",
                              "gene_set_INTERGIN_FOCAL_160",
                              "gene_set_MAPK_signaling_pathway_294",
                              "gene_set_mTOR_signaling_pathway_156",
                              "gene_set_Notch_signaling_pathway_53",
                              "gene_set_Rap1_signaling_pathway_210",
                              "gene_set_TGF_beta_signaling_pathway_94",
                              "gene_set_VEGF_signaling_pathway_59",
                              "gene_set_Wnt_signaling_pathway_160",
                              "gene_set_BMP_190",
                              "gene_set_CRK_CRKL_60",
                              "gene_set_FGF_141",
                              "gene_set_Intergin_and_focal_adhesion_signaling_164",
                              "gene_set_MAPK_76",
                              "gene_set_PDGF_97",
                              "gene_set_SMAD_171",
                              "gene_set_TGF_192",
                              "gene_set_WNT_159",
                              "gene_set_Collagen_73",
                              "gene_set_gevir_percentile_5_967",
                              "gene_set_loeuf_percentile_5_967",
                              "gene_set_virlof_percentile_5_968",
                              "gene_set_HIS_percentile_5_767",
                              "gene_set_GHIS_percentile_5_818",
                              "gene_set_RVIS_EVS_percentile_5_845",
                              "gene_set_RVIS_percentile_ExAC_5_579"]

    synonymous_table = "synonymous_snp"
    if should_log:
        logging.basicConfig(filename="{0}{1}.log".format(sys._getframe().f_code.co_name, job_id),
                            level=logging.DEBUG, format=log_format, filemode="w")
        logging.debug("begin")
        logging.debug("db_file=[{0}]".format(db_file))
        logging.debug("phenotype=[{0}]".format(phenotype))
        logging.debug("sample_table_name=[{0}]".format(sample_table_name))
        logging.debug("sample_restrict=[{0}]".format(sample_restrict))
        logging.debug("variance_restrict=[{0}]".format(variance_restrict))
        logging.debug("synonymous_restrict=[{0}]".format(synonymous_restrict))
        logging.debug("output=[{0}]".format(output))

    if phenotype not in ["CTD"]:
        if should_log:
            logging.error("illegal phenotype [{}]".format(phenotype))
        return

    if should_log:
        logging.debug("begin select control list and case list...")
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()
    gene_region_dict = get_gene_region_dict(cursor)
    format_sample_restrict = "" if not sample_restrict else " AND {}".format(sample_restrict)
    cursor.execute("SELECT s.gen_id FROM {2} AS s WHERE s.{0}='0'{1}".format(phenotype,
                                                                             format_sample_restrict,
                                                                             sample_table_name))
    control_id_list = [i[0] for i in cursor.fetchall()]  # m
    print("control number={}".format(len(control_id_list)))
    cursor.execute("SELECT s.gen_id FROM {2} AS s WHERE s.{0}='1'{1}".format(phenotype,
                                                                             format_sample_restrict,
                                                                             sample_table_name))
    case_id_list = [i[0] for i in cursor.fetchall()]  # n
    print("case number={}".format(len(case_id_list)))
    if should_log:
        logging.debug("begin handle genes")
    with open(output, "w") as fp:
        fp.write("##phenotype:\"{0}\"\n##control number:{1}\n##case number:{2}\n"
                 "##sample restrict:\"{3}\"\n##version: 2 (select variance with gene name instead of gene region)\n"
                 "".format(phenotype, len(control_id_list), len(case_id_list), sample_restrict))
        fp.write("""##                         case     control
##                      ---------------------
##                      |         |         |
##    alt allele number |    A    |    B    |
##                      |         |         |
##                      ---------------------
##                      |         |         |
##    ref allele number |    C    |    D    |
##                      |         |         |
##                      ---------------------
##
##
##                              case     control
##                           ---------------------
##                           |         |         |
##          subjects have alt| A1+cnvA | B1+cnvB |
##                           |         |         |
##                           ---------------------
##                           |         |         |
##  subjects do not have alt | C1+cnvC | D1+cnvD |
##                           |         |         |
##                           ---------------------
""")
        fp.write(
            "#gene_id\tgene_name\tA\tB\tC\tD\tp_value\todds_ratio\tsynonymous_case_alt(C)\tsynonymous_control_alt(D)\todds_ratio_synonymous\tA1\tB1\tC1\tD1\tp_value1\todds_ratio1\t"
            "n_variance_in_gene\tvariance_in_gene_name\tn_variance_control\tvariance_in_control_name\t"
            "n_variance_case\tvariance_in_case_name\tmouse_mean_WT_OFT_TPM_2\tmouse_mean_WT_LV_TPM_3\tmouse_mean_WT_PA_TPM_3\t"
            "Interactions_IntAct\tInteractions_BioGRID\tInteractions_ConsensusPathDB\tHIPred\tGHIS\tgnomAD_pLI\t"
            "GDI_Phred\tGene_damage_prediction__all_disease_causing_genes\tEssential_gene\tEssential_gene_CRISPR\t"
            "Gene_indispensability_pred\tgevir_percentile\tloeuf_percentile\tvirlof_percentile\t"
            "HHE_E14.5\t{}\n".format("\t".join(gene_set_col_name_list)))
        # gene data
        format_gene_restrict = "" if not gene_restrict else " WHERE {}".format(gene_restrict)
        cmd_str = "SELECT g.gene_id, g.gene_name, mouse_mean_WT_OFT_TPM_2, mouse_mean_WT_LV_TPM_3, mouse_mean_WT_PA_TPM_3, " \
                  "Interactions_IntAct_, Interactions_BioGRID_, Interactions_ConsensusPathDB_, " \
                  "HIPred, GHIS, gnomAD_pLI, GDI_Phred, Gene_damage_prediction__all_disease_causing_genes_ , " \
                  "Essential_gene, Essential_gene_CRISPR , Gene_indispensability_pred , gevir_percentile, " \
                  "loeuf_percentile, virlof_percentile, " \
                  "\"HHE_E14.5\", {2} FROM {0} AS g{1}".format(gene_table_name, format_gene_restrict,
                                                               ", ".join(gene_set_col_name_list))
        logging.debug("gene sql = [{}]".format(cmd_str))
        cursor.execute(cmd_str)
        gene_data = cursor.fetchall()
        print("gene number = {0}".format(len(gene_data)))
        # prepare the variance data
        print("build id data ...")
        cmd_str = "select id, {0} from {1} where {2}".format(sub_gene_name(variance_table, variance_restrict),
                                                             variance_table,
                                                             variance_restrict)
        cursor.execute(cmd_str)
        variance_id_data = [[i[0], build_gene_name_list(filter(lambda x: x is not None, i[1:]))] for i in
                            cursor.fetchall()]

        cmd_str = "select id, {0} from {1} where {2}".format(sub_gene_name(synonymous_table, synonymous_restrict),
                                                             synonymous_table,
                                                             synonymous_restrict)
        cursor.execute(cmd_str)
        synonymous_id_data = [[i[0], build_gene_name_list(filter(lambda x: x is not None, i[1:]))] for i in
                              cursor.fetchall()]
        cmd_str = "SELECT vcf_id, chr, pos, "
        for control_id in control_id_list:
            cmd_str = "{0}sample_{1}, ".format(cmd_str, control_id)
        for case_id in case_id_list:
            cmd_str = "{0}sample_{1}, ".format(cmd_str, case_id)
        cmd_str = cmd_str.strip(", ")
        assert variance_restrict != ""
        # format_variance_restrict = " WHERE {0}".format(variance_restrict)
        cmd_synonymous_str = "{0} FROM synonymous_snp as v".format(cmd_str)
        cmd_str = "{0} FROM variance as v".format(cmd_str)

        ret_str = ""
        icounter = 0
        gene_data_len = len(gene_data)
        for gene_data_list in gene_data:
            gene_id = gene_data_list[0]
            gene_name = gene_data_list[1]
            mouse_mean_WT_OFT_TPM_2 = gene_data_list[2]
            mouse_mean_WT_LV_TPM_3 = gene_data_list[3]
            mouse_mean_WT_PA_TPM_3 = gene_data_list[4]
            Interactions_IntAct_ = gene_data_list[5]
            Interactions_BioGRID_ = gene_data_list[6]
            Interactions_ConsensusPathDB_ = gene_data_list[7]
            HIPred = gene_data_list[8]
            GHIS = gene_data_list[9]
            gnomAD_pLI = gene_data_list[10]
            GDI_Phred = gene_data_list[11]
            Gene_damage_prediction__all_disease_causing_genes_ = gene_data_list[12]
            Essential_gene = gene_data_list[13]
            Essential_gene_CRISPR = gene_data_list[14]
            Gene_indispensability_pred = gene_data_list[15]
            gevir_percentile = gene_data_list[16]
            loeuf_percentile = gene_data_list[17]
            virlof_percentile = gene_data_list[18]
            hhe = gene_data_list[19]
            if should_log:
                if icounter % 100 == 0 and icounter > 0:
                    logging.debug("handled {0} / {1} genes".format(icounter, gene_data_len))
            variance_id_list = [str(i[0]) for i in variance_id_data if gene_id is not None and (
                    gene_name in i[1] or gene_id in i[1]) or gene_id is None and gene_name in i[1]]
            synonymous_id_list = [str(i[0]) for i in synonymous_id_data if gene_id is not None and (
                    gene_name in i[1] or gene_id in i[1]) or gene_id is None and gene_name in i[1]]
            if len(variance_id_list) == 0:
                icounter += 1
                continue
            tmp_cmd_str = cmd_str + " where v.id in (" + ", ".join(variance_id_list) + ")"
            # print(tmp_cmd_str)
            cursor.execute(tmp_cmd_str)

            variance_selected_list = [list(i) for i in cursor.fetchall()]

            tmp_cmd_str = cmd_synonymous_str + " where v.id in (" + ", ".join(synonymous_id_list) + ")"
            cursor.execute(tmp_cmd_str)
            synonymous_selected_list = [list(i) for i in cursor.fetchall()]

            all_variance_in_gene_num = len(variance_selected_list)
            all_variance_in_gene_name_list = [i[0] for i in variance_selected_list]  # vcf_id
            control_data = [i[3:3 + len(control_id_list)] for i in variance_selected_list]
            case_data = [i[3 + len(control_id_list):3 + len(control_id_list) + len(case_id_list)] for i in
                         variance_selected_list]

            synonymous_control_data = [i[3:3 + len(control_id_list)] for i in synonymous_selected_list]
            synonumous_case_data = [i[3 + len(control_id_list):3 + len(control_id_list) + len(case_id_list)] for i in
                                    synonymous_selected_list]
            tmp_list = []
            control_variance_in_gene_num = 0
            control_variance_in_gene_name_list = []
            for i in xrange(len(control_data)):
                control_line_data = filter(lambda x: type(x) == int, control_data[i])
                tmp_list.extend(control_line_data)
                if sum(control_line_data) > 0:
                    control_variance_in_gene_num += 1
                    control_variance_in_gene_name_list.append(variance_selected_list[i][0])
            B = sum(tmp_list)
            D = 2 * len(tmp_list) - B

            tmp_list = []
            for i in xrange(len(synonymous_control_data)):
                synonymous_control_line_data = filter(lambda x: type(x) == int, synonymous_control_data[i])
                tmp_list.extend(synonymous_control_line_data)
            DD = sum(tmp_list)

            tmp_list = []
            case_variance_in_gene_num = 0
            case_variance_in_gene_name_list = []
            for i in xrange(len(case_data)):
                case_line_data = filter(lambda x: type(x) == int, case_data[i])
                tmp_list.extend(case_line_data)
                if sum(case_line_data) > 0:
                    case_variance_in_gene_num += 1
                    case_variance_in_gene_name_list.append(variance_selected_list[i][0])
            A = sum(tmp_list)
            C = 2 * len(tmp_list) - A

            tmp_list = []
            for i in xrange(len(synonumous_case_data)):
                synonymous_case_line_data = filter(lambda x: type(x) == int, synonumous_case_data[i])
                tmp_list.extend(synonymous_case_line_data)
            CC = sum(tmp_list)

            control_people_data = [sum(filter(lambda x: type(x) == int, people_line)) for people_line in
                                   zip(*control_data)]
            B1 = len(filter(lambda x: x > 0, control_people_data))
            D1 = len(control_people_data) - B1
            case_people_data = [sum(filter(lambda x: type(x) == int, people_line)) for people_line in zip(*case_data)]
            A1 = len(filter(lambda x: x > 0, case_people_data))
            C1 = len(case_people_data) - A1

            gene_region_line = gene_region_dict['{0}_{1}'.format(gene_id, gene_name)]
            cnv_A = len(Counter(
                [cnv_line[3] for cnv_line in cnv_case if
                 is_cnv_gene_overlap(gene_region_line, cnv_line)]).most_common())
            cnv_B = len(Counter(
                [cnv_line[3] for cnv_line in cnv_control if
                 is_cnv_gene_overlap(gene_region_line, cnv_line)]).most_common())
            cnv_case_num = len(Counter([i[3] for i in cnv_case]).most_common())
            cnv_control_num = len(Counter([i[3] for i in cnv_control]).most_common())
            cnv_C = cnv_case_num - cnv_A
            cnv_D = cnv_control_num - cnv_B
            A1 = A1 + cnv_A
            B1 = B1 + cnv_B
            C1 = C1 + cnv_C
            D1 = D1 + cnv_D
            # if A1 + B1 <= 2:
            #     icounter += 1
            #     continue  # the number of people with variance is less than 2. Go to next gene.

            oddsratio, pvalue = stats.fisher_exact([[A, B], [C, D]])
            oddsratio1, pvalue1 = stats.fisher_exact([[A1, B1], [C1, D1]])
            or_new = "inf" if B == 0 or CC == 0 else float(A * DD) / (B * CC)
            ret_str = "{0}{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t" \
                      "{11}\t{12}\t{13}\t{14}\t{15}\t{16}\t{17}\t{18}\t{19}\t{20}\t" \
                      "{21}\t{22}\t{23}\t{24}\t{25}\t{26}\t{27}\t{28}\t{29}\t{30}\t" \
                      "{31}\t{32}\t{33}\t{34}\t{35}\t{36}\t{37}\t{38}\t{39}\t{40}\t" \
                      "{41}\t{42}\n" \
                      "".format(ret_str, gene_id, gene_name, A, B, C, D, pvalue, oddsratio, CC, DD, or_new,
                                A1, B1, C1, D1,
                                pvalue1, oddsratio1,
                                all_variance_in_gene_num,
                                ";".join(all_variance_in_gene_name_list) if len(
                                    all_variance_in_gene_name_list) > 0 else ".",
                                control_variance_in_gene_num,
                                ";".join(control_variance_in_gene_name_list) if len(
                                    control_variance_in_gene_name_list) > 0 else ".",
                                case_variance_in_gene_num,
                                ";".join(case_variance_in_gene_name_list) if len(
                                    case_variance_in_gene_name_list) > 0 else ".",
                                mouse_mean_WT_OFT_TPM_2, mouse_mean_WT_LV_TPM_3, mouse_mean_WT_PA_TPM_3,
                                Interactions_IntAct_, Interactions_BioGRID_,
                                Interactions_ConsensusPathDB_, HIPred,
                                GHIS, gnomAD_pLI,
                                GDI_Phred, Gene_damage_prediction__all_disease_causing_genes_,
                                Essential_gene, Essential_gene_CRISPR,
                                Gene_indispensability_pred, gevir_percentile,
                                loeuf_percentile, virlof_percentile,
                                hhe, "\t".join([str(i) for i in gene_data_list[20:]]))
            icounter += 1
            print("handled {0} / {1} genes".format(icounter, gene_data_len))

        print("handled {0} / {1} genes in total".format(icounter, gene_data_len))
        fp.write(ret_str)
    if should_log:
        logging.debug("all done")


def analyze_CTD_gene_fisher_with_cnv(database, hg38_cnv_file, path):
    table_name_element_dict = {"": "",
                               " AND (annovar = 1)": "annovar",
                               " AND (bystro = 1)": "bystro",
                               " AND (dmis = 1)": "dmis",
                               " AND (dsplicing = 1)": "dsplicing",
                               " AND (bystro=1 or annovar=1 or vep=1 or dmis=1 or dsplicing=1 or spliceAI=1)": "all.6",
                               " AND (spliceAI = 1)": "spliceAI",
                               " AND (vep = 1)": "vep",
                               " AND (annovar = 1 or bystro = 1 or vep = 1)": "LOF",
                               " AND (dsplicing = 1 or spliceAI = 1)": "dsplicing.all",
                               " AND (annovar = 1 and bystro = 1)": "annovar.bystro",
                               "tof6": "tof",
                               "CTD": "CTD",
                               "heart6": "heart6",
                               "bystro_sampleMaf <= 0.01": "maf.01",
                               "bystro_sampleMaf <= 0.05": "maf.05",
                               "bystro_sampleMaf <= 0.1": "maf.1",
                               "bystro_sampleMaf <= 1": "",
                               # "bystro_sampleMaf <= 0.2": "maf.2",
                               # "bystro_sampleMaf <= 0.3": "maf.3",
                               " AND (bystro_cadd>=10)": "cadd10",
                               " AND (bystro_cadd>=15)": "cadd15",
                               " AND (bystro_cadd>=20)": "cadd20",
                               " AND bystro_phastCons >= 0.4": "Cons.4",
                               " AND bystro_phastCons >= 0.5": "Cons.5",
                               " AND bystro_phastCons >= 0.6": "Cons.6",
                               " AND bystro_phyloP >= -1": "loPn1",
                               " AND bystro_phyloP >= 0": "loP0",
                               " AND bystro_phyloP >= 1": "loP1",
                               " AND (ccrs >= 95)": "ccrs95",
                               " AND (ccrs >= 90)": "ccrs90",
                               " AND (ccrs >= 85)": "ccrs85",
                               " AND (ccrs >= 80)": "ccrs80",
                               " AND (domain_limbr <= {})".format(domain_limbr_quartile_dict[40]): "dlimbr40",
                               " AND (domain_limbr <= {})".format(domain_limbr_quartile_dict[50]): "dlimbr50",
                               " AND (domain_limbr <= {})".format(domain_limbr_quartile_dict[60]): "dlimbr60",
                               " AND (domain_limbr <= {})".format(domain_limbr_quartile_dict[70]): "dlimbr70",
                               " AND (exone_limbr <= {})".format(exone_limbr_quartile_dict[40]): "elimbr40",
                               " AND (exone_limbr <= {})".format(exone_limbr_quartile_dict[50]): "elimbr50",
                               " AND (exone_limbr <= {})".format(exone_limbr_quartile_dict[60]): "elimbr60",
                               " AND (exone_limbr <= {})".format(exone_limbr_quartile_dict[70]): "elimbr70",
                               " AND (is_ccds = 1)": "ccds1"}

    pp = Popen(["mkdir -p {}".format(path)], shell=True)
    pp.wait()
    annotator_list = [" AND (bystro=1 or annovar=1 or vep=1 or dmis=1 or dsplicing=1 or spliceAI=1)",
                      " AND (annovar = 1)",
                      " AND (bystro = 1)",
                      " AND (dmis = 1)",
                      " AND (vep = 1)",
                      " AND (annovar = 1 or bystro = 1 or vep = 1)",
                      " AND (dsplicing = 1 or spliceAI = 1)",
                      " AND (annovar = 1 and bystro = 1)"]
    icounter = 1
    icounter2 = 0
    script_str = ''
    for phenotype in ["tof6", "CTD", "heart6"]:
        for annotator in annotator_list:
            for freq in ["bystro_sampleMaf <= 1"]:
                for cadd in ["", " AND (bystro_cadd>=10)", " AND (bystro_cadd>=15)", " AND (bystro_cadd>=20)"]:
                    for ph in ["",
                               " AND bystro_phastCons >= 0.4",
                               " AND bystro_phastCons >= 0.5",
                               " AND bystro_phastCons >= 0.6",
                               " AND bystro_phyloP >= -1",
                               " AND bystro_phyloP >= 0",
                               " AND bystro_phyloP >= 1"]:
                        for regional_constraint in ["",
                                                    " AND (ccrs >= 95)",
                                                    " AND (ccrs >= 90)",
                                                    " AND (ccrs >= 85)",
                                                    " AND (ccrs >= 80)",
                                                    " AND (domain_limbr <= {})".format(domain_limbr_quartile_dict[40]),
                                                    " AND (domain_limbr <= {})".format(domain_limbr_quartile_dict[50]),
                                                    " AND (domain_limbr <= {})".format(domain_limbr_quartile_dict[60]),
                                                    " AND (domain_limbr <= {})".format(domain_limbr_quartile_dict[70]),
                                                    " AND (exone_limbr <= {})".format(exone_limbr_quartile_dict[40]),
                                                    " AND (exone_limbr <= {})".format(exone_limbr_quartile_dict[50]),
                                                    " AND (exone_limbr <= {})".format(exone_limbr_quartile_dict[60]),
                                                    " AND (exone_limbr <= {})".format(exone_limbr_quartile_dict[70])]:
                            for is_ccds in ["", " AND (is_ccds = 1)"]:
                                variance_restrict = "{0}{1}{2}{3}{4}{5}".format(freq, annotator, cadd, ph,
                                                                                regional_constraint, is_ccds)
                                name_element_list = [table_name_element_dict[phenotype],
                                                     table_name_element_dict[annotator],
                                                     table_name_element_dict[freq],
                                                     table_name_element_dict[cadd],
                                                     table_name_element_dict[ph],
                                                     table_name_element_dict[regional_constraint],
                                                     table_name_element_dict[is_ccds]]
                                output_name = "{0}.cnv.fet" \
                                              "".format("_".join(filter(lambda x: len(x) > 0, name_element_list)))
                                output_name = os.path.join(path, output_name)
                                if script_str == "":
                                    script_str = """#!/bin/bash
#$ -j y
#$ -S /bin/bash
#$ -cwd
#$ -N fgc{3}
#$ -pe smp 1
#$ -l h_vmem=8G
#$ -m bes
# -q all.q
export LC_ALL=C
export MALLOC_ARENA_MAX=4
module load python/2.7.15/gcc.4.4.7
module load sqlite3/3.8.11/gcc.4.4.7
time=`date`
echo "==START $time =="
wgsa.py build_contingency_table_gene_name_cnv variance {0} CTD sampleChdPhenotype \"\" gene_table \"\" \"{1}\" \"vep=1\" {2} {3} {4}
echo {2} is done
time=`date`
echo == $time ==
""".format(database, variance_restrict, output_name, icounter2, hg38_cnv_file)
                                # build_contingency_table_gene_name_cnv(variance_table, db_file, phenotype,
                                #                                       sample_table_name, sample_restrict,
                                #                                       gene_table_name, gene_restrict,
                                #                                       variance_restrict, synonymous_restrict,
                                #                                       output, job_id, hg38_cnv_file,
                                else:
                                    script_str += "wgsa.py build_contingency_table_gene_name_cnv " \
                                                  "variance {0} CTD sampleChdPhenotype \"\" gene_table \"\" \"{1}\" \"vep=1\" {2} {3} {4}" \
                                                  "\necho {2} is done" \
                                                  "\ntime=`date`\necho == $time ==\n" \
                                                  "".format(database, variance_restrict, output_name, icounter2,
                                                            hg38_cnv_file)

                                if icounter == 53:
                                    shell_file_name = os.path.join(path, "qsub{}.sh".format(icounter2))
                                    with open(shell_file_name, "w") as fp:
                                        fp.write(script_str + "\ndate; echo \"==END==\"")
                                    while True:
                                        if os.access(shell_file_name, os.R_OK):
                                            break
                                        time.sleep(1)
                                    print("qsubing {} batch".format(icounter2 + 1))
                                    pp = Popen(["qsub {}".format(shell_file_name)], shell=True)
                                    pp.wait()
                                    script_str = ""
                                    icounter = 0
                                    icounter2 += 1
                                icounter += 1
    if len(script_str) > 0:
        shell_file_name = os.path.join(path, "qsub{}.sh".format(icounter2))
        with open(shell_file_name, "w") as fp:
            fp.write(script_str + "\ndate; echo \"==END==\"")
        while True:
            if os.access(shell_file_name, os.R_OK):
                break
            time.sleep(1)
        print("qsubing {} batch".format(icounter2 + 1))
        pp = Popen(["qsub {}".format(shell_file_name)], shell=True)
        pp.wait()

    print("All the jobs has been submitted. Job number = {}\n"
          "If any job has problem, kill it. And qsub the corresponding sh file".format(icounter2 + 1))


def analyze_variant_in_cnv(hg38_cnv_file, database):
    logging.basicConfig(filename="analyze_variant_in_cnv.log", level=logging.DEBUG, format=log_format, filemode="w")
    annotator_list = [" AND (bystro=1 or annovar=1 or vep=1 or dmis=1 or dsplicing=1 or spliceAI=1)",
                      " AND (annovar = 1)",
                      " AND (bystro = 1)",
                      " AND (dmis = 1)",
                      " AND (vep = 1)",
                      " AND (annovar = 1 or bystro = 1 or vep = 1)",
                      " AND (dsplicing = 1 or spliceAI = 1)",
                      " AND (annovar = 1 and bystro = 1)"]

    for annotator in annotator_list:
        for freq in ["bystro_sampleMaf <= 1"]:
            for cadd in ["", " AND (bystro_cadd>=10)", " AND (bystro_cadd>=15)", " AND (bystro_cadd>=20)"]:
                for ph in ["",
                           " AND bystro_phastCons >= 0.4",
                           " AND bystro_phastCons >= 0.5",
                           " AND bystro_phastCons >= 0.6",
                           " AND bystro_phyloP >= -1",
                           " AND bystro_phyloP >= 0",
                           " AND bystro_phyloP >= 1"]:
                    for regional_constraint in ["",
                                                " AND (ccrs >= 95)",
                                                " AND (ccrs >= 90)",
                                                " AND (ccrs >= 85)",
                                                " AND (ccrs >= 80)",
                                                " AND (domain_limbr <= {})".format(domain_limbr_quartile_dict[40]),
                                                " AND (domain_limbr <= {})".format(domain_limbr_quartile_dict[50]),
                                                " AND (domain_limbr <= {})".format(domain_limbr_quartile_dict[60]),
                                                " AND (domain_limbr <= {})".format(domain_limbr_quartile_dict[70]),
                                                " AND (exone_limbr <= {})".format(exone_limbr_quartile_dict[40]),
                                                " AND (exone_limbr <= {})".format(exone_limbr_quartile_dict[50]),
                                                " AND (exone_limbr <= {})".format(exone_limbr_quartile_dict[60]),
                                                " AND (exone_limbr <= {})".format(exone_limbr_quartile_dict[70])]:
                        for is_ccds in ["", " AND (is_ccds = 1)"]:
                            variance_restrict = "{0}{1}{2}{3}{4}{5}".format(freq, annotator, cadd, ph,
                                                                            regional_constraint, is_ccds)
                            shell_str = './variant_in_cnv.sh {0} {1} \"{2}\"'.format(hg38_cnv_file, database,
                                                                                     variance_restrict)
                            logging.debug("regional_constraint=[{}]".format(regional_constraint))
                            logging.debug("variance_restrict=[{}]".format(variance_restrict))
                            logging.debug("shell str=[{}]".format(shell_str))
                            phandle = Popen([shell_str], shell=True)
                            phandle.wait()


def variance_in_cnv_sample(sample_id_mapping, db_file, variance_restrict, hg38_cnv):
    from scipy.stats import mannwhitneyu
    def load_sample_id_dict(sample_id_mapping):
        with open(sample_id_mapping, 'r') as fp:
            data = [j for j in [i.strip().split("\t") for i in fp.readlines()] if len(j) == 2]
        data2 = [[i[1], i[0]] for i in data]
        return [dict(data), dict(data2)]

    def load_hg38_cnv(hg38_cnv, famid2genid_dict):
        with open(hg38_cnv, 'r') as fp:
            data = [i.strip().split("\t") for i in fp.readlines()]
        ret_dict = {}
        for i in data:
            if i[3] not in famid2genid_dict:
                continue
            igenid = int(famid2genid_dict[i[3]])
            if igenid not in ret_dict:
                ret_dict[igenid] = [[i[0], int(i[1]), int(i[2])]]
            else:
                ret_dict[igenid].append([i[0], int(i[1]), int(i[2])])
        return ret_dict

    variance_restrict = variance_restrict.strip()
    famid2genid_dict, genid2famid = load_sample_id_dict(sample_id_mapping)
    genid2cnv_cict = load_hg38_cnv(hg38_cnv, famid2genid_dict)
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()
    cursor.execute('select gen_id from sampleChdPhenotype where CTD="1"')
    case_id_list = [i[0] for i in cursor.fetchall()]
    cursor.execute('select gen_id from sampleChdPhenotype where CTD="0"')
    control_id_list = [i[0] for i in cursor.fetchall()]
    case_num_list = []
    control_num_list = []
    icounter = 0
    for case_id in case_id_list:
        icounter += 1
        print("case {0} / {1}".format(icounter, len(case_id_list)))
        if case_id not in genid2cnv_cict:
            continue
        cnv_list = genid2cnv_cict[case_id]
        tmp = 0
        for cnv_chr, cnv_start, cnv_end in cnv_list:
            cnv_chr = cnv_chr.strip('chr')
            sql_str = 'select chr, pos from variance where sample_{0} > 0 and sample_{0}!="na" ' \
                      'and chr="{1}" and pos>{2} and pos<{3} {4}' \
                      ''.format(case_id, cnv_chr, cnv_start, cnv_end,
                                "" if len(variance_restrict) == 0 else "and ({})".format(variance_restrict))
            cursor.execute(sql_str)
            variance_data = cursor.fetchall()
            tmp += len(variance_data)
        case_num_list.append(tmp)
    icounter = 0
    for control_id in control_id_list:
        icounter += 1
        print("control {0} / {1}".format(icounter, len(control_id_list)))
        if control_id not in genid2cnv_cict:
            continue
        cnv_list = genid2cnv_cict[control_id]
        tmp = 0
        for cnv_chr, cnv_start, cnv_end in cnv_list:
            # cnv_chr, cnv_start, cnv_end = genid2cnv_cict[control_id]
            cnv_chr = cnv_chr.strip('chr')
            sql_str = 'select chr, pos from variance where sample_{0} > 0 and sample_{0}!="na" ' \
                      'and chr="{1}" and pos>{2} and pos<{3} {4}' \
                      ''.format(control_id, cnv_chr, cnv_start, cnv_end,
                                "" if len(variance_restrict) == 0 else "and ({})".format(variance_restrict))
            cursor.execute(sql_str)
            variance_data = cursor.fetchall()
            tmp += len(variance_data)
        control_num_list.append(tmp)
    print('case_num_list={}'.format(case_num_list))
    print('control_num_list={}'.format(control_num_list))
    ret = mannwhitneyu(case_num_list, control_num_list)
    print(ret.pvalue)
    return ret.pvalue


def cnv_overlap_gene(hg38_cnv_file, db_file, output):
    def get_gene_region_dict(cursor):
        ret_dict = {}
        sql_str = 'SELECT gene_id, gene_name, chr, start_pos, end_pos, chr2, start_pos2, end_pos2 FROM gene_table'
        cursor.execute(sql_str)
        gene_data = cursor.fetchall()
        for i in gene_data:
            ret_dict["{0}_{1}".format(i[0], i[1])] = i
        return ret_dict

    def is_cnv_gene_overlap(gene_line, cnv_line):
        def region_overlap(chr1, start1, end1, chr2, start2, end2):
            if chr1 != chr2:
                return False
            if end1 < start2 or start1 > end2:
                return False
            return True

        chrom, start_pos, end_pos, chrom2, start_pos2, end_pos2 = gene_line[2:]
        cnv_chrom, cnv_start, cnv_end = cnv_line[:3]
        start_pos = int(start_pos)
        end_pos = int(end_pos)
        chrom = str(chrom)
        if chrom2 is not None:
            start_pos2 = int(start_pos2)
            end_pos2 = int(end_pos2)
            chrom2 = str(chrom2)
        cnv_chrom = cnv_chrom.strip('chr')
        cnv_start = int(cnv_start)
        cnv_end = int(cnv_end)
        if region_overlap(chrom, start_pos, end_pos, cnv_chrom, cnv_start, cnv_end):
            return True
        if chrom2:
            if region_overlap(chrom2, start_pos2, end_pos2, cnv_chrom, cnv_start, cnv_end):
                return True
        return False

    with open(hg38_cnv_file, 'r') as fp:
        cnv_data = [re.split("[\t;]", i.strip()) for i in fp.readlines()]
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()
    gene2region_dict = get_gene_region_dict(cursor)
    conn.commit()
    conn.close()
    # total_len = len(gene2region_dict)
    total_len = len(cnv_data)
    icounter = 0
    with open(output, "w") as fp:
        # fp.write("cnv_chr\tcnv_start\tcnv_end\tcnv_famID\tcnv_heart6\tcnv_CTD\tcnv_sampleID\tcnv_alternate_CEL_file\t"
        #          "cnv_cluster\tcnv_width\tcnv_snpn\tcnv_typecode\tcnv_index\tcnv_freq_cnt\tcnv_freq\tcnv_gene\t"
        #          "gene_id\tgene_name\tgene_chr\t"
        #          "gene_start_pos\tgene_end_pos\tgene_chr2\tgene_start_pos2\tgene_end_pos2\n")
        # for gene in gene2region_dict:
        #     gene_line = gene2region_dict[gene]
        #     for overlapped_cnv_line in [cnv_line for cnv_line in cnv_data if is_cnv_gene_overlap(gene_line, cnv_line)]:
        #         fp.write("{0}\t{1}\n".format("\t".join(overlapped_cnv_line),
        #                                    "\t".join([str(i) for i in gene_line])))
        #     icounter += 1
        #     print("{0} / {1}".format(icounter, total_len))
        fp.write("cnv_chr\tcnv_start\tcnv_end\tcnv_famID\tcnv_heart6\tcnv_CTD\tcnv_sampleID\tcnv_alternate_CEL_file\t"
                 "cnv_cluster\tcnv_width\tcnv_snpn\tcnv_typecode\tcnv_index\tcnv_freq_cnt\tcnv_freq\tcnv_gene\t"
                 "genes\n")
        for cnv_line in cnv_data:
            gene_str = ",".join(
                [gene for gene in gene2region_dict if is_cnv_gene_overlap(gene2region_dict[gene], cnv_line)])
            fp.write("{0}\t{1}\n".format("\t".join(cnv_line),
                                         gene_str))
            icounter += 1
            print("{0} / {1}".format(icounter, total_len))


def grep_gatk_variance(gatk_path, variant_sample_info_file, output_path, mode, window=100):
    def load_variance_sample_info_file(variant_sample_info_file):
        variance_dict = {}
        with open(variant_sample_info_file, 'r') as fp:
            data_line = fp.readline()
            while True:
                data_line = fp.readline()
                if not data_line:
                    break
                if not data_line.strip():
                    continue
                target_file, chrom, pos, ref, alt, gen_id = data_line.strip().split('\t')
                if target_file not in variance_dict:
                    variance_dict[target_file] = [[chrom, pos, ref, alt, gen_id]]
                else:
                    variance_dict[target_file].append([chrom, pos, ref, alt, gen_id])
        return variance_dict

    variance_dict = load_variance_sample_info_file(variant_sample_info_file)
    for target_file in variance_dict:
        output_file = os.path.join(output_path, "target_{}".format(target_file))
        target_file_whole = os.path.join(gatk_path, target_file)
        print('handling {}'.format(target_file))
        sys.stdout.flush()
        with open(output_file, 'w') as fp:
            for chrom, pos, ref, alt, gen_id in variance_dict[target_file]:
                if mode == 0:
                    cmd_str = "grep -P '{0}\\t{1}\\t' {2}".format(chrom, pos, target_file_whole)
                else:
                    cmd_str = "awk '{{if($1=={0} && $2<{1} && $2>{2}){{print $0}}}}' {3}".format(chrom,
                                                                                                 int(pos) + int(window),
                                                                                                 int(pos) - int(window),
                                                                                                 target_file_whole)
                # print(cmd_str)
                pp = Popen([cmd_str], shell=True, stdout=PIPE)
                ret = pp.stdout.readlines()
                if not ret:
                    print("ref={0} alt={1} cmdstr {2} get no result".format(ref, alt, cmd_str))
                    continue
                if len(ret) > 1:
                    print('ref={0} alt={1} cmdstr {2} get more than 1 retults'.format(ref, alt, cmd_str))
                for i in range(len(ret)):
                    fp.write(
                        "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}".format(target_file, chrom, pos, ref, alt, gen_id, ret[i]))


def check_PE_GATK_master(db_file, table):
    def return_task_queue():
        global task_queue
        return task_queue

    def return_result_queue():
        global result_queue
        return result_queue

    def return_msgc_queue():
        global msgc_queue
        return msgc_queue

    multiprocessing.freeze_support()
    BaseManager.register('get_task_queue', callable=return_task_queue)
    BaseManager.register('get_result_queue', callable=return_result_queue)
    BaseManager.register('get_msgc_queue', callable=return_msgc_queue)
    manager = BaseManager(address=('', 5000), authkey=b"wgsawgsa")
    # start Queue:
    manager.start()
    # get Queue objects:
    task_queue = manager.get_task_queue()
    result_queue = manager.get_result_queue()
    msgc_queue = manager.get_msgc_queue()
    # print("starting writer...")
    # logging.debug("starting writer...")
    # proc_write_result = multiprocessing.Process(target=cross_genehancer_snp_writer, args=(result_queue, output_file))
    # proc_write_result.daemon = True
    # proc_write_result.start()
    print("loading variance data...")
    logging.debug("loading variance data...")
    # Assign task
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()
    cursor.execute('select * from {}'.format(table))
    variance_data = cursor.fetchall()
    job_num = 1000
    step = len(variance_data) / (job_num - 1)

    for i in range(job_num):
        if i < job_num - 1:
            sub_variance_data = variance_data[i * step:(i * step + step)]
        else:
            sub_variance_data = variance_data[i * step:]
        task_queue.put(sub_variance_data)

    print("waiting for tasks...")
    logging.debug("waiting for tasks...")
    # wait for tasks
    task_queue.join()
    # result_queue.put("STOP", block=True)
    # proc_write_result.join()
    worker_ip_list = []
    while 1:
        if msgc_queue.empty():
            break
        else:
            worker_ip_list.append(msgc_queue.get())
    logging.debug("got {} worker in total".format(len(worker_ip_list)))
    for i in range(len(worker_ip_list)):
        task_queue.put("STOP", block=True)
    task_queue.join()
    time.sleep(10)
    manager.shutdown()
    logging.debug("sent {} done message".format(len(worker_ip_list)))
    print('master exit.')


def check_gatk(chrom, pos_hg19, ref, alt, gatk_file, tolerance=100):
    gatk_pos = ''
    gatk_ref = ''
    gatk_alt = ''
    in_gatk = 0
    if len(ref) == 1 and len(alt) == 1:
        pp = Popen(["grep -P '{0}\\t{1}\\t' {2}".format(chrom, pos_hg19, gatk_file)],
                   shell=True, stdout=PIPE)
        pipe_ret = pp.stdout.readlines()
        if not pipe_ret:
            reason = 'not in corresponding gatk file'
        else:
            in_gatk = 1
            reason = ''
            gatk_list = pipe_ret[0].strip().split('\t')
            gatk_pos = gatk_list[1]
            gatk_ref = gatk_list[3]
            gatk_alt = gatk_list[4]
        return [in_gatk, reason, gatk_pos, gatk_ref, gatk_alt, pipe_ret[0] if pipe_ret else '']
    # indel
    cmd = 'awk -v chr="{0}" -v pos="{1}"  \'{{if($1==chr && $2<pos+{3} && $2>pos-{3} && length($4)!=length($5)){{print $0}}}}\' {2}' \
          ''.format(chrom, pos_hg19, gatk_file, tolerance)
    print(cmd)
    sys.stdout.flush()
    pp = Popen([cmd], shell=True, stdout=PIPE)
    pipe_ret = pp.stdout.readlines()
    if not pipe_ret:
        reason = 'not in corresponding gatk file'
    else:
        in_gatk = 1
        reason = ''
        gatk_list = pipe_ret[0].strip().split('\t')
        gatk_pos = gatk_list[1]
        gatk_ref = gatk_list[3]
        gatk_alt = gatk_list[4]
    return [in_gatk, reason, gatk_pos, gatk_ref, gatk_alt, pipe_ret[0] if pipe_ret else '']


def handle_task(variance_data, converter, sample_col_index_list,
                col_index2gen_id_dict, genid2slid_dict, gatk_path,
                fp, fp_org):
    print('new job {} variants'.format(len(variance_data)))
    icounter = 1
    for variance_line in variance_data:
        print('handling {} variant'.format(icounter))
        sys.stdout.flush()
        icounter += 1
        chrom = variance_line[1]
        pos_hg38 = variance_line[2]
        hg19 = converter[chrom][int(pos_hg38)]
        ref = variance_line[3]
        alt = variance_line[4]
        if len(ref) == len(alt):  # todo 只分析indel
            continue
        if len(hg19) == 0:
            pos_hg19 = ''
        else:
            pos_hg19 = hg19[0][1]

        for col_index in sample_col_index_list:
            genotype = variance_line[col_index]
            if genotype and genotype != 'na' and int(genotype) > 0:
                genid = col_index2gen_id_dict[col_index]
                slid = genid2slid_dict[genid]
                gatk_file = pack_gatk_filename(gatk_path, slid)
                gatk_pos = ''
                gatk_ref = ''
                gatk_alt = ''
                in_gatk = 0
                if pos_hg19 == '':
                    reason = 'can not liftover to hg19'
                elif not os.path.exists(gatk_file):
                    reason = 'no corresponding gatk file'
                else:
                    in_gatk, reason, gatk_pos, gatk_ref, gatk_alt, gatk_data = check_gatk(chrom, pos_hg19, ref,
                                                                                          alt, gatk_file, 50)

                out_str = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\n" \
                          "".format(chrom, pos_hg38, pos_hg19, ref, alt, genid,
                                    slid, in_gatk, reason, gatk_pos, gatk_ref, gatk_alt,
                                    "\t".join([str(i) for i in variance_line[:sample_col_index_list[0]]]),
                                    "\t".join([str(i) for i in variance_line[sample_col_index_list[-1] + 1:]]))
                fp.write(out_str)
                logging.debug('write: {}'.format(out_str))
                if in_gatk == 1:
                    fp_org.write('{0}\t{1}\t{2}'.format(os.path.basename(gatk_file, ), genid, gatk_data))


def check_PE_GATK_worker(db_file, table, slid2genid_file, gatk_path, output_path, worker_id, server_addr):
    output_file = os.path.join(output_path, worker_id) + '.subret'
    outorgdata_file = os.path.join(output_path, worker_id) + '.orgdata'
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()
    cursor.execute('PRAGMA table_info([{}])'.format(table))
    table_info = cursor.fetchall()
    sample_col_index_list = [i[0] for i in table_info if i[1].startswith('sample_')]
    # print('sample_col_index_list={}'.format(sample_col_index_list))
    col_index2gen_id_dict = dict([[i[0], i[1][7:]] for i in table_info if i[1].startswith('sample_')])
    genid2slid_dict = load_genid2slid(slid2genid_file)

    converter = get_lifter('hg38', 'hg19')

    BaseManager.register("get_task_queue")
    BaseManager.register("get_result_queue")
    BaseManager.register('get_msgc_queue')
    print("worker{0} Connect to server {1}...".format(worker_id, server_addr))
    manager = BaseManager(address=(server_addr, 5000), authkey=b"wgsawgsa")
    # wait for server
    icounter = 0
    while 1:
        try:
            manager.connect()
        except:
            time.sleep(1)
            icounter += 1
            print("waiting for server {}s".format(icounter))
            continue
        break
    print("server {} Connected".format(server_addr))
    # get Queue object
    icounter = 0
    while True:
        try:
            task_queue = manager.get_task_queue()
            result_queue = manager.get_result_queue()
            msgc_queue = manager.get_msgc_queue()
            icounter += 1
        except:
            print("try to get queue object {} times".format(icounter))
            continue
        break
    msgc_queue.put("{}".format(get_ip()), block=False)
    with open(output_file, 'w') as fp, open(outorgdata_file, 'w') as fp_org:
        while 1:
            try:
                msg_task = task_queue.get(block=True)
                if msg_task == "STOP":
                    task_queue.task_done()
                    print("I'm done.")
                    break
                # handle the task here
                variance_data = msg_task
                handle_task(variance_data, converter, sample_col_index_list,
                            col_index2gen_id_dict, genid2slid_dict, gatk_path,
                            fp, fp_org)
                task_queue.task_done()
            except Queue.Empty:
                print("task queue is empty.")


def chek_PE_GATK_multiple(num, db_file, table, slid2genid_file, gatk_path, output_path):
    logging.basicConfig(filename="{}.log".format(sys._getframe().f_code.co_name),
                        level=logging.DEBUG, format=log_format, filemode="w")
    if float(sys.version[:3]) != 2.7:
        print("CRITICAL: Python version must be 2.7!\n")
        sys.exit(1)
    curr_dir = os.getcwd()
    server_addr = get_ip()
    print("master ip = {}".format(server_addr))
    with open(curr_dir + "/tmp.sh", "w") as fp:
        fp.write("""#!/bin/bash
#$ -j y
#$ -S /bin/bash
#$ -cwd
#$ -N worker
#$ -pe smp 1
#$ -l h_vmem=16G
#$ -m bes
#$ -q all.q
export LC_ALL=C
export MALLOC_ARENA_MAX=4
export PATH=$PATH:$HOME/wyj/.code
wgsa.py check_PE_GATK_worker {0} {5} {1} {2} {3} $1 {4}
""".format(db_file, slid2genid_file, gatk_path, output_path, server_addr, table))
        fp.flush()
        fp.close()

    while True:
        if os.access(curr_dir + "/tmp.sh", os.R_OK):
            print("{} read OK".format(curr_dir + "/tmp.sh"))
            break
        time.sleep(1)
    for i in range(int(num)):
        cmd_str = "qsub {0}/tmp.sh {1}".format(curr_dir, i + 1)
        logging.debug(cmd_str)
        Popen([cmd_str], shell=True).wait()

    check_PE_GATK_master(db_file, table)
    Popen(["rm {0}/tmp.sh".format(curr_dir)], shell=True).wait()
    Popen(["cat {} > overlap.ret".format(os.path.join(output_path, "*\.subret"))], shell=True).wait()
    Popen(["cat {} > found_gatk".format(os.path.join(output_path, "*\.orgdata"))], shell=True).wait()


def pack_gatk_filename(gatk_path, slid):
    return os.path.join(gatk_path, slid) + '.GATK-4.beta.6.vcf'


def load_genid2slid(slid2genid_file):
    with open(slid2genid_file, 'r') as fp:
        return dict([i.strip().split('\t')[::-1] for i in fp.readlines() if len(i.strip()) > 0])


def check_PE_GATK_overlap(slid2genid_file, gatk_path, db_file, output):
    logging.basicConfig(filename="check_PE_GATK_overlap.log", level=logging.DEBUG, format=log_format, filemode="w")
    genid2slid_dict = load_genid2slid(slid2genid_file)
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()
    cursor.execute('PRAGMA table_info([variance])')
    table_info = cursor.fetchall()
    sample_col_index_list = [i[0] for i in table_info if i[1].startswith('sample_')]
    col_index2gen_id_dict = dict([[i[0], i[1][7:]] for i in table_info if i[1].startswith('sample_')])
    cursor.execute('select * from variance')
    variance_data = cursor.fetchall()
    converter = get_lifter('hg38', 'hg19')
    start_time = time.time()
    icounter = 0
    with open(output, 'w') as fp:
        fp.write('chr\thg38_pos\thg19_pos\tref\talt\tgenid\tslid\tin_gatk\treason\tgatk_pos\tgatk_ref\tgatk_alt\n')
        for variance_line in variance_data:
            chrom = variance_line[1]
            pos_hg38 = variance_line[2]
            hg19 = converter[chrom][int(pos_hg38)]
            ref = variance_line[3]
            alt = variance_line[4]
            curr_time = time.time()
            icounter += 1
            average_time = (curr_time - start_time) / icounter
            est_time = (len(variance_data) - icounter) * average_time
            logging.debug('handling {0}_{1}_{2}_{3}_{4} average time {5:.2}s    est:{6:.2}s'
                          ''.format(variance_line[0], chrom, pos_hg38, ref, alt, average_time, est_time))
            if len(hg19) == 0:
                pos_hg19 = ''
            else:
                pos_hg19 = hg19[0][1]
            for col_index in sample_col_index_list:
                genotype = variance_line[col_index]
                if genotype != 'na' and int(genotype) > 0:
                    genid = col_index2gen_id_dict[col_index]
                    slid = genid2slid_dict[genid]
                    gatk_file = pack_gatk_filename(gatk_path, slid)
                    gatk_pos = ''
                    gatk_ref = ''
                    gatk_alt = ''
                    in_gatk = 0
                    if pos_hg19 == '':
                        reason = 'can not liftover to hg19'
                    elif not os.path.exists(gatk_file):
                        reason = 'no corresponding gatk file'
                    else:
                        pp = Popen(["grep -P '{0}\\t{1}\\t' {2}".format(chrom, pos_hg19, gatk_file)],
                                   shell=True, stdout=PIPE)
                        pipe_ret = pp.stdout.readlines()
                        if not pipe_ret:
                            reason = 'not in corresponding gatk file'
                        else:
                            in_gatk = 1
                            reason = ''
                            gatk_list = pipe_ret[0].strip().split('\t')
                            gatk_pos = gatk_list[1]
                            gatk_ref = gatk_list[3]
                            gatk_alt = gatk_list[4]
                    out_str = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\n" \
                              "".format(chrom, pos_hg38, pos_hg19, ref, alt, genid,
                                        slid, in_gatk, reason, gatk_pos, gatk_ref, gatk_alt)
                    fp.write(out_str)
                    logging.debug('write: {}'.format(out_str))
        logging.debug('all done')


def person_split_PE_variance(PE_variance, output_path, id_pairs):
    genid2slid_dict = load_genid2slid(id_pairs)
    with open(PE_variance, 'r') as fp:
        while True:
            data_line = fp.readline()
            if data_line.startswith('##'):
                continue
            if data_line.startswith('#'):
                head_list = data_line.strip().split("\t")
                break

    for i in range(len(head_list)):
        target_col = i + 1
        if target_col < 10:
            continue
        slid = genid2slid_dict[head_list[i]]
        with open('tmp.sh', 'w') as fp:
            fp.write("""
#$ -j y
#$ -S /bin/bash
#$ -cwd
#$ -N cut
#$ -pe smp 1
#$ -l h_vmem=32G
#$ -m bes
#$ -q all.q
export LC_ALL=C
export MALLOC_ARENA_MAX=4
echo start
echo 'cut -f 1,2,3,4,5,6,7,8,9,{0} {1} | grep -v '0/0' | grep -v '\./\.' > {2}'
cut -f 1,2,3,4,5,6,7,8,9,{0} {1} | grep -v '0/0' | grep -v '\./\.' > {2}
echo end
        """.format(target_col,
                   PE_variance,
                   os.path.join(output_path,
                                os.path.basename(PE_variance))[
                   :-3] + slid + '.vcf'))  # | wgsa.py split_vcf_with_sample | wgsa.py left_normalization
        # cmd_str = "cut -f 1,2,3,4,5,6,7,8,9,{0} {1} | grep -v '0/0' | grep -v '\./\.' > {2}" \
        #           "".format(target_col, PE_variance,
        #                     os.path.join(output_path,
        #                                  os.path.basename(PE_variance))[:-3] + slid + '.vcf')
        Popen(['qsub tmp.sh'], shell=True).wait()


def vcf_reduce_PE_alt():
    target_chr_set = set(['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14',
                          '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y'])
    sys.stdout.write("##fileformat=VCFv4.0\n"
                     "##fileDate=201737\n"
                     "##reference=/gs/gsfs0/users/yizhao/wyj/ref/human_g1k_v37_decoy.fasta\n"
                     "##phasing=none\n"
                     "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">\n"
                     "##FORMAT=<ID=GQ,Number=1,Type=Float,Description=\"Genotype Quality\">\n"
                     "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
                     "##contig=<ID=1,length=249250621>\n"
                     "##contig=<ID=2,length=243199373>\n"
                     "##contig=<ID=3,length=198022430>\n"
                     "##contig=<ID=4,length=191154276>\n"
                     "##contig=<ID=5,length=180915260>\n"
                     "##contig=<ID=6,length=171115067>\n"
                     "##contig=<ID=7,length=159138663>\n"
                     "##contig=<ID=8,length=146364022>\n"
                     "##contig=<ID=9,length=141213431>\n"
                     "##contig=<ID=10,length=135534747>\n"
                     "##contig=<ID=11,length=135006516>\n"
                     "##contig=<ID=12,length=133851895>\n"
                     "##contig=<ID=13,length=115169878>\n"
                     "##contig=<ID=14,length=107349540>\n"
                     "##contig=<ID=15,length=102531392>\n"
                     "##contig=<ID=16,length=90354753>\n"
                     "##contig=<ID=17,length=81195210>\n"
                     "##contig=<ID=18,length=78077248>\n"
                     "##contig=<ID=19,length=59128983>\n"
                     "##contig=<ID=20,length=63025520>\n"
                     "##contig=<ID=21,length=48129895>\n"
                     "##contig=<ID=22,length=51304566>\n"
                     "##contig=<ID=X,length=155270560>\n"
                     "##contig=<ID=Y,length=59373566>\n"
                     "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")
    while True:
        vcf_line = sys.stdin.readline()
        if not vcf_line:
            break
        if vcf_line.startswith("#"):
            # sys.stdout.write(vcf_line)
            continue
        vcf_line = vcf_line.strip()
        if len(vcf_line) == 0:
            continue
        # sys.stderr.write('handling {}\n'.format(vcf_line))
        vcf_list = vcf_line.split('\t')
        vcf_list[7] = 'NS=1'
        if vcf_list[0].startswith('chr'):
            vcf_list[0] = vcf_list[0][3:]
        if vcf_list[0] not in target_chr_set:
            continue
        if ',' in vcf_list[4]:
            sample_info_list = vcf_list[9].split(':')
            gene_type_set = set([int(i) for i in sample_info_list[0].split('/')])
            alt_list = vcf_list[4].split(',')
            genotype2alt_dict = dict(zip(range(1, len(alt_list) + 1, 1), alt_list))
            alt_list = []
            for i in genotype2alt_dict:
                if i in gene_type_set:
                    alt_list.append(genotype2alt_dict[i])
            if 0 in gene_type_set:
                sample_info_list[0] = '0/1'
            elif len(gene_type_set) == 1:
                sample_info_list[0] = '1/1'
            else:
                sample_info_list[0] = '1/2'
            vcf_list[9] = ':'.join(sample_info_list)
            vcf_list[4] = ','.join(alt_list)
        # sys.stderr.write('vcf_list={}'.format(vcf_list))
        ret_str = '\t'.join(vcf_list)

        ret_str = line_left_normalization(ret_str) + '\n'

        sys.stdout.write(ret_str)


def liftover_PE_variance2hg19(PE_variance, output, unmapped):
    converter = get_lifter('hg38', 'hg19')
    with open(PE_variance, 'r') as fp_hg38, open(output, 'w') as fp_out, open(unmapped, 'w') as fp_unmapped:
        while True:
            hg38_line = fp_hg38.readline()
            if not hg38_line:
                break
            if hg38_line.startswith('#'):
                fp_out.write(hg38_line)
                continue
            hg38_list = hg38_line.split('\t')
            if len(hg38_list) <= 2:
                continue
            try:
                hg19 = converter[hg38_list[0]][int(hg38_list[1])]
            except:
                print(traceback.format_exc())
                print('hg38_list = {}'.format(hg38_list))
                exit(0)
            if len(hg19):
                fp_out.write('{0}\t{1}\t{2}'.format(hg19[0][0], hg19[0][1], '\t'.join(hg38_list[2:])))
            else:
                fp_unmapped.write(hg38_line)


def analyze_PE_GATK_overlap(hg19_pe_vcf, hg19_gatk_vcf, outpath, mode):
    pe_only_out = os.path.join(outpath, 'pe_only')
    gatk_only_out = os.path.join(outpath, 'gatk_only')
    pe_in_gatk_out = os.path.join(outpath, 'pe_in_gatk')
    gatk_in_pe_out = os.path.join(outpath, 'gatk_in_pe')
    pe_dict = {}
    with open(hg19_pe_vcf, 'r') as fp:
        while True:
            data_line = fp.readline()
            if not data_line:
                break
            if data_line.startswith('#'):
                continue
            data_list = data_line.strip().split('\t')
            if len(data_list) <= 2:
                continue
            gt = data_list[9].split(':')[0]
            if ',' in data_list[4]:
                alt_list = data_list[4].split(',')
                alt_list.sort()
                data_list[4] = ','.join(alt_list)
            if mode == '0':
                key = '{0}_{1}_{2}_{3}'.format(data_list[0], data_list[1], data_list[3], data_list[4])
            else:
                key = '{0}_{1}_{2}_{3}_{4}'.format(data_list[0], data_list[1], data_list[3], data_list[4], gt)
            pe_dict[key] = data_line
    gatk_dict = {}
    with open(hg19_gatk_vcf, 'r') as fp:
        while True:
            data_line = fp.readline()
            if not data_line:
                break
            if data_line.startswith('#'):
                continue
            data_list = data_line.strip().split('\t')
            if len(data_list) <= 2:
                continue
            gt = data_list[9].split(':')[0]
            if ',' in data_list[4]:
                alt_list = data_list[4].split(',')
                alt_list.sort()
                data_list[4] = ','.join(alt_list)
            if mode == '0':
                key = '{0}_{1}_{2}_{3}'.format(data_list[0], data_list[1], data_list[3], data_list[4])
            else:
                key = '{0}_{1}_{2}_{3}_{4}'.format(data_list[0], data_list[1], data_list[3], data_list[4], gt)
            gatk_dict[key] = data_line
    pe_set = set(pe_dict.keys())
    gatk_set = set(gatk_dict.keys())
    with open(pe_only_out, 'w') as fp:
        for key in pe_set - gatk_set:
            fp.write(pe_dict[key])
    with open(gatk_only_out, 'w') as fp:
        for key in gatk_set - pe_set:
            fp.write(gatk_dict[key])
    with open(pe_in_gatk_out, 'w') as fp:
        for key in gatk_set & pe_set:
            fp.write("{0}\n".format(pe_dict[key].strip()))
    with open(gatk_in_pe_out, 'w') as fp:
        for key in gatk_set & pe_set:
            fp.write('{0}\n'.format(gatk_dict[key].strip()))


def cal_concordance(path, slid):
    gt_path = os.path.join(path, slid, 'consider_gt')
    none_gt_path = os.path.join(path, slid, 'not_consider_gt')

    def handle_path(path):
        def count_num(file):
            all_num = 0
            snp_num = 0
            indel_num = 0
            with open(file, 'r') as fp:
                while True:
                    data_line = fp.readline()
                    if not data_line:
                        break
                    if data_line.startswith('#'):
                        continue
                    data_list = data_line.strip().split('\t')
                    all_num += 1
                    if ',' not in data_list[4]:
                        if len(data_list[3]) != len(data_list[4]):
                            indel_num += 1
                        else:
                            snp_num += 1
                    else:
                        alt_list = data_list[4].split(',')
                        for alt in alt_list:
                            if len(alt) != len(data_list[3]):
                                indel_num += 1
                                break
                        else:
                            snp_num += 1
            return [all_num, snp_num, indel_num]

        if not os.path.exists(os.path.join(path, 'pe_only')):
            exit(0)
        pe_only_all, pe_only_snp, pe_only_indel = count_num(os.path.join(path, 'pe_only'))
        pe_in_gatk_all, pe_in_gatk_snp, pe_in_gatk_indel = count_num(os.path.join(path, 'pe_in_gatk'))
        gatk_only_all, gatk_only_snp, gatk_only_indel = count_num(os.path.join(path, 'gatk_only'))
        con_pe_all = pe_in_gatk_all / float(pe_only_all + pe_in_gatk_all)
        con_gatk_all = pe_in_gatk_all / float(gatk_only_all + pe_in_gatk_all)

        con_pe_snp = pe_in_gatk_snp / float(pe_only_snp + pe_in_gatk_snp)
        con_gatk_snp = pe_in_gatk_snp / float(gatk_only_snp + pe_in_gatk_snp)

        con_pe_indel = pe_in_gatk_indel / float(pe_only_indel + pe_in_gatk_indel)
        con_gatk_indel = pe_in_gatk_indel / float(gatk_only_indel + pe_in_gatk_indel)
        return [con_pe_all, con_gatk_all, con_pe_snp, con_gatk_snp, con_pe_indel, con_gatk_indel]

    con_pe_all, con_gatk_all, con_pe_snp, con_gatk_snp, con_pe_indel, con_gatk_indel = handle_path(gt_path)
    print('{0}\t1\tall\t{1}\t{2}\n'
          '{0}\t1\tsnp\t{3}\t{4}\n'
          '{0}\t1\tindel\t{5}\t{6}'.format(slid, con_pe_all, con_gatk_all,
                                           con_pe_snp, con_gatk_snp,
                                           con_pe_indel, con_gatk_indel))
    con_pe_all, con_gatk_all, con_pe_snp, con_gatk_snp, con_pe_indel, con_gatk_indel = handle_path(none_gt_path)
    print('{0}\t0\tall\t{1}\t{2}\n'
          '{0}\t0\tsnp\t{3}\t{4}\n'
          '{0}\t0\tindel\t{5}\t{6}'.format(slid, con_pe_all, con_gatk_all,
                                           con_pe_snp, con_gatk_snp,
                                           con_pe_indel, con_gatk_indel))


def gene_size_analyse(database, output_table, gene_list_file, inputlist):
    gene_list_dict = {}
    gene_list_name_list = []

    def load_gene_list(gene_list_file, gene_list_dict, gene_list_name_list):
        with open(gene_list_file, 'r') as fp:
            data = [i.strip().split('\t') for i in fp.readlines() if (not i.startswith('#')) and len(i.strip()) > 0]
        for line in data:
            gene_list_dict[line[0]] = set(line[1:])
            gene_list_name_list.append(line[0])

    def load_input_list(input_list):
        with open(input_list, 'r') as fp:
            return set([i for i in fp.readline().strip().split('\t')])

    input_set = load_input_list(inputlist)
    load_gene_list(gene_list_file, gene_list_dict, gene_list_name_list)
    conn = sqlite3.connect(database)
    cursor = conn.cursor()
    cursor.execute("select gene_id, gene_name, start_pos, end_pos, exon_number, total_exon_size  from gene_table")
    db_data = cursor.fetchall()
    with open(output_table, 'w') as fp:
        fp.write('#gene_id\tgene_name\tgene_size\texon_number\ttotal_exon_size')
        for gene_list_name in gene_list_name_list:
            fp.write('\t{}'.format(gene_list_name))
        fp.write('\n')
        for gene_id, gene_name, start_pos, end_pos, exon_number, total_exon_size in db_data:
            key = "{0}_{1}".format(gene_name, gene_id)
            # if not key in input_set:
            #     continue
            gene_size = abs(int(start_pos) - int(end_pos))
            output_str = "{0}\t{1}\t{2}\t{3}\t{4}".format(gene_id, gene_name, gene_size,
                                                          exon_number if exon_number else '',
                                                          total_exon_size if total_exon_size else '')

            for gene_list_name in gene_list_name_list:
                output_str += '\t1' if key in gene_list_dict[gene_list_name] else '\t0'
            fp.write(output_str + '\n')


def grep_human_gencode_gtf3(input_gtf3, output, db_file):
    def merge_regions(regions):
        regions.sort(key=lambda x: x[0])
        i = 1
        while i < len(regions):
            last_start = regions[i - 1][0]
            last_end = regions[i - 1][1]
            curr_start = regions[i][0]
            curr_end = regions[i][1]
            if last_end < curr_start:
                i += 1
                continue
            if curr_start <= last_end <= curr_end:
                regions[i - 1][1] = curr_end
                regions.pop(i)
                continue
            regions.pop(i)
            continue

    print("grep_human_gencode_gtf3 begin", file=sys.stderr)
    data_dict = {}
    print('#chr\tstart\tend\ttype\tsource\tgene_id\tgene_name\ttranscript_name\texon_number\texon_id')
    icounter = 0
    with open(input_gtf3, 'r') as fp:
        print("transforming data", file=sys.stderr)
        while True:
            line = fp.readline()
            if not line:
                break
            if line.startswith("#"):
                continue
            data_list = line.strip().split('\t')
            if data_list[2] != 'exon':
                continue
            info = data_list[8]
            chrom = data_list[0]
            source = data_list[1]
            type = data_list[2]
            start = int(data_list[3])
            end = int(data_list[4])
            gene_id = re.findall('gene_id=(\S+?);', info)[0]
            gene_id = gene_id.split('.')[0]
            gene_name = re.findall('gene_name=(\S+?);', info)[0]
            transcript_name = re.findall('transcript_name=(\S+?);', info)[0]
            exon_num = re.findall('exon_number=(\d+?);', info)[0]
            exon_id = re.findall('exon_id=(\S+?);', info)[0]
            print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}"
                  "".format(chrom, start, end, type, source,
                            gene_id, gene_name, transcript_name, exon_num, exon_id))
            key = "{0}\t{1}".format(gene_name, gene_id)
            if key not in data_dict:
                data_dict[key] = [[min(start, end), max(start, end)]]
            else:
                data_dict[key].append([min(start, end), max(start, end)])
            icounter += 1
            if icounter % 100 == 0:
                print('transforming data {}'.format(icounter), file=sys.stderr)

    # merge regions
    print("merging begin", file=sys.stderr)
    for i in data_dict:
        merge_regions(data_dict[i])
    with open(output, 'w') as fp:
        fp.write('#gene_name\tgene_id\texon_number\ttotal_distance\n')
        for i in data_dict:
            num = len(data_dict[i])
            dist = sum([end - start + 1 for start, end in data_dict[i]])
            fp.write("{0}\t{1}\t{2}\n".format(i, num, dist))
    print("collecting sql", file=sys.stderr)
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()
    db_add_col(cursor, 'gene_table', 'exon_number', 'varchr(20)')
    db_add_col(cursor, 'gene_table', 'total_exon_size', 'varchr(20)')
    cmd = 'select id, gene_name, gene_id from gene_table'
    cursor.execute(cmd)
    table_data = cursor.fetchall()
    sql_list = []
    for id, gene_name, gene_id in table_data:
        key = "{0}\t{1}".format(gene_name, gene_id)
        if key in data_dict:
            num = len(data_dict[key])
            dist = sum([end - start + 1 for start, end in data_dict[key]])
            sql_list.append("UPDATE gene_table SET exon_number = '{0}' WHERE id = {1}".format(num, id))
            sql_list.append("UPDATE gene_table SET total_exon_size = '{0}' WHERE id = {1}".format(dist, id))
    print("executing sql", file=sys.stderr)
    icounter = 0
    for sql in sql_list:
        if icounter % 100 == 0:
            print("executing sql {0}  /  {1}".format(icounter, len(sql_list)), file=sys.stderr)
        cursor.execute(sql)
        icounter += 1
    cursor.close()
    conn.commit()
    conn.close()
    print("all done", file=sys.stderr)


def db_build_PE_GATK_variants_table(db_file, overlap_result):
    def load_overlap_ret(overlap_result):
        ret = {}
        with open(overlap_result, 'r') as fp:
            while True:
                line = fp.readline()
                if not line:
                    break
                if not line.strip('\n'):
                    continue
                data_list = line.strip().split('\t')
                chrom, pos_hg38, pos_hg19, ref, alt, genid, slid, in_gatk, reason, gatk_pos, gatk_ref, gatk_alt = data_list[
                                                                                                                  :12]
                key = "{0}_{1}_{2}_{3}_{4}".format(chrom, pos_hg38, ref, alt, genid)
                if reason == 'not in corresponding gatk file':
                    ret[key] = False
                else:
                    ret[key] = True
        return ret

    new_table = 'PE_GATK_variant'
    print('loading overlap file...')
    sys.stdout.flush()
    overlap_dict = load_overlap_ret(overlap_result)
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()

    # create new table
    cursor.execute("DROP TABLE IF EXISTS {0}".format(new_table))
    cmd = """CREATE TABLE "PE_GATK_variant" ("id"	INTEGER NOT NULL UNIQUE, 	"chr"	varchar(20), 	"pos"	int, 	"ref"	varchr(40), 	"alt"	varchr(40), 	"is001"	int, 	"annovar"	varchar(1), 	"bystro"	varchar(1), 	"dmis"	varchar(1), 	"dsplicing"	varchar(1), 	"spidex"	varchar(1), 	"spliceAI"	varchar(1), 	"vep"	varchar(1), 	"spliceAI_anno"	varchr(500), 	"spliceAI_score"	varchr(20), 	"spidex_score"	varchr(20), 	"vcf_id"	varchr(40), 	"vcf_qual"	varchr(20), 	"vcf_filter"	varchr(20), 	"vcf_info"	varchr(20), 	"vcf_format"	varchr(20), 	"sample_1668"	varchr(2), 	"sample_1262"	varchr(2), 	"sample_901"	varchr(2), 	"sample_1229"	varchr(2), 	"sample_1787"	varchr(2), 	"sample_1859"	varchr(2), 	"sample_1065"	varchr(2), 	"sample_1777"	varchr(2), 	"sample_1855"	varchr(2), 	"sample_1661"	varchr(2), 	"sample_1794"	varchr(2), 	"sample_1075"	varchr(2), 	"sample_674"	varchr(2), 	"sample_1790"	varchr(2), 	"sample_1070"	varchr(2), 	"sample_676"	varchr(2), 	"sample_649"	varchr(2), 	"sample_1245"	varchr(2), 	"sample_1226"	varchr(2), 	"sample_1040"	varchr(2), 	"sample_1778"	varchr(2), 	"sample_1666"	varchr(2), 	"sample_50009"	varchr(2), 	"sample_641"	varchr(2), 	"sample_1244"	varchr(2), 	"sample_1242"	varchr(2), 	"sample_1049"	varchr(2), 	"sample_1856"	varchr(2), 	"sample_686"	varchr(2), 	"sample_1235"	varchr(2), 	"sample_900"	varchr(2), 	"sample_1663"	varchr(2), 	"sample_1064"	varchr(2), 	"sample_1797"	varchr(2), 	"sample_1191"	varchr(2), 	"sample_690"	varchr(2), 	"sample_1799"	varchr(2), 	"sample_1784"	varchr(2), 	"sample_1780"	varchr(2), 	"sample_1247"	varchr(2), 	"sample_1800"	varchr(2), 	"sample_1082"	varchr(2), 	"sample_50028"	varchr(2), 	"sample_1239"	varchr(2), 	"sample_1057"	varchr(2), 	"sample_1196"	varchr(2), 	"sample_661"	varchr(2), 	"sample_1250"	varchr(2), 	"sample_678"	varchr(2), 	"sample_1236"	varchr(2), 	"sample_680"	varchr(2), 	"sample_1665"	varchr(2), 	"sample_1861"	varchr(2), 	"sample_1234"	varchr(2), 	"sample_1232"	varchr(2), 	"sample_1241"	varchr(2), 	"sample_685"	varchr(2), 	"sample_1796"	varchr(2), 	"sample_1246"	varchr(2), 	"sample_1081"	varchr(2), 	"sample_50033"	varchr(2), 	"sample_1857"	varchr(2), 	"sample_1233"	varchr(2), 	"sample_1195"	varchr(2), 	"sample_1035"	varchr(2), 	"sample_1193"	varchr(2), 	"sample_1194"	varchr(2), 	"sample_1793"	varchr(2), 	"sample_1779"	varchr(2), 	"sample_1050"	varchr(2), 	"sample_642"	varchr(2), 	"sample_1243"	varchr(2), 	"sample_692"	varchr(2), 	"sample_1238"	varchr(2), 	"sample_1231"	varchr(2), 	"sample_1237"	varchr(2), 	"sample_1225"	varchr(2), 	"sample_671"	varchr(2), 	"sample_1420"	varchr(2), 	"sample_1083"	varchr(2), 	"sample_1802"	varchr(2), 	"sample_1255"	varchr(2), 	"sample_1192"	varchr(2), 	"sample_675"	varchr(2), 	"sample_693"	varchr(2), 	"sample_1664"	varchr(2), 	"sample_1078"	varchr(2), 	"sample_1786"	varchr(2), 	"sample_1228"	varchr(2), 	"sample_1881"	varchr(2), 	"sample_1074"	varchr(2), 	"sample_697"	varchr(2), 	"sample_1795"	varchr(2), 	"sample_544"	varchr(2), 	"sample_1412"	varchr(2), 	"sample_1436"	varchr(2), 	"sample_1416"	varchr(2), 	"sample_1430"	varchr(2), 	"sample_50038"	varchr(2), 	"sample_50022"	varchr(2), 	"sample_550"	varchr(2), 	"sample_565"	varchr(2), 	"sample_1435"	varchr(2), 	"sample_1425"	varchr(2), 	"sample_1424"	varchr(2), 	"sample_559"	varchr(2), 	"sample_1389"	varchr(2), 	"sample_1427"	varchr(2), 	"sample_564"	varchr(2), 	"sample_567"	varchr(2), 	"sample_1380"	varchr(2), 	"sample_1372"	varchr(2), 	"sample_1406"	varchr(2), 	"sample_548"	varchr(2), 	"sample_549"	varchr(2), 	"sample_1437"	varchr(2), 	"sample_569"	varchr(2), 	"sample_570"	varchr(2), 	"sample_1384"	varchr(2), 	"sample_1431"	varchr(2), 	"sample_552"	varchr(2), 	"sample_1399"	varchr(2), 	"sample_554"	varchr(2), 	"sample_1422"	varchr(2), 	"sample_551"	varchr(2), 	"sample_1421"	varchr(2), 	"sample_541"	varchr(2), 	"sample_540"	varchr(2), 	"sample_546"	varchr(2), 	"sample_1381"	varchr(2), 	"sample_1417"	varchr(2), 	"sample_547"	varchr(2), 	"sample_1379"	varchr(2), 	"sample_50032"	varchr(2), 	"sample_1391"	varchr(2), 	"sample_1433"	varchr(2), 	"sample_558"	varchr(2), 	"sample_1387"	varchr(2), 	"sample_1392"	varchr(2), 	"sample_1374"	varchr(2), 	"sample_1388"	varchr(2), 	"sample_543"	varchr(2), 	"sample_1409"	varchr(2), 	"sample_556"	varchr(2), 	"sample_568"	varchr(2), 	"sample_1393"	varchr(2), 	"sample_553"	varchr(2), 	"sample_1434"	varchr(2), 	"sample_1426"	varchr(2), 	"sample_1398"	varchr(2), 	"sample_1107"	varchr(2), 	"sample_1713"	varchr(2), 	"sample_578"	varchr(2), 	"sample_1095"	varchr(2), 	"sample_582"	varchr(2), 	"sample_464"	varchr(2), 	"sample_1726"	varchr(2), 	"sample_576"	varchr(2), 	"sample_1544"	varchr(2), 	"sample_1559"	varchr(2), 	"sample_1553"	varchr(2), 	"sample_577"	varchr(2), 	"sample_1550"	varchr(2), 	"sample_466"	varchr(2), 	"sample_591"	varchr(2), 	"sample_479"	varchr(2), 	"sample_1112"	varchr(2), 	"sample_1725"	varchr(2), 	"sample_1546"	varchr(2), 	"sample_1105"	varchr(2), 	"sample_586"	varchr(2), 	"sample_1540"	varchr(2), 	"sample_572"	varchr(2), 	"sample_574"	varchr(2), 	"sample_50012"	varchr(2), 	"sample_573"	varchr(2), 	"sample_1554"	varchr(2), 	"sample_1556"	varchr(2), 	"sample_593"	varchr(2), 	"sample_1541"	varchr(2), 	"sample_595"	varchr(2), 	"sample_1100"	varchr(2), 	"sample_1722"	varchr(2), 	"sample_1727"	varchr(2), 	"sample_598"	varchr(2), 	"sample_580"	varchr(2), 	"sample_585"	varchr(2), 	"sample_1563"	varchr(2), 	"sample_1707"	varchr(2), 	"sample_1715"	varchr(2), 	"sample_571"	varchr(2), 	"sample_1545"	varchr(2), 	"sample_1548"	varchr(2), 	"sample_1552"	varchr(2), 	"sample_470"	varchr(2), 	"sample_1110"	varchr(2), 	"sample_588"	varchr(2), 	"sample_1561"	varchr(2), 	"sample_1549"	varchr(2), 	"sample_1728"	varchr(2), 	"sample_1538"	varchr(2), 	"sample_1539"	varchr(2), 	"sample_1101"	varchr(2), 	"sample_1117"	varchr(2), 	"sample_590"	varchr(2), 	"sample_1542"	varchr(2), 	"sample_1723"	varchr(2), 	"sample_1547"	varchr(2), 	"sample_1096"	varchr(2), 	"sample_1716"	varchr(2), 	"sample_1097"	varchr(2), 	"sample_784"	varchr(2), 	"sample_1543"	varchr(2), 	"sample_596"	varchr(2), 	"sample_1492"	varchr(2), 	"sample_1868"	varchr(2), 	"sample_1497"	varchr(2), 	"sample_1477"	varchr(2), 	"sample_1459"	varchr(2), 	"sample_1453"	varchr(2), 	"sample_1491"	varchr(2), 	"sample_1140"	varchr(2), 	"sample_1478"	varchr(2), 	"sample_1443"	varchr(2), 	"sample_1479"	varchr(2), 	"sample_1517"	varchr(2), 	"sample_1513"	varchr(2), 	"sample_1125"	varchr(2), 	"sample_1500"	varchr(2), 	"sample_1148"	varchr(2), 	"sample_1637"	varchr(2), 	"sample_1496"	varchr(2), 	"sample_1498"	varchr(2), 	"sample_1509"	varchr(2), 	"sample_1467"	varchr(2), 	"sample_1865"	varchr(2), 	"sample_1456"	varchr(2), 	"sample_1499"	varchr(2), 	"sample_1130"	varchr(2), 	"sample_1510"	varchr(2), 	"sample_1462"	varchr(2), 	"sample_1142"	varchr(2), 	"sample_1515"	varchr(2), 	"sample_1482"	varchr(2), 	"sample_1636"	varchr(2), 	"sample_1121"	varchr(2), 	"sample_1446"	varchr(2), 	"sample_1493"	varchr(2), 	"sample_1866"	varchr(2), 	"sample_1506"	varchr(2), 	"sample_1146"	varchr(2), 	"sample_1458"	varchr(2), 	"sample_1638"	varchr(2), 	"sample_1138"	varchr(2), 	"sample_533"	varchr(2), 	"sample_1454"	varchr(2), 	"sample_1133"	varchr(2), 	"sample_1476"	varchr(2), 	"sample_1122"	varchr(2), 	"sample_1505"	varchr(2), 	"sample_1448"	varchr(2), 	"sample_1495"	varchr(2), 	"sample_535"	varchr(2), 	"sample_1123"	varchr(2), 	"sample_1135"	varchr(2), 	"sample_1440"	varchr(2), 	"sample_1127"	varchr(2), 	"sample_1514"	varchr(2), 	"sample_1472"	varchr(2), 	"sample_1507"	varchr(2), 	"sample_1447"	varchr(2), 	"sample_1502"	varchr(2), 	"sample_1867"	varchr(2), 	"sample_532"	varchr(2), 	"sample_1134"	varchr(2), 	"sample_1461"	varchr(2), 	"sample_1471"	varchr(2), 	"sample_1501"	varchr(2), 	"sample_1147"	varchr(2), 	"sample_1444"	varchr(2), 	"sample_1141"	varchr(2), 	"sample_1480"	varchr(2), 	"sample_1504"	varchr(2), 	"sample_1129"	varchr(2), 	"sample_1445"	varchr(2), 	"sample_1450"	varchr(2), 	"sample_1451"	varchr(2), 	"sample_1463"	varchr(2), 	"sample_1126"	varchr(2), 	"sample_1449"	varchr(2), 	"sample_1508"	varchr(2), 	"sample_1457"	varchr(2), 	"sample_1442"	varchr(2), 	"sample_1474"	varchr(2), 	"sample_1464"	varchr(2), 	"sample_1512"	varchr(2), 	"sample_1143"	varchr(2), 	"sample_1494"	varchr(2), 	"sample_1460"	varchr(2), 	"sample_1465"	varchr(2), 	"sample_1137"	varchr(2), 	"sample_1441"	varchr(2), 	"sample_1131"	varchr(2), 	"sample_1475"	varchr(2), 	"sample_687"	varchr(2), 	"sample_689"	varchr(2), 	"sample_879"	varchr(2), 	"sample_537"	varchr(2), 	"sample_1614"	varchr(2), 	"sample_1634"	varchr(2), 	"sample_1605"	varchr(2), 	"sample_1598"	varchr(2), 	"sample_1602"	varchr(2), 	"sample_682"	varchr(2), 	"sample_882"	varchr(2), 	"sample_876"	varchr(2), 	"sample_1613"	varchr(2), 	"sample_1585"	varchr(2), 	"sample_1611"	varchr(2), 	"sample_1594"	varchr(2), 	"sample_1587"	varchr(2), 	"sample_897"	varchr(2), 	"sample_1581"	varchr(2), 	"sample_1635"	varchr(2), 	"sample_688"	varchr(2), 	"sample_1568"	varchr(2), 	"sample_1625"	varchr(2), 	"sample_1593"	varchr(2), 	"sample_691"	varchr(2), 	"sample_1592"	varchr(2), 	"sample_895"	varchr(2), 	"sample_1615"	varchr(2), 	"sample_1768"	varchr(2), 	"sample_1586"	varchr(2), 	"sample_1628"	varchr(2), 	"sample_1618"	varchr(2), 	"sample_1567"	varchr(2), 	"sample_1607"	varchr(2), 	"sample_1591"	varchr(2), 	"sample_896"	varchr(2), 	"sample_1705"	varchr(2), 	"sample_669"	varchr(2), 	"sample_1600"	varchr(2), 	"sample_1565"	varchr(2), 	"sample_1580"	varchr(2), 	"sample_874"	varchr(2), 	"sample_881"	varchr(2), 	"sample_890"	varchr(2), 	"sample_1633"	varchr(2), 	"sample_1573"	varchr(2), 	"sample_1570"	varchr(2), 	"sample_870"	varchr(2), 	"sample_869"	varchr(2), 	"sample_884"	varchr(2), 	"sample_1703"	varchr(2), 	"sample_1616"	varchr(2), 	"sample_1767"	varchr(2), 	"sample_1578"	varchr(2), 	"sample_1590"	varchr(2), 	"sample_1617"	varchr(2), 	"sample_1577"	varchr(2), 	"sample_893"	varchr(2), 	"sample_677"	varchr(2), 	"sample_888"	varchr(2), 	"sample_1610"	varchr(2), 	"sample_1596"	varchr(2), 	"sample_891"	varchr(2), 	"sample_1595"	varchr(2), 	"sample_894"	varchr(2), 	"sample_1584"	varchr(2), 	"sample_684"	varchr(2), 	"sample_1627"	varchr(2), 	"sample_1579"	varchr(2), 	"sample_538"	varchr(2), 	"sample_672"	varchr(2), 	"sample_1572"	varchr(2), 	"sample_1608"	varchr(2), 	"sample_1766"	varchr(2), 	"sample_1597"	varchr(2), 	"sample_1606"	varchr(2), 	"sample_1149"	varchr(2), 	"sample_770"	varchr(2), 	"sample_1278"	varchr(2), 	"sample_1085"	varchr(2), 	"sample_1662"	varchr(2), 	"sample_1038"	varchr(2), 	"sample_1188"	varchr(2), 	"sample_1879"	varchr(2), 	"sample_1876"	varchr(2), 	"sample_1524"	varchr(2), 	"sample_1883"	varchr(2), 	"sample_1522"	varchr(2), 	"sample_1169"	varchr(2), 	"sample_1771"	varchr(2), 	"sample_1033"	varchr(2), 	"sample_1186"	varchr(2), 	"sample_1058"	varchr(2), 	"sample_1280"	varchr(2), 	"sample_1525"	varchr(2), 	"sample_1872"	varchr(2), 	"sample_1151"	varchr(2), 	"sample_1173"	varchr(2), 	"sample_1180"	varchr(2), 	"sample_1671"	varchr(2), 	"sample_1281"	varchr(2), 	"sample_1185"	varchr(2), 	"sample_1672"	varchr(2), 	"sample_1521"	varchr(2), 	"sample_1268"	varchr(2), 	"sample_1154"	varchr(2), 	"sample_1878"	varchr(2), 	"sample_1885"	varchr(2), 	"sample_1166"	varchr(2), 	"sample_1120"	varchr(2), 	"sample_1184"	varchr(2), 	"sample_1284"	varchr(2), 	"sample_1882"	varchr(2), 	"sample_1152"	varchr(2), 	"sample_1523"	varchr(2), 	"sample_1659"	varchr(2), 	"sample_1090"	varchr(2), 	"sample_767"	varchr(2), 	"sample_1150"	varchr(2), 	"sample_1052"	varchr(2), 	"sample_1177"	varchr(2), 	"sample_1094"	varchr(2), 	"sample_1886"	varchr(2), 	"sample_50034"	varchr(2), 	"sample_1089"	varchr(2), 	"sample_1267"	varchr(2), 	"sample_1273"	varchr(2), 	"sample_1670"	varchr(2), 	"sample_1656"	varchr(2), 	"sample_1174"	varchr(2), 	"sample_1773"	varchr(2), 	"sample_1772"	varchr(2), 	"sample_768"	varchr(2), 	"sample_1286"	varchr(2), 	"sample_1274"	varchr(2), 	"sample_1172"	varchr(2), 	"sample_1270"	varchr(2), 	"sample_1176"	varchr(2), 	"sample_1041"	varchr(2), 	"sample_1160"	varchr(2), 	"sample_1283"	varchr(2), 	"sample_1161"	varchr(2), 	"sample_1277"	varchr(2), 	"sample_1673"	varchr(2), 	"sample_1275"	varchr(2), 	"sample_772"	varchr(2), 	"sample_1093"	varchr(2), 	"sample_1162"	varchr(2), 	"sample_1875"	varchr(2), 	"sample_1873"	varchr(2), 	"sample_1880"	varchr(2), 	"sample_769"	varchr(2), 	"sample_1182"	varchr(2), 	"sample_637"	varchr(2), 	"sample_773"	varchr(2), 	"sample_1159"	varchr(2), 	"sample_1520"	varchr(2), 	"sample_1657"	varchr(2), 	"sample_1158"	varchr(2), 	"sample_1884"	varchr(2), 	"sample_1658"	varchr(2), 	"sample_1153"	varchr(2), 	"sample_1272"	varchr(2), 	"sample_1770"	varchr(2), 	"sample_1874"	varchr(2), 	"sample_621"	varchr(2), 	"sample_1674"	varchr(2), 	"sample_814"	varchr(2), 	"sample_1682"	varchr(2), 	"sample_626"	varchr(2), 	"sample_623"	varchr(2), 	"sample_1676"	varchr(2), 	"sample_1689"	varchr(2), 	"sample_1678"	varchr(2), 	"sample_1699"	varchr(2), 	"sample_602"	varchr(2), 	"sample_828"	varchr(2), 	"sample_1686"	varchr(2), 	"sample_836"	varchr(2), 	"sample_827"	varchr(2), 	"sample_859"	varchr(2), 	"sample_617"	varchr(2), 	"sample_854"	varchr(2), 	"sample_1687"	varchr(2), 	"sample_849"	varchr(2), 	"sample_50021"	varchr(2), 	"sample_830"	varchr(2), 	"sample_861"	varchr(2), 	"sample_1696"	varchr(2), 	"sample_810"	varchr(2), 	"sample_622"	varchr(2), 	"sample_813"	varchr(2), 	"sample_1693"	varchr(2), 	"sample_625"	varchr(2), 	"sample_610"	varchr(2), 	"sample_863"	varchr(2), 	"sample_816"	varchr(2), 	"sample_852"	varchr(2), 	"sample_50019"	varchr(2), 	"sample_818"	varchr(2), 	"sample_1700"	varchr(2), 	"sample_50006"	varchr(2), 	"sample_1675"	varchr(2), 	"sample_832"	varchr(2), 	"sample_841"	varchr(2), 	"sample_611"	varchr(2), 	"sample_629"	varchr(2), 	"sample_851"	varchr(2), 	"sample_1679"	varchr(2), 	"sample_838"	varchr(2), 	"sample_824"	varchr(2), 	"sample_1695"	varchr(2), 	"sample_826"	varchr(2), 	"sample_856"	varchr(2), 	"sample_605"	varchr(2), 	"sample_855"	varchr(2), 	"sample_1697"	varchr(2), 	"sample_1692"	varchr(2), 	"sample_618"	varchr(2), 	"sample_853"	varchr(2), 	"sample_1774"	varchr(2), 	"sample_865"	varchr(2), 	"sample_1684"	varchr(2), 	"sample_1683"	varchr(2), 	"sample_609"	varchr(2), 	"sample_843"	varchr(2), 	"sample_825"	varchr(2), 	"sample_50027"	varchr(2), 	"sample_1690"	varchr(2), 	"sample_603"	varchr(2), 	"sample_50005"	varchr(2), 	"sample_831"	varchr(2), 	"sample_857"	varchr(2), 	"sample_844"	varchr(2), 	"sample_848"	varchr(2), 	"sample_862"	varchr(2), 	"sample_1688"	varchr(2), 	"sample_1698"	varchr(2), 	"sample_627"	varchr(2), 	"sample_821"	varchr(2), 	"sample_787"	varchr(2), 	"sample_716"	varchr(2), 	"sample_728"	varchr(2), 	"sample_743"	varchr(2), 	"sample_791"	varchr(2), 	"sample_740"	varchr(2), 	"sample_708"	varchr(2), 	"sample_706"	varchr(2), 	"sample_715"	varchr(2), 	"sample_50024"	varchr(2), 	"sample_744"	varchr(2), 	"sample_796"	varchr(2), 	"sample_738"	varchr(2), 	"sample_781"	varchr(2), 	"sample_711"	varchr(2), 	"sample_717"	varchr(2), 	"sample_733"	varchr(2), 	"sample_1224"	varchr(2), 	"sample_783"	varchr(2), 	"sample_632"	varchr(2), 	"sample_1221"	varchr(2), 	"sample_777"	varchr(2), 	"sample_729"	varchr(2), 	"sample_698"	varchr(2), 	"sample_1211"	varchr(2), 	"sample_713"	varchr(2), 	"sample_737"	varchr(2), 	"sample_725"	varchr(2), 	"sample_1205"	varchr(2), 	"sample_797"	varchr(2), 	"sample_727"	varchr(2), 	"sample_714"	varchr(2), 	"sample_709"	varchr(2), 	"sample_723"	varchr(2), 	"sample_1207"	varchr(2), 	"sample_704"	varchr(2), 	"sample_1210"	varchr(2), 	"sample_1208"	varchr(2), 	"sample_1198"	varchr(2), 	"sample_630"	varchr(2), 	"sample_795"	varchr(2), 	"sample_734"	varchr(2), 	"sample_1059"	varchr(2), 	"sample_700"	varchr(2), 	"sample_705"	varchr(2), 	"sample_1206"	varchr(2), 	"sample_779"	varchr(2), 	"sample_1051"	varchr(2), 	"sample_703"	varchr(2), 	"sample_726"	varchr(2), 	"sample_1213"	varchr(2), 	"sample_50023"	varchr(2), 	"sample_786"	varchr(2), 	"sample_1222"	varchr(2), 	"sample_1215"	varchr(2), 	"sample_736"	varchr(2), 	"sample_778"	varchr(2), 	"sample_721"	varchr(2), 	"sample_785"	varchr(2), 	"sample_730"	varchr(2), 	"sample_699"	varchr(2), 	"sample_776"	varchr(2), 	"sample_712"	varchr(2), 	"sample_793"	varchr(2), 	"sample_1204"	varchr(2), 	"sample_782"	varchr(2), 	"sample_799"	varchr(2), 	"sample_724"	varchr(2), 	"sample_731"	varchr(2), 	"sample_735"	varchr(2), 	"sample_741"	varchr(2), 	"sample_631"	varchr(2), 	"sample_395"	varchr(2), 	"sample_400"	varchr(2), 	"sample_154"	varchr(2), 	"sample_155"	varchr(2), 	"sample_156"	varchr(2), 	"sample_157"	varchr(2), 	"sample_158"	varchr(2), 	"sample_159"	varchr(2), 	"sample_160"	varchr(2), 	"sample_452"	varchr(2), 	"sample_161"	varchr(2), 	"sample_162"	varchr(2), 	"sample_163"	varchr(2), 	"sample_164"	varchr(2), 	"sample_165"	varchr(2), 	"sample_168"	varchr(2), 	"sample_170"	varchr(2), 	"sample_171"	varchr(2), 	"sample_172"	varchr(2), 	"sample_123"	varchr(2), 	"sample_125"	varchr(2), 	"sample_129"	varchr(2), 	"sample_131"	varchr(2), 	"sample_137"	varchr(2), 	"sample_140"	varchr(2), 	"sample_141"	varchr(2), 	"sample_143"	varchr(2), 	"sample_444"	varchr(2), 	"sample_138"	varchr(2), 	"sample_124"	varchr(2), 	"sample_134"	varchr(2), 	"sample_132"	varchr(2), 	"sample_133"	varchr(2), 	"sample_128"	varchr(2), 	"sample_87"	varchr(2), 	"sample_406"	varchr(2), 	"sample_50004"	varchr(2), 	"sample_432"	varchr(2), 	"sample_433"	varchr(2), 	"sample_436"	varchr(2), 	"sample_151"	varchr(2), 	"sample_150"	varchr(2), 	"sample_50002"	varchr(2), 	"sample_152"	varchr(2), 	"sample_449"	varchr(2), 	"sample_146"	varchr(2), 	"sample_42"	varchr(2), 	"sample_38"	varchr(2), 	"sample_314"	varchr(2), 	"sample_39"	varchr(2), 	"sample_46"	varchr(2), 	"sample_43"	varchr(2), 	"sample_44"	varchr(2), 	"sample_316"	varchr(2), 	"sample_27"	varchr(2), 	"sample_76"	varchr(2), 	"sample_77"	varchr(2), 	"sample_79"	varchr(2), 	"sample_80"	varchr(2), 	"sample_81"	varchr(2), 	"sample_85"	varchr(2), 	"sample_403"	varchr(2), 	"sample_98"	varchr(2), 	"sample_100"	varchr(2), 	"sample_108"	varchr(2), 	"sample_407"	varchr(2), 	"sample_118"	varchr(2), 	"sample_119"	varchr(2), 	"sample_82"	varchr(2), 	"sample_84"	varchr(2), 	"sample_88"	varchr(2), 	"sample_86"	varchr(2), 	"sample_90"	varchr(2), 	"sample_115"	varchr(2), 	"sample_121"	varchr(2), 	"sample_440"	varchr(2), 	"sample_144"	varchr(2), 	"sample_145"	varchr(2), 	"sample_441"	varchr(2), 	"sample_122"	varchr(2), 	"sample_136"	varchr(2), 	"sample_139"	varchr(2), 	"sample_442"	varchr(2), 	"sample_279"	varchr(2), 	"sample_10"	varchr(2), 	"sample_2"	varchr(2), 	"sample_317"	varchr(2), 	"sample_312"	varchr(2), 	"sample_36"	varchr(2), 	"sample_315"	varchr(2), 	"sample_390"	varchr(2), 	"sample_321"	varchr(2), 	"sample_329"	varchr(2), 	"sample_339"	varchr(2), 	"sample_350"	varchr(2), 	"sample_384"	varchr(2), 	"sample_385"	varchr(2), 	"sample_322"	varchr(2), 	"sample_323"	varchr(2), 	"sample_324"	varchr(2), 	"sample_325"	varchr(2), 	"sample_326"	varchr(2), 	"sample_327"	varchr(2), 	"sample_328"	varchr(2), 	"sample_330"	varchr(2), 	"sample_331"	varchr(2), 	"sample_332"	varchr(2), 	"sample_333"	varchr(2), 	"sample_334"	varchr(2), 	"sample_335"	varchr(2), 	"sample_336"	varchr(2), 	"sample_337"	varchr(2), 	"sample_338"	varchr(2), 	"sample_340"	varchr(2), 	"sample_341"	varchr(2), 	"sample_342"	varchr(2), 	"sample_500"	varchr(2), 	"sample_501"	varchr(2), 	"sample_507"	varchr(2), 	"sample_513"	varchr(2), 	"sample_364"	varchr(2), 	"sample_367"	varchr(2), 	"sample_368"	varchr(2), 	"sample_371"	varchr(2), 	"sample_372"	varchr(2), 	"sample_374"	varchr(2), 	"sample_379"	varchr(2), 	"sample_380"	varchr(2), 	"sample_343"	varchr(2), 	"sample_344"	varchr(2), 	"sample_345"	varchr(2), 	"sample_346"	varchr(2), 	"sample_347"	varchr(2), 	"sample_348"	varchr(2), 	"sample_349"	varchr(2), 	"sample_351"	varchr(2), 	"sample_353"	varchr(2), 	"sample_354"	varchr(2), 	"sample_355"	varchr(2), 	"sample_356"	varchr(2), 	"sample_50020"	varchr(2), 	"sample_358"	varchr(2), 	"sample_359"	varchr(2), 	"sample_360"	varchr(2), 	"sample_362"	varchr(2), 	"sample_363"	varchr(2), 	"sample_365"	varchr(2), 	"sample_373"	varchr(2), 	"sample_377"	varchr(2), 	"sample_381"	varchr(2), 	"sample_383"	varchr(2), 	"sample_516"	varchr(2), 	"sample_50003"	varchr(2), 	"sample_526"	varchr(2), 	"sample_498"	varchr(2), 	"sample_499"	varchr(2), 	"sample_487"	varchr(2), 	"sample_127"	varchr(2), 	"sample_319"	varchr(2), 	"sample_18"	varchr(2), 	"sample_53"	varchr(2), 	"sample_51"	varchr(2), 	"sample_45"	varchr(2), 	"sample_21"	varchr(2), 	"sample_55"	varchr(2), 	"sample_29"	varchr(2), 	"sample_60"	varchr(2), 	"sample_72"	varchr(2), 	"sample_73"	varchr(2), 	"sample_74"	varchr(2), 	"sample_75"	varchr(2), 	"sample_404"	varchr(2), 	"sample_93"	varchr(2), 	"sample_12"	varchr(2), 	"sample_299"	varchr(2), 	"sample_300"	varchr(2), 	"sample_271"	varchr(2), 	"sample_272"	varchr(2), 	"sample_474"	varchr(2), 	"sample_475"	varchr(2), 	"sample_465"	varchr(2), 	"sample_467"	varchr(2), 	"sample_468"	varchr(2), 	"sample_484"	varchr(2), 	"sample_50007"	varchr(2), 	"sample_307"	varchr(2), 	"sample_40"	varchr(2), 	"sample_48"	varchr(2), 	"sample_20"	varchr(2), 	"sample_22"	varchr(2), 	"sample_28"	varchr(2), 	"sample_370"	varchr(2), 	"sample_23"	varchr(2), 	"sample_41"	varchr(2), 	"sample_54"	varchr(2), 	"sample_67"	varchr(2), 	"sample_68"	varchr(2), 	"sample_70"	varchr(2), 	"sample_71"	varchr(2), 	"sample_83"	varchr(2), 	"sample_114"	varchr(2), 	"sample_91"	varchr(2), 	"sample_113"	varchr(2), 	"sample_117"	varchr(2), 	"sample_310"	varchr(2), 	"sample_1"	varchr(2), 	"sample_11"	varchr(2), 	"sample_302"	varchr(2), 	"sample_303"	varchr(2), 	"sample_305"	varchr(2), 	"sample_308"	varchr(2), 	"sample_309"	varchr(2), 	"sample_297"	varchr(2), 	"sample_266"	varchr(2), 	"sample_268"	varchr(2), 	"sample_311"	varchr(2), 	"sample_483"	varchr(2), 	"sample_472"	varchr(2), 	"sample_485"	varchr(2), 	"sample_445"	varchr(2), 	"sample_180"	varchr(2), 	"sample_181"	varchr(2), 	"sample_182"	varchr(2), 	"sample_183"	varchr(2), 	"sample_184"	varchr(2), 	"sample_185"	varchr(2), 	"sample_454"	varchr(2), 	"sample_186"	varchr(2), 	"sample_187"	varchr(2), 	"sample_188"	varchr(2), 	"sample_189"	varchr(2), 	"sample_190"	varchr(2), 	"sample_191"	varchr(2), 	"sample_192"	varchr(2), 	"sample_193"	varchr(2), 	"sample_194"	varchr(2), 	"sample_195"	varchr(2), 	"sample_197"	varchr(2), 	"sample_199"	varchr(2), 	"sample_200"	varchr(2), 	"sample_201"	varchr(2), 	"sample_202"	varchr(2), 	"sample_203"	varchr(2), 	"sample_204"	varchr(2), 	"sample_205"	varchr(2), 	"sample_206"	varchr(2), 	"sample_207"	varchr(2), 	"sample_208"	varchr(2), 	"sample_209"	varchr(2), 	"sample_210"	varchr(2), 	"sample_211"	varchr(2), 	"sample_212"	varchr(2), 	"sample_213"	varchr(2), 	"sample_214"	varchr(2), 	"sample_215"	varchr(2), 	"sample_217"	varchr(2), 	"sample_219"	varchr(2), 	"sample_50001"	varchr(2), 	"sample_221"	varchr(2), 	"sample_222"	varchr(2), 	"sample_223"	varchr(2), 	"sample_224"	varchr(2), 	"sample_225"	varchr(2), 	"sample_226"	varchr(2), 	"sample_227"	varchr(2), 	"sample_228"	varchr(2), 	"sample_229"	varchr(2), 	"sample_230"	varchr(2), 	"sample_457"	varchr(2), 	"sample_232"	varchr(2), 	"sample_233"	varchr(2), 	"sample_234"	varchr(2), 	"sample_235"	varchr(2), 	"sample_236"	varchr(2), 	"sample_237"	varchr(2), 	"sample_238"	varchr(2), 	"sample_239"	varchr(2), 	"sample_240"	varchr(2), 	"sample_242"	varchr(2), 	"sample_243"	varchr(2), 	"sample_244"	varchr(2), 	"sample_246"	varchr(2), 	"sample_247"	varchr(2), 	"sample_248"	varchr(2), 	"sample_249"	varchr(2), 	"sample_250"	varchr(2), 	"sample_251"	varchr(2), 	"sample_252"	varchr(2), 	"sample_254"	varchr(2), 	"sample_173"	varchr(2), 	"sample_174"	varchr(2), 	"sample_175"	varchr(2), 	"sample_176"	varchr(2), 	"sample_178"	varchr(2), 	"sample_256"	varchr(2), 	"sample_50008"	varchr(2), 	"sample_258"	varchr(2), 	"sample_259"	varchr(2), 	"sample_260"	varchr(2), 	"sample_261"	varchr(2), 	"sample_262"	varchr(2), 	"sample_263"	varchr(2), 	"sample_460"	varchr(2), 	"sample_264"	varchr(2), 	"sample_461"	varchr(2), 	"sample_198"	varchr(2), 	"sample_147"	varchr(2), 	"sample_97"	varchr(2), 	"sample_109"	varchr(2), 	"sample_120"	varchr(2), 	"sample_167"	varchr(2), 	"sample_753"	varchr(2), 	"sample_126"	varchr(2), 	"sample_130"	varchr(2), 	"sample_757"	varchr(2), 	"sample_1526"	varchr(2), 	"sample_1534"	varchr(2), 	"sample_759"	varchr(2), 	"sample_764"	varchr(2), 	"sample_755"	varchr(2), 	"sample_529"	varchr(2), 	"sample_528"	varchr(2), 	"sample_531"	varchr(2), 	"sample_530"	varchr(2), 	"sample_439"	varchr(2), 	"sample_418"	varchr(2), 	"sample_420"	varchr(2), 	"sample_422"	varchr(2), 	"sample_409"	varchr(2), 	"sample_414"	varchr(2), 	"sample_415"	varchr(2), 	"sample_386"	varchr(2), 	"sample_1622"	varchr(2), 	"sample_1624"	varchr(2), 	"sample_8"	varchr(2), 	"sample_13"	varchr(2), 	"sample_774"	varchr(2), 	"sample_289"	varchr(2), 	"sample_290"	varchr(2), 	"sample_291"	varchr(2), 	"sample_292"	varchr(2), 	"sample_294"	varchr(2), 	"sample_296"	varchr(2), 	"sample_298"	varchr(2), 	"sample_301"	varchr(2), 	"sample_1640"	varchr(2), 	"sample_1641"	varchr(2), 	"sample_1642"	varchr(2), 	"sample_1643"	varchr(2), 	"sample_1644"	varchr(2), 	"sample_1645"	varchr(2), 	"sample_1646"	varchr(2), 	"sample_1647"	varchr(2), 	"sample_1648"	varchr(2), 	"sample_1649"	varchr(2), 	"sample_1651"	varchr(2), 	"sample_1652"	varchr(2), 	"sample_1654"	varchr(2), 	"sample_1655"	varchr(2), 	"sample_1805"	varchr(2), 	"sample_1806"	varchr(2), 	"sample_1807"	varchr(2), 	"sample_1811"	varchr(2), 	"sample_1813"	varchr(2), 	"sample_1816"	varchr(2), 	"sample_1819"	varchr(2), 	"sample_1820"	varchr(2), 	"sample_1823"	varchr(2), 	"sample_1827"	varchr(2), 	"sample_1828"	varchr(2), 	"sample_1830"	varchr(2), 	"sample_1834"	varchr(2), 	"sample_1836"	varchr(2), 	"sample_1842"	varchr(2), 	"sample_1849"	varchr(2), 	"sample_1852"	varchr(2), 	"sample_1853"	varchr(2), 	"sample_1858"	varchr(2), 	"sample_1860"	varchr(2), 	"sample_1862"	varchr(2), 	"sample_1864"	varchr(2), 	"sample_754"	varchr(2), 	"sample_756"	varchr(2), 	"sample_758"	varchr(2), 	"sample_761"	varchr(2), 	"sample_762"	varchr(2), 	"sample_1527"	varchr(2), 	"sample_1762"	varchr(2), 	"sample_1763"	varchr(2), 	"sample_1764"	varchr(2), 	"sample_1869"	varchr(2), 	"sample_1871"	varchr(2), 	"sample_898"	varchr(2), 	"sample_902"	varchr(2), 	"sample_906"	varchr(2), 	"sample_907"	varchr(2), 	"sample_909"	varchr(2), 	"sample_911"	varchr(2), 	"sample_912"	varchr(2), 	"sample_913"	varchr(2), 	"sample_915"	varchr(2), 	"sample_916"	varchr(2), 	"sample_928"	varchr(2), 	"sample_1824"	varchr(2), 	"sample_1535"	varchr(2), 	"sample_1537"	varchr(2), 	"sample_1758"	varchr(2), 	"sample_1759"	varchr(2), 	"sample_1760"	varchr(2), 	"sample_1532"	varchr(2), 	"sample_50030"	varchr(2), 	"sample_1325"	varchr(2), 	"sample_1364"	varchr(2), 	"sample_1356"	varchr(2), 	"sample_1848"	varchr(2), 	"sample_1340"	varchr(2), 	"sample_1336"	varchr(2), 	"sample_1300"	varchr(2), 	"sample_1313"	varchr(2), 	"sample_1321"	varchr(2), 	"sample_1318"	varchr(2), 	"sample_1309"	varchr(2), 	"sample_1341"	varchr(2), 	"sample_50040"	varchr(2), 	"sample_1333"	varchr(2), 	"sample_1331"	varchr(2), 	"sample_944"	varchr(2), 	"sample_1317"	varchr(2), 	"sample_994"	varchr(2), 	"sample_1002"	varchr(2), 	"sample_1297"	varchr(2), 	"sample_998"	varchr(2), 	"sample_1323"	varchr(2), 	"sample_1838"	varchr(2), 	"sample_997"	varchr(2), 	"sample_1293"	varchr(2), 	"sample_1006"	varchr(2), 	"sample_1012"	varchr(2), 	"sample_1334"	varchr(2), 	"sample_1294"	varchr(2), 	"sample_1345"	varchr(2), 	"sample_1303"	varchr(2), 	"sample_1010"	varchr(2), 	"sample_1298"	varchr(2), 	"sample_1306"	varchr(2), 	"sample_1304"	varchr(2), 	"sample_1302"	varchr(2), 	"sample_1305"	varchr(2), 	"sample_995"	varchr(2), 	"sample_1311"	varchr(2), 	"sample_1009"	varchr(2), 	"sample_1342"	varchr(2), 	"sample_1329"	varchr(2), 	"sample_1362"	varchr(2), 	"sample_936"	varchr(2), 	"sample_1326"	varchr(2), 	"sample_1332"	varchr(2), 	"sample_1338"	varchr(2), 	"sample_1291"	varchr(2), 	"sample_1295"	varchr(2), 	"sample_1296"	varchr(2), 	"sample_1344"	varchr(2), 	"sample_1361"	varchr(2), 	"sample_1324"	varchr(2), 	"sample_1327"	varchr(2), 	"sample_1347"	varchr(2), 	"sample_1343"	varchr(2), 	"sample_1363"	varchr(2), 	"sample_1021"	varchr(2), 	"sample_1020"	varchr(2), 	"sample_1330"	varchr(2), 	"sample_1320"	varchr(2), 	"sample_933"	varchr(2), 	"sample_1353"	varchr(2), 	"sample_1307"	varchr(2), 	"sample_934"	varchr(2), 	"sample_1354"	varchr(2), 	"sample_1022"	varchr(2), 	"sample_1843"	varchr(2), 	"sample_1322"	varchr(2), 	"sample_1029"	varchr(2), 	"sample_993"	varchr(2), 	"sample_50029"	varchr(2), 	"sample_1299"	varchr(2), 	"sample_1335"	varchr(2), 	"sample_1355"	varchr(2), 	"sample_1290"	varchr(2), 	"sample_1352"	varchr(2), 	"sample_1358"	varchr(2), 	"sample_965"	varchr(2), 	"sample_1001"	varchr(2), 	"sample_1014"	varchr(2), 	"sample_957"	varchr(2), 	"sample_989"	varchr(2), 	"sample_979"	varchr(2), 	"sample_1011"	varchr(2), 	"sample_953"	varchr(2), 	"sample_948"	varchr(2), 	"sample_926"	varchr(2), 	"sample_1017"	varchr(2), 	"sample_984"	varchr(2), 	"sample_50014"	varchr(2), 	"sample_50031"	varchr(2), 	"sample_50025"	varchr(2), 	"sample_976"	varchr(2), 	"sample_958"	varchr(2), 	"sample_992"	varchr(2), 	"sample_988"	varchr(2), 	"sample_1831"	varchr(2), 	"sample_1025"	varchr(2), 	"sample_991"	varchr(2), 	"sample_980"	varchr(2), 	"sample_1026"	varchr(2), 	"sample_977"	varchr(2), 	"sample_960"	varchr(2), 	"sample_1031"	varchr(2), 	"sample_983"	varchr(2), 	"sample_963"	varchr(2), 	"sample_1829"	varchr(2), 	"sample_1013"	varchr(2), 	"sample_922"	varchr(2), 	"sample_925"	varchr(2), 	"sample_966"	varchr(2), 	"sample_50037"	varchr(2), 	"sample_949"	varchr(2), 	"sample_1027"	varchr(2), 	"sample_908"	varchr(2), 	"sample_962"	varchr(2), 	"sample_964"	varchr(2), 	"sample_919"	varchr(2), 	"sample_1008"	varchr(2), 	"sample_950"	varchr(2), 	"sample_1809"	varchr(2), 	"sample_943"	varchr(2), 	"sample_921"	varchr(2), 	"sample_985"	varchr(2), 	"sample_978"	varchr(2), 	"sample_923"	varchr(2), 	"sample_1018"	varchr(2), 	"sample_981"	varchr(2), 	"sample_50036"	varchr(2), 	"sample_931"	varchr(2), 	"sample_932"	varchr(2), 	"sample_990"	varchr(2), 	"sample_969"	varchr(2), 	"sample_946"	varchr(2), 	"sample_971"	varchr(2), 	"sample_1016"	varchr(2), 	"sample_940"	varchr(2), 	"sample_905"	varchr(2), 	"sample_1000"	varchr(2), 	"sample_961"	varchr(2), 	"sample_951"	varchr(2), 	"sample_1028"	varchr(2), 	"sample_982"	varchr(2), 	"sample_941"	varchr(2), 	"sample_967"	varchr(2), 	"sample_1629"	varchr(2), 	"sample_50026"	varchr(2), 	"sample_903"	varchr(2), 	"sample_987"	varchr(2), 	"sample_1007"	varchr(2), 	"sample_959"	varchr(2), 	"sample_918"	varchr(2), 	"sample_50015"	varchr(2), 	"sample_1832"	varchr(2), 	"bystro_sampleMaf"	varchr(20), 	"bystro_phastCons"	varchr(20), 	"bystro_phyloP"	varchr(20), 	"bystro_cadd"	varchr(20), 	"ccrs"	varchr(20), 	"domain_limbr"	varchr(20), 	"exone_limbr"	varchr(20), 	"is_ccds"	varchr(1), 	"bystro_final_function"	varchr(256), 	"bystro_gene_name"	varchr(256), 	"annovar_gene_name"	varchr(256), 	"annovar_function"	varchr(256), 	"annovar_function_detail"	varchr(256), 	"vep_gene_id"	varchr(256), 	"vep_feature"	varchr(256), 	"vep_consequence"	varchr(256), 	"spliceAI_gene_name"	varchr(256), 	"dmis_gene_name"	varchr(256), 	"dmis_HGVSc"	varchr(256), 	"dmis_HGVSp"	varchr(256), 	"dsplicing_gene_name"	varchr(256), 	"dsplicing_region"	varchr(256), 	"dsplicing_detailed_consequence"	varchr(256), 	"dsplicing_ada_score"	varchr(256), 	"dsplicing_rf_score"	varchr(256), 	PRIMARY KEY("id" AUTOINCREMENT))"""
    cursor.execute(cmd)

    # get data structure of variance table
    cursor.execute('PRAGMA table_info([variance])')
    table_info = cursor.fetchall()
    sample_col_index_list = [i[0] for i in table_info if i[1].startswith('sample_')]
    col_index2gen_id_dict = dict([[i[0], i[1][7:]] for i in table_info if i[1].startswith('sample_')])

    # get the data to handle
    cursor.execute('select * from variance')
    data = cursor.fetchall()
    data = [list(i) for i in data]
    filter_list = []
    # check overlap
    total_len = len(data)
    process_counter = 0
    for idx_row in range(len(data)):
        process_counter += 1
        if process_counter % 100 == 0:
            print('check overlap {0}\t/\t{1}'.format(process_counter, total_len))
            sys.stdout.flush()

        icounter = 0
        for idx_col in sample_col_index_list:
            value = data[idx_row][idx_col]
            if value == 0 or value == 'na':
                continue

            genid = col_index2gen_id_dict[idx_col]
            chrom = data[idx_row][1]
            pos_hg38 = data[idx_row][2]
            ref = data[idx_row][3]
            alt = data[idx_row][4]
            key = "{0}_{1}_{2}_{3}_{4}".format(chrom, pos_hg38, ref, alt, genid)
            assert key in overlap_dict
            if overlap_dict[key]:
                icounter += 1
                continue
            data[idx_row][idx_col] = 0
        filter_list.append(icounter > 0)
    new_data = list(compress(data, filter_list))
    total_len = len(new_data)
    process_counter = 0
    for new_line in new_data:
        process_counter += 1
        if process_counter % 100 == 0:
            print('executing sql {0}\t/\t{1}'.format(process_counter, total_len))
            sys.stdout.flush()
        value_tmp = ''
        for i in new_line:
            value_tmp += ('"' + str(i) + '",') if i is not None else 'NULL,'

        cmd = "insert into {0} ({1}) values ({2})" \
              "".format(new_table,
                        ', '.join([i[1] for i in table_info]),
                        value_tmp.strip(","))
        # print(cmd)
        cursor.execute(cmd)
    cursor.close()
    conn.commit()
    conn.close()
    print("all done")


def variant_table_info(db_file, table_name):
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()
    cursor.execute('PRAGMA table_info([variance])')
    table_info = cursor.fetchall()
    sample_col_index_list = [i[0] for i in table_info if i[1].startswith('sample_')]
    cursor.execute('select * from {0}'.format(table_name))
    data = cursor.fetchall()
    snp_sample_num = 0
    indel_sample_num = 0
    snp_num = 0
    indel_num = 0
    for data_line in data:
        is_indel = (len(data_line[3]) != len(data_line[4]))
        tmp = 0
        for idx in sample_col_index_list:
            value = data_line[idx]
            if value == 1 or value == 2:
                tmp += 1
        if is_indel:
            indel_sample_num += tmp
            indel_num += 1
        else:
            snp_sample_num += tmp
            snp_num += 1
    print('db file: {}'.format(db_file))
    print("table: {}".format(table_name))
    print('snp number: {}'.format(snp_num))
    print('indel number: {}'.format(indel_num))
    print('snp sample pair number: {}'.format(snp_sample_num))
    print('indel sample pair number: {}'.format(indel_sample_num))


def select_candidate_case_from_binom_pileup_sort(pileup_file, binom_case_c, topN):
    case2control_dict = {"A2.B0": "A0.B2",
                         "A3.2.B1.0": "A1.0.B3.2",
                         "A3.B0": "A0.B3",
                         "A3.B1": "A1.B3"}
    ret_list = []

    def get_criteria(pCHD, pCTD, pTof, pControl):
        vCHD = -10 * math.log10(pCHD)
        vCTD = -10 * math.log10(pCTD)
        vTof = -10 * math.log10(pTof)
        vControl = -10 * math.log10(pControl)
        ave = (vCHD + vCTD + vTof) / 3
        if ave == 0:
            return -1
        min_v = min(vCHD, vCTD, vTof)
        if vControl > min_v:
            return -1
        return ave - vControl

    if binom_case_c not in case2control_dict:
        print("{0} not in case2control dict. Can not convert to control condition".format(binom_case_c),
              file=sys.stderr)
        exit(0)
    binom_control_c = case2control_dict[binom_case_c]
    print("#pileup file: {}".format(pileup_file))
    print("#binom case condition: {}".format(binom_case_c))
    print("#binom control condition: {}".format(binom_control_c))
    print("#criteria\tpCHD\tpCTD\tpTof\tpControl\tgene_set_name\tfet_c")
    print('topN={}'.format(topN), file=sys.stderr)
    # total_line_num = file_line_num(pileup_file, '-')
    prekey = ''
    icounter = 0
    pCHD = -1
    pControl = -1
    pCTD = -1
    pTof = -1
    with open(pileup_file, 'r') as fp:
        while True:
            line = fp.readline()
            icounter += 1
            if icounter % 10000 == 0:
                print('handling {0} '.format(icounter), file=sys.stderr)
                sys.stderr.flush()
            if not line:
                break
            if line.startswith('#'):
                continue
            data_list = line.strip().split('\t')
            gene_set_name = data_list[0]
            fet_c = data_list[1]
            phenotype = data_list[2]
            binom_c = data_list[3]
            p2 = data_list[9]
            key = '{0}\t{1}'.format(gene_set_name, fet_c)
            if key != prekey and prekey:
                assert pCHD != -1
                assert pControl != -1
                assert pCTD != -1
                assert pTof != -1
                criteria = get_criteria(pCHD, pCTD, pTof, pControl)
                if len(ret_list) < topN:
                    ret_list.append([criteria, pCHD, pCTD, pTof, pControl, prekey])
                    ret_list = sorted(ret_list, key=lambda x: x[0], reverse=True)
                else:
                    ret_list.append([criteria, pCHD, pCTD, pTof, pControl, prekey])
                    ret_list = sorted(ret_list, key=lambda x: x[0], reverse=True)
                    ret_list.pop()
                pCHD = -1
                pControl = -1
                pCTD = -1
                pTof = -1
                prekey = key
            else:
                prekey = key
            if phenotype == 'heart6' and binom_c == binom_case_c:
                pCHD = float(p2)
                continue
            if phenotype == 'heart6' and binom_c == binom_control_c:
                pControl = float(p2)
                continue
            if phenotype == 'CTD' and binom_c == binom_case_c:
                pCTD = float(p2)
                continue
            if phenotype == 'tof' and binom_c == binom_case_c:
                pTof = float(p2)
                continue
    assert pCHD != -1
    assert pControl != -1
    assert pCTD != -1
    assert pTof != -1
    criteria = get_criteria(pCHD, pCTD, pTof, pControl)
    if len(ret_list) < topN:
        ret_list.append([criteria, pCHD, pCTD, pTof, pControl, prekey])
        ret_list = sorted(ret_list, key=lambda x: x[0], reverse=True)
    else:
        ret_list.append([criteria, pCHD, pCTD, pTof, pControl, prekey])
        ret_list = sorted(ret_list, key=lambda x: x[0], reverse=True)
        ret_list.pop()
    print('all done', file=sys.stderr)
    for ret_ele in ret_list:
        print("\t".join([str(i) for i in ret_ele]), file=sys.stdout)


def get_file_list(gatk_path):
    g = os.walk(gatk_path)
    ret = []
    for path, dir_list, file_list in g:
        for file in file_list:
            if not file.endswith('.vcf'):
                continue
            ret.append(os.path.join(path, file))
        break
    return ret


def parse_sample(format_list, sample_list, target):
    assert target in format_list
    return sample_list[format_list.index(target)]


def combine_gatk(gatk_path, output, target_chr):
    gatk_file_list = get_file_list(gatk_path)
    data_dict = {}
    slid_list = []
    counter2 = 0
    for gatk_file in gatk_file_list:
        counter2 += 1
        slid = re.findall('\S+/(\S+?)\.GATK', gatk_file)[0]
        slid_list.append(slid)
        data_loaded = False
        print('loading data ... {1} / {2}\t{0}\ttarget_chr={3}'.format(slid, counter2,
                                                                       len(gatk_file_list), target_chr))
        with open(gatk_file, 'r') as fp:
            while True:
                line = fp.readline()
                if not line:
                    break
                line = line.strip()
                if not line:
                    continue
                if line.startswith('#'):
                    continue
                data_list = line.split('\t')
                if data_list[6] != 'PASS':
                    continue
                if data_list[0] != target_chr:
                    if not data_loaded:
                        continue
                    else:
                        break
                if not data_loaded:
                    data_loaded = True
                keyword = "{0}\t{1}\t{2}\t{3}".format(data_list[0], data_list[1], data_list[3], data_list[4])
                format_list = data_list[8].split(':')
                sample_list = data_list[9].split(':')
                gt = parse_sample(format_list, sample_list, 'GT')
                gq = parse_sample(format_list, sample_list, 'GQ')
                if keyword in data_dict:
                    data_dict[keyword][slid] = '{0}:{1}'.format(gt, gq)
                else:
                    data_dict[keyword] = {slid: '{0}:{1}'.format(gt, gq)}
                # curr_percentage = math.floor(counter / float(variants_num) * 100)


    print('writing result ...')
    percentage = -1
    counter = 0
    with open(output, 'w') as fp:
        fp.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")
        for slid in slid_list:
            fp.write('\t' + slid)
        fp.write('\n')
        variants_num = len(data_dict)
        for keyword in data_dict:
            genotype_dict = data_dict[keyword]
            chrom, pos, ref, alt = keyword.split('\t')
            ret = '{0}\t{1}\t.\t{2}\t{3}\t.\tPASS\tNS={4}\tGT:GQ' \
                  ''.format(chrom, pos, ref, alt, len(slid_list))
            for slid in slid_list:
                if slid in genotype_dict:
                    ret += "\t" + genotype_dict[slid]
                else:
                    ret += '\t0/0:99'
            ret += '\n'
            fp.write(ret)
            counter += 1
            curr_percentage = math.floor(counter / float(variants_num) * 100)
            if curr_percentage != percentage:
                percentage = curr_percentage
                print('writing result ... {}%'.format(percentage))
    print('all done')


def combine_gatk2(gatk_path, output):
    def my_print(msg):
        sys.stdout.write("{0} {1}\n".format(time.asctime(time.localtime(time.time())), msg))
        sys.stdout.flush()

    gatk_file_list = get_file_list(gatk_path)
    # print(gatk_file_list)
    slid_list = []
    conn = sqlite3.connect('gatk.db')
    cursor = conn.cursor()
    cursor.execute("DROP TABLE IF EXISTS variant")
    sql_cmd = 'create table variant (id int primary key, chrom varchr(64), pos varchr(64),' \
              'vcf_id varchr(64), ref varchr(512), alt varchr(512), qual varchr(8),' \
              'filter varchr(64), info varchr(64), format varchr(64)'
    for gatk_file in gatk_file_list:
        slid = re.findall('\S+/(\S+?)\.GATK', gatk_file)[0]
        slid_list.append(slid)
        sql_cmd += ', {} varchr(64)'.format(slid)
    sql_cmd += ')'
    # print(sql_cmd)
    cursor.execute(sql_cmd)
    cursor.execute('CREATE INDEX idx_chrom ON variant(chrom);')
    cursor.execute('CREATE INDEX idx_pos ON variant(pos);')
    cursor.execute('CREATE INDEX idx_vcf_id ON variant(vcf_id);')
    cursor.execute('CREATE INDEX idx_ref ON variant(ref);')
    cursor.execute('CREATE INDEX idx_alt ON variant(alt);')
    counter2 = 0
    key_set = set([])
    for gatk_file in gatk_file_list:
        sql_list = []
        counter2 += 1
        variants_num = file_line_num(gatk_file, "#")
        slid = re.findall('\S+/(\S+?)\.GATK', gatk_file)[0]
        percentage = -1
        counter = 0
        with open(gatk_file, 'r') as fp:
            while True:
                line = fp.readline()
                if not line:
                    break
                line = line.strip()
                if not line:
                    continue
                if line.startswith('#'):
                    continue
                counter += 1
                # if 'SL98377' == slid and (counter % 10000 == 0 or counter == 1):
                #     print("{0} / {1}".format(counter, variants_num))
                data_list = line.split('\t')
                if data_list[6] != 'PASS':
                    continue
                format_list = data_list[8].split(':')
                sample_list = data_list[9].split(':')
                gt = parse_sample(format_list, sample_list, 'GT')
                gq = parse_sample(format_list, sample_list, 'GQ')
                condition = "chrom='{0}' and pos='{1}' and vcf_id='{2}' and ref='{3}' and alt='{4}'" \
                            "".format(data_list[0], data_list[1], data_list[2], data_list[3], data_list[4])
                keyword = "{0}\t{1}\t{2}\t{3}".format(data_list[0], data_list[1], data_list[3], data_list[4])
                # Determine whether there is a data row
                if keyword in key_set:
                    sql_cmd = "UPDATE variant SET {0}='{1}' WHERE {2}" \
                              "".format(slid, '{0}:{1}'.format(gt, gq), condition)
                else:
                    sql_cmd = "insert into variant (chrom, pos, vcf_id, ref, alt, {0}) " \
                              "values ('{1}', '{2}','{3}','{4}','{5}','{6}')" \
                              "".format(slid, data_list[0], data_list[1], data_list[2], data_list[3], data_list[4],
                                        '{0}:{1}'.format(gt, gq))
                    key_set.add(keyword)
                # print(sql_cmd)
                # if data_list[1] == '10409':
                #     print('keyword=[{0}] sql_cmd={1}'.format(keyword, sql_cmd))
                #     print(key_set)
                sql_list.append(sql_cmd)
                curr_percentage = math.floor(counter / float(variants_num) * 100)
                if curr_percentage != percentage:
                    percentage = curr_percentage
                    my_print('collecting sql {0} / {1}\t{2} ... {3}%'.format(counter2, len(gatk_file_list), slid,
                                                                             percentage))
        counter = 0
        percentage = -1
        for sql_cmd in sql_list:
            cursor.execute(sql_cmd)
            counter += 1
            if slid == 'SL98377' and counter % 100 == 0:
                print("{0} / {1}".format(counter, len(sql_list)))
            curr_percentage = math.floor(counter / float(len(sql_list)) * 100)
            if curr_percentage != percentage:
                percentage = curr_percentage
                my_print('executing {0} {1}%'.format(slid, percentage))

    cursor.close()
    conn.commit()
    conn.close()


def grep_data(file, chrom, pos, ref, alt):
    cmd = "grep -P '{0}\\t{1}\\t.\\t{2}\\t{3}' {4}".format(chrom, pos, ref, alt, file)
    #print(cmd)
    pp = Popen([cmd], shell=True, stdout=PIPE)
    pipe_ret = pp.stdout.readlines()

    if len(pipe_ret) > 0:
        assert len(pipe_ret) == 1
        data_list = pipe_ret[0].strip().split('\t')
        format_list = data_list[8].split(':')
        sample_list = data_list[9].split(':')
        gt = parse_sample(format_list, sample_list, 'GT')
        gq = parse_sample(format_list, sample_list, 'GQ')
        ret = '{0}:{1}'.format(gt, gq)
    else:
        ret = '0/0:99'
    #print(ret)
    return ret


def combine_gatk3(gatk_path, output):
    gatk_file_list = get_file_list(gatk_path)
    slid_list = []
    key_set = set([])
    counter2 = 0
    for gatk_file in gatk_file_list:
        counter2 += 1
        variants_num = file_line_num(gatk_file, "#")
        slid = re.findall('\S+/(\S+?)\.GATK', gatk_file)[0]
        slid_list.append(slid)
        counter = 0
        percentage = -1
        with open(gatk_file, 'r') as fp:
            while True:
                line = fp.readline()
                if not line:
                    break
                line = line.strip()
                if not line:
                    continue
                if line.startswith('#'):
                    continue
                counter += 1
                data_list = line.split('\t')
                if data_list[6] != 'PASS':
                    continue
                keyword = "{0}\t{1}\t{2}\t{3}".format(data_list[0], data_list[1], data_list[3], data_list[4])
                key_set.add(keyword)

                curr_percentage = math.floor(counter / float(variants_num) * 100)
                if curr_percentage != percentage:
                    percentage = curr_percentage
                    print('loading data ... {2} / {3}\t{0}\t{1}%'.format(slid, percentage,
                                                                         counter2, len(gatk_file_list)))

    print('writing result ...')
    percentage = -1
    counter = 0
    with open(output, 'w') as fp:
        fp.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")
        for slid in slid_list:
            fp.write('\t' + slid)
        fp.write('\n')
        variants_num = len(key_set)
        for keyword in key_set:
            chrom, pos, ref, alt = keyword.split('\t')
            ret = '{0}\t{1}\t.\t{2}\t{3}\t.\tPASS\tNS={4}\tGT:GQ' \
                  ''.format(chrom, pos, ref, alt, len(slid_list))
            for gatk_file in gatk_file_list:
                ret += '\t' + grep_data(gatk_file, chrom, pos, ref, alt)
            ret += '\n'
            fp.write(ret)
            counter += 1
            curr_percentage = math.floor(counter / float(variants_num) * 100)
            if curr_percentage != percentage:
                percentage = curr_percentage
                print('writing result ... {}%'.format(percentage))
                sys.stdout.flush()


def combine_gatk_get_key_list(gatk_path, output):
    gatk_file_list = get_file_list(gatk_path)
    key_set = set([])
    with open(output, 'w') as fp_out:
        for gatk_file in gatk_file_list:
            with open(gatk_file, 'r') as fp:
                while True:
                    line = fp.readline()
                    if not line:
                        break
                    line = line.strip()
                    if not line:
                        continue
                    if line.startswith('#'):
                        continue
                    data_list = line.split('\t')
                    if data_list[6] != 'PASS':
                        continue

                    keyword = "{0}\t{1}\t{2}\t{3}".format(data_list[0], data_list[1], data_list[3], data_list[4])
                    if keyword not in key_set:
                        key_set.add(keyword)
                        fp_out.write("{0}\t{1}\t{2}\t{3}\t{4}\n"
                                     "".format(data_list[0], data_list[1], data_list[2], data_list[3], data_list[4]))


def gatk_split():
    def remove_redundancy(vcf_line):
        if vcf_line.startswith("#"):
            return vcf_line
        if not vcf_line.strip():
            return ""
        vcf_list = vcf_line.strip("\r").strip("\n").strip("\r").split("\t")
        if len(vcf_list) < 2:
            print('vcf line is [{}]'.format(vcf_line))
            sys.stdout.flush()
        pos = int(vcf_list[1])
        ref = vcf_list[3]
        alt = vcf_list[4]
        iright_delete = 0
        for index in range(min(len(alt), len(ref)) - 1):
            if ref[-(index + 1)] == alt[-(index + 1)]:
                iright_delete += 1
            else:
                break
        if iright_delete > 0:
            ref = ref[:-iright_delete]
            alt = alt[:-iright_delete]

        vcf_list[1] = str(pos)
        vcf_list[3] = ref
        vcf_list[4] = alt
        return "\t".join(vcf_list)

    while True:
        vcf_line = sys.stdin.readline()
        if not vcf_line:
            break
        if vcf_line.startswith("#"):
            sys.stdout.write(vcf_line)
            continue
        if not vcf_line.strip():
            continue
        data_list = vcf_line.strip().split('\t')
        if data_list[6] != 'PASS':
            continue
        format_list = data_list[8].split(':')
        sample_list = data_list[9].split(':')
        gq = parse_sample(format_list, sample_list, 'GQ')
        for alt in data_list[4].split(','):
            tmp_list = copy.copy(data_list)
            tmp_list[4] = alt
            tmp_list[2] = '.'
            tmp_list[7] = 'NS=1'
            tmp_list[8] = 'GT:GQ'
            tmp_list[9] = '0/1:{}'.format(gq)
            sys.stdout.write(remove_redundancy('\t'.join(tmp_list)) + '\n')

def export_geneset(db_file, geneset_out):
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()
    gene_table_info = cursor.execute("PRAGMA table_info([gene_table])")
    gene_set_col_name_list = [i[1] for i in gene_table_info if i[1].startswith('gene_set_')]
    sql_cmd = 'select ' + ', '.join(gene_set_col_name_list) + ' from gene_table'
    cursor.execute(sql_cmd)
    gene_set_data = cursor.fetchall()
    # gene_set_data_T = list(zip(*gene_set_data))
    cursor.execute('select gene_id from gene_table')
    gene_name_list = [i[0] for i in cursor.fetchall()]
    with open(geneset_out, 'w') as fp:
        for i in range(len(gene_set_col_name_list)):
            ret_str = gene_set_col_name_list[i][9:]
            for j in range(len(gene_name_list)):
                cell_data = gene_set_data[j][i]
                if cell_data == 1:
                    ret_str += '\t' + gene_name_list[j]
            fp.write(ret_str + '\n')
    cursor.close()
    conn.commit()
    conn.close()

def main():
    if len(sys.argv) == 1:
        print("""   Usage:
        1 filter1:
            According to the result file tsv of bystro，Create a new tsv file，
            To filter the original vcf file. The filter condition is the frequency of 
            occurrence in the population <= rate
            python wgsa.py filter1 <rate> <bystro_tsv> <output>
        2 data_reduce:
            Perform data pruning on certain columns specified in the tsv result of bystro. 
            The pruning rule is to keep only the largest value in the current column. 
            col separated by ','
            python wgsa.py data_reduce <file_in> <col> <output>
            eg:python wgsa.py data_reduce in.tsv 1,2,3 out.tsv
        3 get_words:
            Perform lexical analysis on the specified column of tsv to analyze 
            how many words appear
            python wgsa.py get_words <file_in> <col_list>
            eg: python wgsa.py get_words xxx.tsv 1/2/3
        4 reduce_words:
            According to word_order, perform vocabulary pruning on the specified column in 
            tsv, and keep only the top word in the word list.
            Which columns and word lists are analyzed are recorded in a config file 
            in the following format:
                col_num1        word1_1 word1_2 word1_3
                col_num2        word2_1 word2_2
            ps: col_num从1开始计算
            python wgsa.py reduce_words <input_file> <config_file> <output>
        5 filter_by_word:
            Filter tsv files according to word, filter conditions:
                mode==0: Keep the row as long as one of the selected columns hits
                mode==1: All selected columns must be hit to keep the row
            python wgsa.py filter_by_word <file_in> <word_file> <mode> <output>
            word file The format is the same as the word_order format in function 4
        6 reclassify_ALT:
            Process the lof part of the annovar result, and treat the alt that is not lof as ref
            python wgsa.py reclassify_ALT <lof_file> <output>
        7 filter_alt:
            Process vcf and delete genotypes that do not appear in the sample. The result is 
            saved in output. If there is a change in the snp, record a separate copy in output2.
            python wgsa.py filter_alt <vcf_file> <output> <output2>    
        8 left_normalization:
            To process vcf, process the alt, ref and pos in indels with the method of 
            left normalization. support pipeline
            python wgsa.py left_normalization <vcf_file> <output>
            cat <vcf_file> | python wgsa.py left_normalization
        9 reduce_duplicate:
            Process the annovar results and de-duplicate the same mutations 
            according to the given priority list
            The priority list format is the same as the word_order format in function 4
            python wgsa.py reduce_duplicate <annovar_result> <col_ref> <col_alt> <col_chr> <col_pos> <list_file> <output>  
        10 rebuild_vcf:
            According to the filtered annovar results, restore the vcf. Treat the alt 
            that does not appear as a ref in the genotype。（group non-damaging alt 
            alleles back to reference）
            python wgsa.py rebuild_vcf <file_in> <output>    
        11 rebuild_vcf_bystro:
            Based on the filtered bystro results, rebuild the vcf. Ignore the 
            sample and only process the first 7 rows.
            python wgsa.py rebuild_vcf_bystro <tsv_in> <vcf_in> <output>    
        12 rebuild_vcf_vep:
            Based on the LOF results annotated by VEP, VCF was initially screened for ID.
            python wgsa.py rebuild_vcf_vep <lof_in> <vcf_in> <need_head> <output>    
        13 rebuild_vcf_vep2:
            According to the lof result of the vep comment, refactor the vcf, and 
            treat the alt that is not lof as a ref. Delete the same id but not lof.
            python wgsa.py rebuild_vcf_vep2 <lof_in> <vcf_in> <need_head> <output>    
        14 check_id:
            Check whether the id in vcf is repeated, and print the repeated id on the screen
            python wgsa.py check_id <vcf_in> <id_col>
        15 UI:
            Find the intersection or union of multiple files
            python wgsa.py UI <file_list> <col_list> <output> <mode>
            file_list: file names separated by ","
            col_list: column numbers separated by "," It is used to specify the key columns 
            of the union of intersection sets, and only these columns are kept in 
            the output result.
            mode: "union" or "intersection"
            eg: python wgsa.py UI file1,file2,file3 1,2,3,4 oo union
        16 split_vcf:
            Split each row of data in the original vcf into multiple rows of data 
            according to the principle of one alt per row (only the first 8 columns are kept). 
            Support for pipes.
            python wgsa.py split_vcf <vcf_in> <output>
            cat <vcf_in> | python wgsa.py split_vcf
        17 add_colums:
            According to the key column, add the column pointed to by lib file 
            add_colums to the left or right of data_file.
            python wgsa.py add_colums <data_file> <data_col> <lib_file> <lib_col> <add_cols> <output> <mode>
            mode: left Data is padded to the left. right Data is filled to the right
            eg python wgsa.py add_colums in.vcf 1,2,4,5 lib_file 1,2,4,5 8 output left
        18 splice_ai_filter:
            Screen spliceAI genome-wide data (all SNPs) annotated vcf, support pipeline
            python wgsa.py splice_ai_filter <input> <score> <output>
            eg python wgsa.py splice_ai_filter input_vcf 0.8 output_file
               cat input_vcf | python wgsa.py splice_ai_filter 0.8 | less -S
            ps:DS_AG, DS_AL, DS_DG, As long as one of the four scores of DS_DL is 
            greater than or equal to the input score, it will be selected
        19 cut:
            Divide the data into equal-width bins and support pipelines
            python wgsa.py cut <data_in> <num> <max_in> <min_in> <output>
            cat data_in | python wgsa.py cut num max_in min_in | less
        20 add_vcf_head:
            Add a head to the vcf without a head, and the contig is generated 
            from the index fai file of fa.
            python wgsa.py add_vcf_head <ver> <fai_in> <vcf_in> <output> <org_vcf>
            ver is hg19 or hg38
            org_vcf the original vcf file to get the table head
        21 select_vcf:
            Extract vcf data, extract a row every num-1 data
            python wgsa.py select_vcf <vcf_in> <num> <output>
            cat vcf_in | python wgsa.py select_vcf num | less -S
        22 component:
            Split a vcf into four files: one snp, multiple snp, one indel, and multiple indels.
            python wgsa.py component <vcf_in> <output_path> <need_head>
        23 splice_ai_filter_index:
            Filter the indel data of spliceAI
            python wgsa.py splice_ai_filter_indel <vcf_in> <score> <output>
            cat vcf_in | python wgsa.py splice_ai_filter_indel score | less
        24 select_by_cols:
            Filter vcf according to some columns of list file
            Both vcf_cols and list_cols are separated by ","
            python wgsa.py select_by_cols <vcf_in> <vcf_cols> <list_file> <list_cols> <output> <not_found>
            eg: python wgsa.py select_by_cols in 1,2,4,5 list 1,2,3,4 output not_found
        25 splice_ai_merge_snp_indel:
            Merge the filtered spliceAI snp anno and spliceAI indel anno into one file
            The comment column and comment format of snp are different from indel 
            and need to be converted
            The annotation of snp is in the first column, format:
                SYMBOL=TUBB8;STRAND=-;TYPE=E;DIST=-53;DS_AG=0.0000;DS_AL=0.0000;DS_DG=0.0000;DS_DL=0.0000;DP_AG=-26;DP_AL=-10;DP_DG=3;DP_DL=35
            The comment of indel is in the eighth column, the format is:
                ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL
                SpliceAI=C|SAMD11|0.00|0.01|0.00|0.00|-126|-131|396|-53
            The anno format of the result:
                SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL
            python wgsa.py splice_ai_merge_snp_indel <snp_in> <indel_in> <output>
        26 build_vcf_index:
            Create an index for the original vcf, through which the offset of the data 
            row where the original vcf is located can be queried from the vcf data after
            split and left normalization
            The index file name is vcfi
            python wgsa.py build_vcf_index <vcf_in>
        27 create_db:
            Create a table with the data of union0.01
            python wgsa.py create_db <union01_file> <union001_file> <annovar01> <bystro01> <dmis01> <dsplicing01> <spidex01> <spliceai01> <vep01> <db_file>
        28 db_add_sample2table_variance:
            Add the sample data to the table variance
            python wgsa.py db_add_sample2table_variance <vcf_file> <db_file> <table_name>
        29 check_vcf_missing:
            Check whether half of the missing genotypes such as ./1 ./0 exist in the sample of vcf
            python wgsa.py check_vcf_missing <vcf_in> 
        30 export_vcf_from_db:
            Export the vcf data from the variance table (the chr pos in the database
            is the hg38 version)
            python wgsa.py export_vcf_from_db <db_file> <table_name> <variance_restrict> <org_vcf> <output> <fai_in>
            variance_restrict eg: v.is001=1 AND v.vep=1
        31 db_add_spliceAI_anno:
            Add the annotation of spliceAI to the variance table.
            python wgsa.py db_add_spliceAI_anno <db_file> <annotation>
            ps: annotation is the merged file of snp and indel
        32 db_add_spidex_anno
            Add the score of spidex to the database (to be completed)
            python wgsa.py db_add_spidex_anno <db_file> <spidex_file>
        33 db_add_vcf_pos
            Add the pos of the original vcf to db
            python wgsa.py db_add_vcf_pos <db_file> <vcf_file>
        34 db_add_bystro_anno
            Save bystro's comments (47 phastCons, 48 ​​phyloP, 49 cadd columns) to db
            python wgsa.py db_add_bystro_anno <db_file> <tsv_file> <table_name>
        35 check_dbNSFP
            Check whether the dbNSFP annotation results are empty after column 
            411 when there are multiple genenames.
            If not empty, it will be printed in the log.
            python wgsa.py check_dbNSFP <anno_in>
        36 check_protein_coding_gene
            Check whether the protein_coding_gene file has the same gene id and gene 
            name but different positions.
            python wgsa.py check_protein_coding_gene <file_in>
        37 check_protein_coding_transcript
            Examine whether different transcripts of the same gene are distributed
            in multiple isolated regions
            python wgsa.py check_protein_coding_transcript <file_in>
            python wgsa.py check_dbNSFP anno_in
        38 check_dms
            Check the dms file downloaded on ucsc
            python wgsa.py check_dms <dms_in>
        39 split_gene_centric_annotation (problematic)
            Split the gene centric annotation file according to Ensembl_geneid, and the corresponding annotation column should also be split. After splitting, deduplication
        40 db_build_synonymous_snp_table
            Build the synonymous_snp table
            python wgsa.py db_build_synonymous_snp_table <db_file> <table_name> <annovar01> <bystro01> <vep01> <union001_file> <union01_file>
        41 db_build_gene_table
            Build the gene table
            python wgsa.py db_build_gene_table <protein_coding_gene_transcript> <db_file> <table_name>
        42 db_copy_table
            Import the data table saved in a file into the database
            python wgsa.py db_copy_table <db_file> <table_name> <table_file>
        43 db_add_gene_anno
            Add the annotation of the gene to the gene_table.
            python wgsa.py db_add_gene_anno <db_file> <table_name> <gene_file> <gene_id_col>[spliter]
            spliter: According to gene_file, it can be tab or comma, the default value is tab
            gene_file: Comments are in the gene_file file, starting from the third column
            gene_id_col: gene The column where the id is located, can be 1 or 2
        44 db_reduce_col_duplicate
            Deduplicate a column of a table in the database
            python wgsa.py db_reduce_col_duplicate <db_file> <table_name> <col_name>
        45 collapse
            For each gene, select the snp within its range from the snp table
            python wgsa.py collapse <db_file> <variance_table_name> <gene_table_name> <variance_restrict> <gene_restrict> <output>
            variance_restrict eg: v.is001=1 AND v.bystro=1
            gene_restrict eg: g.gene_name IN ('SOCS4', 'RAB22A')
        46 build_contingency_table
            According to the collapse result and phenptype, a four-table table is generated
            for the sample
            python wgsa.py build_contingency_table <db_file> <collapse_file> <phenotype> <output> <sample_restrict>
            
            build_contingency_table_new <db_file> <phenotype> <sample_table_name> <sample_restrict> <gene_table_name> 
                                        <gene_restrict> <variance_table> <variance_restrict> <fai_in> <output> <job_id>
            
            build_contingency_table_new_with_gene_list_file <db_file> <phenotype>
                                                    <sample_table_name> <sample_restrict>
                                                    <gene_table_name> <gene_list_file>
                                                    <variance_table> <variance_restrict>
                                                    <fai_in>
            gene_list_file:
            # comment
            name1   gene1_1   gene1_2   gene1_3
            name2   gene2_1   gene2_2   gene2_3
            
            The new version uses gene name to select variance, see 114
            
        47 filter_collapse_and_vcf_with_gene_list
            Filter collapse and vcf according to gene_list
            python wgsa.py filter_collapse_and_vcf_with_gene_list <org_vcf> <org_collapse> <gene_list_file> <vcf_out> <collapse_out>
        48 get_variance_num_in_sample
            Calculate the number of alleles with variance for each person
            python wgsa.py get_variance_num_in_sample <db_file> <variance_table> <variance_restrict> <sample_table> <sample_restrict> <output>
        49 get_variance_function_number
            Preparing data tables for R plotting
            python wgsa.py get_variance_function_number <file_all> <file01> <file001> <col> <output>
        50 t_test
            T test
            python wgsa.py t_test <file_in> <value_col> <group_cols>
            group_cols: separated by commas
        51 split_vcf_with_sample
            Split vcf with sample
            cat <vcf_in> | python wgsa.py split_vcf_with_sample > output
        52 cross_genehancer_snp
            Merge the genehancer and snp files together by location
            python wgsa.py cross_genehancer_snp <hancer_file> <snp_file> <output>
        53 cross_genehancer_snp_worker
            Add a cross_genehancer_snp_worker for use by cross_genehancer_snp_multiple.
            You can also use this function to manually add the number of workers
            python wgsa.py cross_genehancer_snp_worker <server_addr> <snp_file> <worker_id>
        54 cross_genehancer_snp_multiple
            cross_genehancer_snp的分布式启动
            python wgsa.py cross_genehancer_snp_multiple <num> <hancer_file> <snp_file> <output_file>
        55 test_ks
            Computes the Kolmogorov-Smirnov statistic on 2 samples.
            test_ks <file_in> <value_col> <group_cols> <out_path> <p>
            Use ks to analyze the data and save p significantly
            value_col: observation column
            group_cols: Grouping variable columns, can be multiple columns, 
            separated by commas. Example: 2,3,4

        56 test_MWW
            test_MWW <file_in> <value_col> <group_cols> <out_path> <p>
            Compute the Mann-Whitney rank test on samples x and y.

            In statistics, the Mann–Whitney U test (also called the 
            Mann–Whitney–Wilcoxon (MWW), Wilcoxon rank-sum test, or 
            Wilcoxon–Mann–Whitney test) is a nonparametric test of the 
            null hypothesis that it is equally likely that a randomly 
            selected value from one population will be less than or 
            greater than a randomly selected value from a second population.
            https://en.wikipedia.org/wiki/Mann–Whitney_U_test
            Analyze the data with MWW and save p significant
            Parameters: Reference 50 t_test 55 test_ks
        57 get_geneset_variance_num_in_sample
            get_geneset_variance_num_in_sample <db_file> <variance_table> <variance_restrict> <sample_table>
                                       <sample_restrict> <output> <gene_id_list>
            Similar to 48, a constraint on the gene is added to the variance constraint.
            Genes are represented by gene ids, separated by commas: 
            ENSG00000074755, ENSG00000159840
        58 test_GLM
            Analyze the data with a generalized linear model and save p significantly
            References: https://blog.csdn.net/qq_41518277/article/details/85100652#Generalized_Linear_92
                    https://www.statsmodels.org/stable/glm.html
            usage:
            test_GLM <file_in> <value_col> <group_cols> [Binomial|Gamma|Gaussian|InverseGaussian|NegativeBinomial|Poisson|Tweedie] <out_path> <p>
            # Use the default Poisson distribution
            python wgsa.py test_GLM file 1 2,3,4 . 0.05
            # Use Poisson distribution
            python wgsa.py test_GLM file 1 2,3,4 Poisson . 0.05
            ps:
            Requires the statsmodels package
        59 Check the .table file. If the minimum value of p_value or p_value1 is 
           less than the input threshold, print the file name.
            check_pvalue_in_table <p_threshold_in> <table_file>
        60 bar_hist_plot <number_file> <phenotype>
            Draw the error histogram and distribution histogram according 
            to the number file and fit the distribution curve
        61 analyze_fisher_test_variance <database> <path>
            Perform a fisher test with various combinations of conditions 
            on the variance table, and store the results in the path
            p.s.: The number of combinations is as high as 2^7*4*3*14*7=150528 combinations! ! !
        62 whole_genome_variance_test <db_file> <variance_table> <sample_table> <sample_restrict> <outputpath> <p_threshold>
            Perform ks test, MWW test, and GLM test on variance. The results of the test 
            are screened by the p threshold and saved as a file
        63 analyze_fisher_test_variance_new <database> <variance_table> <sample_restrict> <fai_in> <path>
            Perform a fisher test with various combinations of conditions on the variance 
            table, and store the results in the path
            p.s.: The number of variance restrict combinations is as high as 
            2^7*4*3*14*7=150528 combinations! ! !
        64 map_vcf_sample_id <vcf_in> <id_pair_in> <vcf_out>
            Convert the sample id of vcf from the first column of id_pair to the 
            second column of id_pair according to id_pair
            Other data remains unchanged
        65 rerun_selected_combination <db_file> <variance_table_name> <selected_o_file> <path> <org_vcf> 
                                      <fai_in> <sample_restrict> <with_ccrs>
            Re-run the filtered table according to the output file of the filtered p value
        
        66 export_vcf_with_gene_list <db_file> <variance_table> <fai_in> <gene_list_in> <output>
            According to the genelist, filter out all the variances falling in the 
            gene from the variance table
            The first column of the genelist file is gene id, no header
        67 db_add_ccrs_limbr <db_file> <ccrs_file> <domain_file> <exone_file> <fai_in>
            According to ccrs, domain, exone files, add ccrs, domain_limbr, domain_limbr_quartile,
            exone_limbr、exone_limbr_quartile列
        68 analyze_fisher_test_variance_2 <database> <variance_table> <sample_restrict> <fai_in> <path>
            Added ccrs, domain_limbr, exone_limbr in different situations based on the 
            original fisher test.
            The number of combinations of variance restrict is increased to: 
            2^7*4*3*14*7*5*5*5=18,816,000 Refer to 69
        69 permutation_burden <file_in> <out_put> <times>
            Perform a permutation on the last column of the given file, the last 
            column consists of "1", "2", "NA".
        70 permutation_fisher <db_file> <phenotype> <sample_table_name> <sample_restrict> <gene_table_name> 
                              <gene_restrict> <variance_table> <variance_restrict> <fai_in> <path> <times>
            Perform permtation on the fisher test, and randomly exchange samples between 
            the case group and the control group while keeping the number unchanged. 
            Generate the results of times group permutation.
            path: The path to store the results needs to be established in advance
        71 analyze_fisher_test_variance_3 <database> <variance_table> <sample_restrict> <fai_in> <path>
            Conditional simplification of 66, the number of combinations of variance restrict is 1750
        72 modify_gene_name2id <db_file> <file_in> <gene_table_name> <hg19_gene_info> <output>
            According to the mapping between gene id and gene name in the database and hg19 version file (hg19_gene_info), replace the gene name in the input file (gene set) with
            gene id, Unable to convert typed in the log
            gene set file format:
            gene set name\tgene name1\tgene name2\n
            
            hg19_gene_info: /gs/gsfs0/users/yizhao/annotated.wgs/hg19.gene.name.txt
        73 modify_vcfid <list_in> <vcf_in> <output>
            According to the mapping relationship in the list file, map the id in vcf, 
            if it does not exist in the list file, keep it as it is
            
        74 add_gene_name <file_in> <db_file> <gene_table_name> <id_col> <output>
            Specify the gene id column in the file, and find the gene name according 
            to the gene id. Insert the gene name into the last column of the file
        75 modify_duplicateID <vcf_in> <output>
            Check the id of the vcf. Duplicates are found, and suffixes are added to 
            avoid repetitions.
            Note: vcf needs to be sorted by id in advance
        76 manhattan_direct <file_in> <name_col> <chr_col> <start_col> <end_col> <p_col> <fai_in> <output>
            Draw the Manhattan map directly from the input file
        77 get_gene_set_info <db_file> <gene_list_file> <output> <fai_in> <phenos> <variance_table_name>
            According to the gene list file, find the SNPs in each gene from the variance table, and classify the SNPs, and calculate the number of cases carried by each SNP in different phenotypes,
            And control the number of people to carry.
            Multiple phenotypes separated by commas
            The gene list contains only one column, which can be gene id or gene name or 
            a mixture of both
        78 copy_vcf2table <vcf_in> <table_name> <db_file>
            Import the vcf into the database.
            Use with caution: the table_name table in the database will be deleted first
        79 db_add_ccds <db_file> <ccds_bed> <fai_in>
            Add the is_ccds column to the variance table and synonymous_snp table in the database
        80 firth_logistic_regression <db_file> <sample_table> <phenotype_in> <sample_restrict>
                              <variance_table> <variance_restrict>
            Genome-wide firth_logistic_regression analysis, analyzing data in all variance tables
            return value [coef, lci, uci, p_value]
        81 firth_logistic_regression_with_gene <db_file> <sample_table> <phenotype_in> <sample_restrict>
                                               <variance_table> <variance_restrict> <gene_table> <gene_restrict>
                                               <fai_in>
            The firth_logistic_regression analysis of the limited gene, the gene is limited 
            by gene_table, gene_restrict
            return value [coef, lci, uci, p_value]
        82 firth_logistic_regression_genelist <db_file> <sample_table> <phenotype_in> <sample_restrict>
                                              <variance_table> <variance_restrict> <gene_table> <gene_restrict>
                                              <fai_in> <gene_list_file_in> <id_rows>
            Define the firth_logistic_regression analysis of the gene, and jointly define the gene through gene_table, gene_restrict, gene_list_file_in, id_rows
            return value [coef, lci, uci, p_value]
        83 firth_logistic_regression_plot <script_list_file_in>
            Use the plot script to draw the firth_logistic plot. Scripts are separated into multiple subplots by blank lines. The first line of the script is a fixed parameter. The first row in the data of each subgraph is title and xlable.
            Starting from the second line are the parameters of firth_logistic_regression_with_gene.
            For example:
            
            y1	y2	yticks ytickposes fontsize lablesize marker markersize	linewidth	db_file	fai_in left right   top bottom  hspace
            title1	xlable1
            sample_table11	phenotype11	sample_restrict11	variance_table11	variance_restrict11	gene_table11	gene_restrict11
            sample_table12	phenotype12	sample_restrict12	variance_table12	variance_restrict12	gene_table12	gene_restrict12
    
            y1，y2：ylim
            yticks：y-axis labels, separated by ","
            ytickposes：y-axis label position, separated by ","
            Marker: "." Dot See https://matplotlib.org/2.2.5/api/markers_api.html#module-matplotlib.markers for the rest
            db_file：database file
            fai_in：hg38's ref index
            left, right, top, bottom, hspace: Margin see https://blog.csdn.net/qq_33039859/article/details/79424858
        84 variance_distribution_plot <db_file> <phenotype> <sample_table> <sample_restrict> <variance_table> 
                                <variance_restrict> <gene_table> <gene_restrict> <fai_in> <output> <h>
            Use gene_restrict to select the gene region, and then select the variance that satisfies the variance_restrict and is located in the gene region in the variance table
            Count the number of people and alleles with variance in cases and controls, take the number of cases as the abscissa, and the number of controls as the ordinate, and draw the variance data on the graph.
            gene_restrict example:
            genelist:gene list file name  # Use a list of genes to select genes
            restrict:g.chr=1  # Select genes with constraints on one gene
            directly  # Ignore the gene, directly use the variance in all variance tables as the candidate
        
        85 select_enriched_variance <db_file> <sample_table> <sample_restrict> <phenotype_in> <variance_table_name>
                             <variance_restrict> <gene_table> <gene_restrict> <fai_in> <bigger_fraction>
                             <smaller_fraction> <enrich_in_case_out> <enrich_in_control_out>
            Define the enrichment degree of a variance in a case by dividing the number of people who appear in a variance in a case by the total number of people in a case.
            Similarly, the number of people who appear in the variance in the control divided by the total number of people in the control is the enrichment of the variance in the control.
            enrich_in_case_out: The file is saved in the case where the enrichment is greater than bigger_fraction, and in the control where the enrichment is less than
                                The variance of smaller_fraction. And the attributes of the variance such as gene, ccrs, limbr, etc.
            enrich_in_control_out: This file is saved in the control, the enrichment degree is greater than bigger_fraction, enriched in the case
                                    variance with degree less than smaller_fraction. And the attributes of the variance such as gene, ccrs, limbr, etc.
            
            
            select_enriched_variance_gene_name <db_file> <sample_table> <sample_restrict> <phenotype_in>
                                       <variance_table_name> <variance_restrict> <gene_table> <gene_restrict>
                                       <bigger_fraction> <smaller_fraction> <output> <group>
            Enrichment is defined as above.
            output: Record the variance that is enriched in the specified group and not enriched in another group, and the gene, ccrs, limbr and other attributes of the variance
            group: case or control. Specifies the group to analyze
        86 gene_distribution_plot <db_file> <phenotype_in> <sample_table> <sample_restrict> <variance_table>
                           <variance_restrict> <gene_table> <gene_restrict> <fai_in> <output> <h>
            Investigate in units of gene. There is a variance in the gene called the gene being hit. Use the average number of hits in the case and the average number of hits in the control as coordinates to draw the gene graph
        87 qqplot_fisher_permutation <path> <fai_in> <output>
           Used in conjunction with 70 permutation_fisher. The path here is the output path of 70.
           output为png文件
           Internally, R is automatically called to execute the qqplot_permutation.R script. Dependencies, see R script internals for details.
        88 gene_set_fisher_test <gene_set_file> <db_file> <sample_table> <sample_restrict> <phenotype_in>
                         <variance_table> <variance_restrict> <gene_table> <fai_in> <output>
           Perform fisher exact analysis on a gene set basis.
           gene set 文件格式如下：
           set1    gene11    gene12    gene13
           set2    gene21    gene22
        89 gene_set_fisher_test2 <db_file> <sample_table> <sample_restrict> <phenotype_in> <variance_table>
                          <variance_restrict> <gene_table> <gene_restrict> <gene_set_name> <fai_in> <output>
           Same as 88. With gene_restrict as the constraint, select a gene set from the gene_table and conduct fisher exact test analysis.
        90 analyze_fisher_test_variance_May19 <database> <sample_restrict> <fai_in> <path>
           Multiple combinations of fisher exact test
        91 analyze_geneset_fisher_test <gene_set_file> <database> <sample_restrict> <fai_in> <path>
           Take the gene set as the unit, and carry out fisher exact test according to various constraint combinations of variance
        92 fisher_test_synonymous <database> <sample_restrict> <fai_in> <path>
           Do gene based fisher test on synonymous table
        93 multiple_geneset_variance_num_in_sample <db_file> <variance_table> <variance_restrict> <sample_table> 
                                            <sample_restrict> <output> <gene_set_file>
            Refer to 57, 48. For the geneset in the geneset file, count the number of variance alelle in the gene set for each person
        94 prepare_liftover_bed <org_bed_in> <value_col_num> <bed_out>
            Before doing the lift over, it is necessary to split the area into positions for the bed of hg19. It is wrong to directly lift over without splitting.
                   hg19     ------------>        hg38
            chr1  100  103                  chr1  200  205
            
            It is wrong to convert directly like above.
            The right way to conver should be:
            
                            prepare bed                   lift over
                  hg19     -------------->    hg19     --------------->        hg38
            chr1  100  103               chr1  100  101                  chr1  200  201
                                         chr1  101  102                  chr1  204  205
        95 gether_bed_region <input> <output>
            Classify the position with the same value as an interval, which is the reverse process of 94
        
        96 check_gene_overlap <database> <gene_table> <fai_in> <output>
            检查基因的重叠情况
        97 db_add_bystro_gene_name <db_file> <bystro_tsv> <ref_in>
            Add bystro annotations to the database RefSeq gene name 和 RefSeq transcript ID
        98 select_gene_name_in_tsv <tsv_in> <tsv_out>
            For the original bystro tsv file, delete the gene name according to refSeq.siteType and refSeq.exonicAlleleFunction
        99 select_gene_name_in_annovar_result <file_in> <file_out>
            For the annovar annotation file, calculate the pos of the variance in the database, and merge the function and detail
        100 reduce_dsplicing_result <file_in> <output>
            Prune the result of dsplicing, keep only the data with ada_score>=0.6, and delete the columns that are not written into the database
        101 reduce_vep_result <file_in> <output> <ref_in>
            Delete the results of vep, merge the same variance rows, and delete the columns that are not written to the database
        102 reduce_splice_ai <file_in> <output>
            Delete the results of spliceAI and delete the information that is not written into the database
        103 db_add_annovar_gene_name <db_file> <selected_annovar_file> <ref_in>
            Add annovar-annotated function and detail to the database
        104 db_add_vep_gene_name <db_file> <selected_vep_file>
            Add the Gene Feature Consequence of vep annotation to the database
        105 db_add_spliceAI_gene_name <db_file> <reduced_spliceAI_file>
            Add the gene name of the spliceAI annotation to the database
        106 db_add_dmis_gene_name <db_file> <dmis_file>
            Add dmis annotation gene name, HGVSc_ANNOVAR, HGVSp_ANNOVAR information to the database
        107 db_add_dsplicing_geneinfo <db_file> <dsplicing_file>
            Add dsplicing to the database region, gene_name, detailed_consequence, ada_score, rf_score 
        108 db_add_synonymous_annovar_gene_name <db_file> <file_in>
            Add the gene name and detail information of synonymous annovar to the database
        109 db_add_synonymous_vep_gene_name <db_file> <selected_vep_file>
            Add information such as the gene name of synonymous vep to the database
        110 select_gene_name_in_synonymous_tsv <tsv_in> <tsv_out>
            For the synonymous original bystro tsv file, delete the gene name according to refSeq.siteType and refSeq.exonicAlleleFunction
        111 db_add_synonymous_bystro_gene_name <db_file> <bystro_tsv, ref_in>
            Add information such as the gene name of synonymous bystro to the database
        112 check_gene_name <db_file>
            Check the variance table, whether the gene name added in the synonymous_snp table exists in the gene table. If it doesn't exist, print it out
        113 modify_gene_name <db_file> <search_list> <variance_not_in_gene_table>
            Process each row in variance_not_in_gene_table.
            1 If the gene name has multiple records in the gene table, no processing will be done
            2 If the gene name is in the specified reserved set, do not process
            3 If the gene name is not in the search_list. If there is no comment, you need to set the corresponding annotator to 0
            4 If the gene name is in the search_list, check it again with the id of the search_list
              If found, update the corresponding gene name in the variance table and synonymous table to be consistent with the gene name in the gene_table;
              If not found, print out
        114 build_contingency_table_gene_name <db_file> <phenotype> <sample_table_name> <sample_restrict>
                                      <gene_table_name> <gene_restrict> <variance_restrict> 
                                      <synonymous_variance_restrict>
                                      <output> <job_id>
            Refer to 46 build_contingency_table, select the variance method of the gene, select from the range box, and change it to the gene name in the comment
            db_file: Only applicable to the database version with the gene name added
            variance_restrict: When selecting annotator, please avoid applying "!="
            job_id: The log file name used for parallel processing. Non-parallel processing, can be set at will
            
            
        115 get_ensembl_mouse_gene_name <file_in> <output>
            Find the gene name corresponding to the gene id from Ensembl Genome Browser, and save it in the last column of the file
        116 db_update_mouse_data <db_file> <gene_table> <mouse_data>
            Adding gene expression data in mouse heart tissue to the database
        117 test_gene_binom_test <gene_set_file> <fet_result> <mode> 
            Use the proportion of the gene set appearing in the gene table as the probability, the number of genes in the gene list selected by FET as n, and the number of n in the gene set as m，
            Perform binom test, print p value
            gene_set_file: It needs to be a gene set with id. See 136
            mode:
            A1: Use B1<1 and A1 >1 in the FET result as a condition to select the gene id
            B1: Use A1<1 and B1 >1 in the FET result as a condition to select the gene id
            p1_control: select gene id with p_value1<=0.05 and odds_ratio1<1 in the FET result
            p1_case: Use p_value1<=0.05 and odds_ratio1>1 in the FET result to select the gene id
        118 fisher_gene_name_binom_test <database> <sample_restrict> <gene_set> <path> <script>
            Use a variety of variance constraint combinations, perform FET, and then perform binomial test
            Submit a task on hpc for every 53 constraint combinations of variance
        119 binom_pileup <path> <output>
            Organize information on all .binom files in the path, pileup into one file
        120 flat_pileup <pileup_in> <out_put> <heart6_binom_c1> <ctd_binom_c> <tof_binom_c> <heart6_binom_c2>
            For the binom pileup of 119, according to the order of heart6_binom1, CTD, tof, heart6_binom2, roll out according to the constraints of the specified binom, and print to the screen. If the binom condition is not met, the data is ignored.
            binom pileup should be sorted in advance
        121 build_contingency_table_gene_name_synonymous <db_file> <phenotype> <sample_table_name> <sample_restrict>
                                      <gene_table_name> <gene_restrict> <variance_restrict> <synonymous_restrict>
                                      <output> <job_id>
        122 fisher_gene_name_binom_test_synonymous <database> <sample_restrict> <gene_set> <path>
            See 118, use 121 for FET
        123 gene_set2gene_set_with_name <db_file> <gene_set_file>
            Convert the gene id in the gene set to a combination of gene name and gene id (name_id)
        124 db_update_gene_table <db_file> <file_in> <key_col>
            Add mouse gene expression data to gene table
            key_col can be gene_id or gene_name
        125 modify_mouse_data <db_file>
            Divide mouse_mean_WT_LV_TPM_3 in the gene table by 3
        126 db_add_gene_info <db_file> <file_in>
            Add the gene comment in the file in to the gene table according to the gene name.
            Requirement: the first column of file in is gene name
        127 selected_flat2variance_info <phenotype> <selected_flat> <fet_path> <db_file> <output>
            Process the filtered flat file and convert it into gene and variance related information
            
        128 selected_flat2sample_info <phenotype> <selected_flat> <fet_path> <db_file> <output>
            Process the filtered flat files and convert them into information about gene, variance and sample
        
        129 db_update_mouse_data2 <db_file> <gene_table> <mouse_data> <target_mouse_gene_file>
            the gene name in db is human gene. human gene name ----> target mouse gene name ----> target mouse data
            Update the MGI_mouse_gene in the database, mouse_mean_WT_LV_TPM_3, mouse_mean_WT_OFT_TPM_2, mouse_mean_WT_PA_TPM_3
            Biomart homology mapping from human gene name to mouse gene name. Then use the homologous gene name of the mouse to search the gene name found on the ensemble and the original gene name.
            
        130 build_contingency_table_gene_name_with_gene_list_file<db_file> <phenotype> 
                                                                 <sample_restrict> <gene_set_file> <variance_restrict>
                                                                 <output>
            Analyze in units of gene set, use annotation to determine whether variance belongs to gene
            gene_set_file: It needs to be the version with gene name, the format is
            set name1    name1_id1    name2_id2
            set name2    name3_id3    name4_id4
        
        131 analyze_geneset_fisher_test_v3 <database> <sample_restrict> <gene_set_file>  <path>
            A variety of variance constraints and phenotype combinations are used, and the gene set is used as the analysis unit. Annotation is used to judge whether the variance belongs to a gene, and a FET test is performed.
            
        132 db_update_mouse_data3 <db_file> <gene_table> <mouse_data> <mouse2human_file>
            Find the homologous human gene name corresponding to the mouse gene name through ncbi, and then add the mouse data to the database
            mouse_data: mouse data
            mouse2human_file: The gene mapping relationship between mice and humans detected by ncbi
        133 db_add_human_gene_set2gene_table <db_file> <human_gene_set_file>
            Add columns to the gene table according to the human gene set. Each gene set is a column. The data represents the inclusion of the gene in this gene set况。
            human_gene_set_file: It needs to be a gene set with all gene symbols or all gene ids
        134 db_add_mouse_gene_set2gene_table <db_file> <mouse_set_file> <mouse2human>
            mouse2human: A mapping file of mouse to human gene names. The format is as follows:
            mouse_gene_name    human_official_gene_name    human_also_known_as
            Prrx1\tPRRX1\tPMX1;PRX1;AGOTC;PHOX1;PRX-1
        135 transform_mouse2human_gene_set <mouse_gene_set_file> <mapping_file> <output>
            Convert mouse gene set to human gene set. Because of the former name, there will be many more genes in each gene set after conversion.
        
        136 gene_set2gene_set_with_id <db_file> <gene_set_file>
            Convert mouse gene set to human gene set. Because of the former name, there will be many more genes in each gene set after conversion.
        
        137 fet_pileup <selected_list_file> <fet_path> <threshold> <mode>
            Pileup the qualified lines in the selected FET file together
            
        138 (discarded) variance_restrict2variance_info_sample_info <db_file> <variance_restrict> <variance_info_out> <sample_info_out>
            According to the constraints of variance, export the variance information and the information of the sample with these variances
            
        139 variance_restrict2variance_info <db_file> <variance_restrict> <variance_info_out> [variance_table]
            According to the constraint of variance, derive the relevant information of gene and variance
            variance_table: variance   or   PE_GATK_variant (default: variance)
            
        140 variance_restrict2sample_info <db_file> <variance_restrict> <sample_info_out> [variance_table]
            According to the constraints of variance, export the relevant information of gene, variance, and sample
            variance_table: variance   or   PE_GATK_variant (default: variance)
            
        141 add_exon_count2file <raw_file> <gene_symbol_col> <excon_count_col> <output>
            By going to the ncbi website to query, grab the website information, and update the excon count in the file
            raw_file file containing gene symbol
            gene_symbol_col gene The column where the symbol is located. Symbol separated by ","
            excon_count_col excon The column where count is located. count is separated by ",", consistent with the symbol order
            output output file
            
        142 add_gene_size2file <raw_file> <gene_symbol_col> <gene_size_col> <db_file> <output>
            Update the gene size in the file. The gene size is defined as: the absolute value of the difference between start_pos and end_pos in the database.
        143 build_contingency_table_gene_name_cnv <variance_table> <db_file> <phenotype>
                                          <sample_table_name> <sample_restrict>
                                          <gene_table_name> <gene_restrict>
                                          <variance_restrict> <synonymous_restrict>
                                          <output> <job_id> <hg38_cnv_file>
            hg38_cnv_file: chrom start end sample_id heart6 CTD
            Add the A1 B1 C1 D1 in the original 114 to the ABCD in the cnv analysis, and then calculate the fisher test 
        144 analyze_CTD_gene_fisher_with_cnv <database> <hg38_cnv_file> <path>
            With a variety of variance combinations, the gene is used as the analysis unit, and 143 FETs with cnv are used for analysis.
            hg38_cnv_file：See 143
            path：output path
        145 analyze_variant_in_cnv <hg38_cnv_file> <database>
            Analyze the number of variance in cnv with different variance constraints
        
        146 cnv_overlap_gene <hg38_cnv_file> <db_file> <output>
            Find the intersection of the cnv of the hg38 version and the gene in the database
        147 export_variance_info_with_gene_list_gene_name <db_file> <gene_list_file> <output_variance> <output_sample>
           Extract variance information and sample information according to the gene list.
           Extract the variance in all genes in the gene list, and all sample information containing these variances.
           Only extract the variance table, use annotation to judge whether the variance belongs to a gene, and add the number of people with mutations in the ctd case group and the number of people with mutations in the heart6 case control group after the variance information
           The gene list can be a row or a column, it can be a gene id or a gene name, and it can be without or with a header.
           
        148 grep_gatk_variance <gatk_path> <variant_sample_info_file> <output_path> <mode> [window]
            Find the specified variance from gatk's whole genome sequencing variance (hg19)
            gatk_path stores the path of gatk's vcf
            variant_sample_info_file 指定variance，格式如下：
                SLID.vcf.file   chr     pos_hg19        ref     alt     gen_id
                SL204704.GATK-4.beta.6.vcf      2       215865477       A       G       1678
            output_path Path to store results
            mode 0 Use grep to match chrom and pos
                 1 Use awk to match chrom and +— window pos
            window The default value is 100
        
        149 check_PE_GATK_overlap <slid2genid_file> <gatk_path> <db_file> <output>
            Extract all the variance of hg38 from the variance table of the database, convert it to hg19, then query the genid of the variance, convert it to slid,
            Go to the gatk folder to find whether the corresponding vcf file contains the hg19 version of the variance
            The results are saved in the following format:
            chr\thg38_pos\thg19_pos\tref\talt\tgenid\tslid\tin_gatk\treason\tgatk_pos\tgatk_ref\tgatk_alt\n
        150 chek_PE_GATK_multiple <num> <db_file> <table> <slid2genid_file> <gatk_path> <output_path>
            The parallel version of 149 needs to run on hpc
        
        151 liftover_PE_variance2hg19 <PE_variance> <output> <unmapped>
            Liftover the PE variance of hg38 to the version of hg1
            
        152 person_split_PE_variance <PE_variance> <output_path> <id_pairs>
            Split PE variance into one file per person, and save it with slid as the file suffix
        
        153 analyze_PE_GATK_overlap <hg19_pe_vcf> <hg19_gatk_vcf> <outpath> <mode>
            Analyze the overlap between the variance of PE and the variance of GATK
            GATK's variance does not need split, leftnorm
            The variance of PE has already reduced the redundant alt during the person split, adjusted the number of gt, and also done leftnorm (see 154)
            mode 0 Don't consider gt, only chr pos ref alt
                 1 Consider chr pos ref alt gt 5 factors
            output
                pe_only_out     unique to pe
                gatk_only_out   unique to gatk
                pe_in_gatk_out  The variance shared by pe and gatk, the data of pe
                gatk_in_pe_out  The variance shared by pe and gatk, the data of gatk
        154 vcf_reduce_PE_alt
            cat pe.vcf | wgsa.py vcf_reduce_PE_alt > reduced.pe.vcf
        
        155 cal_concordance <path> <slid>
        156 gene_size_analyse <database> <output_table> <gene_list> <input_list>
            Divide all genes in the database into two groups, the group in gene_list and the group not in gene_list. Save in form. Table Format Compatibility 50 t_test.
            output_table: The table format is as follows:
            #gene_id    gene_name    distance    gene_list
            
            gene_list: The format is as follows
            gene_list_name1    name1_id1    name2_id2    name3_id3
            gene_list_name2    name4_id4    name5_id5    name6_id6
            After the run is complete, use 50 for t test
        157 grep_human_gencode_gtf3 <input_gtf3> <output> <db_file>
            Organize the gtf3 format file downloaded by gencode into the following format and print it out:
            chr    start    end    type    source    gene_id    gene_name    transcript_name    exon_number    exon_id
            
            output: The output format is
            gene_name    exon_number    total_distence
            The exon number here is the number of exon merged
        158 db_build_PE_GATK_variants_table <db_file> <overlap_result>
            Create a new PE_GATK_variant table based on the variance table. The data structure of the new table is consistent with the variance table. Whether each sample contains the variance, and compare it with the GATK data.
            If GATK also has the variance, keep the PE data. If GATK does not contain the variance, set the PE data to 0
            If it is an indel, consider whether there is an indel within the range of 100bp above and below the corresponding position of GATK. If an indel also exists, it is considered that the indel also exists in GATK.
            overlap_result: Analysis results of 150
        159 variant_table_info <db_file> <table_name>
            Analyze the indel number, snp number, indel sample pair number, snp sample pair number in a variant table
            
        160 select_candidate_case_from_binom_pileup_sort <pileup_file> <binom_case_c> <topN>
            Select good results from already sorted binom pileup files. Given a binom case condition (such as A2.B0) CTD, heart6, p2 in the A2.B0 row of tof
            It is best to be sunny. And the p2 in row A0B2 of heart6 is preferably not positive. Choose this result..
            binom_case_c can be："A2.B0" "A3.2.B1.0" "A3.B0" "A3.B1"
            topN：The number of best results
        
        161 fet_binomc2sample_variance_info <db_file> <fet_file> <binom_c> <variance_table> <sample_info_out> <variance_info_out>
            Extract sample and variance information according to the FET file and binom conditions (such as A2.B0).
        
        162 fet_genelist2sample_variance_info <db_file> <fet_file> <gene_list> <variance_table> <sample_info_out> <variance_info_out>
            Extract sample and variance information according to the FET file and gene list.
            
        163 combine_gatk <gatk_path> <output> <target_chr>
            Combine gatk results into one file
        
        164 export_geneset <db_file> <geneset_out>
            Export all gene sets (only gene name) from the gene_table of the database
            
            """.format(sys.argv[0]))
        exit(0)
    if sys.argv[1] == "filter1":
        filter1(float(sys.argv[2]), sys.argv[3], sys.argv[4])
    elif sys.argv[1] == "data_reduce":
        data_reduce(sys.argv[2], sys.argv[3], sys.argv[4])
    elif sys.argv[1] == "get_words":
        get_words(sys.argv[2], sys.argv[3])
    elif sys.argv[1] == "reduce_words":
        reduce_words(sys.argv[2], sys.argv[3], sys.argv[4])
    elif sys.argv[1] == "filter_by_word":
        filter_by_word(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
    elif sys.argv[1] == "reclassify_ALT":
        reclassify_ALT(sys.argv[2], sys.argv[3])
    elif sys.argv[1] == "filter_alt":
        filter_alt(sys.argv[2], sys.argv[3], sys.argv[4])
    elif sys.argv[1] == "left_normalization" and len(sys.argv) == 4:
        left_normalization(sys.argv[2], sys.argv[3])
    elif sys.argv[1] == "reduce_duplicate":
        reduce_duplicate(sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]),
                         int(sys.argv[6]), sys.argv[7], sys.argv[8])
    elif sys.argv[1] == "rebuild_vcf":
        rebuild_vcf(sys.argv[2], sys.argv[3])
    elif sys.argv[1] == "rebuild_vcf_bystro":
        rebuild_vcf_bystro(sys.argv[2], sys.argv[3], sys.argv[4], False)
    elif sys.argv[1] == "rebuild_vcf_vep":
        rebuild_vcf_vep(sys.argv[2], sys.argv[3], sys.argv[4] in ["true", "True", "TRUE"], sys.argv[5])
    elif sys.argv[1] == "rebuild_vcf_vep2":
        rebuild_vcf_vep2(sys.argv[2], sys.argv[3], sys.argv[4] in ["true", "True", "TRUE"], sys.argv[5])
    elif sys.argv[1] == "check_id":
        check_id(sys.argv[2], sys.argv[3])
    elif sys.argv[1] == "UI":
        union_intersection(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
    elif sys.argv[1] == "split_vcf" and len(sys.argv) == 4:
        split_vcf(sys.argv[2], sys.argv[3])
    elif sys.argv[1] == "split_vcf" and len(sys.argv) == 2:
        split_vcf(".", ".", 1)
    elif sys.argv[1] == "left_normalization" and len(sys.argv) == 2:
        left_normalization(".", ".", 1)
    elif sys.argv[1] == "add_colums":
        add_colums(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7], sys.argv[8])
    elif sys.argv[1] == "splice_ai_filter" and len(sys.argv) == 5:
        splice_ai_filter(sys.argv[2], sys.argv[3], sys.argv[4])
    elif sys.argv[1] == "splice_ai_filter" and len(sys.argv) == 3:
        splice_ai_filter(".", sys.argv[2], ".", 1)
    elif sys.argv[1] == "cut" and len(sys.argv) == 7:
        cut(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6])
    elif sys.argv[1] == "cut" and len(sys.argv) == 5:
        cut(".", sys.argv[2], sys.argv[3], sys.argv[4], ".", 1)
    elif sys.argv[1] == "add_vcf_head":
        add_vcf_head(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6])
    elif sys.argv[1] == "select_vcf" and len(sys.argv) == 5:
        select_vcf(sys.argv[2], sys.argv[3], sys.argv[4])
    elif sys.argv[1] == "select_vcf" and len(sys.argv) == 3:
        select_vcf(".", sys.argv[2], ".", 1)
    elif sys.argv[1] == "component":
        component(sys.argv[2], sys.argv[3], sys.argv[4] in ["true", "True", "TRUE"])
    elif sys.argv[1] == "splice_ai_filter_indel" and len(sys.argv) == 5:
        splice_ai_filter_indel(sys.argv[2], sys.argv[3], sys.argv[4], 0)
    elif sys.argv[1] == "splice_ai_filter_indel" and len(sys.argv) == 3:
        splice_ai_filter_indel(".", sys.argv[2], ".", 1)
    elif sys.argv[1] == "select_by_cols":
        if len(sys.argv) == 8:
            select_by_cols(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7])
        elif len(sys.argv) == 7:
            select_by_cols(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], None)
    elif sys.argv[1] == "splice_ai_merge_snp_indel":
        splice_ai_merge_snp_indel(sys.argv[2], sys.argv[3], sys.argv[4])
    elif sys.argv[1] == "build_vcf_index":
        build_vcf_index(sys.argv[2])
    elif sys.argv[1] == "create_db":
        create_db(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7], sys.argv[8],
                  sys.argv[9], sys.argv[10], sys.argv[11])
    elif sys.argv[1] == "db_add_sample2table_variance":
        db_add_sample2table_variance(sys.argv[2], sys.argv[3], sys.argv[4])
    elif sys.argv[1] == "check_vcf_missing":
        check_vcf_missing(sys.argv[2])
    elif sys.argv[1] == "export_vcf_from_db":
        # print sys.argv
        export_vcf_from_db(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7])
    elif sys.argv[1] == "db_add_spliceAI_anno":
        db_add_spliceAI_anno(sys.argv[2], sys.argv[3])
    elif sys.argv[1] == "db_add_spidex_anno":
        db_add_spidex_anno(sys.argv[2], sys.argv[3])
    elif sys.argv[1] == "db_add_vcf_pos":
        db_add_vcf_pos(sys.argv[2], sys.argv[3])
    elif sys.argv[1] == "db_add_bystro_anno":
        db_add_bystro_anno(sys.argv[2], sys.argv[3], sys.argv[4])
    elif sys.argv[1] == "check_dbNSFP":
        check_dbNSFP(sys.argv[2])
    elif sys.argv[1] == "check_protein_coding_gene":
        check_protein_coding_gene(sys.argv[2])
    elif sys.argv[1] == "check_protein_coding_transcript":
        check_protein_coding_transcript(sys.argv[2])
    elif sys.argv[1] == "check_dms":
        check_dms(sys.argv[2])
    elif sys.argv[1] == "split_gene_centric_annotation":
        split_gene_centric_annotation(sys.argv[2], sys.argv[3])
    elif sys.argv[1] == "db_build_synonymous_snp_table":
        db_build_synonymous_snp_table(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7],
                                      sys.argv[8])
    elif sys.argv[1] == "db_build_gene_table":
        db_build_gene_table(sys.argv[2], sys.argv[3], sys.argv[4])
    elif sys.argv[1] == "db_copy_table":
        db_copy_table(sys.argv[2], sys.argv[3], sys.argv[4])
    elif sys.argv[1] == "db_add_gene_anno":
        if len(sys.argv) == 6:
            db_add_gene_anno(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], "tab")
        elif len(sys.argv) == 7:
            db_add_gene_anno(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6])
        else:
            print("wrong argument")
    elif sys.argv[1] == "db_reduce_col_duplicate":
        db_reduce_col_duplicate(sys.argv[2], sys.argv[3], sys.argv[4])
    elif sys.argv[1] == "collapse":
        collapse(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7])
    elif sys.argv[1] == "collapse_new":
        collapse_new(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7], sys.argv[8])
    elif sys.argv[1] == "build_contingency_table":
        # print sys.argv
        build_contingency_table(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6])
    elif sys.argv[1] == "build_contingency_table_new":
        build_contingency_table_new(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6],
                                    sys.argv[7], sys.argv[8], sys.argv[9], sys.argv[10], sys.argv[11],
                                    sys.argv[12])
    elif sys.argv[1] == "filter_collapse_and_vcf_with_gene_list":
        filter_collapse_and_vcf_with_gene_list(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6])
    elif sys.argv[1] == "get_variance_num_in_sample":
        get_variance_num_in_sample(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7])
    elif sys.argv[1] == "get_variance_function_number":
        get_variance_function_number(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6])
    elif sys.argv[1] == "t_test":
        test_t(sys.argv[2], sys.argv[3], sys.argv[4])
    elif sys.argv[1] == "split_vcf_with_sample":
        split_vcf_with_sample()
    elif sys.argv[1] == "cross_genehancer_snp":
        cross_genehancer_snp(sys.argv[2], sys.argv[3], sys.argv[4])
    elif sys.argv[1] == "cross_genehancer_snp_worker":
        cross_genehancer_snp_worker(sys.argv[2], sys.argv[3], sys.argv[4])
    elif sys.argv[1] == "cross_genehancer_snp_multiple":
        cross_genehancer_snp_multiple(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
    elif sys.argv[1] == "check_log":
        check_log(sys.argv[2])
    elif sys.argv[1] == "test_ks":
        test_ks(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6])
    elif sys.argv[1] == "test_MWW":
        test_MWW(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6])
    elif sys.argv[1] == "get_geneset_variance_num_in_sample":
        get_geneset_variance_num_in_sample(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7],
                                           sys.argv[8])
    elif sys.argv[1] == "test_GLM":
        if len(sys.argv) == 8:
            test_GLM(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7])
        elif len(sys.argv) == 7:
            test_GLM(sys.argv[2], sys.argv[3], sys.argv[4], "Poisson", sys.argv[5], sys.argv[6])
        else:
            print("wrong argument")
    elif sys.argv[1] == "check_pvalue_in_table":
        check_pvalue_in_table(sys.argv[2], sys.argv[3])
        # check_pvalue_in_table_new(sys.argv[2], sys.argv[3])
    elif sys.argv[1] == "bar_hist_plot":
        bar_hist_plot(sys.argv[2], sys.argv[3])
    elif sys.argv[1] == "analyze_fisher_test_variance":
        analyze_fisher_test_variance(sys.argv[2], sys.argv[3])
    elif sys.argv[1] == "whole_genome_variance_test":
        whole_genome_variance_test(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7])
    elif sys.argv[1] == "analyze_fisher_test_variance_new":
        analyze_fisher_test_variance_new(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6])
    elif sys.argv[1] == "cut_org_exported_vcf_samples":
        cut_org_exported_vcf_samples_new(sys.argv[2], sys.argv[3])
    elif sys.argv[1] == "map_vcf_sample_id":
        map_vcf_sample_id(sys.argv[2], sys.argv[3], sys.argv[4])
    elif sys.argv[1] == "cut_org_exported_vcf_samples_multiple":
        cut_org_exported_vcf_samples_multiple(sys.argv[2], sys.argv[3], sys.argv[4])
    elif sys.argv[1] == "rerun_selected_combination":
        if len(sys.argv) == 9:
            rerun_selected_combination(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7],
                                       sys.argv[8], 0)
        elif len(sys.argv) == 10:
            rerun_selected_combination(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7],
                                       sys.argv[8], sys.argv[9])
    elif sys.argv[1] == "export_vcf_with_gene_list":
        export_vcf_with_gene_list(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6])
    elif sys.argv[1] == "db_add_ccrs_limbr":
        db_add_ccrs_limbr(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6])
    elif sys.argv[1] == "analyze_fisher_test_variance_2":
        analyze_fisher_test_variance_2(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6])
    elif sys.argv[1] == "permutation_burden":
        permutation_burden(sys.argv[2], sys.argv[3], sys.argv[4])
    elif sys.argv[1] == "permutation_fisher":
        permutation_fisher(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6],
                           sys.argv[7], sys.argv[8], sys.argv[9], sys.argv[10], sys.argv[11],
                           sys.argv[12])
    elif sys.argv[1] == "analyze_fisher_test_variance_3":
        analyze_fisher_test_variance_3(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6])
    elif sys.argv[1] == "modify_gene_name2id":
        modify_gene_name2id(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6])
    elif sys.argv[1] == "modify_vcfid":
        modify_vcfid(sys.argv[2], sys.argv[3], sys.argv[4])
    elif sys.argv[1] == "add_gene_name":
        add_gene_name(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6])
    elif sys.argv[1] == "modify_duplicateID":
        modify_duplicateID(sys.argv[2], sys.argv[3])
    elif sys.argv[1] == "manhattan_direct":
        manhattan_direct(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7], sys.argv[8],
                         sys.argv[9])
    elif sys.argv[1] == "get_gene_set_info":
        get_gene_set_info(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7])
    elif sys.argv[1] == "copy_vcf2table":
        copy_vcf2table(sys.argv[2], sys.argv[3], sys.argv[4])
    elif sys.argv[1] == "db_add_ccds":
        db_add_ccds(sys.argv[2], sys.argv[3], sys.argv[4])
    elif sys.argv[1] == "firth_logistic_regression":
        firth_logistic_regression(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7])
    elif sys.argv[1] == "firth_logistic_regression_with_gene":
        firth_logistic_regression_with_gene(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6],
                                            sys.argv[7], sys.argv[8], sys.argv[9], sys.argv[10])
    elif sys.argv[1] == "firth_logistic_regression_genelist":
        firth_logistic_regression_genelist(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6],
                                           sys.argv[7], sys.argv[8], sys.argv[9], sys.argv[10], sys.argv[11],
                                           sys.argv[12])
    elif sys.argv[1] == "firth_logistic_regression_plot":
        firth_logistic_regression_plot(sys.argv[2])
    elif sys.argv[1] == "build_contingency_table_new_with_gene_list_file":
        build_contingency_table_new_with_gene_list_file(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6],
                                                        sys.argv[7], sys.argv[8], sys.argv[9], sys.argv[10])
    elif sys.argv[1] == "variance_distribution_plot":
        variance_distribution_plot(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7],
                                   sys.argv[8], sys.argv[9], sys.argv[10], sys.argv[11], sys.argv[12])
    elif sys.argv[1] == "select_enriched_variance":
        select_enriched_variance(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7],
                                 sys.argv[8], sys.argv[9], sys.argv[10], sys.argv[11], sys.argv[12], sys.argv[13],
                                 sys.argv[14])
    elif sys.argv[1] == "gene_distribution_plot":
        gene_distribution_plot(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7],
                               sys.argv[8], sys.argv[9], sys.argv[10], sys.argv[11], sys.argv[12])
    elif sys.argv[1] == "qqplot_fisher_permutation":
        qqplot_fisher_permutation(sys.argv[2], sys.argv[3], sys.argv[4])
    elif sys.argv[1] == "gene_set_fisher_test":
        gene_set_fisher_test(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6],
                             sys.argv[7], sys.argv[8], sys.argv[9], sys.argv[10], sys.argv[11])
    elif sys.argv[1] == "gene_set_fisher_test2":
        gene_set_fisher_test2(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6],
                              sys.argv[7], sys.argv[8], sys.argv[9], sys.argv[10], sys.argv[11],
                              sys.argv[12])
    elif sys.argv[1] == "analyze_fisher_test_variance_May19":
        analyze_fisher_test_variance_May19(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
    elif sys.argv[1] == "analyze_geneset_fisher_test":
        analyze_geneset_fisher_test(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6])
    elif sys.argv[1] == "fisher_test_synonymous":
        fisher_test_synonymous(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
    elif sys.argv[1] == "multiple_geneset_variance_num_in_sample":
        multiple_geneset_variance_num_in_sample(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6],
                                                sys.argv[7], sys.argv[8])
    elif sys.argv[1] == "prepare_liftover_bed":
        prepare_liftover_bed(sys.argv[2], sys.argv[3], sys.argv[4])
    elif sys.argv[1] == "gether_bed_region":
        gether_bed_region(sys.argv[2], sys.argv[3])
    elif sys.argv[1] == "check_gene_overlap":
        check_gene_overlap(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
    elif sys.argv[1] == "db_add_bystro_gene_name":
        db_add_bystro_gene_name(sys.argv[2], sys.argv[3], sys.argv[4])
    elif sys.argv[1] == "select_gene_name_in_tsv":
        select_gene_name_in_tsv(sys.argv[2], sys.argv[3])
    elif sys.argv[1] == "select_gene_name_in_annovar_result":
        select_gene_name_in_annovar_result(sys.argv[2], sys.argv[3])
    elif sys.argv[1] == "reduce_dsplicing_result":
        reduce_dsplicing_result(sys.argv[2], sys.argv[3])
    elif sys.argv[1] == "reduce_vep_result":
        reduce_vep_result(sys.argv[2], sys.argv[3], sys.argv[4])
    elif sys.argv[1] == "reduce_splice_ai":
        reduce_splice_ai(sys.argv[2], sys.argv[3])
    elif sys.argv[1] == "db_add_annovar_gene_name":
        db_add_annovar_gene_name(sys.argv[2], sys.argv[3], sys.argv[4])
    elif sys.argv[1] == "db_add_vep_gene_name":
        db_add_vep_gene_name(sys.argv[2], sys.argv[3])
    elif sys.argv[1] == "db_add_spliceAI_gene_name":
        db_add_spliceAI_gene_name(sys.argv[2], sys.argv[3])
    elif sys.argv[1] == "db_add_dmis_gene_name":
        db_add_dmis_gene_name(sys.argv[2], sys.argv[3])
    elif sys.argv[1] == "db_add_dsplicing_geneinfo":
        db_add_dsplicing_geneinfo(sys.argv[2], sys.argv[3])
    elif sys.argv[1] == "db_add_synonymous_annovar_gene_name":
        db_add_synonymous_annovar_gene_name(sys.argv[2], sys.argv[3])
    elif sys.argv[1] == "db_add_synonymous_vep_gene_name":
        db_add_synonymous_vep_gene_name(sys.argv[2], sys.argv[3])
    elif sys.argv[1] == "db_add_synonymous_bystro_gene_name":
        db_add_synonymous_bystro_gene_name(sys.argv[2], sys.argv[3], sys.argv[4])
    elif sys.argv[1] == "select_gene_name_in_synonymous_tsv":
        select_gene_name_in_synonymous_tsv(sys.argv[2], sys.argv[3])
    elif sys.argv[1] == "check_gene_name":
        check_gene_name(sys.argv[2])
    elif sys.argv[1] == "modify_gene_name":
        modify_gene_name(sys.argv[2], sys.argv[3], sys.argv[4])
    elif sys.argv[1] == "build_contingency_table_gene_name":
        build_contingency_table_gene_name(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6],
                                          sys.argv[7], sys.argv[8], sys.argv[9], sys.argv[10], sys.argv[11])
    elif sys.argv[1] == "select_enriched_variance_gene_name":
        select_enriched_variance_gene_name(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5],
                                           sys.argv[6], sys.argv[7], sys.argv[8], sys.argv[9],
                                           sys.argv[10], sys.argv[11], sys.argv[12], sys.argv[13])
    elif sys.argv[1] == "get_ensembl_mouse_gene_name":
        get_ensembl_mouse_gene_name(sys.argv[2], sys.argv[3])
    elif sys.argv[1] == "db_update_mouse_data":
        db_update_mouse_data(sys.argv[2], sys.argv[3], sys.argv[4])
    elif sys.argv[1] == "test_gene_binom_test":
        test_gene_binom_test(sys.argv[2], sys.argv[3], sys.argv[4])
    elif sys.argv[1] == "fisher_gene_name_binom_test":
        fisher_gene_name_binom_test_slurm(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6])
    elif sys.argv[1] == "binom_pileup":
        binom_pileup(sys.argv[2], sys.argv[3])
    elif sys.argv[1] == "flat_pileup":
        flat_pileup(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7])
    elif sys.argv[1] == "build_contingency_table_gene_name_synonymous":
        build_contingency_table_gene_name_synonymous(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6],
                                                     sys.argv[7], sys.argv[8], sys.argv[9], sys.argv[10], sys.argv[11])
    elif sys.argv[1] == "fisher_gene_name_binom_test_synonymous":
        fisher_gene_name_binom_test_synonymous(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
    elif sys.argv[1] == "gene_set2gene_set_with_name":
        gene_set2gene_set_with_name(sys.argv[2], sys.argv[3])
    elif sys.argv[1] == "db_update_gene_table":
        db_update_gene_table(sys.argv[2], sys.argv[3], sys.argv[4])
    elif sys.argv[1] == "modify_mouse_data":
        modify_mouse_data(sys.argv[2])
    elif sys.argv[1] == "db_add_gene_info":
        db_add_gene_info(sys.argv[2], sys.argv[3])
    elif sys.argv[1] == "selected_flat2variance_info":
        selected_flat2variance_info(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6])
    elif sys.argv[1] == "selected_flat2sample_info":
        selected_flat2sample_info(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6])
    elif sys.argv[1] == "db_update_mouse_data2":
        db_update_mouse_data2(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
    elif sys.argv[1] == "build_contingency_table_gene_name_with_gene_list_file":
        build_contingency_table_gene_name_with_gene_list_file(sys.argv[2], sys.argv[3], sys.argv[4],
                                                              sys.argv[5], sys.argv[6], sys.argv[7])
    elif sys.argv[1] == "analyze_geneset_fisher_test_v3":
        analyze_geneset_fisher_test_v3(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
    elif sys.argv[1] == "db_update_mouse_data3":
        db_update_mouse_data3(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
    elif sys.argv[1] == "db_add_human_gene_set2gene_table":
        db_add_human_gene_set2gene_table(sys.argv[2], sys.argv[3])
    elif sys.argv[1] == "db_add_mouse_gene_set2gene_table":
        db_add_mouse_gene_set2gene_table(sys.argv[2], sys.argv[3], sys.argv[4])
    elif sys.argv[1] == "transform_mouse2human_gene_set":
        transform_mouse2human_gene_set(sys.argv[2], sys.argv[3], sys.argv[4])
    elif sys.argv[1] == "gene_set2gene_set_with_id":
        gene_set2gene_set_with_id(sys.argv[2], sys.argv[3])
    elif sys.argv[1] == "fet_pileup":
        fet_pileup(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
    elif sys.argv[1] == 'variance_restrict2variance_info':
        if len(sys.argv) == 5:
            variance_restrict2variance_info(sys.argv[2], sys.argv[3], sys.argv[4])
        elif len(sys.argv) == 6:
            variance_restrict2variance_info(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
        else:
            print('wrong argument!')
    elif sys.argv[1] == 'variance_restrict2sample_info':
        if len(sys.argv) == 5:
            variance_restrict2sample_info(sys.argv[2], sys.argv[3], sys.argv[4], "variance")
        elif len(sys.argv) == 6:
            variance_restrict2sample_info(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
    elif sys.argv[1] == 'add_exon_count2file':
        add_exon_count2file(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
    elif sys.argv[1] == 'add_gene_size2file':
        add_gene_size2file(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6])
    elif sys.argv[1] == 'build_contingency_table_gene_name_cnv':
        build_contingency_table_gene_name_cnv(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6],
                                              sys.argv[7], sys.argv[8], sys.argv[9], sys.argv[10], sys.argv[11],
                                              sys.argv[12], sys.argv[13])
    elif sys.argv[1] == 'analyze_CTD_gene_fisher_with_cnv':
        analyze_CTD_gene_fisher_with_cnv(sys.argv[2], sys.argv[3], sys.argv[4])
    elif sys.argv[1] == "analyze_variant_in_cnv":
        analyze_variant_in_cnv(sys.argv[2], sys.argv[3])
    elif sys.argv[1] == 'variance_in_cnv_sample':
        variance_in_cnv_sample(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
    elif sys.argv[1] == "cnv_overlap_gene":
        cnv_overlap_gene(sys.argv[2], sys.argv[3], sys.argv[4])
    elif sys.argv[1] == 'export_variance_info_with_gene_list_gene_name':
        export_variance_info_with_gene_list_gene_name(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
    elif sys.argv[1] == 'grep_gatk_variance':
        if len(sys.argv) == 6:
            grep_gatk_variance(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
        elif len(sys.argv) == 7:
            grep_gatk_variance(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6])
    elif sys.argv[1] == 'check_PE_GATK_overlap':
        check_PE_GATK_overlap(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
    elif sys.argv[1] == 'check_PE_GATK_worker':
        check_PE_GATK_worker(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7], sys.argv[8])
    elif sys.argv[1] == 'chek_PE_GATK_multiple':
        chek_PE_GATK_multiple(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7])
    elif sys.argv[1] == 'liftover_PE_variance2hg19':
        liftover_PE_variance2hg19(sys.argv[2], sys.argv[3], sys.argv[4])
    elif sys.argv[1] == 'person_split_PE_variance':
        person_split_PE_variance(sys.argv[2], sys.argv[3], sys.argv[4])
    elif sys.argv[1] == 'analyze_PE_GATK_overlap':
        analyze_PE_GATK_overlap(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
    elif sys.argv[1] == 'vcf_reduce_PE_alt':
        vcf_reduce_PE_alt()
    elif sys.argv[1] == 'cal_concordance':
        cal_concordance(sys.argv[2], sys.argv[3])
    elif sys.argv[1] == 'gene_size_analyse':
        gene_size_analyse(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
    elif sys.argv[1] == 'grep_human_gencode_gtf3':
        grep_human_gencode_gtf3(sys.argv[2], sys.argv[3], sys.argv[4])
    elif sys.argv[1] == 'db_build_PE_GATK_variants_table':
        db_build_PE_GATK_variants_table(sys.argv[2], sys.argv[3])
    elif sys.argv[1] == 'variant_table_info':
        variant_table_info(sys.argv[2], sys.argv[3])
    elif sys.argv[1] == "select_candidate_case_from_binom_pileup_sort":
        select_candidate_case_from_binom_pileup_sort(sys.argv[2], sys.argv[3], int(sys.argv[4]))
    elif sys.argv[1] == 'fet_binomc2sample_variance_info':
        fet_binomc2sample_variance_info(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7])
    elif sys.argv[1] == 'fet_genelist2sample_variance_info':
        fet_genelist2sample_variance_info(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7])
    elif sys.argv[1] == 'combine_gatk':
        combine_gatk(sys.argv[2], sys.argv[3], sys.argv[4])
    elif sys.argv[1] == 'gatk_split':
        gatk_split()
    elif sys.argv[1] == 'grep_data':
        grep_data(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6])
    elif sys.argv[1] == 'combine_gatk3':
        combine_gatk3(sys.argv[2], sys.argv[3])
    elif sys.argv[1] == 'export_geneset':
        export_geneset(sys.argv[2], sys.argv[3])
    elif sys.argv[1] == "test":
        test()
    else:
        print("wrong argument")
        print(sys.argv)


def test():
    grep_data('/Users/yinjie/PycharmProjects/yingjie/wgsa/rr.vcf', '14', '76736078', 'A', 'AAGAAGAGAAGAAG')
    # combine_gatk2('/Users/yinjie/PycharmProjects/yingjie/wgsa/gatk_test', 'output')
    # combine_gatk_get_key_list('/Users/yinjie/PycharmProjects/yingjie/wgsa/gatk_test', 'output')
    # conn = sqlite3.connect('test.db')
    # cursor = conn.cursor()
    # # cursor.execute('create table test (id int primary key, gene_id varchar(40), gene_name varchr(20), '
    # #                'chr varchr(20), start_pos int, end_pos int, '
    # #                'chr2 varchr(20), start_pos2 int, end_pos2 int)')
    # sql_cmd = "select 1 FROM test  WHERE gene_id='2' and gene_name='2' AND chr='3' LIMIT 1"
    # cursor.execute(sql_cmd)
    # data = cursor.fetchall()
    # cursor.close()
    # conn.commit()
    # conn.close()
    # print(data)


if __name__ == "__main__":
    main()
