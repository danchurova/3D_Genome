#!/usr/bin/env python3
# -*- coding:utf-8 -*-

import os
import sys
import pandas as pd
import numpy as np

directory = 'data/K562_interchromosomal/100kb_resolution_interchromosomal'
frequency = {}
total = 0
percent = sys.argv[1]

def main():
    global total
    for chrom_pair in os.listdir(directory):
        if chrom_pair != '.DS_Store':
            for file in os.listdir(directory + '/' + chrom_pair + '/MAPQGE30/'):
                if file.endswith(".RAWobserved"):
                    with open(directory+'/'+chrom_pair+'/MAPQGE30/' + file) as raw:
                        for line in raw.readlines():
                            _, _, value = line.strip().split('\t')
                            value = float(value)
                            frequency[value] = frequency.get(value,0) + 1
                            total += 1
    percentage = round(total*(int(percent)/100))
    freq_sum = 0
    prev_element = 0
    for key in reversed(sorted(frequency)):
        freq_sum += frequency[key]
        if freq_sum >= percentage:
            threshold = prev_element
            break
        prev_element = key

    concatenated_data = pd.DataFrame(columns = ['chrom1', 'coord1', 'chrom2', 'coord2', 'value'])

    for chrom_pair in os.listdir(directory):
        if chrom_pair != '.DS_Store':
            for file in os.listdir(directory + '/' + chrom_pair + '/MAPQGE30/'):
                if file.endswith(".RAWobserved"):
                    chr1, chr2 = chrom_pair.split('_')
                    chr_id1 = chr1.replace('chr', '')
                    chr_id2 = chr2.replace('chr', '')
                    if chr_id1 == "X":
                        chr_id1 = 23
                    else:
                        chr_id1 = int(chr_id1)
                    if chr_id2 == "X":
                        chr_id2 = 23
                    else:
                        chr_id2 = int(chr_id2)
                    with open(directory+'/'+chrom_pair+'/MAPQGE30/' + file) as raw:
                        data = pd.read_table(raw, header=None)
                        data.columns = ['coord1', 'coord2', 'value']
                        data.insert(0, 'chrom1', chr_id1)
                        data.insert(2, 'chrom2', chr_id2)
                        trimmed_data = data[data.value > threshold]
                        concatenated_data = pd.concat([concatenated_data, trimmed_data])
    print(concatenated_data)

    #freq_data = pd.DataFrame(frequency, columns = ['value', 'frequency'])
 

      #  chr1, chr2 = chrom_pair.split('_')
      #  chr_id1 = int(chr1.replace('chr', ''))
      #  chr_id2 = int(chr2.replace('chr', ''))
if __name__ == '__main__':
    main()
