#!/usr/bin/env python3
# -*- coding:utf-8 -*-

import os
import sys
import pandas as pd
import numpy as np

directory = 'data/K562_interchromosomal/100kb_resolution_interchromosomal'
frequency = {}
total = 0

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
    #freq_data = pd.DataFrame(frequency, columns = ['value', 'frequency'])
 

      #  chr1, chr2 = chrom_pair.split('_')
      #  chr_id1 = int(chr1.replace('chr', ''))
      #  chr_id2 = int(chr2.replace('chr', ''))
if __name__ == '__main__':
    main()
