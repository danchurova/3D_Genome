#!/usr/bin/env python3
# -*- coding:utf-8 -*-

import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pylab as plt

directory = 'data/K562_interchromosomal/100kb_resolution_interchromosomal'
frequency = {}
total = 0
data_type = sys.argv[1] # raw, KRnorm, VCnorm or SQRTVCnorm
percent = sys.argv[2]


def main():
    global total
    for chrom_pair in os.listdir(directory):
        if chrom_pair != '.DS_Store':
            def select(name):
                return directory+'/'+chrom_pair+'/MAPQGE/'+name

            chr1, chr2 = chrom_pair.split('_')
            chr_id1 = chr1.replace('chr', '')
            chr_id2 = chr2.replace('chr', '')
            if chr_id1 != "X":
                chr_id1 = int(chr_id1)
            if chr_id2 != "X":
                chr_id2 = int(chr_id2)

            if data_type == 'raw':
                with open(select('chr'+str(chr_id1)+'_'+str(chr_id2)+'_100kb.RAWobserved')) as data:

            elif data_type == 'KRnorm':
                with open(select('chr'+str(chr_id1)+'_'+str(chr_id2)+'_100kb_KRnormalized')) as data:

            elif data_type == 'VCnorm':
                with open(select('chr'+str(chr_id1)+'_'+str(chr_id2)+'_100kb_VCnormalized')) as data:

            elif data_type == 'SQRTVCnorm':
                with open(select('chr'+str(chr_id1)+'_'+str(chr_id2)+'_100kb_SQRTVCnormalized')) as data:

                    for line in data.readlines():
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

    # draw a plot with genomic interchromosomal contact frequency distribution for each datatype
    lists = sorted(frequency.items())
    x, y = zip(*lists)
    plt.plot(x, y)
    plt.title(str(data_type)+'data interchromosomal contact frequency distribution')
    plt.ylabel('Frequency')
    plt.xlabel('Contact value')
    plt.savefig('data/K562_interchromosomal/100kb_resolution_interchromosomal/'+data_type+'plot.png')

    # draft for concatenated data of top interactions
    concatenated_data = pd.DataFrame(columns = ['chrom1', 'coord1', 'chrom2', 'coord2', 'value'])

    for chrom_pair in os.listdir(directory):
        if chrom_pair != '.DS_Store':
            for file in os.listdir(directory + '/' + chrom_pair + '/MAPQGE30/'):
                if file.endswith(".RAWobserved"):
                    with open(directory+'/'+chrom_pair+'/MAPQGE30/' + file) as raw:
                        data = pd.read_table(raw, header=None)
                        data.columns = ['coord1', 'coord2', 'value']
                        data.insert(0, 'chrom1', chr_id1)
                        data.insert(2, 'chrom2', chr_id2)
                        trimmed_data = data[data.value > threshold]
                        concatenated_data = pd.concat([concatenated_data, trimmed_data])
    concatenated_data.to_csv('data/K562_interchromosomal/K562_100kb./100kb_resolution_interchromosomal/'+data_type+'top_interactions.csv')

    #freq_data = pd.DataFrame(frequency, columns = ['value', 'frequency'])


      #  chr1, chr2 = chrom_pair.split('_')
      #  chr_id1 = int(chr1.replace('chr', ''))
      #  chr_id2 = int(chr2.replace('chr', ''))
if __name__ == '__main__':
    main()
