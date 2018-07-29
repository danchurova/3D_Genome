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
percent = sys.argv[2] # 1, 5 and so on


def main():
    global total
    for chrom_pair in os.listdir(directory):
        if chrom_pair != '.DS_Store':
            def select(name):
                return directory+'/'+chrom_pair+'/MAPQGE30/'+name
            print(chrom_pair)
            chr1, chr2 = chrom_pair.split('_')
            chr_id1 = chr1.replace('chr', '')
            chr_id2 = chr2.replace('chr', '')
            if chr_id1 != "X":
                chr_id1 = int(chr_id1)
            if chr_id2 != "X":
                chr_id2 = int(chr_id2)

            if data_type == 'raw':
                with open(select('chr'+str(chr_id1)+'_'+str(chr_id2)+'_100kb.RAWobserved')) as raw:
                    for line in raw.readlines():
                        _, _, value = line.strip().split('\t')
                        value = float(value)
                        frequency[value] = frequency.get(value,0) + 1
                        total += 1

            elif data_type == 'KRnorm':
                with open(select('chr'+str(chr_id1)+'_'+str(chr_id2)+'_100kb_KRnormalized')) as KRnorm:
                    for line in KRnorm.readlines():
                        if 'coord' in line:
                            continue
                        _, _, _, value = line.strip().split('\t')
                        value = float(value)
                        frequency[value] = frequency.get(value,0) + 1
                        total += 1

            elif data_type == 'VCnorm':
                with open(select('chr'+str(chr_id1)+'_'+str(chr_id2)+'_100kb_VCnormalized')) as VCnorm:
                    for line in VCnorm.readlines():
                        if 'coord' in line:
                            continue
                        _, _, _, value = line.strip().split('\t')
                        if value != 'value':
                            value = float(value)
                        frequency[value] = frequency.get(value,0) + 1
                        total += 1

            elif data_type == 'SQRTVCnorm':
                with open(select('chr'+str(chr_id1)+'_'+str(chr_id2)+'_100kb_SQRTVCnormalized')) as SQRTVCnorm:
                    for line in SQRTVCnorm.readlines():
                        if 'coord' in line:
                            continue
                        _, _, _, value = line.strip().split('\t')
                        if value != 'value':
                            value = float(value)
                        frequency[value] = frequency.get(value,0) + 1
                        total += 1
    print("I counted all pairs")
    percentage = round(total*(int(percent)/100))
    print("One percentage is " + str(percentage))
    freq_sum = 0
    prev_element = 0
    for key in reversed(sorted(frequency)):
        freq_sum += frequency[key]
        if freq_sum >= percentage:
            print('This is frequency sum that is bigger than 1%: ' + str(freq_sum))
            threshold = prev_element
            print('This is threshold value: ' + str(threshold))
            break
        prev_element = key

    with open("data/K562_interchromosomal/frequency_stats_100kb_"+data_type, 'w') as freqs:
        for key in reversed(sorted(frequency)):
            freqs.write(str(key)+'\t'+str(frequency[key])+'\n')
        freqs.close()

    print("I finished frequency dict writing")

    # draw a plot with genomic interchromosomal contact frequency distribution for each datatype
    #lists = sorted(frequency.items())
    #x, y = zip(lists)
    #print(x,y)
    #plt.plot(x, y)
    #plt.title(str(data_type)+'data interchromosomal contact frequency distribution')
    #plt.ylabel('Frequency')
    #plt.xlabel('Contact value')
    #plt.savefig('data/K562_interchromosomal/plots/'+data_type+'plot_100kb.png')
    #print("I saved plot with frequency dictribution!")
    # draft for concatenated data of top interactions
    concatenated_data = pd.DataFrame(columns = ['chrom1', 'coord1', 'chrom2', 'coord2', 'value'])

    for chrom_pair in os.listdir(directory):
        if chrom_pair != '.DS_Store':
            print(chrom_pair)
            chr1, chr2 = chrom_pair.split('_')
            chr_id1 = chr1.replace('chr', '')
            chr_id2 = chr2.replace('chr', '')
            if chr_id1 != "X":
                chr_id1 = int(chr_id1)
            if chr_id2 != "X":
                chr_id2 = int(chr_id2)

            if data_type == 'raw':
                with open(select('chr'+str(chr_id1)+'_'+str(chr_id2)+'_100kb.RAWobserved')) as raw:
                    data = pd.read_table(raw, header=None)
                    data.columns = ['coord1', 'coord2', 'value']
                    data.insert(0, 'chrom1', chr_id1)
                    data.insert(2, 'chrom2', chr_id2)
                    trimmed_data = data[data.value > threshold]
                    concatenated_data = pd.concat([concatenated_data, trimmed_data])

            elif data_type == 'KRnorm':
                with open(select('chr'+str(chr_id1)+'_'+str(chr_id2)+'_100kb_KRnormalized')) as KRnorm:
                    data = pd.read_table(KRnorm)
                    data.drop('Unnamed: 0', 1, inplace=True)
                    data.insert(0, 'chrom1', chr_id1)
                    data.insert(2, 'chrom2', chr_id2)
                    trimmed_data = data[data.value > threshold]
                    concatenated_data = pd.concat([concatenated_data, trimmed_data])

            elif data_type == 'VCnorm':
                with open(select('chr'+str(chr_id1)+'_'+str(chr_id2)+'_100kb_VCnormalized')) as VCnorm:
                    data = pd.read_table(VCnorm)
                    data.drop('Unnamed: 0', 1, inplace=True)
                    data.insert(0, 'chrom1', chr_id1)
                    data.insert(2, 'chrom2', chr_id2)
                    trimmed_data = data[data.value > threshold]
                    concatenated_data = pd.concat([concatenated_data, trimmed_data])

            elif data_type == 'SQRTVCnorm':
                with open(select('chr'+str(chr_id1)+'_'+str(chr_id2)+'_100kb_SQRTVCnormalized')) as SQRTVCnorm:
                    data = pd.read_table(SQRTVCnorm)
                    data.drop('Unnamed: 0', 1, inplace=True)
                    data.insert(0, 'chrom1', chr_id1)
                    data.insert(2, 'chrom2', chr_id2)
                    trimmed_data = data[data.value > threshold]
                    concatenated_data = pd.concat([concatenated_data, trimmed_data])

    concatenated_data.to_csv('data/K562_interchromosomal/top_interactions/'+data_type+'_top_interactions_100kb.csv')
    print("Concatenation is complete!")
    #freq_data = pd.DataFrame(frequency, columns = ['value', 'frequency'])


      #  chr1, chr2 = chrom_pair.split('_')
      #  chr_id1 = int(chr1.replace('chr', ''))
      #  chr_id2 = int(chr2.replace('chr', ''))
if __name__ == '__main__':
    main()
