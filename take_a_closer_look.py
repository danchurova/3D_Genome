#!/usr/bin/env python3
# -*- coding:utf-8 -*-

import os
import sys
import pandas as pd

directory = 'data/K562_interchromosomal'
initial_resolution = sys.argv[1] # almost always here 100
desired_resolution = sys.argv[2] # in kb 5 -> for 5kb, 10 -> for 10kb and so on
data_type = sys.argv[3] # raw, KRnorm, VCnorm, SQRTVCnorm

def main():
    for top_inter in os.listdir(directory+'/top_interactions'):
        ### getting access to top_interations of 100kb data ###
        if top_inter != '.DS_Store':
            if data_type == 'raw':
                with open(directory+'/top_interactions/raw_top_interactions_100kb.csv') as raw_top_inter:
                    prev_chrom1, prev_chrom2 = 0, 0
                    current_chrom1, current_chrom2 = 0, 0

                    for line in raw_top_inter.readlines():

                        if 'coord' not in line:
                            id_number, chrom1, coord1, chrom2, coord2, value= line.strip().split(',')
                            current_chrom1, current_chrom2 = chrom1, chrom2

                            if prev_chrom1 == 0:
                                prev_chrom1, prev_chrom2 = current_chrom1, current_chrom2

                            while current_chrom1 == prev_chrom1 and current_chrom2 == current_chrom2:
                                ### while we work with the same chromosome_pair ###
                                start1 = int(coord1)
                                stop1 = int(coord1)+1*(int(initial_resolution)*1000)
                                start2 = int(coord2)
                                stop2 = int(coord2)+1*(int(initial_resolution)*1000)
                                print(start1)
                                print(stop1)
                                chrom1_list = [] # we use these lists for collecting lines from higher resolution data files
                                chrom2_list = [] # that are located within top interactions regions
                                if data_type == 'raw':

                                    def select(chr1, chr2, resolution, name):
                                        return directory+'/'+resolution+'kb_resolution_interchromosomal/chr'+ \
                                                chr1+'_chr'+chr2+'/MAPQGE30/chr'+chr1+'_'+chr2+'_'+resolution+'kb'+name

                                    for chrom_pair in os.listdir(directory+'/'+desired_resolution+'kb_resolution_interchromosomal/chr'+ \
                                            current_chrom1+'_chr'+current_chrom2+'/MAPQGE30/'):

                                        with open(select(current_chrom1, current_chrom2, desired_resolution,'.RAWobserved')) as raw_higher_resolution:
                                            ### we open ech chrom_pair_file and get only those lines which are within ###
                                            ### top interactions regions
                                            for line in raw_higher_resolution.readlines():
                                                coord1, coord2, value= line.strip().split('\t')

                                                if int(coord1) >= start1 and int(coord1) <= stop1:
                                                    chrom1_list.append(line)

                                                if int(coord2) >=start2 and int(coord2) <= stop2:
                                                    chrom2_list.append(line)

#                                                print(chrom1_list)
#                                                print(chrom2_list)
                                                common_lines = list(set(chrom1_list).intersection(set(chrom2_list)))
                                                print(common_lines)


                            prev_chrom1, prev_chrom2 = current_chrom1, current_chrom2




if __name__ == '__main__':
    main()
