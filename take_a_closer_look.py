#!/usr/bin/env python3
# -*- coding:utf-8 -*-

import os
import sys
import pandas as pd

directory = 'data/K562_interchromosomal'
resolution = sys.argv[1] # 5kb, 10kb and so on
data_type = sys.argv[2] # raw, KRnorm, VCnorm, SQRTVCnorm

def main():
    for top_inter in os.listdir(directory+'/top_interactions'):
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
                                if data_type == 'raw':

                                    def select(chr1, chr2, name):
                                        return directory+'/'+resolution+'kb_resolution_interchromosomal/chr'+ \
                                                chr1+'_chr'+chr2+'/MAPQGE30/chr'+chr1+'_'+chr2+'_'+resolution+'kb'+name

                                    for chrom_pair in os.listdir(directory+'/'+resolution+'kb_resolution_interchromosomal/chr'+ \
                                            current_chrom1+'_chr'+current_chrom2+'/MAPQGE30/'):

                                        with open(select(current_chrom1, current_chrom2,'.RAWobserved')) as raw_higher_resolution:


                            prev_chrom1, prev_chrom2 = current_chrom1, current_chrom2




if __name__ == '__main__':
    main()
