#!/usr/bin/env python3
# -*- coding:utf-8 -*-

import os
import sys
import pandas as pd
import numpy as np

normalization = sys.argv[1] #KR, VC or SQRTVC
resolution = sys.argv[2] #100, 50, 25, 10, 5
directory = 'data/K562_interchromosomal/'+resolution+'kb_resolution_interchromosomal'

def normalize(raw_data, norm_vector_1, norm_vector_2, resolution):
    # raw_data format = ['coord1', 'coord2', 'value']
    # norm_vector_i format = ['coef']
    # resolution = integer, ex. 100kb -> 100000, 5kb -> 5000
    norm_data_temp = raw_data.copy() # temporary dataset for normalization
    norm_data_temp['coord1'] = (raw_data['coord1']/resolution)+1
    norm_data_temp['coord2'] = (raw_data['coord2']/resolution)+1

    # access to normalization coefficients
    norm_data_temp.insert(1, 'coord1_coef', [norm_vector_1.coef.iloc[int(i)] for i in norm_data_temp.coord1])
    norm_data_temp.insert(3, 'coord2_coef', [norm_vector_2.coef.iloc[int(i)] for i in norm_data_temp.coord2])

    norm_data_temp.insert(5, 'norm_value', round(norm_data_temp.value/(norm_data_temp.coord1_coef*norm_data_temp.coord2_coef), 3))
    print(norm_data_temp[norm_data_temp['value'] == float('Inf')])
    print(len(norm_data_temp))
    with pd.option_context('mode.use_inf_as_null', True):
        res = norm_data_temp.dropna()
    print(len(res))
    norm_data = raw_data[['coord1', 'coord2']]
    norm_data.insert(2, 'value', res.norm_value)
    return norm_data

def main():
        for chrom_pair in os.listdir(directory):
            if chrom_pair != '.DS_Store':
                def select(name):
                    return directory + '/' + chrom_pair + '/MAPQGE30/' + name

                chr1, chr2 = chrom_pair.split('_')
                chr_id1 = chr1.replace('chr', '')
                chr_id2 = chr2.replace('chr', '')
                if chr_id1 != "X":
                    chr_id1 = int(chr_id1)
                print('This is first chromosome #' + str(chr_id1))
                if chr_id2 != "X":
                    chr_id2 = int(chr_id2)
                print('This is second chromosome #' + str(chr_id2))

                with open(select('chr'+str(chr_id1)+'_'+str(chr_id2)+'_'+resolution+'kb.RAWobserved')) as raw:
                    data = pd.read_table(raw, header=None)
                    data.columns = ['coord1', 'coord2', 'value']


                if normalization == 'KR':
                    with open(select('chr'+str(chr_id1)+'_'+resolution+'kb.KRnorm')) as KRnorm1:
                        KRnorm_vector1 = pd.read_table(KRnorm1, header=None)
                        KRnorm_vector1.rename(columns={0:'coef'}, inplace=True)

                    with open(select('chr'+str(chr_id2)+'_'+resolution+'kb.KRnorm')) as KRnorm2:
                        KRnorm_vector2 = pd.read_table(KRnorm2, header=None)
                        KRnorm_vector2.rename(columns={0:'coef'}, inplace=True)

                    KRnorm_data = normalize(data, KRnorm_vector1, KRnorm_vector2, 100000)
                    KRnorm_data.to_csv(select('chr'+str(chr_id1)+'_'+str(chr_id2)+'_'+resolution+'kb_KRnormalized'), sep='\t', na_rep=0)


                elif normalization == 'VC':
                    with open(select('chr'+str(chr_id1)+'_'+resolution+'kb.VCnorm')) as VCnorm1:
                        VCnorm_vector1 = pd.read_table(VCnorm1, header=None)
                        VCnorm_vector1.rename(columns={0:'coef'}, inplace=True)

                    with open(select('chr'+str(chr_id2)+'_'+resolution+'kb.VCnorm')) as VCnorm2:
                        VCnorm_vector2 = pd.read_table(VCnorm2, header=None)
                        VCnorm_vector2.rename(columns={0:'coef'}, inplace=True)

                    VCnorm_data = normalize(data, VCnorm_vector1, VCnorm_vector2, 100000)
                    VCnorm_data.to_csv(select('chr'+str(chr_id1)+'_'+str(chr_id2)+'_'+resolution+'kb_VCnormalized'), sep='\t', na_rep=0)


                if normalization == 'SQRTVC':
                    with open(select('chr'+str(chr_id1)+'_'+resolution+'kb.SQRTVCnorm')) as SQRTVCnorm1:
                        SQRTVCnorm_vector1 = pd.read_table(SQRTVCnorm1, header=None)
                        SQRTVCnorm_vector1.rename(columns={0:'coef'}, inplace=True)

                    with open(select('chr'+str(chr_id2)+'_'+resolution+'kb.SQRTVCnorm')) as SQRTVCnorm2:
                        SQRTVCnorm_vector2 = pd.read_table(SQRTVCnorm2, header=None)
                        SQRTVCnorm_vector2.rename(columns={0:'coef'}, inplace=True)

                    SQRTVCnorm_data = normalize(data, SQRTVCnorm_vector1, SQRTVCnorm_vector2, 100000)
                    SQRTVCnorm_data.to_csv(select('chr'+str(chr_id1)+'_'+str(chr_id2)+'_'+resolution+'kb_SQRTVCnormalized'), sep='\t', na_rep=0)
        print(str(normalization) + ' normalization is over!')

if __name__ == '__main__':
    main()



