#!/usr/bin/env python3
# -*- coding:utf-8 -*-

import sys
import pandas as pd

dataframe1 = sys.argv[1] # files in .csv format with columns
dataframe2 = sys.argv[2] # 'chrom1', 'coord1', 'chrom2', 'coord2'

def main():
    df1 = pd.read_csv(dataframe1)
    df2 = pd.read_csv(dataframe2)

    df1_for_jaccard = df1[['chrom1','coord1', 'chrom2', 'coord2']]
    df2_for_jaccard = df2[['chrom1','coord1', 'chrom2', 'coord2']]

    intersect = pd.merge(df1_for_jaccard, df2_for_jaccard, how='inner',on=['chrom1', 'coord1', 'chrom2', 'coord2'])
    join = pd.concat([df1_for_jaccard,df2_for_jaccard]).drop_duplicates().reset_index(drop=True)
    jaccard_coef = len(intersect)/len(join)
    print('Jaccard coefficient for these datasets is: ' + str(jaccard_coef))

if __name__ == '__main__':
    main()

