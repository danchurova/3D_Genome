#!/usr/bin/env python3
# -*- coding:utf-8 -*-

import os
import sys
import pandas as pd
import numpy as np

directory = 'data/K562_interchromosomal/100kb_resolution_interchromosomal'
normalization = sys.arg[1]


def main():
        for chrom_pair in os.listdir(directory):
            if chrom_pair != '.DS_Store':
                for file in os.listdir(directory + '/' + chrom_pair + '/MAPQGE30/'):
                    if file.endswith(".RAWobserved"):

