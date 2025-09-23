#!/bin/bash
# Database searching

source activate openms_36

DIR='../data/mzML-converted/'

python calc_tic.py "amundsen_prot_normalization_factors.csv" $(ls ../data/mzML-converted/210420_0977_097_SM*.mzML) &




