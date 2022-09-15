# Import Libraries
import numpy as np
import scipy
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import allel; print('scikit-allel', allel.__version__)
import hmmlearn

narw_vcf = allel.read_vcf('narw_unrelated_filtered_aug5.vcf.gz', fields=['variants/*','calldata/*'])


narw_roh = pd.DataFrame()
for SCAFFOLD in np.unique(narw_vcf['variants/CHROM']):
    narw_gt = allel.GenotypeArray(narw_vcf['calldata/GT'][narw_vcf['variants/CHROM'] == SCAFFOLD])
    for IND in range(0,9):
        gv = narw_gt[:,IND]
        pos = narw_vcf['variants/POS'][narw_vcf['variants/CHROM'] == SCAFFOLD]
        scaffold_size = pos[len(pos)-1]
        roh = allel.roh_mhmm(gv, pos,min_roh=10000, contig_size=scaffold_size)
        roh[0]['SAMPLE']=IND
        roh[0]['CHROMOSOME']=SCAFFOLD
        narw_roh=narw_roh.append(roh[0])

# Save ROH data per sample per scaffold to a .csv file
narw_roh.to_csv('narw_roh_details.csv')

