# Import Libraries
import numpy as np
import scipy
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import allel; print('scikit-allel', allel.__version__)
import hmmlearn

srw_vcf = allel.read_vcf('srw_unrelated_filtered_aug5.vcf.gz', fields=['variants/*','calldata/*'])


srw_roh = pd.DataFrame()
for SCAFFOLD in np.unique(srw_vcf['variants/CHROM']):
    srw_gt = allel.GenotypeArray(srw_vcf['calldata/GT'][srw_vcf['variants/CHROM'] == SCAFFOLD])
    for IND in range(0,9):
        gv = srw_gt[:,IND]
        pos = srw_vcf['variants/POS'][srw_vcf['variants/CHROM'] == SCAFFOLD]
        scaffold_size = pos[len(pos)-1]
        roh = allel.roh_mhmm(gv, pos,min_roh=10000, contig_size=scaffold_size, phet_roh=0)
        roh[0]['SAMPLE']=IND
        roh[0]['CHROMOSOME']=SCAFFOLD
        srw_roh=srw_roh.append(roh[0])

# Save ROH data per sample per scaffold to a .csv file
srw_roh.to_csv('srw_roh_details_nov28.csv')

