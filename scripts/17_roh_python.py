# Import Libraries
import numpy as np
import scipy
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import allel; print('scikit-allel', allel.__version__)
import hmmlearn

narw_vcf = allel.read_vcf('narw_unrelated_filtered_aug5.vcf.gz', fields=['variants/*','calldata/*'])

sample_no = []
scaff = []
narw_percent_roh = []
narw_num_roh = []
for SCAFFOLD in np.unique(narw_vcf['variants/CHROM']):
    narw_gt = allel.GenotypeArray(narw_vcf['calldata/GT'][narw_vcf['variants/CHROM'] == SCAFFOLD])
    for IND in range(0,9):
        gv = narw_gt[:,IND]
        pos = narw_vcf['variants/POS'][narw_vcf['variants/CHROM'] == SCAFFOLD]
        scaffold_size = pos[len(pos)-1]
        roh = allel.roh_mhmm(gv, pos,min_roh=1000000, contig_size=scaffold_size)
        sample_no.append(IND)
        scaff.append(SCAFFOLD)
        narw_percent_roh.append(roh[1])
        narw_num_roh.append(len(roh[0]))

# Save ROH data per sample per scaffold to a .csv file
narw_roh_wide = pd.DataFrame((sample_no, scaff, narw_num_roh, narw_percent_roh))
narw_roh_persamp = narw_roh_wide.transpose()
narw_roh_persamp.columns = ['Sample_No', 'Scaffold', 'Number_ROHs', 'Percent_Scaffold_in_ROH']
pd.DataFrame(narw_roh_persamp).to_csv('narw_roh_persamp_1Mb.csv')

