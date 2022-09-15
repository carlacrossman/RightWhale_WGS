#bin/bash

module load poplddecay/3.41

PopLDdecay    -InVCF  ~/projects/def-frasiert/RW_WGS/vcf/locked/srw_unrelated_filtered_aug5.vcf.gz  -OutStat LDdecay
PopLDdecay    -InVCF  ~/projects/def-frasiert/RW_WGS/vcf/locked/narw_unrelated_filtered_aug5.vcf.gz  -OutStat NARW_LDdecay
PopLDdecay    -InVCF  ~/projects/def-frasiert/RW_WGS/vcf/locked/all_on_narw_unrelated_filtered_aug5.vcf.gz  -OutStat ALL_LDdecay
