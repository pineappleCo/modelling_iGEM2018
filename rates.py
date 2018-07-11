from math import log

copy_num = [1.]

#relative strength of anderson promoters
# http://parts.igem.org/Promoters/Catalog/Anderson
anderson_str = [1., 0.7, 0.86, 0.01, 0.72, 0.24, 0.47, 0.36, 0.51, 0.04, 0.33, 0.58, 0.01, 0.1, 0.15, 0.16, 0.06, 0.56]

#modelled as in salis lab's rbs calculator - https://salislab.net/software/ - a proportional scale is used from 0 to 100,000+ in which 1 is the translation initiation rate of unaltered shine-dalgarno sequence
rbs_affinity = [1.]

#speed of RNA polymerase (60nt/s) / length (nt), imm not expressed in maxicells
imm_transcription_rate = 60./258.
col_transcription_rate = 60./1743.
imm_transcription_rate_maxicell = 0.

# average mRNA half-life in e. coli - 5min - http://book.bionumbers.org/
# exponential decay constant = ln(2)/half-life
mRNA_deg_rate = log(2)/300.

# movement speed of ribosome (20 nt/s - http://book.bionumbers.org/) / num
# amino acids
imm_translation_rate = 20./83.
col_translation_rate = 20./581.

# effective halflife of proteins with halflife greater than generation length
# in micro-organisms is the generation length - 40min in e coli
# exponential decay constant = ln(2)/half-life
imm_deg_rate = log(2)/((40./60.)*3600.)
col_deg_rate = log(2)/((40./60.)*3600.)

# expasy protparam - https://web.expasy.org/protparam/ - estimates half-life
# >10hrs in E. coli and 30hrs in reticulocyte, average at 20hrs
# exponential decay constant = ln(2)/half-life
imm_deg_rate_maxicell = log(2)/(20.*3600.)
col_deg_rate_maxicell = log(2)/(20.*3600.)
