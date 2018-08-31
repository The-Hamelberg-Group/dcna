
############ PARAMETERS YOU CAN ADJUST ################################
# Filename of the first contact statistics file
contfile1 = '../example/free'

# Filename of the second contact statistics file
contfile2 = '../example/cis'

# Filename of the pdb file
pdbfile = '../example/cypa_wt_nowater.pdb'

# threshold to define "significant" contact difference between two states
# (0.1 means if contact frequency difference between the two states is less than 10%, ignore them)
fcutoff = 0.1

############# DON'T CHANGE FOLLOWING LINES (UNLESS YOU KNOW WHAT IT MEANS) #########
library(bio3d)
source('vmd.cylinder.R')

freqs <- lapply(list(f1=contfile1, f2=contfile2), read.table)
pdb <- read.pdb(pdbfile)

f1 <- freqs$f1
f2 <- freqs$f2
f1[, 3] <- f1[, 3]/max(f1[, 3])
f2[, 3] <- f2[, 3]/max(f2[, 3])
df <- f1
df[, 3] <- f2[, 3] - f1[, 3]
vmd.cylinder(df, pdb=pdb, file='bar', cutoff=fcutoff)

