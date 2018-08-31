
############ PARAMETERS YOU CAN ADJUST ################################
# Filename of the first contact statistics file
contfile1 = '../example/free'

# Filename of the second contact statistics file
contfile2 = '../example/cis'

# Filename of the pdb file
pdbfile = '../example/cypa_wt_nowater.pdb'

# threshold to define "stable" contact (e.g. 90% or 0.9)
qcutoff = 0.9

# threshold to define "significant" contact difference between two states
# (0.1 means if contact frequency difference between the two states is less than 10%, ignore them)
fcutoff = 0.1

# number of communities to show; If set `NULL`, use the automatically determined number.
num.comm = NULL


############# DON'T CHANGE FOLLOWING LINES (UNLESS YOU KNOW WHAT IT MEANS) ##########
library(bio3d)
library(igraph)
source('funs.R')

freqs <- lapply(list(f1=contfile1, f2=contfile2), read.table)
pdb <- read.pdb(pdbfile)
nets <- buildnet(freqs, pdb, qcut=qcutoff, fcut=fcutoff, ncom=num.comm)
plot_community(nets[[1,2]], pdb=pdb, fcut=fcutoff)

## Save manually adjusted layout and replot
# layout <- tk_coords(1)  ## The number "1" is the number shown on the top of the plot window
# plot_community(nets[[1,2]], pdb=pdb, layout=layout, fcut=fcutoff, interactive=FALSE)

## Adjust vertex and edge sizes
# plot_community(nets[[1,2]], pdb=pdb, layout=layout, fcut=fcutoff, interactive=FALSE, mag.edge=5, mag.vertex=0.8)

## check modularity scores for alternative "num.comm" values
# rem <- community.tree(nets[[1,2]])
# plot(rem$num.of.comm, rem$modularity, typ='b', cex=0.6, xlab='Number of Communities', ylab='Modularity')


