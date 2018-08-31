vmd.cylinder <- function(edges, pdb, file="cylinder", cutoff=0.1, scale=1.0) {
	if(is.pdb(pdb)) {
      ca.inds <- atom.select(pdb, 'calpha', verbose=FALSE)
      xyz <- pdb$xyz[ca.inds$xyz]
    } else {
      stop("PDB must be provided")
    }

    pdbfile <- paste(file, ".pdb", sep="")
    vmdfile <- paste(file, ".vmd", sep="")
    
    write.pdb(pdb, file=pdbfile)

    ## Generate the script for drawing cylinders between C-alpha atoms.
    cat("mol new ", basename(pdbfile), " type pdb\n", 
    	"mol delrep 0 top\n",
    	"mol representation NewCartoon\n", 
    	"mol color colorID 8\n", 
    	"mol addrep top\n", file=vmdfile)

    cat("display update off\n", file=vmdfile, append=TRUE)

    if(any(edges[, 3]<0)) { 
        cat("draw color blue\n", file=vmdfile, append=TRUE)
        inds <- which(edges[, 3] > cutoff)
        if(length(inds)>0) {
          for(i in 1:length(inds)) {
          	a <- edges[inds[i], 1]
          	b <- edges[inds[i], 2]
          	val <- edges[inds[i], 3]
          	cat("draw cylinder {", paste(xyz[atom2xyz(a)], collapse=" "), 
          	  "} {", paste(xyz[atom2xyz(b)], collapse=" "), "} radius ",
          	  abs(val)*scale, " resolution 20 filled 0\n", file=vmdfile, append=TRUE)
          }
        }
  
        cat("draw color red\n", file=vmdfile, append=TRUE)
        inds <- which(edges[, 3] < -cutoff)
        if(length(inds)>0) {
          for(i in 1:length(inds)) {
          	a <- edges[inds[i], 1]
          	b <- edges[inds[i], 2]
          	val <- edges[inds[i], 3]
          	cat("draw cylinder {", paste(xyz[atom2xyz(a)], collapse=" "), 
          	  "} {", paste(xyz[atom2xyz(b)], collapse=" "), "} radius ",
          	  abs(val)*scale, " resolution 20 filled 0\n", file=vmdfile, append=TRUE)
          }
        }
    } else {
        cat("draw color blue\n", file=vmdfile, append=TRUE)
        inds <- which(edges[, 3] >= cutoff)
        if(length(inds)>0) {
          for(i in 1:length(inds)) {
              a <- edges[inds[i], 1]
              b <- edges[inds[i], 2]
              val <- edges[inds[i], 3]
              cat("draw cylinder {", paste(xyz[atom2xyz(a)], collapse=" "), 
                "} {", paste(xyz[atom2xyz(b)], collapse=" "), "} radius ",
                abs(val)*scale, " resolution 20 filled 0\n", file=vmdfile, append=TRUE)
          }
        }

        cat("draw color red\n", file=vmdfile, append=TRUE)
        inds <- which(edges[, 3] < cutoff)
        if(length(inds)>0) {
          for(i in 1:length(inds)) {
              a <- edges[inds[i], 1]
              b <- edges[inds[i], 2]
              val <- edges[inds[i], 3]
              cat("draw cylinder {", paste(xyz[atom2xyz(a)], collapse=" "), 
                "} {", paste(xyz[atom2xyz(b)], collapse=" "), "} radius ",
                abs(val)*scale, " resolution 20 filled 0\n", file=vmdfile, append=TRUE)
          }
        }
        # cols <- colorRamp(c("white", "blue"))
        # inds <- which(edges[, 3] > cutoff)
        # cat("set color_start [colorinfo num]\n", file=vmdfile, append=TRUE)
        # for(i in 1:length(inds)) {
        #     a <- edges[inds[i], 1]
        #     b <- edges[inds[i], 2]
        #     val <- edges[inds[i], 3]
        #     col <- cols(val)[1:3] / 255
        #     cat("color change rgb [expr ", i, " + $color_start] ", 
        #         paste(col, collapse=" "), "\n", sep="", file=vmdfile, append=TRUE)
        #     cat("graphics top color [expr ", i, " + $color_start]\n", sep="", 
        #         file=vmdfile, append=TRUE)
        #     cat("draw cylinder {", paste(xyz[atom2xyz(a)], collapse=" "), 
        #       "} {", paste(xyz[atom2xyz(b)], collapse=" "), "} radius ",
        #       abs(val)*scale, " resolution 6 filled 0\n", file=vmdfile, append=TRUE)
        # }
    }
  
    # turn on display update
    cat("display update on\n", file=vmdfile, append=TRUE)

}
