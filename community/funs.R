buildnet <- function(freqs, pdb, qcut=0.9, fcut=0.1, mcut=0.005, scut=3, 
  prune.net=FALSE, ncom=NULL, ncore=NULL) {
  if(!is.list(freqs) && length(freqs)<2) 
     stop('Input `freqs` must be a list and has at least two components.')

  require(bio3d)
  ncore <- setup.ncore(ncore)

  freqs <- lapply(freqs, function(f) {
     f[, 3] <- f[, 3] / max(f[, 3])
     f
  })
  n <- length(freqs)
  np <- pairwise(n)

#  ft <- sapply(freqs, '[', 3)
#  if(is.matrix(ft) && ncol(ft)>1) favg <- rowMeans(ft)
#  else favg <- ft

  ## find consensus stable contacts.
  fdiff <- matrix(vector('list', n*n), n, n)
  fdiff[lower.tri(fdiff)] <- parallel::mclapply(1:nrow(np), function(i) {
     ii <- np[i, 1]; jj <- np[i, 2]
     freqs[[jj]][, 3] - freqs[[ii]][, 3]
  }, mc.cores=ncore)
  fdiff <- t(fdiff)
  fdiff[lower.tri(fdiff)] <- parallel::mclapply(1:nrow(np), function(i) {
     ii <- np[i, 1]; jj <- np[i, 2]
     freqs[[ii]][, 3] - freqs[[jj]][, 3]
  }, mc.cores=ncore)

  flags <- rep(TRUE, nrow(freqs[[1]]))
  for(i in 1:n) {
    flags <- flags & freqs[[i]][, 3]>=qcut
  }
  for(i in 1:nrow(np)) {
     flags <- flags & abs(fdiff[[np[i, 1], np[i, 2]]])<fcut
  }
#  flags <- flags | (freqs[[1]][, 2] - freqs[[1]][, 1]) <= 2  # force 1-3 contact
#  flags <- flags & !((freqs[[1]][, 2] - freqs[[1]][, 1]) == 2)  # force 1-3 contact out
#  flags <- flags & (freqs[[1]][, 2] - freqs[[1]][, 1]) >= scut # non-local contact only

  ## network shell
  ca.inds <- atom.select(pdb, elety='CA', verbose=FALSE)
#  ca.inds <- combine.select(ca.inds, atom.select(pdb, elety="C4'", verbose=FALSE), operator='+', verbose=FALSE)
#  ca.inds <- combine.select(ca.inds, atom.select(pdb, resid=c("FEM", "WAT"), verbose=FALSE), operator='+', verbose=FALSE)
  cij <- matrix(0, length(ca.inds$atom), length(ca.inds$atom))
  cij[as.matrix(apply(freqs[[1]][, 1:2], 2, match, pdb$atom$resno[ca.inds$atom]))] <- 
      as.numeric(flags)
  cij[lower.tri(cij)] <- t(cij)[lower.tri(cij)]

  net <- cna(cij, cutoff.cij=0.0)
  if(!is.null(ncom)) {
     net <- prunenet(net, ncom=ncom)
  } else if(prune.net==TRUE) {
     net <- prunenet(net, cutoff=mcut)
  }

  ## update community networks
  nets <- matrix(vector('list', n*n), n, n)  
  nets[upper.tri(nets)] <- parallel::mclapply(fdiff[upper.tri(fdiff)], function(fdiff) {
     dcij <- matrix(0, length(ca.inds$atom), length(ca.inds$atom))
     dcij[as.matrix(apply(freqs[[1]][, 1:2], 2, match, pdb$atom$resno[ca.inds$atom]))] <- fdiff
     dcij[diag.ind(dcij, n=scut)] <- 0
     dcij[lower.tri(dcij)] <- t(dcij)[lower.tri(dcij)]
     remodel.cna(net, dcij, fcut=fcut)
#     dnet.free_cis <- cna(dcij, cutoff.cij=0.0)
  }, mc.cores=ncore)
  nets[lower.tri(nets)] <- parallel::mclapply(fdiff[lower.tri(fdiff)], function(fdiff) {
     dcij <- matrix(0, length(ca.inds$atom), length(ca.inds$atom))
     dcij[as.matrix(apply(freqs[[1]][, 1:2], 2, match, pdb$atom$resno[ca.inds$atom]))] <- fdiff
     dcij[diag.ind(dcij, n=scut)] <- 0
     dcij[lower.tri(dcij)] <- t(dcij)[lower.tri(dcij)]
     remodel.cna(net, dcij, fcut=fcut)
#     dnet.free_cis <- cna(dcij, cutoff.cij=0.0)
  }, mc.cores=ncore)
  
  return(nets)
#  return(list(free_cis=net.free_cis, dfree_cis=dnet.free_cis, cis_ts=net.cis_ts, dcis_ts=dnet.cis_ts))
}

# To update network with new cij values.
remodel.cna <- function(net, new.cij, method=c('sum', 'max'), edge.col=NULL, fcut=0.1) {
  require(igraph)
  method <- match.arg(method)
  
  cg.cij <- matrix(0, nrow(net$community.cij), ncol(net$community.cij))
  for(i in 1:nrow(cg.cij)) {
  	for(j in 1:nrow(cg.cij)) {
  		cij.sub <- 
  		   new.cij[net$communities$membership==i, 
  		           net$communities$membership==j]
  		cg.cij[i, j] <- switch(method,  
  		  "sum" = sum(cij.sub),
  		  "max" = max(cij.sub) )
  	}
  }
  net$community.cij <- cg.cij
  net$community.network <- graph.adjacency(net$community.cij,
              mode = "undirected", weighted = TRUE, diag = FALSE)
  if(is.null(edge.col)) {
    	E(net$community.network)$color <- rep('blue', length(E(net$community.network)$weight))
    	E(net$community.network)$color[E(net$community.network)$weight<0] <- 'red'
    	E(net$community.network)$color[abs(E(net$community.network)$weight)<fcut] <- 'gray70'
  } 
  else {
    	E(net$community.network)$color <- rep(edge.col, length(E(net$community.network)$weight))
  }
  V(net$community.network)$color <- vmd_colors()[1:max(net$communities$membership)]
  V(net$community.network)$size <- table(net$communities$membership)
  net
}

prunenet <- function(net, cutoff=0.005, ncom=NULL) {
  ## prune network communities by checking modularities.

  remodel <- community.tree(net, rescale=TRUE)
  #plot(remodel$num.of.comms, remodel$modularity, typ="p")

  # choose network to plot
  if(!is.null(ncom)) {
     nncomm=ncom
  } else {
     n.max = length(unique(net$communities$membership))
     ind.max = which(remodel$num.of.comms == n.max)
     v = remodel$modularity[length(remodel$modularity):ind.max]
     v = rev(diff(v))
     fa = which(v>=cutoff)[1] - 1
     nncomm = ifelse(is.na(fa), min(remodel$num.of.comms), n.max - fa)
  }
  ind <- which(remodel$num.of.comms == nncomm)
  network.amendment(net, remodel$tree[ind, ])
}

plot_community <- function(net, pdb, layout=NULL, fcut=0.1, interactive=TRUE, 
   mag.edge=10, mag.vertex=1) {
  vmd(net, pdb, vmdfile='network.vmd', pdbfile='network.pdb')
  wt <- E(net$community.network)$weight
  elabel <- round(wt, 1)
  elabel[abs(wt) < fcut] <- ''
  if(is.null(layout)) {
     plot(net, pdb=pdb, weights=abs(wt)*mag.edge, vertex.size=V(net$community.network)$size*mag.vertex, 
        interactive=interactive, edge.label=elabel)
  } else {
    plot(net, layout=layout, weights=abs(wt)*mag.edge, vertex.size=V(net$community.network)$size*mag.vertex, 
        interactive=interactive, edge.label=elabel)
  } 
}
