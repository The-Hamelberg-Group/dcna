#  ## TEST: logodd-ratio of contact probabilities as weight
#  ## logodd: log(p/(1-p))
buildnet_logodd_ratio <- function(f1, f2=NULL, f3=NULL, f4=NULL, directed=TRUE, omit.zero.change=FALSE, pc=0.2, nc=2, ec=0.0, ddf.only=FALSE) {

   kb <- 0.001985875
   Temp <- 300.0
   penalty <- 2.0  # kcal/mol

   # pc: contact cutoff; all dg < pc are set to 2 kcal/mol and >(1-pc) set to -2 kcal/mol
   # significant pair: not all zeros (over corrected f1 and f2)
   # covalent neighbors dg will be finite; but it is okay because we will skip by nc anyway
   maxn1 <- max(f1[, 3])
   f1[, 3] <- f1[, 3]/maxn1
   f1[f1[, 3]<pc, 3] <- 0
   f1[f1[, 3]>(1-pc), 3] <- 1
   signif.inds <- which((f1[, 3]>0 & f1[, 3]<1))
   mlogit1 <- -kb*Temp*log(f1[, 3]/(1-f1[, 3]))
   mlogit1[mlogit1==Inf] <- penalty
   mlogit1[mlogit1==-Inf] <- -penalty
   ddf <- f1
   ddf[, 3] <- mlogit1

   if(!is.null(f2)) {
      maxn2 <- max(f2[, 3])
      f2[, 3] <- f2[, 3]/maxn2
      f2[f2[, 3]<pc, 3] <- 0
      f2[f2[, 3]>(1-pc), 3] <- 1
      signif.inds <- which(!((f1[, 3]==0 & f2[, 3]==0) | (f1[, 3]==1 & f2[, 3]==1)))
   
      mlogit2 <- -kb*Temp*log(f2[, 3]/(1-f2[, 3]))
      mlogit2[mlogit2==Inf] <- penalty
      mlogit2[mlogit2==-Inf] <- -penalty
   
      df1 <- f1
      df1[, 3] <- mlogit2 - mlogit1
      ddf <- df1
   }

   if(!is.null(f2) && !is.null(f3) && !is.null(f4)) {
      maxn3 <- max(f3[, 3])
      maxn4 <- max(f4[, 3])
      f3[, 3] <- f3[, 3]/maxn3
      f4[, 3] <- f4[, 3]/maxn4
      f3[f3[, 3]<pc, 3] <- 0
      f3[f3[, 3]>(1-pc), 3] <- 1
      f4[f4[, 3]<pc, 3] <- 0
      f4[f4[, 3]>(1-pc), 3] <- 1
      signif.inds <- which(!((f1[, 3]==0 & f2[, 3]==0 & f3[, 3]==0 & f4[, 3]==0) | 
                             (f1[, 3]==1 & f2[, 3]==1 & f3[, 3]==1 & f4[, 3]==1)))

      mlogit3 <- -kb*Temp*log(f3[, 3]/(1-f3[, 3]))
      mlogit4 <- -kb*Temp*log(f4[, 3]/(1-f4[, 3]))
      mlogit3[mlogit3==Inf] <- penalty
      mlogit3[mlogit3==-Inf] <- -penalty
      mlogit4[mlogit4==Inf] <- penalty
      mlogit4[mlogit4==-Inf] <- -penalty

      df2 <- f3
      df2[, 3] <- mlogit4 - mlogit3

      ddf <- f1
      ddf[, 3] <- df2[, 3] - df1[, 3]
   }

   meaningful.inds <- signif.inds
   if(omit.zero.change) {
      meaningful.inds <- intersect(signif.inds, which(abs(ddf[, 3])>0))
   }
   if(ddf.only) {
      ddf[setdiff(1:nrow(ddf), meaningful.inds), 3] <- NA
      return(ddf)
   }
   else if(is.null(f2)) {
      stop("Single f provide. ddf.only must be TRUE")
   }
   if(directed) {
      maxv <- max(ddf[meaningful.inds, 3])
      minv <- min(ddf[meaningful.inds, 3])
      ddf[meaningful.inds, 3] <- ddf[meaningful.inds, 3] - minv - log(0.9999)
   } else {
      maxv <- max(abs(ddf[meaningful.inds, 3]))
      ddf[meaningful.inds, 3] <- maxv - abs(ddf[meaningful.inds, 3]) - log(0.9999)
   }
   nres <- max(ddf[, 1])
   cij <- matrix(0, nres, nres)

   cij[as.matrix(ddf[, 1:2])] <- ddf[, 3]
   # mask neighboring residues, i - i+n, n<nc
   cij[diag.ind(cij, n=nc)] <- 0
   cij[lower.tri(cij)] <- t(cij)[lower.tri(cij)]

   cna(cij, cutoff.cij=ec, network.only=TRUE, minus.log=FALSE)
}

## -log(|df|)
buildnet <- function(f1, f2, pc=0.0, nc=2) {
   f1[, 3] <- f1[, 3] / max(f1[, 3])
   f2[, 3] <- f2[, 3] / max(f2[, 3])
   df <- f1
   df[, 3] <- f2[, 3] - f1[, 3]
   nres <- max(df[, 1])
   cij <- matrix(0, nres, nres)
   cij[as.matrix(df[, 1:2])] <- df[, 3]
   # mask neighboring residues, i - i+n, n<nc
   cij[diag.ind(cij, n=nc)] <- 0
   cij[lower.tri(cij)] <- t(cij)[lower.tri(cij)]

   cna(cij, cutoff.cij=pc, network.only=TRUE)
}

## Dynamical network and protein contact network
buildnet_dn_pcn <- function(f1, cij=NULL, pc=0.75, nc=2) {
   f1[, 3] <- f1[, 3] / max(f1[, 3])
   f1[, 3] <- as.numeric( f1[, 3] >= pc)
   nres <- max(f1[, 1])
   cm <- matrix(0, nres, nres)
   cm[as.matrix(f1[, 1:2])] <- f1[, 3]
   # mask neighboring residues, i - i+n, n<nc
   cm[diag.ind(cm, n=nc)] <- 0
   cm[lower.tri(cm)] <- t(cm)[lower.tri(cm)]

   if(is.null(cij)) {
      cna(cm, cutoff.cij=0.0, network.only=TRUE)
   }
   else {
      cna(cij, cutoff.cij=0.0, cm=cm, network.only=TRUE)
   }
}

