# Title: Modified Rao's quadratic entropy
# The original source function was extracted from the package "SYNCSA: Analysis of Functional and Phylogenetic Patterns in Metacommunities"
# In the original rao.diversity() function, the Gower's dissimilarity is used to calculate the dissimilarity matrix for the traits, hence with Euclidean properties.
# As the trait matrice extracted from metagenomic data has specific properties (sparse, i.e. zero-inflated, discrete data), I chose the Bray-Curtis dissimilarity metric which is recognized to be appropriate for such datasets.
# I hence modified line 88 of the present code and replaced method="gower" by method="bray 
# original: dist.functional <- sqrt(as.matrix(FD::gowdis(x=traits, asym.bin = NULL, ord = ord, w = weights, ...)))
# modified: dist.functional <- as.matrix(vegan::vegdist(x=traits, method="bray", binary=FALSE, diag=FALSE, upper=FALSE,na.rm = FALSE, ...))


rao.diversity.mod

function(comm, traits = NULL, phylodist = NULL, checkdata = TRUE, ord = "metric",
         put.together = NULL, standardize = TRUE, ...)
{
  diver.internal <- function(community, distance){
    if(any(is.na(distance))){
      distance.na <- ifelse(is.na(distance), 0, 1)
      inter.na <- community%*%distance.na
      adjustment <- rowSums(sweep(community, 1, inter.na, "*", check.margin = FALSE))
      distance[is.na(distance)] <- 0
      inter <- community%*%distance
      res <- rowSums(sweep(community, 1, inter, "*", check.margin = FALSE))
      res <- ifelse(adjustment>0, res/adjustment, res)
    } else{
      inter <- community%*%distance
      res <- rowSums(sweep(community, 1, inter, "*", check.margin = FALSE))
    }
    return(res)
  }
  res <- list(call = match.call())
  if (inherits(comm, "metacommunity.data")) {
    if (!is.null(traits) | !is.null(phylodist) | !is.null(put.together)) {
      stop("\n When you use an object of class metacommunity.data the arguments traits, phylodist and put.together must be null. \n")
    }
    traits <- comm$traits
    phylodist <- comm$phylodist
    put.together <- comm$put.together
    comm <- comm$community
  }
  list.warning <- list()
  if(checkdata){
    organize.temp <- organize.syncsa(comm, traits = traits, phylodist = phylodist, check.comm = TRUE)
    if(!is.null(organize.temp$stop)){
      organize.temp$call <- match.call()
      return(organize.temp)
    }
    list.warning <- organize.temp$list.warning
    comm <- organize.temp$community
    traits <- organize.temp$traits
    phylodist <- organize.temp$phylodist
  }
  if(length(list.warning)>0){
    res$list.warning <- list.warning
  }
  if(any(is.na(comm))){
    stop("\n community data with NA\n")
  }
  comm <- as.matrix(comm)
  N <- nrow(comm)
  S <- ncol(comm)
  dist.1 <- 1 - diag(x = rep(1, S))
  if (!is.null(traits)) {
    traits <- as.data.frame(traits)
    m <- ncol(traits)
    weights <- rep(1, m)
    make.names <- is.null(colnames(traits))
    colnames(traits) <- colnames(traits, do.NULL = FALSE, prefix = "T")
    names(weights) <- colnames(traits)
    if(!is.null(put.together)){
      if(!inherits(put.together, "list")){
        stop("\n put.together must be a object of class list\n")
      }
      if(make.names){
        for(k in 1:length(put.together)){
          put.together[[k]] <- paste("T", put.together[[k]], sep = "")
        }
      }
      if(max(table(unlist(put.together)))>1){
        stop("\n The same trait appears more than once in put.together\n")
      }
      if(length(setdiff(unlist(put.together), colnames(traits)))>0){
        stop("\n Check traits names in put.together\n")
      }
      for(k in 1:length(put.together)){
        weights[put.together[[k]]] <- 1/length(put.together[[k]])
      }
    }
    dist.functional <- as.matrix(vegan::vegdist(x=traits, method="bray", binary=FALSE, diag=FALSE, upper=FALSE,na.rm = FALSE, ...))
    if (checkdata) {
      if(any(is.na(dist.functional))){
        # stop("\n traits with too much NA \n")
        warning("Warning: NA in distance between species", call. = FALSE)
      }
    }
  }
  if (!is.null(phylodist)) {
    dist.phylogenetic <- as.matrix(phylodist)
    if (checkdata) {
      if(any(is.na(dist.phylogenetic))){
        # stop("\n phylodist with NA \n")
        warning("Warning: NA in phylodist", call. = FALSE)
      }
    }
    if(standardize){
      dist.phylogenetic <- dist.phylogenetic/max(dist.phylogenetic, na.rm = TRUE)
    }
  }
  comm <- sweep(comm, 1, rowSums(comm, na.rm = TRUE), "/")
  SD <- diver.internal(comm, dist.1)
  res$Simpson <- SD
  if (!is.null(traits)){
    FD <- diver.internal(comm, dist.functional)
    res$FunRao <- FD
    res$FunRedundancy <- SD-FD
  }
  if (!is.null(phylodist)){
    PD <- diver.internal(comm, dist.phylogenetic)
    res$PhyRao <- PD
    res$PhyRedundancy <- SD-PD
  }
  return(res)
}
<bytecode: 0x0000014942a466c8>
