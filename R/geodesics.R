##########################################################################
#
# Geodesic routines for ElectroGraph (sigh).
# Andrew C. Thomas
# January 19, 2011
#
# Edge-based methods are implemented.
#
##########################################################################


#This routine assumes a 1:n notation for the node IDs.
geodesic.selected.from.edgelist <- function(edgelist, sourcepoint=1, destpoint=NULL, node.ids=NULL) {
  if (dim(edgelist)[2] != 2 & dim(edgelist)[2] != 3) {
    stop ("edgelist needs to be k-by-2 or k-by-3 in dimension.")
  } else if (dim(edgelist)[2] == 2) edgelist <- cbind(edgelist, 1)

  #if the edgelist is unprocessed?
  edgelist <- as.matrix(edgelist); if (!is.numeric(edgelist)) stop ("geodesic.selected.from.edgelist error: edgelist input must be numeric, and should correspond to the vector (1:n).")

  nn.temp <- max(c(edgelist))
    
  if (is.null(node.ids)) {
    node.ids <- 1:nn.temp
    nn <- nn.temp
  } else {
    nn <- length(node.ids)
    if (nn.temp > nn) stop("geodesic-selected error: Length of node.ids is shorter than the apparent number of nodes in edgelist.")
  }

  if (is.null(destpoint)) destpoint <- nn+1;  #no destpoint.
  output <- cbind(.C("dijkstra_single",
              edges=as.integer(edgelist[,1:2]-1),
              edgevalues=as.double(edgelist[,3]),
              output=as.double(rep(0, nn)),
              dim1=as.integer(nn),
              nnedges=as.integer(dim(edgelist)[1]),
              sourcey=as.integer(sourcepoint-1),
              destpoint=as.integer(destpoint-1))$output)
  
  shortest <- min(output[output>0]); output[output>1e42*shortest] <- Inf
  rownames(output) <- node.ids
  return (output)
}

#This routine assumes a 1:n notation for the node IDs.
geodesic.matrix.from.edgelist <- function(edgelist, node.ids=NULL) {
  
  if (dim(edgelist)[2] != 2 & dim(edgelist)[2] != 3) {
    stop ("edgelist needs to be k-by-2 or k-by-3 in dimension.")
  } else if (dim(edgelist)[2] == 2) edgelist <- cbind(edgelist, 1)

  #if the edgelist is unprocessed?
  edgelist <- as.matrix(edgelist); if (!is.numeric(edgelist)) stop ("geodesic matrix error: edgelist input must be numeric, and should correspond to the vector (1:n).")

  nn.temp <- max(c(edgelist))
    
  if (is.null(node.ids)) {
    node.ids <- 1:nn.temp
    nn <- nn.temp
  } else {
    nn <- length(node.ids)
    if (nn.temp > nn) stop("geodesic matrix error: Length of node.ids is shorter than the apparent number of nodes in edgelist.")
  }
  
  output <- array(.C("dijkstra_all",
                     edges=as.integer(edgelist[,1:2]-1),
                     strengths=as.double(edgelist[,3]),
                     output=as.double(rep(0, nn^2)),
                     nn=as.integer(nn),
                     pedges=as.integer(dim(edgelist)[1]))$output, rep(nn,2))
  
  shortest <- min(output[output>0]); output[output>1e42*shortest] <- Inf
  colnames(output) <- rownames(output) <- node.ids   
  return(output)  
}


single.source.betweenness <- function(edgelist, sourcepoint=1, destpoint=NULL, node.ids=NULL) {
  if (dim(edgelist)[2] != 2 & dim(edgelist)[2] != 3) {
    stop ("edgelist needs to be k-by-2 or k-by-3 in dimension.")
  } else if (dim(edgelist)[2] == 2) edgelist <- cbind(edgelist, 1)
  nn <- max(edgelist[,1:2]); edgecount <- dim(edgelist)[1]
  if (is.null(node.ids)) node.ids <- 1:nn else if (length(node.ids)<nn) stop("Length of node.ids does not match number of nodes.")

  if (is.null(destpoint)) destpoint <- nn+1;
  pathcount.proto <- .C("betweenness_geodesic_single", edges=as.integer(edgelist[,1:2]-1),
                        strength=as.double(edgelist[,3]),
                        output=as.double(rep(0, nn*edgecount)), path.lengths=as.double(rep(0, nn)),
                        dim1=as.integer(nn), edgedim=as.integer(edgecount),
                        sourcey=as.integer(sourcepoint-1), desty=as.integer(destpoint-1))
  pathcount <- array(pathcount.proto$output, c(edgecount, nn))

  rownames(pathcount) <- paste(node.ids[edgelist[,1]],"->",node.ids[edgelist[,2]],", ",sep="")
  colnames(pathcount) <- paste(node.ids[sourcepoint],":",node.ids,sep="")

  return (pathcount)

}





betweenness.centralities <- function(edgelist, path.weight=c("constant","closeness"), verbose=TRUE, node.ids=NULL) {
  #sociomatrix = eg$grand.socio; path.weight="constant"
  if (dim(edgelist)[2] != 2 & dim(edgelist)[2] != 3) {
    stop ("edgelist needs to be k-by-2 or k-by-3 in dimension.")
  } else if (dim(edgelist)[2] == 2) edgelist <- cbind(edgelist, 1)
  if (all(path.weight == c("constant","closeness"))) path.weight <- "constant"
  nn <- max(edgelist[,1:2])
  if (is.null(node.ids)) node.ids <- 1:nn else if (length(node.ids)<nn) stop("Length of node.ids does not match number of nodes.")
  
  fullout <- .C("betweenness_geodesic_full", edgelist=as.integer(edgelist[,1:2]-1), strengths=as.double(edgelist[,3]),
                output=as.double(0*edgelist[,3]), nodeoutput=as.double(rep(0,nn)),
                nn=as.integer(nn), edges=as.integer(dim(edgelist)[1]),
                distweight=as.integer(1*(path.weight=="closeness")))

  edgewise <- array(fullout$output, c(dim(edgelist)[1],1))
  rownames(edgewise) <- paste(node.ids[edgelist[,1]],"->",node.ids[edgelist[,2]],", ",sep="")
  
  nodal <-array(fullout$nodeoutput, c(nn, 1)) 
  rownames(nodal) <- node.ids

  return(list(edgewise=edgewise, nodal=nodal)) #paths=paths, dists=dists, 
  
}


recourse.betweenness.one <- function(edgelist, sourcepoint=1, destpoint=NULL, penalty=20, node.ids=NULL) {

  if (dim(edgelist)[2] != 2 & dim(edgelist)[2] != 3) {
    stop ("edgelist needs to be k-by-2 or k-by-3 in dimension.")
  } else if (dim(edgelist)[2] == 2) edgelist <- cbind(edgelist, 1)

  nn <- max(edgelist[,1:2]); edgecount <- dim(edgelist)[1]
  if (is.null(destpoint)) destpoint <- nn
  if (is.null(node.ids)) node.ids <- 1:nn else if (length(node.ids)<nn) stop("Length of node.ids does not match number of nodes.")
  
  pathcount.proto <- .C("recourse_betweenness_single_sd",
                        edges=as.integer(edgelist[,1:2]-1), strength=as.double(abs(edgelist[,3])),
                        output=as.double(rep(0, edgecount)), path.lengths=as.double(rep(0, nn)),
                        
                        dim1=as.integer(nn), edgedim=as.integer(edgecount), sourcey=as.integer(sourcepoint-1),
                        desty=as.integer(destpoint-1), ppenalty=as.double(penalty))

  pathcount <- array(pathcount.proto$output, c(edgecount, 1))

  rownames(pathcount) <- paste(node.ids[edgelist[,1]],"->",node.ids[edgelist[,2]],", ",sep="")
  colnames(pathcount) <- paste(node.ids[sourcepoint],":",node.ids[destpoint],sep="")

  return (pathcount)
}


recourse.betweenness.source <- function(edgelist, sourcepoint=1, penalty=20, node.ids=NULL) {

  if (dim(edgelist)[2] != 2 & dim(edgelist)[2] != 3) {
    stop ("edgelist needs to be k-by-2 or k-by-3 in dimension.")
  } else if (dim(edgelist)[2] == 2) edgelist <- cbind(edgelist, 1)

  nn <- max(edgelist[,1:2]); edgecount <- dim(edgelist)[1]
  if (is.null(node.ids)) node.ids <- 1:nn else if (length(node.ids)<nn) stop("Length of node.ids does not match number of nodes.")
  #if (is.null(destpoint)) destpoint <- nn
  
  pathcount.proto <- .C("recourse_betweenness_one_source", edges=as.integer(edgelist[,1:2]-1),
                        strength=as.double(abs(edgelist[,3])), output=as.double(rep(0, edgecount*nn)),
                        path.lengths=as.double(rep(0, nn)),
                        dim1=as.integer(nn), edgedim=as.integer(edgecount),
                        sourcey=as.integer(sourcepoint-1), ppenalty=as.double(penalty))

  pathcount <- array(pathcount.proto$output, c(edgecount, nn))

  rownames(pathcount) <- paste(node.ids[edgelist[,1]],"->",node.ids[edgelist[,2]],", ",sep="")
  colnames(pathcount) <- paste(node.ids[sourcepoint],":",node.ids,sep="")

  return (pathcount)
}

recourse.betweenness.full <- function(edgelist, penalty=20, path.weight=c("constant","closeness"), node.ids=NULL) {

  if (dim(edgelist)[2] != 2 & dim(edgelist)[2] != 3) {
    stop ("edgelist needs to be k-by-2 or k-by-3 in dimension.")
  } else if (dim(edgelist)[2] == 2) edgelist <- cbind(edgelist, 1)

  if (all(path.weight == c("constant","closeness"))) path.weight <- "constant"
  nn <- max(edgelist[,1:2]); edgecount <- dim(edgelist)[1]
  if (is.null(node.ids)) node.ids <- 1:nn else if (length(node.ids)<nn) stop("Length of node.ids does not match number of nodes.")
  #if (is.null(destpoint)) destpoint <- nn
  
  pathcount.proto <- .C("recourse_betweenness_full",
                        edges=as.integer(edgelist[,1:2]-1),
                        strength=as.double(abs(edgelist[,3])),
                        output=as.double(rep(0, edgecount)),
                       
                        dim1=as.integer(nn),
                        edgedim=as.integer(edgecount),
                        
                        distweight=as.integer(1*(path.weight=="closeness")),
                        ppenalty=as.double(penalty))

  pathcount <- array(pathcount.proto$output, c(edgecount, 1))

  rownames(pathcount) <- paste(node.ids[edgelist[,1]],"->",node.ids[edgelist[,2]],", ",sep="")
  colnames(pathcount) <- path.weight

  return (pathcount)
}


##############################################################################
#
# Legacy code.
#
##############################################################################

#now in C! It's actually Floyd-Warshall... with absolute values.
geodesic.mat.socio <- function(sociomatrix) {
  n.pts <- dim(sociomatrix)[1]
  distance <- 1/abs(sociomatrix)

  if (is.infinite(min(distance))) {
    out <- distance
    diag(out) <- 0
    
  } else {
    
    diag(distance) <- 0
    distance <- as.double(distance)
    maxval <- sum(distance[is.finite(distance) & !is.na(distance)])
    distance[is.infinite(distance)] <- maxval
    if (any(is.na(distance))) {
      warning("When solving for geodesic path lengths, NA terms detected in sociomatrix. Substituting.")
      distance[is.na(distance)] <- maxval
    }
    
    dothis <- .C("floyd_warshall",
                 distance=as.double(distance),
                 nn = as.integer(n.pts))
  
    out <- array(dothis$distance,rep(n.pts,2))
    out[out >= maxval] <- Inf
  }
  
  return(out)
}



geodesic.selected.from.edgelist.socio <- function(sociomatrix, sourcepoint=1) {
  dists <- .C("dijkstra_single_socio", input=as.double(sociomatrix),
              output=as.double(rep(0, dim(sociomatrix)[1])),
              dim1=as.integer(dim(sociomatrix)[1]),
              sourcey=as.integer(sourcepoint-1))$output
  return (dists)
}


betweenness.centralities.in.r <- function(sociomatrix, path.weight=c("constant","closeness"), verbose=TRUE) {
  #sociomatrix = eg$grand.socio; path.weight="constant"
  if (all(path.weight == c("constant","closeness"))) path.weight <- "constant"
  
  ccc <- dim(sociomatrix)
  pathcounts <- sapply(1:ccc[1], FUN=function(kk) {
    #print(kk);
    .C("betweenness_geodesic_single_socio", input=as.double(sociomatrix),
       output=as.double(rep(0, ccc[1]^3)),
       path.lengths=as.double(rep(0, ccc[1])),
       dim1=as.integer(dim(sociomatrix)[1]),
       sourcey=as.integer(kk-1))
                     })
  paths <- array(NA, c(ccc[1]^2, ccc[1], ccc[1]))
  dists <- array(NA, ccc)

  for (kk in 1:ccc[1]) {
    paths[,,kk] <- pathcounts[[2+5*(kk-1)]]
    dists[kk,] <- pathcounts[[3+5*(kk-1)]]
  }
  weights <- array(1, ccc); if (path.weight == "closeness") weights <- 1/dists
  diag(weights) <- 0
    
  #nodal <- rep(0, ccc[1])
  edgewise <- array(0, ccc);
  edgesaver <- sapply(1:ccc[1]^2, FUN=function(ii) {
    ss <- floor((ii-1)/ccc[1])+1; dd <- ((ii-1) %% ccc[1])+1
    if (ss != dd) {
      total.paths <- sum(paths[ccc[1]*(dd-1)+1:ccc[1], dd, ss])+1*(sum(paths[ccc[1]*(dd-1)+1:ccc[1], dd, ss])==0)
      out <- paths[,dd,ss]/total.paths/weights[dd,ss]
    } else out <- 0*sociomatrix
    return(out)
  })
  for (kk in 1:length(edgesaver)) edgewise <- edgewise + edgesaver[[kk]];
  
  #for (ss in 1:ccc[1]) for (dd in (1:ccc[2])[-ss]) {
  #  total.paths <- sum(paths[ccc[1]*(dd-1)+1:ccc[1], dd, ss])+1*(sum(paths[ccc[1]*(dd-1)+1:ccc[1], dd, ss])==0)
  #  edgewise <- edgewise + paths[,dd,ss]/total.paths/weights[dd,ss]
  #}
  edgewise <- edgewise/sum(weights)
  nodal <- apply(edgewise, 2, sum) #weight of incoming edges.
    
  return(list(edgewise=edgewise, nodal=nodal)) #paths=paths, dists=dists, 
  
}


betweenness.centralities.socio <- function(sociomatrix, path.weight=c("constant","closeness"), verbose=TRUE) {
  #sociomatrix = eg$grand.socio; path.weight="constant"
  if (all(path.weight == c("constant","closeness"))) path.weight <- "constant"

  fullout <- .C("betweenness_geodesic_full_socio", socio=as.double(sociomatrix),
                          output=as.double(0*sociomatrix),
                          nodeoutput=as.double(rep(0,dim(sociomatrix)[1])),
                          nn=as.integer(dim(sociomatrix)[1]),
                          distweight=as.integer(1*(path.weight=="closeness"))
                          )
    
  return(list(edgewise=fullout$output, nodal=fullout$nodeoutput)) #paths=paths, dists=dists, 
  
}
single.source.betweenness.socio <- function(sociomatrix, sourcepoint=1) {

  ccc <- dim(sociomatrix)
  pathcount.proto <- .C("betweenness_geodesic_single_socio", input=as.double(sociomatrix),
                        output=as.double(rep(0, ccc[1]^3)),
                        path.lengths=as.double(rep(0, ccc[1])),
                        dim1=as.integer(dim(sociomatrix)[1]),
                        sourcey=as.integer(sourcepoint-1))

  pathcount <- array(pathcount.proto$output, c(ccc[1]^2, ccc[1]))

  rownames(pathcount) <- paste(rep(1:ccc[1],ccc[1])-1,"->",sort(rep(1:ccc[1],ccc[1]))-1,", ",sep="")
  colnames(pathcount) <- 1:ccc[1]-1

  return (pathcount)

}
