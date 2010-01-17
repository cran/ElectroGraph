#Ohmic/Electro-social routines for R
#How we get importance of ties, and to a point individuals,
#to information communication.
#Andrew C. Thomas
#Last updated: January 20, 2010

debug.mode <- 1

#produces all ordered pairs between 1 and nn.
pair.sequence <- function(nn) {
  out <- array(NA,c(choose(nn,2),2))
  count <- 0
  for (kk in 1:(nn-1)) { 
    out[count+1:(nn-kk),] <- cbind(kk,(kk+1):nn)
    count <- count+(nn-kk)
  }
  return(out)
}

#inputs an edgelist, a valued edgelist or a valued dyad list.
make.sociomatrix.from.edges <- function(inputmat, symmetric=FALSE, fidelities=NULL) {

  if (dim(inputmat)[2]>4) stop("Incorrect dimensionality for edge matrix.")
  fidelity.dim <- 2-1*(dim(inputmat)[2]<4)
  if (is.null(fidelities)) {
    fidelities <- array(1,c(dim(inputmat)[1],fidelity.dim))
  } else {
    fidelities <- cbind(fidelities)
    if (dim(fidelities)[1]!=dim(inputmat)[1]) stop(paste("Number of rows in the fidelity term,",dim(fidelities)[1]," does not match the rows in the tie matrix,", dim(inputmat)[1], "."))
    #if (fidelity.dim != dim(fidelities)[2]) stop(paste("The number of columns in the fidelity term, does not match the tie matrix.")
    if (fidelity.dim == 2 & dim(fidelities)[2] == 1) fidelities <- cbind(fidelities, fidelities)
    if (fidelity.dim == 1 & dim(fidelities)[2] == 2) stop(paste("There are two columns in the fidelity term when only one is possible."))
  }
  
  inputmat <- as.matrix(inputmat)
  id.names <- sort(unique(as.vector(inputmat[,1:2])))
  rows <- length(id.names)

  outmat <- array(0,rep(rows,2))
  out.fidelity <- array(1,rep(rows,2))
  cols <- dim(inputmat)[2]

  #relabel.
  new.pairs <- cbind(match(inputmat[,1],id.names),
                     match(inputmat[,2],id.names))
  old.labels <- inputmat[,1:2]
  inputmat[,1:2] <- new.pairs
  inputmat <- array(as.numeric(inputmat),dim(inputmat))
  
  if (cols==2 & symmetric) for (kk in 1:dim(inputmat)[1]) {
    outmat[inputmat[kk,1],inputmat[kk,2]] <- outmat[inputmat[kk,2],inputmat[kk,1]] <- 1
    out.fidelity[inputmat[kk,1],inputmat[kk,2]] <- out.fidelity[inputmat[kk,2],inputmat[kk,1]] <- fidelities[kk,1]
  }
  if (cols==2 & !symmetric) for (kk in 1:dim(inputmat)[1]) {
    outmat[inputmat[kk,1],inputmat[kk,2]] <- 1
    out.fidelity[inputmat[kk,1],inputmat[kk,2]] <- fidelities[kk,1]
  }

  if (cols==3 & symmetric) for (kk in 1:dim(inputmat)[1]) {
    outmat[inputmat[kk,1],inputmat[kk,2]] <- outmat[inputmat[kk,2],inputmat[kk,1]] <- inputmat[kk,3]
    out.fidelity[inputmat[kk,1],inputmat[kk,2]] <- out.fidelity[inputmat[kk,2],inputmat[kk,1]] <- fidelities[kk,1]
  }
  if (cols==3 & !symmetric) for (kk in 1:dim(inputmat)[1]) {
    outmat[inputmat[kk,1],inputmat[kk,2]] <- inputmat[kk,3]
    out.fidelity[inputmat[kk,1],inputmat[kk,2]] <- fidelities[kk,1]
  }
  
  if (cols==4) for (kk in 1:dim(inputmat)[1]) {
    outmat[inputmat[kk,1],inputmat[kk,2]] <- inputmat[kk,3]
    outmat[inputmat[kk,2],inputmat[kk,1]] <- inputmat[kk,4]
    out.fidelity[inputmat[kk,1],inputmat[kk,2]] <- fidelities[kk,1]
    out.fidelity[inputmat[kk,2],inputmat[kk,1]] <- fidelities[kk,2]
  }
  
  rownames(outmat) <- colnames(outmat) <- id.names
  rownames(out.fidelity) <- colnames(out.fidelity) <- id.names

  output <- list(sociomatrix=outmat, fidelities=out.fidelity)
  return(output)
}


#create a sociomatrix based on distance-1 adjacency in a lattice.
make.sociomatrix.from.lattice <- function(pts.nby2) {
  n.pts <- dim(pts.nby2)[1]
  sociomatrix <- sapply(1:n.pts,FUN=function(ii)
                        sapply(1:n.pts,FUN=function(kk)
                               1*(sum((pts.nby2[ii,]-pts.nby2[kk,])^2)==1)
                               )
                        )
  fidelities <- array(1,dim(sociomatrix))
  return(list(sociomatrix=sociomatrix, fidelities=fidelities))
  
}

#now in C! It's actually Floyd-Warshall... with absolute values.
geodesic.mat <- function(sociomatrix) {
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
    
    dothis <- .C("dijkstra_geodesic_R",
                 distance=as.double(distance),
                 nn = as.integer(n.pts))
  
    out <- array(dothis$distance,rep(n.pts,2))
    out[out >= maxval] <- Inf
  }
  
  return(out)
}


network.components <- function(sociomatrix,
                                    minimum.relative.strength.for.tie=1e-8) {

  n.pts <- dim(sociomatrix)[1]
  maxtie <- max(sociomatrix)
  sociomatrix[sociomatrix < minimum.relative.strength.for.tie*maxtie] <- 0

  #cascade.
  components <- 1:n.pts

  dothis <- .C("network_components",
               sociomatrix=as.double(sociomatrix),
               nn = as.integer(n.pts),
               comps = as.integer(components))
  outputcomps <- dothis$comps

  comps <- sort(unique(outputcomps))
  components <- list(NA)
  component.vector <- rep(NA, n.pts)

  for (cc in 1:length(comps)) {
    component.vector[outputcomps==comps[cc]] <- cc
    components[[cc]] <- which(outputcomps==comps[cc])
  }

  component.vector <- cbind(component.vector)
  rownames(component.vector) <- rownames(sociomatrix)
  
  out <- list(components=components,
              component.vector=component.vector)
  return(out)
}


#assembles weakly connected components.
network.components.old <- function(sociomatrix,
                               pseudo.diameter.bridge=2,
                               minimum.relative.strength.for.tie=1e-8) {
  
  n.pts <- dim(sociomatrix)[1]
  maxtie <- max(sociomatrix)
  sociomatrix[sociomatrix < minimum.relative.strength.for.tie*maxtie] <- 0
  true.geodesic <- geodesic.mat(sociomatrix)

  sociomatrix <- sociomatrix + t(sociomatrix)
  geodesic <- geodesic.mat(sociomatrix)

  col.check <- 1:n.pts
  components <- list(NA)
  n.components <- 0
  component.vector <- rep(NA,n.pts)
  
  while (length(col.check)>0) {
    checker <- col.check[1]
    in.set <- which(is.finite(geodesic[,checker]))
    n.components <- n.components + 1
    components[[n.components]] <- in.set
    component.vector[in.set] <- n.components
    col.check <- sort(setdiff(col.check,in.set))
  }

  component.vector <- cbind(component.vector)
  rownames(component.vector) <- rownames(sociomatrix)
  
  diams <- sapply(1:n.components,FUN=function(kk) max(geodesic[components[[kk]],components[[kk]]]))
  global.pseudo.diameter <- sum(diams)+pseudo.diameter.bridge*(length(diams)-1)
  out <- list(components=components,
              diameters=diams,
              global.pseudo.diameter=global.pseudo.diameter,
              component.vector=component.vector)

  out$geodesic <- true.geodesic
  return(out)
}


electrograph <- function(input, symmetric=TRUE,
                         solve.for.shortest.paths=TRUE,
                         ohmic.properties=TRUE,
                         fidelities=NULL, verbose=FALSE,
                         substitute.names=NULL,...) {

  input <- as.matrix(input)
  
  if (dim(input)[1] == dim(input)[2] & is.numeric(input)) {
    sociomatrix <- input
    if (all(is.null(rownames(sociomatrix)))) {
      rownames(sociomatrix) <- colnames(sociomatrix) <- 1:(dim(sociomatrix)[1])
    }
    if (!is.null(fidelities)) {
      if (!all(dim(fidelities)==dim(sociomatrix))) stop ("Dimensions of the fidelity matrix do not match the sociomatrix.")
    }
  } else {
    sociomatrix.maker <- make.sociomatrix.from.edges(input, symmetric, fidelities)
    sociomatrix <- sociomatrix.maker$sociomatrix
    fidelities <- sociomatrix.maker$fidelities
  }

  if (is.null(fidelities)) fidelities <- array(1,dim(sociomatrix))
  if (max(abs(fidelities))>1) stop("At least one fidelity value is outside the range (-1,1).")
  
  #reset symmetric now.
  symmetric <- all(sociomatrix==t(sociomatrix))
  nonnegative <- all(sociomatrix >= 0)
  if (!nonnegative) {
    message ("electrograph says: input contains negative tie strengths. Interpreting these as antagonistic connections with fidelity -1.") #To include \"enemies\" ties, use the \"fidelities\" option.")
    fidelities[sociomatrix<0] <- -1
    sociomatrix <- abs(sociomatrix)
  }

  #do we have a vector of names to substitute?
  if (!is.null(substitute.names)) {
    substitute.names <- as.matrix(substitute.names)
    if (dim(substitute.names)[2]!=2) stop (paste("The substituted node name matrix should be k-by-2; it is reading as dimension",dim(substitute.names)[2],"."))

    sub.rows <- match(substitute.names[,1],rownames(sociomatrix))
    rownames(sociomatrix)[sub.rows] <- colnames(sociomatrix)[sub.rows] <- substitute.names[,2]
    rownames(fidelities)[sub.rows] <- colnames(fidelities)[sub.rows] <- substitute.names[,2]
  }


  
  pieces <- network.components(sociomatrix)
  sociomatrices <- list(NA)
  fidelity.pieces <- list(NA)
  
  for (kk in 1:length(pieces$components)) {
    sociomatrices[[kk]] <- array(sociomatrix[pieces$components[[kk]],pieces$components[[kk]]],
                                 rep(length(pieces$components[[kk]]),2))
    rownames(sociomatrices[[kk]]) <- colnames(sociomatrices[[kk]]) <-
      rownames(sociomatrix)[pieces$components[[kk]]]

    fidelity.pieces[[kk]] <- array(fidelities[pieces$components[[kk]],pieces$components[[kk]]],
                                 rep(length(pieces$components[[kk]]),2))
    rownames(fidelity.pieces[[kk]]) <- colnames(fidelity.pieces[[kk]]) <-
      rownames(sociomatrix)[pieces$components[[kk]]]

    if (debugmode) print(fidelity.pieces[[kk]])
  }
  
  out <- list(grand.sociomatrix=sociomatrix,
              sociomatrices=sociomatrices,
              grand.fidelity=fidelities,
              fidelities=fidelity.pieces,
              components=pieces$components,
              component.vector=pieces$component.vector,
              symmetric=symmetric)

  #geodesic=pieces$geodesic, diameters=pieces$diameters,
  #global.pseudo.diameter=pieces$global.pseudo.diameter,
              
  class(out) <- "electrograph"

  if (verbose) message("Finished electrograph object loading.")

  if (solve.for.shortest.paths) {
    out$geodesic <- geodesic.mat(sociomatrix)
    if (verbose) message("Finished geodesic path length calculation.")
  }
  
  if (ohmic.properties) {
    out <- electrograph.exam(out, ...)
    if (debug.mode) message("Finished Ohmic properties.")
  }
  
  return (out)
}

clustering.statistics <- function(sociomatrix) {
  nn <- dim(sociomatrix)[1]
  clusts <- .C("clustering_statistics_c",
               socio=as.double(sociomatrix),
               nn=as.integer(nn),
               transitives=as.double(rep(0,nn)),
               cycles=as.double(rep(0,nn)))
  out <- data.frame(transitives=clusts$transitives,
                    cycles=clusts$cycles)
  return(out)
}

summary.electrograph <- function(object, ...) {
  #quantities:
  #node level: electro-distance centrality, standard centrality, average current centrality
  #noted: for one extremely tight, isolated pair, this is a bad measure.

  #weighted in/out degrees.

  #writeLines ("Wake up, Neo.")
  
  nn <- dim(object$grand.sociomatrix)[1]

  out.deg <- apply(object$grand.sociomatrix,1,sum)
  in.deg <- apply(object$grand.sociomatrix,2,sum)
  #construct clustering measures:
  clusts <- .C("clustering_statistics_c",
               socio=as.double(object$grand.sociomatrix),
               nn=as.integer(nn),
               transitives=as.double(rep(0,nn)),
               cycles=as.double(rep(0,nn)))
  transits <- clusts$transitives
  cycles <- clusts$cycles

  out <- data.frame(out.deg, in.deg, transits, cycles)

  
#  out <- data.frame(out,transits,cycles)

  if (!is.null(object$geodesic)) {
    hold.geo <- 1/object$geodesic
    diag(hold.geo) <- 0
    geo.close.out <- apply(hold.geo,1,mean)
    geo.close.in <- apply(hold.geo,2,mean)
    out <- data.frame(out, geo.close.out, geo.close.in)
  }
  
  
  if (!is.null(object$distance.mat)) {
    distance.hold <- 1/object$distance.mat
    diag(distance.hold) <- 0
    ohm.close.out <- apply(distance.hold,1,mean)
    ohm.close.in <- apply(distance.hold,2,mean)
    out <- data.frame(out, ohm.close.out, ohm.close.in)
  }

  if (!is.null(object$avg.current.black.a)) {
#    current.cent <- apply(object$currents.node,1,mean)
    between.a <- apply(object$avg.current.black.a,1,sum)
    between.v <- apply(object$avg.current.black.v,1,sum)
    between.p <- apply(object$avg.current.black.p,1,sum)
    out <- data.frame(out, between.a, between.v, between.p)
  }
  
  return(out)
  
}

print.electrograph <- function(x, ...) {

  for (kk in 1:length(x$sociomatrices)) {
    writeLines(paste("head(sociomatrices[[",kk,"]]):"))
    dd <- min(dim(x$sociomatrices[[kk]]),6)
    print(x$sociomatrices[[kk]][1:dd, 1:dd])
  }

  writeLines(paste("component.vector"))
  print(t(x$component.vector))

  #writeLines(paste("component.diameters"))
  #print(t(x$diameters))
  
  writeLines(paste("Symmetric:",x$symmetric))

  if (!is.null(x$source.sink)) {
    writeLines(paste("source.sink pairs"))
    print(t(x$source.sink))
  }

}

