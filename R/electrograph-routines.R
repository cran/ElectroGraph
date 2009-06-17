#Electro-social routines for R
#How we get importance of ties, and to a point individuals,
#to information communication.
#Andrew C. Thomas
#May 4, 2009

debug.mode <- 0

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
make.sociomatrix.from.edges <- function(inputmat,symmetric=FALSE) {

  if (dim(inputmat)[2]>4) stop("Incorrect dimensionality for edge matrix.")
  
  inputmat <- as.matrix(inputmat)
  id.names <- sort(unique(as.vector(inputmat[,1:2])))
  rows <- length(id.names)

  outmat <- array(0,rep(rows,2))
  cols <- dim(inputmat)[2]

  #relabel.
  new.pairs <- cbind(match(inputmat[,1],id.names),
                     match(inputmat[,2],id.names))
  old.labels <- inputmat[,1:2]
  inputmat[,1:2] <- new.pairs
  
  if (cols==2 & symmetric) for (kk in 1:dim(inputmat)[1]) {
    outmat[inputmat[kk,1],inputmat[kk,2]] <- outmat[inputmat[kk,2],inputmat[kk,1]] <- 1 }
  if (cols==2 & !symmetric) for (kk in 1:dim(inputmat)[1]) {
    outmat[inputmat[kk,1],inputmat[kk,2]] <- 1 }

  if (cols==3 & symmetric) for (kk in 1:dim(inputmat)[1]) {
    outmat[inputmat[kk,1],inputmat[kk,2]] <- outmat[inputmat[kk,2],inputmat[kk,1]] <- inputmat[kk,3] }
  if (cols==3 & !symmetric) for (kk in 1:dim(inputmat)[1]) {
    outmat[inputmat[kk,1],inputmat[kk,2]] <- inputmat[kk,3] }
  
  if (cols==4) for (kk in 1:dim(inputmat)[1]) {
    outmat[inputmat[kk,1],inputmat[kk,2]] <- inputmat[kk,3]
    outmat[inputmat[kk,2],inputmat[kk,1]] <- inputmat[kk,4]
  }
  
  rownames(outmat) <- colnames(outmat) <- id.names
  return(outmat)
}


#create a sociomatrix based on distance-1 adjacency in a lattice.
make.sociomatrix.from.lattice <- function(pts.nby2) {
  n.pts <- dim(pts.nby2)[1]
  return(sapply(1:n.pts,FUN=function(ii)
                sapply(1:n.pts,FUN=function(kk)
                       1*(sum((pts.nby2[ii,]-pts.nby2[kk,])^2)==1)
                       )
                )
         )
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


#assembles weakly connected components.
network.components <- function(sociomatrix,
                       pseudo.diameter.bridge=2) {
  n.pts <- dim(sociomatrix)[1]
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

electrograph <- function(input,symmetric=TRUE,perform.exam=TRUE,
                         enemies.allowed=FALSE, ...) {
  input <- as.matrix(input)
  
  if (dim(input)[1] == dim(input)[2]) {
    sociomatrix <- input
    if (all(is.null(rownames(sociomatrix)))) {
      rownames(sociomatrix) <- colnames(sociomatrix) <- 1:(dim(sociomatrix)[1])
    }
  } else {
    sociomatrix <- make.sociomatrix.from.edges(input,symmetric)
  }

  #reset symmetric now.
  symmetric <- all(sociomatrix==t(sociomatrix))
  nonnegative <- all(sociomatrix >= 0)
  if (!nonnegative & !enemies.allowed) stop ("input contains negative tie strengths, which are disallowed if enemies.allowed is not selected.")

  #if (!nonnegative & enemies.allowed) {
  #  message("For the time being, `enemies' mode is disabled. Setting entries to be nonnegative.")
  #  sociomatrix <- abs(sociomatrix)}
  
  pieces <- network.components(sociomatrix, pseudo.diameter.bridge=2)
  sociomatrices <- list(NA)

  for (kk in 1:length(pieces$components)) {
    sociomatrices[[kk]] <- array(sociomatrix[pieces$components[[kk]],pieces$components[[kk]]],
                                 rep(length(pieces$components[[kk]]),2))
    rownames(sociomatrices[[kk]]) <- colnames(sociomatrices[[kk]]) <-
      rownames(sociomatrix)[pieces$components[[kk]]]               
  }

  out <- list(grand.sociomatrix=sociomatrix, sociomatrices=sociomatrices,
              component.vector=pieces$component.vector,
              geodesic=pieces$geodesic, diameters=pieces$diameters,
              global.pseudo.diameter=pieces$global.pseudo.diameter,
              symmetric=symmetric)
  class(out) <- "electrograph"

  if (debug.mode) message("Finished electrograph step 1.")
  if (perform.exam) {
    out <- electrograph.exam(out, ...)
    if (debug.mode) message("Finished electrograph exam.")
  }
  
  return (out)
}
                         
summary.electrograph <- function(object, ...) {
  #quantities:
  #node level: electro-distance centrality, standard centrality, average current centrality
  #noted: for one extremely tight, isolated pair, this is a bad measure.

  #weighted in/out degrees.

  out.deg <- apply(object$grand.sociomatrix,1,sum)
  in.deg <- apply(object$grand.sociomatrix,2,sum)
  
  hold.geo <- 1/object$geodesic
  diag(hold.geo) <- 0
  out <- NULL
  sp.close <- apply(hold.geo,1,mean)
  
  out <- data.frame(out.deg, in.deg, sp.close)
  nn <- dim(object$grand.sociomatrix)[1]
  
  #construct clustering measures:
  clusts <- .C("clustering_statistics_c",
               socio=as.double(object$grand.sociomatrix),
               nn=as.integer(nn),
               transitives=as.double(rep(0,nn)),
               cycles=as.double(rep(0,nn)))
  transits <- clusts$transitives
  cycles <- clusts$cycles

  out <- data.frame(out,transits,cycles)

  
  if (!is.null(object$distance.mat)) {
    distance.hold <- 1/object$distance.mat
    diag(distance.hold) <- 0
    eg.close <- apply(distance.hold,1,mean)
    out <- data.frame(out,eg.close)
  }

  if (!is.null(object$blk.curr.a)) {
#    current.cent <- apply(object$currents.node,1,mean)
    between.a <- apply(object$blk.curr.a,1,sum)
    between.v <- apply(object$blk.curr.v,1,sum)
    between.p <- apply(object$blk.curr.p,1,sum)
    out <- data.frame(out, between.a, between.v, between.p)
  }
  
  return(out)
  
}

print.electrograph <- function(x, ...) {

  for (kk in 1:length(x$sociomatrices)) {
    writeLines(paste("head(sociomatrices[[",kk,"]]):"))
    print(head(x$sociomatrices[[kk]]))
  }  

  writeLines(paste("component.vector"))
  print(t(x$component.vector))

  writeLines(paste("component.diameters"))
  print(t(x$diameters))
  
  writeLines(paste("Symmetric:",x$symmetric))

  if (!is.null(x$source.sink)) {
    writeLines(paste("source.sink pairs"))
    print(t(x$source.sink))
  }

}

