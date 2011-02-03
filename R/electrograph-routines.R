#Ohmic/Electro-social routines for R
#How we get importance of ties, and to a point individuals,
#to information communication.
#Andrew C. Thomas
#Last updated: January 30, 2011


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

#is an edge matrix representative of a symmetric system?
sym.test <- function(nby3) {
  #nby3 = eg$nby3
  nn <- length(unique(c(nby3[,1:2])))
  ee <- dim(nby3)[1]
  sym <- .C("sym_test", edges=as.integer(nby3[,1:2]),
            strengths=as.double(nby3[,3]),
            nn=as.integer(nn), ee=as.integer(ee),
            sym=as.integer(0))$sym
  return(sym==0)
}

#reprocess the edge matrix to the current format.
edge.processor <- function(inputmat, symmetric, fidelities) {
  #returns edge values as well.
  
  inputmat <- as.matrix(inputmat)

  if (dim(inputmat)[2]>4) stop("Incorrect dimensionality for edge matrix.")
  fidelity.dim <- 2-1*(dim(inputmat)[2]<4)
  if (is.null(fidelities)) {
    fidelities <- array(1,c(dim(inputmat)[1],fidelity.dim))
  } else {
    fidelities <- cbind(fidelities)
    
    if (dim(fidelities)[1]!=dim(inputmat)[1]) stop(paste("Number of rows in the fidelity term,",dim(fidelities)[1]," does not match the rows in the edge matrix,", dim(inputmat)[1], "."))

    if (fidelity.dim == 2 & dim(fidelities)[2] == 1) fidelities <- cbind(fidelities, fidelities)
    if (fidelity.dim == 1 & dim(fidelities)[2] == 2) stop(paste("There are two columns in the fidelity term when only one is possible."))
  }
  
  id.names <- sort(unique(as.vector(inputmat[,1:2])))
  nn <- length(id.names)
  cols <- dim(inputmat)[2]

  #relabel.
  new.pairs <- cbind(match(inputmat[,1],id.names),
                     match(inputmat[,2],id.names))
  old.labels <- inputmat[,1:2]
  if (!any(cols==2:4)) stop("Edge matrix must have 2, 3 or 4 columns.")
  if (cols==4) {new.pairs <- rbind(new.pairs, new.pairs[,2:1]); values <- c(inputmat[,3:4])} else {
    if (!symmetric) {
      if (cols==2) values <- rep(1, length(fidelities))
      if (cols==3) values <- inputmat[,3]
    } else {
      fidelities <- c(fidelities, fidelities)
      if (cols==3) values <- c(inputmat[,3], inputmat[,3]) else values <- rep(1, length(fidelities))
      new.pairs <- rbind(new.pairs, new.pairs[,2:1])
    }
  }
  fidelities[which(values<0)] <- -1*fidelities[which(values<0)]
  values <- abs(values)
  nby3 <- cbind(new.pairs, values)
                
  colnames(nby3) <- c("src","dest","value")
  return(list(nby3=nby3, fidelities=fidelities, id.names=id.names))
  
}



sociomatrix.to.nby3.edges <- function(sociomatrix) {
  if (dim(sociomatrix)[1]!=dim(sociomatrix)[2]) stop ("Input sociomatrix is not square.")
  nn <- dim(sociomatrix)[1]
  nonzeros <- which(sociomatrix > 0)
  return(cbind(floor((nonzeros-1)/nn)+1, (nonzeros-1) %% nn + 1, sociomatrix[nonzeros]))
}


#inputs an edgelist, a valued edgelist or a valued dyad list. Outputs the sociomatrix/fidelity matrix.
make.sociomatrix.from.edges <- function(inputmat, symmetric=FALSE, fidelities=NULL) {
  #inputmat=e.graph$nby3[edge.subset,];  fidelities=e.graph$fidelities[edge.subset]; symmetric=FALSE;
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
  cols <- dim(inputmat)[2]

  #relabel.
  new.pairs <- cbind(match(inputmat[,1],id.names),
                     match(inputmat[,2],id.names))
  old.labels <- inputmat[,1:2]
  inputmat[,1:2] <- new.pairs
  inputmat <- array(as.numeric(inputmat),dim(inputmat))
  outmat <- array(0, rep(rows,2)); out.fidelity <- array(1, rep(rows,2))
  
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




network.components.socio <- function(sociomatrix,
                                    minimum.relative.strength.for.tie=1e-8) {

  n.pts <- dim(sociomatrix)[1]
  maxtie <- max(sociomatrix)
  sociomatrix[sociomatrix < minimum.relative.strength.for.tie*maxtie] <- 0

  #cascade.
  components <- 1:n.pts

  dothis <- .C("network_components_socio",
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


network.components.edges <- function(nby3,
                                     minimum.relative.strength.for.tie=1e-8) {

  pointset <- unique(c(nby3[,1:2]))
  maxtie <- max(nby3[,3])
  edgeset <- which(nby3[,3] > minimum.relative.strength.for.tie*maxtie)
  nby3 <- nby3[edgeset,]
  
  #cascade.
  components <- 1:length(pointset)

  outputcomps <- .C("network_components_edges",
               edgelist=as.integer(as.vector(c(nby3[,1:2]-1))),
               comps = as.integer(components),
               nn = as.integer(length(pointset)),
               edgecount = as.integer(dim(nby3)[1]))$comps

  comps <- sort(unique(outputcomps))
  components <- list(NA)
  component.vector <- rep(NA, length(pointset))

  for (cc in 1:length(comps)) {
    component.vector[outputcomps==comps[cc]] <- cc
    components[[cc]] <- which(outputcomps==comps[cc])
  }

  component.vector <- cbind(component.vector)
  
  out <- list(components=components,
              component.vector=component.vector)
  return(out)
}




#now, edge-dominant version.
electrograph <- function(input, make.edgelist.symmetric=TRUE, solve.for.shortest.paths=TRUE,
                         solve.for.betweenness=FALSE, solve.for.canadian.betweenness=FALSE,
                         ohmic.properties=FALSE, fidelities=NULL, verbose=FALSE,
                         substitute.names=NULL, ...) {
  #input=cbind(c("ee","aa","bb","cc","dd"),c("ff","bb","cc","dd","aa")); symmetric=TRUE; solve.for.shortest.paths=TRUE; ohmic.properties=TRUE; fidelities=NULL; verbose=FALSE; substitute.names=NULL
  
  input <- as.matrix(input)
  
  if (dim(input)[1] == dim(input)[2] & is.numeric(input)) {
    #then it's a sociomatrix.
    nby3 <- sociomatrix.to.nby3.edges(input);  pop.size <- length(unique(c(nby3[,1:2]))) #dim(sociomatrix)[1]
    if (all(is.null(rownames(input)))) {node.ids <- 1:pop.size} else {node.ids <- rownames(input)}
    
    if (!is.null(fidelities)) {
      if (!all(dim(fidelities)==dim(input))) stop ("Dimensions of the fidelity matrix do not match the sociomatrix.")
      fidelities <- fidelities[which(input != 0)]
    }
  } else {
    holder <- edge.processor (input, make.edgelist.symmetric, fidelities)
    nby3 <- holder$nby3;    node.ids <- holder$id.names;
    fidelities <- holder$fidelities
  }

  if (is.null(fidelities)) fidelities <- rep(1, dim(nby3)[1])
  if (max(abs(fidelities))>1) stop("At least one fidelity value is outside the range (-1,1).")
  
  #reset symmetric now.
  symmetric <- sym.test(nby3)

  nonnegative <- all(nby3[,3] >= 0)
  if (!nonnegative) {
    message ("electrograph says: input contains negative tie strengths. Interpreting these as antagonistic connections with fidelity -1.") #To include \"enemies\" ties, use the \"fidelities\" option.")
    fidelities[which(nby3[,3]<0)] <- -1*fidelities[which(nby3[,3]<0)]
    nby3[,3] <- abs(nby3[,3])
  }

                                        #do we have a vector of names to substitute?
  if (!is.null(substitute.names)) {
    substitute.names <- as.matrix(substitute.names)
    if (dim(substitute.names)[2]!=2) stop (paste("The substituted node name matrix should be k-by-2; it is reading as dimension",dim(substitute.names)[2],"."))
 
    sub.vals <- match(substitute.names[,1], node.ids)
    node.ids[sub.vals] <- substitute.names[,2]
  }

  #out <- list(components=components, component.vector=component.vector)
  pieces <- network.components.edges(nby3)
  component.vector <- pieces$component.vector; rownames(component.vector) <- node.ids

  out <- list(nby3=nby3, fidelities=fidelities,
              component.vector=component.vector,
              symmetric=symmetric,
              node.ids=node.ids)

  #outputs: nby3 edgelist, fidelities, component.vector, symmetric?, node.ids

  
              #grand.sociomatrix=sociomatrix,
              #sociomatrices=sociomatrices,
              #grand.fidelity=fidelities,
              #fidelities=fidelity.pieces,
              #components=pieces$components,
              #component.vector=pieces$component.vector,
              #symmetric=symmetric)

  #geodesic=pieces$geodesic, diameters=pieces$diameters,
  #global.pseudo.diameter=pieces$global.pseudo.diameter,
              
  class(out) <- "electrograph"

  #if (verbose) message("Finished electrograph object loading.")

  if (solve.for.shortest.paths) {
    out$geodesic <- geodesic.mat(nby3, node.ids=node.ids)
    if (verbose) message("Finished geodesic path length calculation.")
  }

  if (solve.for.betweenness) {
    bet.hold <- betweenness.centralities(nby3, node.ids=node.ids)
    out$edge.betweenness <- bet.hold$edgewise
    out$node.betweenness <- bet.hold$nodal
    if (verbose) message("Finished standard betweenness centrality calculation.")
  }

  if (solve.for.canadian.betweenness) {
    out$canadian.betweenness <- canadian.betweenness.full(nby3)
    if (verbose) message("Finished Canadian betweenness centrality calculation.")
  }
  
  if (ohmic.properties) {
    out2 <- ohmic.properties.all(out, ...)
    if (verbose) message("Finished Ohmic properties.")
  }
  
  return (out)
}


#gets transitivity/cycles from a 
clustering.statistics.socio <- function(sociomatrix) {
  nn <- dim(sociomatrix)[1]
  clusts <- .C("clustering_statistics_socio",
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
  
  nn <- length(object$node.ids)          #dim(object$grand.sociomatrix)[1]

  intermediates <- .C("short_length_statistics",
                      edges=as.integer(c(object$nby3[,1:2])-1),
                      values=as.double(object$nby3[,3]),
                      pnn=as.integer(nn),
                      pedge=as.integer(dim(object$nby3)[1]),
                      indegree=as.double(rep(0,nn)),
                      outdegree=as.double(rep(0,nn)),
                      transitives=as.double(rep(0,nn)),
                      cycles=as.double(rep(0,nn)))
  
  out <- data.frame(out.deg=intermediates$outdegree,
                    in.deg=intermediates$indegree,
                    transits=intermediates$transitives,
                    cycles=intermediates$cycles)
  
#  out <- data.frame(out,transits,cycles)

  if (!is.null(object$geodesic)) {
    hold.geo <- 1/object$geodesic
    diag(hold.geo) <- 0
    geo.close.out <- apply(hold.geo,1,mean)
    geo.close.in <- apply(hold.geo,2,mean)
    out <- data.frame(out, geo.close.out, geo.close.in)
  }

  if (!is.null(object$node.bet)) {
    out <- data.frame(out, between=object$node.bet)
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

  message("Head of edge list nby3:")
  print(head(x$nby3))

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

