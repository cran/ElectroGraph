#electro-social-specific commands.
#Andrew C. Thomas
#January 20, 2010

#notation:
#voltages, currents in each node suppressed.
#avg.current, red.current obsolete.
#now: blk/red.curr.a/v
#
#transforms: avg.current.black/red.a/v/p
#

debugmode <- 0

solve.twosome <- function(sociomatrix, fidelity, src, snk) {
  if (any(dim(sociomatrix)!=2))
    stop("solve.twosome should only be used on two-by-two sociomatrices.")
  resist <- rep(Inf,length(src))
  voltages <- currents <- array(0,c(2,length(src)))
  avg.current.black.a <- avg.current.black.v <- avg.current.black.p <-
    avg.current.red.a <- avg.current.red.v <- avg.current.red.p <- array(0,rep(2,2))
  fidelity.out <- rep(1,length(src))
  
  for (kk in 1:length(src)) if (src[kk]!=snk[kk]) {
    resist[kk] <- 1/sociomatrix[src[kk],snk[kk]]
    fidelity.out[kk] <- fidelity[src[kk],snk[kk]]
    
    if (src[kk]==1) {
      voltages[,kk] <- c(1,0)
      currents[,kk] <- c(sociomatrix[src[kk],snk[kk]],0)
    } else {
      voltages[,kk] <- c(0,1)
      currents[,kk] <- c(0,sociomatrix[src[kk],snk[kk]])
    }
    
    avg.current.black.a[src[kk],snk[kk]] <- sociomatrix[src[kk],snk[kk]]*(1+fidelity[src[kk],snk[kk]])/2
    avg.current.red.a[src[kk],snk[kk]] <- sociomatrix[src[kk],snk[kk]]*(1-fidelity[src[kk],snk[kk]])/2

    avg.current.black.v[src[kk],snk[kk]] <- avg.current.black.a[src[kk],snk[kk]]/resist[kk]
    avg.current.red.v[src[kk],snk[kk]] <- avg.current.red.a[src[kk],snk[kk]]/resist[kk]
    avg.current.black.p[src[kk],snk[kk]] <- avg.current.black.a[src[kk],snk[kk]]/sqrt(resist[kk])
    avg.current.red.p[src[kk],snk[kk]] <- avg.current.red.a[src[kk],snk[kk]]/sqrt(resist[kk])
  }
  
  out <- list(resist=resist,
              voltages=voltages,
              currents=currents,

              avg.current.black.a=avg.current.black.a,
              avg.current.red.a=avg.current.red.a,
              avg.current.black.v=avg.current.black.v,
              avg.current.red.v=avg.current.red.v,        
              avg.current.black.p=avg.current.black.p,
              avg.current.red.p=avg.current.red.p,        

              fidelity=fidelity.out)
  
  return(out)
}

#as yet to be converted.
#one source-sink pair, mostly for testing.
compute.volts.symmetric <- function(sociomatrix, src, snk, fidelities=NULL) {

  if (is.null(fidelities)) fidelities = array(1,dim(sociomatrix))
  if (any(sociomatrix!=t(sociomatrix))) stop ("Symmetric matrix routine given nonsymmetric input.")
  n.pts <- dim(sociomatrix)[1]
  if (!any(src==1:n.pts))
    stop ("In compute.volts.symmetric, source isn't one of the nodes indicated.")
  if (!any(snk==1:n.pts))
    stop ("In compute.volts.symmetric, sink isn't one of the nodes indicated.")
  if (src==snk) stop ("In compute.volts.symmetric, source and sink can't be identical.")

  voltages <- currents <- rep(0,n.pts)
  resist.eq <- 1234567891; fidelity <- 1

  if (n.pts>2) {
    res <- .C("solve_volts_symmetric_c",
              as.double(sociomatrix), as.double(fidelities),
              as.integer(n.pts),
              as.integer(src-1), as.integer(snk-1),
              resist.eq=as.double(resist.eq),
              fidelity=as.double(fidelities),
              voltages=as.double(voltages), currents=as.double(currents))
  } else {
    res <- solve.twosome(sociomatrix, fidelities,
                         src, snk)
  }
  return(list(resist.eq=res$resist, voltages=res$volt,
              current.out=res$curr, fidelity=res$fidelity))
}


compute.volts.asymmetric <- function(sociomatrix, src, snk, fidelities = NULL) {

  if (is.null(fidelities)) fidelities = array(1,dim(sociomatrix))
  n.pts <- dim(sociomatrix)[1]
  if (!any(src==1:n.pts)) stop ("In compute.volts.asymmetric, source isn't one of the nodes indicated.")
  if (!any(snk==1:n.pts)) stop ("In compute.volts.asymmetric, sink isn't one of the nodes indicated.")
  if (src==snk) stop ("In compute.volts.asymmetric, source and sink can't be identical.")

  voltages <- currents <- rep(0,n.pts)
  resist.eq <- 1234567890; fidelity <- 1
  
  if (n.pts>2) {
    res <- .C("solve_volts_notsymmetric_c",
              as.double(sociomatrix), as.double(fidelities),
              as.integer(n.pts),
              as.integer(src-1), as.integer(snk-1),
              resist.eq=as.double(resist.eq), fidelity=as.double(fidelity),
              voltages=as.double(voltages), currents=as.double(currents))
  } else {
    res <- solve.twosome(sociomatrix, fidelities,
                         src, snk)
  }
  return(list(resist.eq=res$resist, voltages=res$volt,
              current.out=res$curr, fidelity=res$fidelity))
}


#assumes nondirectional for now, but can do enemies... sort of.
#worry about directional in v2.0.
#currently not used.
get.resistance.really.fast <- function(sociomatrix, fidelities,
                                       sources, sinks) {
  #sociomatrix <- test.e$grand.socio; s.s <- pair.sequence(dim(sociomatrix)[1])
  #sources <- s.s[,1]; sinks <- s.s[,2]
  #outcomes: resist, fidelity, volt, curr, avg.current, red.current

  nn <- dim(sociomatrix)[1]
  s.s <- cbind(sources,sinks)
  
  out.resist <- rep(NA,dim(s.s)[1])
  out.fid <- rep(1,dim(s.s)[1])
  avg.current.black.a <- avg.current.black.v <- avg.current.red.a <- avg.current.red.v <- array(0,rep(nn,2))
  
  big.laplacian <- -sociomatrix
  diag(big.laplacian) <- apply(sociomatrix, 1, sum)

  #fixes jams with disconnected pieces.
  diag(big.laplacian)[diag(big.laplacian)==0] <- 1  

  #all non-1, 1-non-2, 1-2-non-3.
  ll.1 <- solve(big.laplacian[-1,-1])
  ll.2 <- solve(big.laplacian[-2,-2])
  ll.3 <- solve(big.laplacian[-3,-3])

  writeLines("Inversions done.")
  
  enemies.present <- any(sociomatrix<0)
  
  #process.
  for (kk in 1:dim(s.s)[1]) {
    voltages <- rep(0, nn)
    if (!any(s.s[kk,]==1)) {
      voltages[-1] <- ll.1[,s.s[kk,1]-1]-ll.1[,s.s[kk,2]-1]
    } else {
      if (any(s.s[kk,]==2)) {
        voltages[-3] <- ll.3[,s.s[kk,1]]-ll.3[,s.s[kk,2]]
      } else {
        ss.hold <- s.s[kk,]-1
        ss.hold[ss.hold==0] <- 1
        voltages[-2] <- ll.2[,ss.hold[1]]-ll.2[,ss.hold[2]]
      }
    }

    out.resist[kk] <- max(voltages)-min(voltages)

    if (enemies.present) {
      #feeding the 1A version.
      enemy.back <- .C("process_for_enemies_c",
                       sociomatrix=as.double(sociomatrix),
                       fidelities=as.double(fidelities),
                       voltages=as.double(voltages),
                       avg.current=as.double(rep(0,nn^2)),
                       red.current=as.double(rep(0,nn^2)),
                       tot.res=as.double(out.resist[kk]),
                       tot.fid=as.double(1),
                       nn=as.integer(nn))
      
      out.fid[kk] <- enemy.back$tot.fid
      avg.current.black.a <- avg.current.black.a + enemy.back$avg.current
      avg.current.red.a <- avg.current.red.a + enemy.back$red.current
      avg.current.black.v <- avg.current.black.v + enemy.back$avg.current/out.resist[kk]
      avg.current.red.v <- avg.current.red.v + enemy.back$red.current/out.resist[kk]
      
    } else {
      
      vv.hold <- array(voltages,c(nn,nn))-t(array(voltages,c(nn,nn)))
      avg.current <- vv.hold*sociomatrix*(vv.hold>0)
      avg.current.black.a <- avg.current.black.a + avg.current
      avg.current.black.v <- avg.current.black.v + avg.current/out.resist[kk]
    }

    if (kk %% 1000 == 0) print(kk)
  }

  #avg.current.black.a <- avg.current.black.a/length(sources)
  #avg.current.red.a <- avg.current.red.a/length(sources)
  #avg.current.black.v <- avg.current.black.v/length(sources)
  #avg.current.red.v <- avg.current.red.v/length(sources)
  
  #outcomes: resist, fidelity, volt, curr, avg.current, red.current
  return(list(resist=out.resist,
              fidelity=out.fid,
              volt=0, curr=0,
              
              avg.current.black.a=avg.current.black.a,
              avg.current.red.a=avg.current.red.a,
              avg.current.black.v=avg.current.black.v,
              avg.current.red.v=avg.current.red.v,

              ll.1=ll.1,
              ll.2=ll.2,
              ll.3=ll.3,
              
              source.sink=s.s))
}





get.resistance <- function(e.graph, sources, sinks) {
  #e.graph <- electrograph(array(c(0,0,0,1,0,0,0,1,0),rep(3,2))); sources=1; sinks=3

  if (debugmode) message("get.resistance")
  
  if (class(e.graph)!="electrograph")
    stop ("get.resistance takes an object of class electrograph.")
  if (length(sources)!=length(sinks))
    stop ("Sources and sinks have different lengths.")
  if (any(sources==sinks))
    stop ("At least one source node is identical to its sink node.")

  n.pts <- dim(e.graph$grand.socio)[1]

  #are the sources and sinks in the same component? 
  source.sink.convert <- cbind(match(sources,rownames(e.graph$component.vector)),
                               match(sinks,rownames(e.graph$component.vector)))
  
  if (any(is.na(source.sink.convert))) {
    warning("Some specified nodes do not appear to exist in the electrograph object.")
    thing <- apply(source.sink.convert,1,any,is.na)
    source.sink.convert <- source.sink.convert[!thing,]
    sources <- sources[!thing]; sinks <- sinks[!thing]
  }
  if (length(as.vector(source.sink.convert))==0)
    stop("In fact, none of ths specified arcs appear to exist in the electrograph object.")

  source.components <- e.graph$component.vector[source.sink.convert[,1]]*
    (e.graph$component.vector[source.sink.convert[,1]]==
     e.graph$component.vector[source.sink.convert[,2]])
  groups <- unique(source.components)
  #identified: (se,sk) with same component.

  if (debugmode) print(source.components)

  #networks over 150: too big to store extra info.
  
  out.conducts <- out.resists <- out.fidelity <- rep(0,length(sources))
  out.avg.current.black.a <- out.avg.current.red.a <-
    out.avg.current.black.v <- out.avg.current.red.v <-
      out.avg.current.black.p <- out.avg.current.red.p <- array(0,rep(n.pts,2))
  
  out.primers <- array(0,c(n.pts,n.pts,3))
  
  for (kk in groups) if (kk>0) {
    if (debugmode) message(paste("in ",kk))
    subset <- which(source.components==kk)
    node.subset <- match(rownames(e.graph$sociomatrices[[kk]]),
                         rownames(e.graph$grand.sociomatrix))
    if (debugmode) {print(subset); print(node.subset)}
    sources.true <- sources[subset]
    sinks.true <- sinks[subset]
    n.pts.piece <- dim(e.graph$sociomatrices[[kk]])[1]
    src <- match(sources.true,rownames(e.graph$sociomatrices[[kk]]))
    snk <- match(sinks.true,rownames(e.graph$sociomatrices[[kk]]))

    this.one.symmetric <- all(e.graph$sociomatrices[[kk]] == t(e.graph$sociomatrices[[kk]]))
    
    if (length(node.subset)==2) {
      #src snk must be (1,2), (2,1) or both.
      result <- solve.twosome(e.graph$sociomatrices[[kk]], e.graph$fidelities[[kk]],
                              src, snk)
      out.conducts[subset] <- 1/result$resist
      out.resists[subset] <- result$resist
      out.fidelity[subset] <- result$fidelity

      out.avg.current.black.a[node.subset,node.subset] <- result$avg.current.black.a
      out.avg.current.red.a[node.subset,node.subset] <- result$avg.current.red.a
      out.avg.current.black.v[node.subset,node.subset] <- result$avg.current.black.v
      out.avg.current.red.v[node.subset,node.subset] <- result$avg.current.red.v
      out.avg.current.black.p[node.subset,node.subset] <- result$avg.current.black.p
      out.avg.current.red.p[node.subset,node.subset] <- result$avg.current.red.p

    } else {
    
      #out.avg.current.piece <- out.red.current.piece <- rep(-7,n.pts.piece^2)

      #result <- get.resistance.really.fast (e.graph$sociomatrices[[kk]], src, snk)

      #don't bother doing disconnected.
      disco <- sapply(1:length(src), FUN=function(kk) {
        kk*is.finite(e.graph$geodesic[node.subset,node.subset][src[kk],snk[kk]])
      })
      disco <- disco[disco>0]
      newsnk <- snk[disco]; newsrc <- src[disco]
      #print(rbind(newsrc,newsnk))
      
      out.conducts.piece <- out.resists.piece <-
        out.fidelity.piece <- rep(-3,length(newsrc))

      if (debugmode) print(e.graph$fidelities[[kk]])
      
      if (this.one.symmetric) {
        result <- .C("get_resistances_fast_symmetric_c",
                     as.double(e.graph$sociomatrices[[kk]]),
                     as.double(e.graph$fidelities[[kk]]),

                     as.integer(n.pts.piece),
                     as.integer(newsrc-1), as.integer(newsnk-1),
                     as.integer(length(newsrc)),
                     
                     resist=as.double(out.resists.piece),
                     fidelity=as.double(out.fidelity.piece),
                   
                     avg.current.black.a=as.double(rep(-7,n.pts.piece^2)),
                     avg.current.red.a=as.double(rep(-7,n.pts.piece^2)),
                     avg.current.black.v=as.double(rep(-7,n.pts.piece^2)),
                     avg.current.red.v=as.double(rep(-7,n.pts.piece^2)),
                     avg.current.black.p=as.double(rep(-7,n.pts.piece^2)),
                     avg.current.red.p=as.double(rep(-7,n.pts.piece^2)),
                     
                     ll.1=as.double(rep(-9,(n.pts.piece-1)^2)),
                     ll.2=as.double(rep(-9,(n.pts.piece-1)^2)),
                     ll.3=as.double(rep(-9,(n.pts.piece-1)^2))               
                     )
      } else {
        result <- .C("get_resistances_c",
                     as.double(e.graph$sociomatrices[[kk]]),
                     as.double(e.graph$fidelities[[kk]]),

                     as.integer(n.pts.piece),
                     as.integer(newsrc-1), as.integer(newsnk-1),
                     as.integer(length(newsrc)),

                     resist=as.double(out.resists.piece),
                     fidelity=as.double(out.fidelity.piece),
                   
                     avg.current.black.a=as.double(rep(-7,n.pts.piece^2)),
                     avg.current.red.a=as.double(rep(-7,n.pts.piece^2)),
                     avg.current.black.v=as.double(rep(-7,n.pts.piece^2)),
                     avg.current.red.v=as.double(rep(-7,n.pts.piece^2)),
                     avg.current.black.p=as.double(rep(-7,n.pts.piece^2)),
                     avg.current.red.p=as.double(rep(-7,n.pts.piece^2)))
        result$ll.1 <- result$ll.2 <- result$ll.3 <- 0
      }
      
      #writeLines("out")

      #print(c(length(result$resist),length(disco)))
      
      out.conducts[subset[disco]] <- 1/result$resist                   
      out.resists[subset[disco]] <- result$resist
      out.fidelity[subset[disco]] <- result$fidelity
      out.conducts[subset[-disco]] <- 0                   
      out.resists[subset[-disco]] <- Inf
      out.fidelity[subset[-disco]] <- 0
      
    
      #if (debugmode) {cat("thru 4",kk,"\n"); print(dim(out.voltages));
      #  print(class(out.voltages));}

      out.avg.current.black.a[node.subset,node.subset] <- result$avg.current.black.a
      out.avg.current.red.a[node.subset,node.subset] <- result$avg.current.red.a
      out.avg.current.black.v[node.subset,node.subset] <- result$avg.current.black.v
      out.avg.current.red.v[node.subset,node.subset] <- result$avg.current.red.v
      out.avg.current.black.p[node.subset,node.subset] <- result$avg.current.black.p
      out.avg.current.red.p[node.subset,node.subset] <- result$avg.current.red.p

      out.primers[node.subset[-1],node.subset[-1],1] <- result$ll.1
      out.primers[node.subset[-2],node.subset[-2],2] <- result$ll.2
      out.primers[node.subset[-3],node.subset[-3],3] <- result$ll.3
      
      
    }
  } else {
    
    subset <- which(source.components==kk)
    sources.true <- source.sink.convert[source.components==kk,1]
    sinks.true <- source.sink.convert[source.components==kk,2]
    out.resists[subset] <- Inf
    
  }

  #fix to averages, zeroes.
  fixer <- function(block) {
    block[block<0] <- 0
    block <- block/length(sources)
    return(block)
  }
  out.avg.current.black.a <- fixer(out.avg.current.black.a)
  out.avg.current.red.a <- fixer(out.avg.current.red.a)
  out.avg.current.black.v <- fixer(out.avg.current.black.v)
  out.avg.current.red.v <- fixer(out.avg.current.red.v)
  out.avg.current.black.p <- fixer(out.avg.current.black.p)
  out.avg.current.red.p <- fixer(out.avg.current.red.p)
  
  return(list(conductances=out.conducts,
              resistances=out.resists,
              fidelities=out.fidelity,
              voltages=0,#out.voltages,
              currents.node=0,#out.currents,
              source.sink=cbind(sources,sinks),

              avg.current.black.a=out.avg.current.black.a,
              avg.current.red.a=out.avg.current.red.a,
              avg.current.black.v=out.avg.current.black.v,
              avg.current.red.v=out.avg.current.red.v,
              avg.current.black.p=out.avg.current.black.p,
              avg.current.red.p=out.avg.current.red.p,

              primers=out.primers
              ))
}



#Now done by default when calling the constructor ``electrograph''.
electrograph.exam <- function(e.graph, sample.edges=FALSE, sample.per=NULL) {

  if (debugmode) message ("electrograph.exam")
  
  #eg.hold <- e.graph
  if (any(e.graph$grand.socio != t(e.graph$grand.socio))) {
    message("Note: Directed/asymmetric graphs are solved at considerably slower speed than symmetric graphs.")
    #e.graph$grand.sociomatrix <- (e.graph$grand.sociomatrix + t(e.graph$grand.sociomatrix))/2
    #for (kk in 1:length(e.graph$sociomat)) e.graph$sociomatrices[[kk]] <- (e.graph$sociomatrices[[kk]] + t(e.graph$sociomatrices[[kk]]))/2
  }

  
  n.pts <- dim(e.graph$grand.socio)[1]
  if (is.null(sample.per)) sample.per <- round(sqrt(n.pts))
  if (sample.edges) {
    e1 <- sort(rep(1:n.pts,sample.per))  
    e2 <- rep(NA,n.pts*sample.per)
    for (kk in 1:n.pts) {
      e2[(kk-1)*sample.per+1:sample.per] <- sample((1:n.pts)[-kk],sample.per)
    }
    edges <- cbind(e1,e2)
  } else {
    edges <- pair.sequence(n.pts)
    if (any(e.graph$grand.socio != t(e.graph$grand.socio))) edges <- rbind(edges,edges[,2:1])
  }
  
  sources <- rownames(e.graph$grand.socio)[edges[,1]]
  sinks <- rownames(e.graph$grand.socio)[edges[,2]]

  
  out <- get.resistance(e.graph,sources,sinks)

  
  distance.mat <- array(-1,rep(n.pts,2))
  diag(distance.mat) <- Inf
  rownames(distance.mat) <- rownames(e.graph$grand.socio)
  colnames(distance.mat) <- colnames(e.graph$grand.socio)
  
  for (kk in 1:dim(out$source.sink)[1]) {
    distance.mat[edges[kk,1],edges[kk,2]] <-
      out$resistances[kk]/abs(out$fidelities[kk])
    
    #if (e.graph$symmetric)
    distance.mat[edges[kk,2],edges[kk,1]] <-
      out$resistances[kk]/abs(out$fidelities[kk])
  }
    
  
  e.graph$resistances <- out$resistances
  e.graph$fidelities <- out$fidelities  
  e.graph$distance.mat <- distance.mat
  e.graph$conductances <- out$conductances
  #e.graph$voltages <- out$voltages
  #e.graph$currents.node <- out$currents.node
  e.graph$source.sink <- out$source.sink

  e.graph$primers <- out$primers
  
  e.graph$avg.current.black.a <- out$avg.current.black.a
  e.graph$avg.current.red.a <- out$avg.current.red.a
  e.graph$avg.current.black.v <- out$avg.current.black.v
  e.graph$avg.current.red.v <- out$avg.current.red.v
  e.graph$avg.current.black.p <- out$avg.current.black.p
  e.graph$avg.current.red.p <- out$avg.current.red.p
  
  if (e.graph$symmetric) {
    e.graph$avg.current.black.a <- e.graph$avg.current.black.a+t(e.graph$avg.current.black.a)
    e.graph$avg.current.red.a <- e.graph$avg.current.red.a+t(e.graph$avg.current.red.a)
    e.graph$avg.current.black.v <- e.graph$avg.current.black.v+t(e.graph$avg.current.black.v)
    e.graph$avg.current.red.v <- e.graph$avg.current.red.v+t(e.graph$avg.current.red.v)
    e.graph$avg.current.black.p <- e.graph$avg.current.black.p+t(e.graph$avg.current.black.p)
    e.graph$avg.current.red.p <- e.graph$avg.current.red.p+t(e.graph$avg.current.red.p)
  }
  
  return(e.graph)
}

process.for.enemies <- function(sociomatrix, fidelities,
                                voltages, current.matrix) {
  #sociomatrix <- socmat2; voltages <- thing$voltage; current.matrix <- thing$average.currents
  
  sociomatrix <- as.matrix(sociomatrix)
  nn <- dim(sociomatrix)[1]
  res <- .C("process_for_enemies_c",
            as.double(sociomatrix),
            as.double(fidelities),
            as.double(voltages),
            pos.curr=as.double(current.matrix), 
            neg.curr=as.double(0*current.matrix),
            con.str=as.double(0),
            con.fid=as.double(0),
            nn=as.integer(nn))

  return(list(pos.curr=array(res$pos.curr,rep(nn,2)),
              neg.curr=array(res$neg.curr,rep(nn,2)),
              str.fid=c(res$con.str, res$con.fid)))

}
