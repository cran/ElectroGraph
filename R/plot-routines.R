#ElectroGraph: plotting routines
#Andrew C. Thomas
#Last updated: February 1, 2011

#refinements: needs sociomatrices replaced with edges.

#already vectorized.

shades <- function(howdark,signs=1*howdark) {
  out <- rep("#000000",length(howdark))
  out[signs>=0] <- rgb(1-howdark[signs>=0],1-howdark[signs>=0],1-howdark[signs>=0])
  out[signs<0] <- rgb(howdark[signs<0],0,0)
  
  return(out)
}

make.edge.colors <- function (relative.thickness,
                              fractional, #from -1 to 1.
                              edge.color.fr.ne.em) {
  #relative.thickness = c(1,1,0.5,0.5); fractional=c(1,0.2,-0.2,-1); edge.color.fr.ne.em=c("black","blue","red")
  #if all were perfect:
  edge.fracs <- col2rgb(edge.color.fr.ne.em)
  if (any(fractional != 1)) {
    rrggbb <- sapply(1:3, FUN=function(kk) {
      (fractional*edge.fracs[kk,1]+(1-fractional)*edge.fracs[kk,2])*
        (fractional>0) +
          (abs(fractional)*edge.fracs[kk,3]+(1-abs(fractional))*edge.fracs[kk,2])*
            (fractional<=0)
    })
    rrggbb <- (255-(255-rrggbb)*rep(relative.thickness,3))/255
    cols <- rgb(rrggbb[,1], rrggbb[,2], rrggbb[,3])
    #cols <- sapply(1:length(fractional), FUN=function(kk) {
    #  c1 <- (fractional[kk]*edge.fracs[,1]+(1-fractional[kk])*edge.fracs[,2])*
    #    (fractional[kk]>0) +
    #      (abs(fractional[kk])*edge.fracs[,3]+(1-abs(fractional[kk]))*edge.fracs[,2])*
    #        (fractional[kk]<=0)
    #  c2 <- (255-(255-c1)*relative.thickness[kk])/255
    #  return(rgb(c2[1],c2[2],c2[3]))
    #})
  } else {
    cols <- rgb(1-relative.thickness, 1-relative.thickness, 1-relative.thickness)
  }
  return(cols) 

}


distance.matrix <- function(coordinates) {
  dist.one <- array(0, rep(dim(coordinates)[1],2))
  for (kk in 1:dim(coordinates)[2]) dist.one <- dist.one + outer(coordinates[,kk], coordinates[,kk], "-")^2
  return(sqrt(dist.one))
}


total.distance <- function(cc1,cc2) sum(sqrt(apply((cc1-cc2)^2,1,sum)))

coordinate.match <- function(cc.move, cc.stat, verbose=1) {

  if (dim(cc.move)[2]!=2) stop ("Moving coordinates aren't two-dimensional.")
  if (dim(cc.stat)[2]!=2) stop ("Stationary coordinates aren't two-dimensional.")
  cc.move.hold <- cc.move
  
  rows <- intersect(rownames(cc.move),rownames(cc.stat))
  if (length(rows)>0) {
    cc.move <- rbind(cc.move[match(rows,rownames(cc.move)),])
    cc.stat <- rbind(cc.stat[match(rows,rownames(cc.stat)),])
    
    cc.stat.cent <- apply(cc.stat,2,mean)
    cc.move.cent <- apply(cc.move,2,mean)
    cc.stat <- t(t(cc.stat)-apply(cc.stat,2,mean))
    cc.move <- t(t(cc.move)-apply(cc.move,2,mean))
    
    angle.stat <- acos(cc.stat[,1]/sqrt(cc.stat[,1]^2+cc.stat[,2]^2))
    angle.move <- acos(cc.move[,1]/sqrt(cc.move[,1]^2+cc.move[,2]^2))
    angles <- (angle.stat-angle.move) %% (2*pi)
    angles <- c(angles,seq(0,pi*2,length=4*dim(cc.move)[1]))
    
    distances <- sapply(1:dim(cc.move)[1],FUN=function(kk) {
      cc2 <- cc.move %*% array(c(cos(angles[kk]),sin(angles[kk]),-sin(angles[kk]),cos(angles[kk])), rep(2,2))
      return(total.distance(cc.stat,cc2))
    })
    distances.flip <- sapply(1:dim(cc.move)[1],FUN=function(kk) {
      cc2 <- cc.move %*% array(c(cos(angles[kk]),-sin(angles[kk]),-sin(angles[kk]),-cos(angles[kk])), rep(2,2))
      return(total.distance(cc.stat,cc2))
    })
    
    pick <- which(distances==min(distances))[1]
    pick2 <- which(distances.flip==min(distances.flip))[1]
    
    if (distances[pick]<distances.flip[pick2]) {
      cc.move.hold <- cc.move.hold %*% array(c(cos(angles[pick]),sin(angles[pick]),
                                               -sin(angles[pick]),cos(angles[pick])), rep(2,2))
    } else {
      cc.move.hold <- cc.move.hold %*% array(c(cos(angles[pick2]),-sin(angles[pick2]),
                                               -sin(angles[pick2]),-cos(angles[pick2])), rep(2,2))
    }
    #cc.move.hold <- t(t(cc.move.hold)+cc.stat.cent)
    cc.move.hold <- t(t(cc.move.hold)-apply(cc.move.hold,2,min))
    
  } else {
    if (verbose) message ("No coordinates to match.")
  }

  return(cc.move.hold)
}



network.projection <- function(connectivity, 
                           force.mode=c("fruchterman.reingold","kamada.kawai"),
                           layout.dimension=2,
                           ego.focus=NULL,
                           initial.coordinates=NULL,
                           verbose=1,
                           projection.iterations=5000) {
  #force.mode="fruchterman.reingold"; layout.dimension=2; ego.focus=NULL; initial.coordinates=NULL; verbose=1; niter=500

  
  force.modes <- c("fruchterman.reingold", "kamada.kawai")
  force.mode <- force.mode[1]; if (!any(force.mode==force.modes))
    stop (paste("Unknown force-directed placement mode:",force.mode))

  nn <- dim(connectivity)[1]
  #if (!any(layout.dimension==1:3)) stop("Layout dimension must be 1, 2 or 3.")

  #symmetrize for calculation.
  connectivity <- (connectivity+t(connectivity))/2
  #fix NAs to zeroes.
  connectivity[is.infinite(connectivity)] <- 0
  connectivity[is.na(connectivity)] <- 0
  
  
  cent.meas <- apply(connectivity,1,sum)
  cent.meas <- cent.meas-0.99*min(cent.meas)
  cent.meas <- cent.meas/max(cent.meas) * 2/nn
  if (verbose) print(head(cent.meas))

  #currently not in use.
  egocenters <- rep(1,nn)
  if (!is.null(ego.focus)) egocenters[ego.focus] <- nn^(1/length(ego.focus))

  
  coords <- array(rnorm(layout.dimension*nn, 0, sqrt(nn)),
                  c(nn,layout.dimension))

  #match row names in initial.coordinates to sociomatrix. For better or worse.
  if (!is.null(initial.coordinates)) {
    if (is.null(rownames(initial.coordinates)))
      rownames(initial.coordinates) <- 1:dim(initial.coordinates)[1]
    
    row.match <- match(rownames(initial.coordinates), rownames(connectivity))
    rows <- rownames(initial.coordinates)[!is.na(row.match)]
    coords[row.match[!is.na(row.match)],] <- initial.coordinates[rows,]

    if (sum(!is.na(row.match))==0)
      stop ("No coordinates were matched for initial locations. Check row names for initial.coordinates.")
    if (sum(!is.na(row.match))< dim(initial.coordinates)[1])
      warning ("Not all coordinates were matched for initial locations.")
  }
  rownames(coords) <- rownames(connectivity)
  force.energy.mode <- 1*(force.mode=="fruchterman.reingold")

  if (verbose>1) {
    #writeLines("Connectivity matrix:"); print(connectivity[1,])
    writeLines("nn, iterations:"); print(c(nn, projection.iterations))
    writeLines("egocenters:"); print(egocenters[1])
    writeLines("dimension:"); print(layout.dimension)
    writeLines("force-energy mode:"); print(force.energy.mode)
    writeLines("coords:"); print(coords[1,])
  }

  pos <- .C("network_projection_by_connectivity",
            connectivity=as.double(connectivity),
            nn=as.integer(nn),
            niter=as.integer(projection.iterations),
            layout.d=as.integer(layout.dimension),
            frucht=as.integer(force.energy.mode),
            egocenters=as.double(egocenters),
            #mincloseness=as.double(0.199),
            coords=as.double(coords),
            verbose=as.integer(1*verbose),
            trialcount=as.integer(0))

  if (verbose) message(paste("Iterations:",projection.iterations - pos$trialcount))
  coords <- array(pos$coords,c(nn,layout.dimension))
  
  coords <- rotate.adapt(coords, verbose=verbose)
  rownames(coords) <- rownames(connectivity)

  
  if (!is.null(initial.coordinates)) if (layout.dimension==2) {
    if (verbose) print ("matching coords")
    coords <- coordinate.match(coords, initial.coordinates)
  }
  
  return(coords)
  
}




#take an n-by-2 coordinate table.
maxrotate <- function(block,tries=16) {
  angles <- seq(0,2*pi,length=(tries+1))[1:tries]
  max.y <- sapply(angles,FUN=function(theta) {
    blocky <- block %*% array(c(cos(theta),-sin(theta),sin(theta),cos(theta)), rep(2,2))
    return(max(blocky[,2])-min(blocky[,2]))
  })
  rot <- angles[which(max.y==max(max.y))][1]
  out <- block %*% array(c(cos(rot),-sin(rot),sin(rot),cos(rot)), rep(2,2))
  out <- t(t(out)-apply(out,2,min))
  return(out)
}


rotate.adapt <- function(coords,first.point=pi, verbose=1) {
  if (verbose>1) writeLines(paste(coords[1:2]))
  if (dim(coords)[2] == 2) {
    coords <- t(t(coords)-apply(coords,2,mean))
    
    #rotate this matrix so that point 1 is on the negative x axis.
    angle <- acos(coords[1,1]/sqrt(coords[1,1]^2+coords[1,2]^2))
    if (coords[1,2]<0) angle <- angle + 2*(pi-angle)
    rotation <- first.point-angle
    coords <- coords%*%array(c(cos(rotation),-sin(rotation),sin(rotation),cos(rotation)),rep(2,2))

    #keep same chirality?
    if (coords[1,2]<0) coords[,2] <- -coords[,2]
  }
  coords <- t(t(coords)-apply(coords,2,min))
  return(coords)
}


block.layout.internal <- function(blockdims, width.height.ratio=1.5) {
  #blockdims <- cbind(c(3,1,2,1,4,1,1),c(2,1,2,2,3,1,1)); width.height.ratio=1.5

  grid.summary <- function(pos, rows) {
    xx <- (pos-1) %% rows + 1
    yy <- trunc((pos-1) / rows) + 1
    return (xx*1.5+yy)
  }

  #x is y, y is x.
  #blockdims <- blockdims[,2:1]
  
  #minimum grid unit?
  min.g <- min(blockdims)/5
  blockdims.original <- blockdims
  blockdims <- ceiling(blockdims/min.g)
  
  #sort by size.
  sizes <- apply(blockdims,1,prod)
  trial <- rev(order(sizes))
  left.bottom <- 0*blockdims
  turns <- rep(0,dim(blockdims)[1])

  out.block <- array(0,c(max(blockdims),
                         max(blockdims)/width.height.ratio))
  ok.out.block <- function(xx,yy, block.x,block.y) {
    if ((xx+block.x-1 > dim(out.block)[1]) | (yy+block.y-1 > dim(out.block)[2])) 
      return(1) else return (sum(out.block[xx-1 + 1:block.x, yy-1 + 1:block.y]))
  }
  
  for (kk in 1:dim(blockdims)[1]) {

    pick <- trial[kk]

    #find the lowest clear point.
    empties <- which(out.block==0)
    if (length(empties)==0) {
      out.block <- rbind(out.block,rep(0,dim(out.block)[2]))
      out.block <- cbind(out.block,rep(0,dim(out.block)[1]))
      empties <- which(out.block==0)
    }
    empties <- empties[order(grid.summary(empties,dim(out.block)[1]))]
    index <- 1
    sx <- (empties[1]-1) %% dim(out.block)[1] + 1
    sy <- trunc((empties[1]-1) / dim(out.block)[1]) + 1
    
    while (ok.out.block(sx,sy,blockdims[pick,1],blockdims[pick,2])>0) {
      index <- index+1
      if (index > length(empties)) {
        if (turns[pick]) {
          turns[pick] <- 0
          blockdims[pick,] <- rev(blockdims[pick,])
          out.block <- rbind(out.block,rep(0,dim(out.block)[2]))
          out.block <- cbind(out.block,rep(0,dim(out.block)[1]))
        } else {
          turns[pick] <- 1
          blockdims[pick,] <- rev(blockdims[pick,])          
        }
        empties <- which(out.block==0)
        empties <- empties[order(grid.summary(empties,dim(out.block)[1]))]
        index <- 1
      }
      sx <- (empties[index]-1) %% dim(out.block)[1] + 1
      sy <- trunc((empties[index]-1) / dim(out.block)[1]) + 1

    }

    out.block[sx-1 + 1:blockdims[pick,1], sy-1 + 1:blockdims[pick,2]] <- pick
    left.bottom[pick,] <- c(sx,sy)
  }

  left.bottom <- (left.bottom-1)*min.g
  #left.bottom[,2] <- max(left.bottom[,2])-left.bottom[,2]  #flip y.
  colnames(out.block) <- rep("C",dim(out.block)[2])
  #out.block <- out.block[,dim(out.block)[2]:1]
  return(list(positions=left.bottom, rotate=turns)) #, master.block=out.block))
}


#first, rotate all coordinates to be longest on the y axis.
#then, place blocks in order of y size.
#finally, give block coordinates.
#Needs an overhaul to maximize floor space.
block.layout <- function(listcoords, bound.factor=1, width.height.factor=1) {
  #listcoords=coords.list
  #estimated, mind you.
  #bound.factor <- 0.5
  new.coords <- listcoords #lapply(listcoords,rotate.adapt) #maxrotate)

  block.size <- sapply(new.coords,FUN=function(block) apply(block,2,max))
  #2 by n.

  orr <- rev(order(block.size[2,]))
  lc.new <- list(NA)
  for (kk in 1:length(orr)) lc.new[[kk]] <- new.coords[[ orr[kk] ]]
  #reorder
  block.size <- cbind(block.size[,orr])
  bound.size <- block.size[1,1]/40*bound.factor+0.001
  for (kk in 1:length(orr)) {
    lc.new[[kk]] <- t(t(lc.new[[kk]]) + rep(bound.size,2))
    block.size[,kk] <- block.size[,kk] + rep(2*bound.size,2)
  }

  holder <- block.layout.internal (t(block.size),width.height.factor)
  origins <- t(holder$positions)
  
  for (kk in 1:length(orr)) if (holder$rotate[kk]) {  ## flip coordinates around.
    lc.new[[kk]] <- lc.new[[kk]][,2:1]
    block.size[,kk] <- block.size[2:1,kk]
  }

  #reflect, so that bigger ones are on top.
  max.y <- max(origins[2,] + block.size[2,])
  origins[2,] <- max.y-origins[2,]-block.size[2,]
    
  old.runs <- function() {
    origins <- array(NA,dim(block.size))
    origins[,1] <- c(0,0)
    max.height <- block.size[2,1]
    new.width <- 0
    current.y <- 0
    current.x <- block.size[1,1]
    
    if (length(orr)>1) for (kk in 2:length(orr)) {
      if (current.y + block.size[2,kk] >= max.height) {
        current.x <- current.x + new.width
        new.width <- 0
        current.y <- 0
      }
      origins[,kk] <- c(current.x,current.y)
      new.width <- max(new.width,block.size[1,kk])
      current.y <- current.y+block.size[2,kk]
    }
  }
  

  out <- list(coords=lc.new,
              block.size=block.size,
              origins=origins)
  return(out)

}

plot.electrograph.core <- function(all.coords, sociomatrix=NULL, fidelity.matrix=NULL,
                                   line.thickness=NULL, new.coord.hold=NULL,
                                   line.colors=NULL,
                                   node.colors=NULL,
                                   label.colors=NULL,

                                   node.cex=1.5, 
                                   label.cex=1.5,
                                   component.border.col=8,
                                   
                                   edges.relative.to.minimum=FALSE,
                                   source.sink.pair=NULL,

                                   max.x=NULL,max.y=NULL,
                                   tick.marks=TRUE,
#                                  edge.colors.specified=NULL,
                                   edge.color.fr.ne.em=NULL,
                                   verbose=1,
                                   set.margins=TRUE,
                                   
                                   ...) {
  
  if (verbose>1) message(paste("Start"))#, paste(proc.time(), collapse="   ")))

  if (is.null(sociomatrix)) sociomatrix <- array(0,rep(dim(all.coords)[1],2))
  #node.colors=NULL; label.colors=NULL; pts.cex=1.5; label.cex=1.5; component.border.col=5;
  arrowheads <- !is.null(source.sink.pair)
  #if (verbose>1) print(paste("line.thickness: ",line.thickness))

  if (is.null(fidelity.matrix)) fidelity.matrix <- array(1,dim(sociomatrix))
  if (is.null(edge.color.fr.ne.em)) edge.color.fr.ne.em <- c("black","blue","red")
  #default value for line.thickness.
  if (is.null(line.thickness)) {
    if (max(sociomatrix)>min(sociomatrix)) line.thickness <- 3*sociomatrix/(max(sociomatrix)-min(sociomatrix)) else line.thickness <- 3*sociomatrix
  }
  
  if (verbose>1) print(paste("line.thickness max: ",max(line.thickness)))
  if (verbose>1) message(paste("Thick"))#, paste(proc.time(), collapse="   ")))

  #redefine max.x, max.y based on blocks
  
  
  if (set.margins) par(mar=c(1, 1, 1, 0.5))
  
  
  #print outline blocks.
  if (!is.null(new.coord.hold)) {
    segs <- array(NA,c(length(new.coord.hold$coords)*4,4))
    for (kk in 1:length(new.coord.hold$coords)) {
      xy <- c(new.coord.hold$origins[,kk],
                  new.coord.hold$origins[,kk]+new.coord.hold$block.size[,kk])
      segs[4*(kk-1)+1:4,] <- xy[c(1,1,1,3,
                                  2,2,4,2,
                                  1,3,3,3,
                                  4,2,4,4)]
    }
    if (is.null(max.x)) max.x <- max(segs[,c(1,3)])
    if (is.null(max.y)) max.y <- max(segs[,c(2,4)])
    
  } else {
    if (is.null(max.x)) max.x <- max(all.coords[,1])*1.05
    if (is.null(max.y)) max.y <- max(all.coords[,2])*1.05
  }
  if (verbose>1) message(paste("Thick1"))#, paste(proc.time(), collapse="   ")))

  plot(c(0,max.x), c(0,max.y), ty="n",axes=FALSE,xlab="",ylab="")#,...)

  if (!is.null(new.coord.hold)) {
    segments(segs[,1], segs[,2], segs[,3], segs[,4], col=component.border.col,lty=2)
  }
  
  #sort edges from dark to light
  n.pts <- dim(all.coords)[1]
  pts.array <- cbind(rep(1:n.pts,n.pts),sort(rep(1:n.pts,n.pts)))

  vals <- abs(as.vector(line.thickness)) #/max(relative.thickness))
  if (verbose>1) message(paste("Thick2"))#, paste(proc.time(), collapse="   ")))

  if (is.null(line.colors)) {
    line.colors <- make.edge.colors(line.thickness/max(line.thickness), fidelity.matrix, edge.color.fr.ne.em)
    #if (is.null(fractional.current)) {
    #  line.colors <- rgb((1-as.vector(fidelity))/2,0,0)
    #} else {
    #  line.colors <- rgb((1-as.vector(fractional.current))/2,0,0)
    #}
  }
  
  #print(rbind(vals,signs))
  
  if (any(vals != 0)) {
    min.v <- min(vals[vals>0])
    max.v <- max(vals[vals>0])
    
    if (edges.relative.to.minimum & min.v<max.v) vals <- (vals-min.v)*(vals>=min.v)/(max.v-min.v)*9/10+1/10*(vals != 0)
  }
  if (verbose>1) message(paste("Thick3"))#, paste(proc.time(), collapse="   ")))

  if (max(abs(vals))>0) { #make sure the resulting graph isn't empty.
    vals <- vals/max(abs(vals))
    pts.array <- pts.array[order(vals),]
    line.colors <- line.colors [order(vals)]
    vals <- vals[order(vals)]
    
    start <- min(which(vals>0),length(vals))
    picks <- start:length(vals)
    
    xx.src <- all.coords[pts.array[picks,1],1]
    xx.snk <- all.coords[pts.array[picks,2],1]
    yy.src <- all.coords[pts.array[picks,1],2]
    yy.snk <- all.coords[pts.array[picks,2],2]

    #note: vals is line.thickness scaled to zero-one.
    segments (xx.src, yy.src, xx.snk, yy.snk,
              col=line.colors[picks], #shades(vals[picks],signs[picks]),
              lwd=max(line.thickness)*vals[picks])

    if (arrowheads) arrows(xx.src, yy.src,
                           (xx.src+xx.snk)/2,
                           (yy.src+yy.snk)/2,
                           length=0.1,
                           lwd=max(line.thickness)*vals[picks],
                           col=line.colors[picks]) #shades(vals[picks],signs[picks]))
    
  }
  if (verbose>1) message(paste("Thick4"))#, paste(proc.time(), collapse="   ")))

  #tick marks.
  if (tick.marks) {
    max.coords <- apply(all.coords,2,max)
    interval <- min(max.coords)/6
    x.seq <- seq(0,max.coords[1],by=interval)
    for (kk in 1:length(x.seq)) lines(rep(x.seq[kk],2),interval*c(-1,1)/20)
    y.seq <- seq(0,max.coords[2],by=interval)
    for (kk in 1:length(y.seq)) lines(interval*c(-1,1)/20,rep(y.seq[kk],2))
  }
  
  if (is.null(node.colors)) node.colors <- rep(4,n.pts)
  if (is.null(label.colors)) label.colors <- rep(1,n.pts)
  
  points (all.coords[,1],all.coords[,2],col=node.colors,pch=19, cex=node.cex)
  text (all.coords[,1],all.coords[,2],rownames(all.coords),col=label.colors,cex=label.cex)

}

plot.electrograph <-
  function(x, connectivity.mode=c("visible.ties", "shortest.path",
                "ohmic", "ohmic.on.ties"),
           force.mode=c("fruchterman.reingold", "kamada.kawai"),
           ego.focus=NULL, manual.coords=NULL,
           redo.coordinates=FALSE,
           just.coordinates=FALSE,
           
           node.colors=NULL, label.colors=NULL, node.cex=1.5, 
           label.cex=1.5, component.border.col=8,

           edge.thickness=c("standard", "geodesic.betweenness", "ohmic.betweenness"),
           max.thick=2,
           source.sink.pair=NULL,
           previous.electrograph.plot.object=NULL,

           tick.marks=FALSE,
           bound.size=1,
           width.height.factor=1,
           edge.colors.specified=NULL,
           edge.color.fr.ne.em=c("black","blue","red"),
           
           verbose=0,

           ...) {

  #connectivity.mode="ohmic"; force.mode="fruchterman.reingold"; ego.focus=NULL; manual.coords=NULL; redo.coordinates=FALSE; just.coordinates=FALSE; max.thick=2; node.colors=NULL; label.colors=NULL; node.cex=1.5; label.cex=1.5; component.border.col=5; bound.size=1;  width.height.factor=1;  tick.marks=FALSE; previous.electrograph.plot.object=NULL; edge.thickness=c("standard", "geodesic.betweenness", "ohmic.betweenness"); source.sink.pair=NULL; verbose=2; edge.colors.specified=NULL; edge.color.fr.ne.em=c("black","blue","red"); x=moxie

  #if (verbose>1) paste(proc.time(), collapse="   ")
  if (class(x)!="electrograph")
    stop("plot.electrograph has (somehow) been given an object of incorrect class.")

  force.modes <- c("fruchterman.reingold", "kamada.kawai")
  connectivity.modes <- c("visible.ties", "shortest.path", "ohmic", "ohmic.on.ties")
  edge.thicknesses <- c("standard", "geodesic.betweenness", "ohmic.betweenness")
  if (!is.null(source.sink.pair)) edge.thickness <- "ohmic.betweenness"
  
  force.mode <- force.mode[1]; if (!any(force.mode==force.modes))
    stop (paste("Unknown force-directed placement mode:", force.mode))
  if (verbose) message(paste("Force mode:",force.mode))
  connectivity.mode <- connectivity.mode[1]; if (!any(connectivity.mode==connectivity.modes))
    stop (paste("Unknown plot distance type. plot.electrograph accepts", connectivity.modes, collapse="; "))
  if (verbose) message(paste("Distance mode:",connectivity.mode))
  edge.thickness <- edge.thickness[1]; if (!any(edge.thickness==edge.thicknesses))
    stop (paste("Unknown edge-thickness mode: ", edge.thickness))
  if (verbose) message(paste("Edge thickness:",edge.thickness))

  if (force.mode == "kamada.kawai" & connectivity.mode == "visible.ties") {
    message ("There is no mode for computing the Kamada-Kawai algorithm based on unprocessed tie information. Setting to shortest.path.")
    connectivity.mode <- "shortest.path"
  }
  if (force.mode == "kamada.kawai" & connectivity.mode == "ohmic.on.ties") {
    message ("There is no mode for computing the Kamada-Kawai algorithm based on tie strengths alone. Setting to ohmic.")
    connectivity.mode <- "ohmic"
  }

  if (!is.null(edge.colors.specified)) {
    if (length(edge.colors.specified) != dim(x$nby3)[1]) stop("In plot.electrograph, edge.colors.specified must have the same number of elements as the number of edges.")   
  }


  grand.stuff <- make.sociomatrix.from.edges (x$nby3, size=length(x$node.ids), fidelities=x$fidel)
  grand.socio <- grand.stuff$socio; grand.fidel <- grand.stuff$fidel
  
  #if (verbose>1) message("plot.e 1")
  #fractional.current <- NULL

  
  if (edge.thickness=="ohmic.betweenness") {
    if (!is.null(source.sink.pair)) {
      if (length(source.sink.pair)!=2) stop("Source-sink pair for edge thickness is not of length 2.") else {
#        if (verbose>1) message("In electro-betweenness thickness routine.")
        
        true.source <- which(source.sink.pair[1]== x$node.ids)
        if (length(true.source)==0) stop ("Specified source node does not appear in the system.")
        true.sink <- which(source.sink.pair[2]== x$node.ids)
        if (length(true.sink)==0) stop ("Specified sink node does not appear in the system.")

        #rerun analysis to make sure it will work.
        if (x$symmetric) {test <- compute.volts.symmetric.edges (x$nby3, true.source, true.sink, x$fidelities)
        } else {test <- compute.volts.asymmetric.edges (x$nby3, true.source, true.sink, x$fidelities)}
        
        volt.hold <- as.vector(test$voltages)
        if (verbose>1) print(volt.hold)
        
        untouched.current <- (array(volt.hold, dim(grand.socio))-
                              t(array(volt.hold,dim(grand.socio))))*grand.socio
        holds <- process.for.enemies (grand.socio, grand.fidel, volt.hold, untouched.current)
        relative.thickness <- holds$pos.curr+holds$neg.curr
        
        frenem <- (holds$pos.curr-holds$neg.curr)/(holds$pos.curr+holds$neg.curr)
        frenem[is.na(frenem)] <- 0

        if (is.null(edge.colors.specified))
          edge.colors.specified <-
            make.edge.colors(relative.thickness/max(relative.thickness), frenem, edge.color.fr.ne.em)  
      }
    } else {
      if (is.null(x$avg.current.pro.p)) x <- ohmic.properties(x)
      relative.thickness <- x$avg.current.pro.p+x$avg.current.con.p+t(x$avg.current.pro.p)+t(x$avg.current.con.p)
      frenem <- (x$avg.current.pro.p+t(x$avg.current.pro.p) - x$avg.current.con.p-t(x$avg.current.con.p))/(x$avg.current.pro.p+t(x$avg.current.pro.p) + x$avg.current.con.p+t(x$avg.current.con.p))
      frenem[is.na(frenem)] <- 0
      
      if (is.null(edge.colors.specified))
        edge.colors.specified <-
          make.edge.colors(relative.thickness, frenem, edge.color.fr.ne.em)
    }
  } else relative.thickness <- grand.socio+t(grand.socio)





  
  #if (verbose) message("plot.e 2")
  
  #first, get coordinates for each connected subset.
  #number of subsets?
  #if (verbose>1) message("Getting coordinates")
  subsets <- length(unique(x$component.vector))
  coords.list <- list(NA)

  if (!is.null(manual.coords)) {
    if (dim(manual.coords)[2]!=2) stop("Coordinates should have two dimensions.")
    if (dim(manual.coords)[1]!=dim(grand.socio)[1])
      stop("The number of coordinates does not match the number of nodes.")
    coords.list[[1]] <- manual.coords
  } else {
  
    if (connectivity.mode == "ohmic" | connectivity.mode == "ohmic.on.ties") {
      if (is.null(x$ohmic.distance.mat) | redo.coordinates) {
        message ("Obtaining Ohmic distances.")
        x <- ohmic.properties(x)
      }
      grand.connectivity <- 1/x$ohmic.distance.mat
      if (connectivity.mode == "ohmic.on.ties")  grand.connectivity <- grand.connectivity*(grand.socio>0)
    } else if (connectivity.mode == "shortest.path") {
      if (is.null(x$geodesic)) {
        message ("Obtaining geodesic distances.")
        x$geodesic <- geodesic.mat.socio(grand.socio)
      }
      grand.connectivity <- 1/x$geodesic
    } else if (connectivity.mode == "visible.ties") {
      grand.connectivity <- grand.socio
    }

    
    
    if (verbose>1) message(paste("grand.connectivity dims: "))#, dim(grand.connectivity)[1], dim(grand.connectivity)[2], paste(proc.time(), collapse="   ")))

    
    for (kk in 1:subsets) if (sum(x$component.vector==kk)>1) {
      if (verbose>1) message(paste("subset",kk, paste(proc.time(), collapse="   ")))

      little.connectivity <- grand.connectivity[which(x$component.vector==kk), which(x$component.vector==kk)]
      if (!is.null(previous.electrograph.plot.object)) {
        last.coords <- previous.electrograph.plot.object$coordinates
      } else last.coords <- NULL

      #if (verbose>1) print(little.connectivity)
      
      #writeLines(paste("Verbose:", print(verbose)))
      t.coords <- network.projection(little.connectivity,
                                 force.mode,
                                 ego.focus=ego.focus,
                                 initial.coordinates=last.coords,
                                 verbose=verbose)
      
      coords.list[[kk]] <- t.coords
      rownames(coords.list[[kk]]) <- x$node.ids[which(x$component==kk)]
      
    } else {
      coords.list[[kk]] <- rbind(c(0,0))
      rownames(coords.list[[kk]]) <- x$node.ids[which(x$component==kk)]
    }

  }
  if (verbose>1) message(paste("Finished layouts."))#, paste(proc.time(), collapse="   ")))


  
  #now, fix coordinates. new-coord-hold.
  #if (verbose>1) message("Getting block layout")
  nc.hld <- block.layout(coords.list, bound.size, width.height.factor)
  sociomatrix <- grand.socio
  
  all.coords <- NULL
  names <- NULL
  for (kk in 1:length(nc.hld$coords)) {
    all.coords <- cbind(all.coords,t(nc.hld$coords[[kk]])+nc.hld$origins[,kk])
    names <- c(names,rownames(nc.hld$coords[[kk]]))
  }
  all.coords <- t(all.coords) #re-sort
  all.coords <- all.coords[match(x$node.ids,rownames(all.coords)),]

  out <- list(coordinates=all.coords,
              new.coord.hold=nc.hld,

              node.colors=node.colors,
              label.colors=label.colors,
              node.cex=node.cex, 
              label.cex=label.cex,
              component.border.col=component.border.col,
              source.sink.pair=source.sink.pair,
              
              line.thickness=max.thick*relative.thickness/max(relative.thickness),
              tick.marks=tick.marks,
              verbose=verbose,
              sociomatrix.hold=sociomatrix,
              grand.fidelity=grand.fidel,

              edge.colors.specified=edge.colors.specified,
              edge.color.fr.ne.em=edge.color.fr.ne.em)

  
  
  class(out) <- "electrograph.plot"            

  if (verbose) message("Finished plot.electrograph preparatory steps.")
  #if (verbose>1) message(paste(proc.time(), collapse="   "))

  if (!just.coordinates)
    plot(out, sociomatrix=sociomatrix, verbose=verbose, ...)

  if (verbose>1) message(paste("Finished"))#, paste(proc.time(), collapse="   ")))
  
  return(invisible(out))
  
}



#the natural move here: make "plot.electrograph.plot", which just takes an electrograph.plot object and sends it to electrograph.core.
#build in overrides.
plot.electrograph.plot <- function(x, sociomatrix=NULL, edgelist=NULL, verbose=1,
                                   node.colors=x$node.colors,
                                   label.colors=x$label.colors,
                                   node.cex=x$node.cex,
                                   label.cex=x$label.cex,
                                   component.border.col=x$component.border.col,
                                   source.sink.pair=x$source.sink.pair,
                                   line.thickness=x$line.thickness,
                                   tick.marks=x$tick.marks,
                                   edge.colors.specified=x$edge.colors.specified,
                                   edge.color.fr.ne.em=x$edge.color.fr.ne.em,
                                   ...) {

  #Add options.
  #if (class(x)!="electrograph.plot") stop ("Object should be of class electrograph.plot.")
  if (is.null(sociomatrix)) {
    sociomatrix <- x$sociomatrix.hold
  }
  if (!is.null(edgelist)) {
    if (verbose) message ("electrograph.plot: Making sociomatrix.")
    socio.hold <- make.sociomatrix.from.edges(edgelist)
    sociomatrix <- socio.hold$sociomatrix
  }
  
  plot.electrograph.core (x$coordinates, sociomatrix=sociomatrix,
                          fidelity.matrix=x$grand.fidelity,
                          new.coord.hold=x$new.coord.hold,
                          node.colors=node.colors, label.colors=label.colors,
                          node.cex=node.cex, label.cex=label.cex,
                          component.border.col=component.border.col,
                          source.sink.pair=source.sink.pair,
                          line.thickness=line.thickness,
                          tick.marks=tick.marks,
                          line.colors=edge.colors.specified,
                          edge.color.fr.ne.em=edge.color.fr.ne.em,
                          verbose=verbose,
                          ...)

}

 



wedding.cake.plot <- function(x, connectivity.mode="visible.ties",
                              force.mode="fruchterman.reingold",
                              ego.focus=NULL,
                              filebase="electrograph-cake", #type="png",
                              lower.bound=NULL, upper.bound=NULL, main.title=TRUE,
                              plot.width=600, plot.height=600,verbose=1, 
                              ...) {
  
  #x=graphy; connectivity.mode="shortest.path"; force.mode="kamada.kawai"; ego.focus=NULL; filebase="electrograph-cake"; lower.bound=NULL; upper.bound=NULL; verbose=1
  
  if (class(x)!="electrograph") stop ("plot.wedding.cake requires an object of class electrograph.")

  nonzero.edges <- x$nby3[,3]  #x$grand.sociomatrix[x$grand.sociomatrix>0]

  #density <- mean(x$grand.sociomatrix>0)
  if (is.null(lower.bound) & is.null(upper.bound))
    lower.bound <- unique(quantile(nonzero.edges,seq(0,1-1/4/length(x$node.ids),length=50)))

  if (!is.null(lower.bound) & !is.null(upper.bound))
    if (length(lower.bound) != length(upper.bound))
      stop ("Upper and lower bounds for edge weights should be of the same length.")

  if (is.null(lower.bound)) lower.bound <- rep(0,length(upper.bound))
  if (is.null(upper.bound)) upper.bound <- rep(Inf,length(lower.bound))
  
  #get coordinates.
  coords <- plot.electrograph(x, connectivity.mode, force.mode, ego.focus, just.coordinates=TRUE)$coord
  #plot(x, connectivity.mode, force.mode, ego.focus, just.coordinates=TRUE)$coord
  #make directory.
  direc <- filebase
  ii <- -1
  works <- suppressWarnings(dir.create(direc))

  ii=0; while (!works) {
    direc <- paste(filebase,ii,sep="-")
    works <- suppressWarnings(dir.create(direc))
    ii <- ii+1
  }
  if (ii >= 0) message (paste("Output directory: ",direc))

  sociomatrix <- make.sociomatrix.from.edges(x$nby3)$sociomatrix

  filename <- rep(NA,length(lower.bound))
  for (kk in 1:length(lower.bound)) {
    socio.hold <- (sociomatrix+t(sociomatrix))/2
    socio.hold <- socio.hold*(socio.hold>lower.bound[kk])*(socio.hold<upper.bound[kk])
    isolates <- which(apply(socio.hold,1,sum)==0)
    
    node.cex <- rep(1.5, dim(socio.hold)[1]); node.cex[isolates] <- 0.5
    node.colors <- rep(5, dim(socio.hold)[1]); node.colors[isolates] <- 8
    label.cex <- rep(1.5, dim(socio.hold)[1]); label.cex[isolates] <- 1
    label.colors <- rep(2, dim(socio.hold)[1]); label.colors[isolates] <- 1

    filename[kk] <- paste(direc,"/",filebase,"-",kk,".png",sep="")
    
    png(filename[kk], width=plot.width, height=plot.height)
    mainy <- ""; if (main.title) {
      if (lower.bound[kk]>0) mainy <- paste(mainy, "Lower bound: ",signif(lower.bound[kk],3))
      if (is.finite(upper.bound[kk]))
        mainy <- paste(mainy, "Upper bound: ",signif(upper.bound[kk],3))
    }
    
    plot.electrograph.core (coords, socio.hold,
                            node.colors=node.colors, label.colors=label.colors,
                            node.cex=node.cex, label.cex=label.cex,
                            component.border.col=0,

                            edges.relative.to.minimum=TRUE,
                            
                            main=mainy)
    
    dev.off()
  }

  #output the "build animation" thing too.

  out.instruct <- paste(direc,"-movie.gif",sep="")
  f.hold <- "";
  for (kk in (length(lower.bound)-1):2)
    f.hold <- paste(filename[kk],f.hold)
  out.instruct <- paste("\\( -delay 10",f.hold,"\\)",out.instruct)
  out.instruct <- paste("\\( -delay 100 ",filename[1],"\\)",out.instruct)
  f.hold <- "";
  for (kk in (2:length(lower.bound)-1))
    f.hold <- paste(filename[kk],f.hold)
  out.instruct <- paste("\\( -delay 10",f.hold,"\\)",out.instruct)
  out.instruct <- paste("\\( -delay 100 ",filename[length(lower.bound)],"\\)",out.instruct)
  
  out.instruct <- paste("convert ",out.instruct)

  sink(paste(direc,"-compile.sh",sep=""))
  writeLines (out.instruct)
  sink()
  message(paste("Movie-making information in file",direc,"-compile.sh.",sep=""))
  
  return (invisible(filename))
  
}

dist.calc <- function(init.pos) { #, fin.pos) {
  init.pos <- cbind(init.pos) #;fin.pos <- cbind(fin.pos)
  comps <- array(NA,c(dim(init.pos)[1],dim(init.pos)))
  ds <- rep(dim(init.pos)[1],2)
  for (kk in 1:dim(init.pos)[2]) {
    comps[,,kk] <- array(init.pos[,kk],ds)-t(array(init.pos[,kk],ds))
  }
  return(sqrt(apply(comps^2,c(1,2),sum)))
}

dist.to.pair <- function(others, pair) {
  if (length(pair)!=dim(others)[2]) stop ("dist.to.pair")
  out <- rep(0,length(pair))
  for (kk in 1:length(pair)) out <- out+(others[,kk]-pair[kk])^2
  return(sqrt(out))
}
#from coordinates, n by 2, to distances, n by n.
dist.coords.old <- function(coords) {
  dists <- array(0,rep(dim(coords)[1],2))
  for (kk in 1:dim(coords)[2]) {
    dists <- dists+(array(coords[,kk],dim(dists))-t(array(coords[,kk],dim(dists))))^2
  }
  dists <- sqrt(dists)
  return(dists)
}
dist.coords <- function(coords, coords2=NULL) {
  if (is.null(coords2)) coords2 <- coords
  coords <- as.matrix(coords)
  coords2 <- as.matrix(coords2)
  t1 <- array(.C("distance_from_coordinates",
                 pos1=as.double(coords),
                 pos2=as.double(coords2),
                 output=as.double(rep(0, dim(coords)[1]*dim(coords2)[1])),
                 dp1=as.integer(dim(coords)[1]),
                 dp2=as.integer(dim(coords2)[1]),
                 dimspace=as.integer(dim(coords)[2]))$output, c(dim(coords)[1],dim(coords2)[1]))
  return(t1)
  
}
  


#Given the plot information for a series of graphs in time, prepare an animation of transitions.
animate.plot.series <- function(plots.to.animate, intermediates=10, filebase="electro-animate",
                                plot.width=600, plot.height=600, verbose=1, ...) {
  #plots.to.animate <- newcomb.plot; intermediates<-10; filebase="electro.animate"; plot.width=600; plot.height=600

  
  #Checking input.
  if (class(plots.to.animate)!="list") stop ("Object plots.to.animate must be of class ``list''.")
  if (length(plots.to.animate)==1) stop ("Object plots.to.animate must be greater than length 1.")
  
  for (kk in 1:length(plots.to.animate)) if (class(plots.to.animate[[kk]])!="electrograph.plot")
    stop (paste("Object",kk,"is not an output from plot.electrograph."))

  for (kk in 1:(length(plots.to.animate)-1)) {
    if (dim(plots.to.animate[[kk]]$coordinates)[1] != dim(plots.to.animate[[kk+1]]$coordinates)[1])
      stop (paste("Plots",kk,"and",kk+1,"do not have the same number of nodes."))

    if (dim(plots.to.animate[[kk]]$coordinates)[2] != dim(plots.to.animate[[kk+1]]$coordinates)[2])
      stop (paste("Plots",kk,"and",kk+1,"do not have the same number of plotting dimensions"))
  }

  weights <- 1-pbeta(seq(0,1,length=intermediates+2),3,3)
  filename <- rep(NA,(length(weights)-1)*(length(plots.to.animate)-1)+1)

  #make directory.
  direc <- filebase
  ii <- -1
  works <- suppressWarnings(dir.create(direc))
  if (!works) ii=0; while (!works) {
    direc <- paste(filebase,ii,sep="-")
    works <- suppressWarnings(dir.create(direc))
    ii <- ii+1
  }
  if (ii >= 0) message (paste("Output directory: ",direc))

  
  #determine coordinate space.
  max.x <- 0; max.y <- 0
  for (kk in 1:length(plots.to.animate)) {
    max.x <- max(max.x,plots.to.animate[[kk]]$coordinates[,1])
    max.y <- max(max.y,plots.to.animate[[kk]]$coordinates[,2])
  }
  
  #All working?
  count <- 0
  for (kk in 1:(length(plots.to.animate)-1)) {
    for (ww in 1:(length(weights)-1)) {
      count <- count+1
      filename[count] <- paste(direc,"/",filebase,"-",kk,"-",ww,".png",sep="")
      png(filename[count], width=plot.width, height=plot.height)

      coords.hold <- weights[ww]*plots.to.animate[[kk]]$coordinates +
                             (1-weights[ww])*plots.to.animate[[kk+1]]$coordinates
      #center image in plot.
      x.range <- range(coords.hold[,1]); 
      coords.hold[,1] <- coords.hold[,1] + (max.x-x.range[2]-x.range[1])/2
      y.range <- range(coords.hold[,2]); 
      coords.hold[,2] <- coords.hold[,2] + (max.y-y.range[2]-y.range[1])/2

      
      #coords.hold <- t(t(coords.hold)-apply(coords.hold,2,mean)+c(max.x,max.y)/2)
      #readjust.
      #coords.hold[,1] <- coords.hold[,1]-min(0,coords.hold[,1])
      #coords.hold[,2] <- coords.hold[,2]-min(0,coords.hold[,2])
      
      
      plot.electrograph.core(coords.hold,

#                             sociomatrix=NULL,
#                             new.coord.hold=NULL,

                             node.cex=weights[ww]*plots.to.animate[[kk]]$node.cex +
                               (1-weights[ww])*plots.to.animate[[kk+1]]$node.cex,
                             label.cex=weights[ww]*plots.to.animate[[kk]]$label.cex +
                               (1-weights[ww])*plots.to.animate[[kk+1]]$label.cex,
                             line.thickness=weights[ww]*plots.to.animate[[kk]]$line.thickness +
                               (1-weights[ww])*plots.to.animate[[kk+1]]$line.thickness,
                             
                             node.colors=plots.to.animate[[kk]]$node.colors,
                             label.colors=plots.to.animate[[kk]]$label.colors,
                             source.sink.pair=plots.to.animate[[kk]]$source.sink.pair,
                             component.border.col=plots.to.animate[[kk]]$component.border.col,

                             max.x=max.x, max.y=max.y,
                             main=paste("Time Point:",kk),...)
      dev.off()            
    }
    if (verbose>1) message(paste("Through",kk))
  }
  count <- count+1
  filename[count] <- paste(direc,"/",filebase,"-",length(plots.to.animate),"-",1,".png",sep="")

  coords.hold <- plots.to.animate[[length(plots.to.animate)]]$coordinates
  x.range <- range(coords.hold[,1]); 
  coords.hold[,1] <- coords.hold[,1] + (max.x-x.range[2]-x.range[1])/2
  y.range <- range(coords.hold[,2]); 
  coords.hold[,2] <- coords.hold[,2] + (max.y-y.range[2]-y.range[1])/2

  plots.to.animate[[length(plots.to.animate)]]$coordinates <- coords.hold
  
  png(filename[count], width=plot.width, height=plot.height)
  plot(plots.to.animate[[length(plots.to.animate)]], max.x=max.x, max.y=max.y, main="Time Point: End")
  dev.off()

  #prepare the file for movie making.
  file.string <- ""
  count <- 0
  for (kk in 1:(length(plots.to.animate)-1)) {
    count <- count+1
    file.string <- paste(file.string,"\\( -delay 50 ",filename[count], "\\)  \\( -delay 10") #,"-delay >100")
    for (ww in 2:(length(weights)-1)) {
      count <- count+1
      file.string <- paste(file.string,filename[count])
    }
    file.string <- paste(file.string," \\)")
  }
  count <- count+1
  file.string <- paste(file.string,"\\( -delay 200 ",filename[count], "\\)") #,"-delay >100")

    
  out.instruct <- paste("convert ",file.string," ",direc,"-movie.gif",sep="")
  
  sink(paste(direc,"-compile.sh",sep=""))
  writeLines (out.instruct)
  sink()
  message(paste("Movie-making information in file ",direc,"-compile.sh.",sep=""))

  return(invisible(filename))
  
}



#energy.force.test <- function (coords, true.dist,
#                               ego.focus=NULL, choice) {
#  efk <- rep(0,dim(coords)[2]+1)
#  if (is.null(ego.focus)) ego.focus <- rep(1,dim(coords)[1])
#  res <- .C("energy_force_kk_c",
#            as.double(coords), as.double(true.dist),
#            as.double(ego.focus), as.integer(choice-1), efk=as.double(efk),
#            dimension=as.integer(dim(coords)[2]), nn=as.integer(dim(coords)[1]))
#  return(res$efk)
#}
