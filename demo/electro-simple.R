#ElectroGraph package
#Reproduce the figures in the accompanying article.
#Andrew C. Thomas
#May 21, 2009

sources <- c(1,1,2,2,2,3,3,4,4,6,6,7,7,8,8,9,10)
sinks <- c(2,3,4,5,6,4,5,5,7,7,8,8,9,9,10,10,11)

edges <- cbind(sources, sinks)
electro.test <- electrograph.exam(electrograph(edges))


par(mfrow=c(4,2))
plot(electro.test,plot.mode="fruchterman.reingold",main="Network Plot, F-R")
plot(electro.test,plot.mode="fruchterman.reingold",distance.mode="electro.social",
     main="Network Plot, F-R, Electro-Social")

plot(electro.test,main="Network Plot, K-K")
plot(electro.test,distance.mode="electro.social",source.sink.pair=c(1,10),
     max.thick=5,main="Network Plot, K-K, Electro-Social")

plot(electro.test,distance.mode="electro.social",max.thick=3,
     edge.thick="electro.betweenness",
     main="Network Plot, K-K, Electro-Social")


#knock out (2,6).
edges <- edges[-5,]
electro.2 <- electrograph.exam(electrograph(edges))

plot(electro.2, edge.thick="electro.betweenness", max.thick=3)

plot(electro.2, distance.mode="electro.social",
     edge.thick="electro.betweenness", max.thick=3)

plot(electro.2, distance.mode="electro.social",source.sink.pair=c(1,10),
     edge.thick="electro.betweenness", max.thick=3)

