#Produce a time-occupying piece of code to test the speed of the system.
#Andrew C. Thomas
#Last updated: May 25, 2009

n.test <- 10
net.obj <- array(rgamma(n.test^2,3,3),rep(n.test,2))
diag(net.obj) <- 0
system.time(electrograph(net.obj))


sparse.obj <- array(rbinom(n.test^2,1,2/n.test),rep(n.test,2))
diag(sparse.obj) <- 0
system.time(electrograph(sparse.obj))


