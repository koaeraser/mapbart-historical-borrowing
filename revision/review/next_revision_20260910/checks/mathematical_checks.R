# Reproducible deterministic checks for the next manuscript revision.
# Uses the local C++ source directly; no simulation/application files are changed.
args <- commandArgs(trailingOnly=TRUE)
root <- normalizePath(args[1]); out <- normalizePath(args[2])
sink(file.path(out, "mathematical_checks.txt")); on.exit(sink())
Rcpp::sourceCpp(file.path(root,"01-code/lrcbart/src/lrcbart.cpp"), cacheDir=file.path(out,"rcpp-cache"))
x <- matrix(c(-1,.5,1,2),ncol=1)
# Nodes encode (variable, cutpoint, left, right, mean, spike).
nodes <- rbind(c(0,0,1,2,0,0),c(-1,0,-1,-1,0,1),c(-1,0,-1,-1,0,0),c(-1,0,-1,-1,0,1))
forest <- list(list(nodes=nodes,roots=as.integer(c(0,3))))
expected <- .25^2*.01 + .75^2*1 + .01
observed <- forest_gvar_cpp(forest,x,.01,1)
stopifnot(abs(observed-expected)<1e-12)
cat("Actual C++ fixed-partition weighted variance:",observed,"expected",expected,"\n")
# A single tree with two independent leaf indicators: exhaustive state sum.
a <- c(.25,.75); w <- .6; V <- .1; t0 <- .01; t1 <- 1
z <- as.matrix(expand.grid(0:1,0:1))
prob <- apply(z,1,function(zz)prod(w^zz*(1-w)^(1-zz)))
vg <- apply(z,1,function(zz)sum(a^2*(zz*t0+(1-zz)*t1)))
exact_proxy <- sum(prob/(V+vg))
one_tree_proxy <- w/(V+sum(a^2)*t0)+(1-w)/(V+sum(a^2)*t1)
stopifnot(abs(exact_proxy-one_tree_proxy)>.1)
cat("Independent leaf-state precision average:",exact_proxy,"one-indicator shortcut:",one_tree_proxy,"\n")
# Conditional Gaussian information vs marginal mixture information.
weights<-c(.5,.5); variances<-c(.01,1)
f<-function(t){d<-vapply(variances,function(v)dnorm(t,0,sqrt(v)),numeric(length(t)));drop(d %*% weights)}
score_num<-function(t){d<-vapply(variances,function(v)-t/v*dnorm(t,0,sqrt(v)),numeric(length(t)));drop(d %*% weights)}
info<-integrate(function(t){p<-f(t); u<-score_num(t); ifelse(p>0,u^2/p,0)},-12,12,subdivisions=1000,rel.tol=1e-10)$value
vr<-1/sum(weights*variances); ci<-sum(weights/variances)
stopifnot(vr<info,info<ci)
cat("Normal scale-mixture: inverse marginal variance",vr,"exact marginal ELIR (sigma^2=1)",info,"average conditional information",ci,"\n")
# Plug-in mean partition variance vs average partitionwise precision.
A <- c(.5,1); plugin<-1/(.1+mean(A)*.2); averaged<-mean(1/(.1+A*.2))
stopifnot(plugin<averaged)
cat("Mean-partition variance plug-in",plugin,"partitionwise average precision",averaged,"\n")
# Shared-leaf constrained quadratic cost: retain spike contributions.
H<-5; t0<-.01; t1<-1; w<-.8; delta<-1; penalty<-log(w/(1-w))+.5*log(t1/t0)
for(m in 0:H){
 v<-c(rep(t1,m),rep(t0,H-m)); alloc<-delta*v/sum(v)
 stopifnot(abs(sum(alloc)-delta)<1e-12)
 calc<-m*penalty+sum(alloc^2/(2*v)); closed<-m*penalty+delta^2/(2*sum(v))
 stopifnot(abs(calc-closed)<1e-12)
 cat("Prior density cost m=",m,":",closed,"\n",sep="")
}
# Pooling loss after optimizing common mean, not n1+n2 scaling.
n1<-50;n2<-200;s1sq<-2.25;s2sq<-3.24; d<-2
A1<-n1/s1sq; A2<-n2/s2sq; common<-A2*d/(A1+A2)
loss<-(A1*common^2+A2*(d-common)^2)/2
stopifnot(abs(loss-d^2/(2*(s1sq/n1+s2sq/n2)))<1e-12)
cat("Optimized Gaussian pooling likelihood loss:",loss,"\n")
# Finite-prior anchor counterexample.
pi<-c(.9,.1); shifted<-c(5/6,1/6); leafvar<-1/(10+100*pi)
v_match<-sum(pi^2*leafvar);v_shift<-sum(shifted^2*leafvar)
stopifnot(v_shift<v_match)
cat("Finite leaf prior variance matched:",v_match,"shifted:",v_shift,"(ceiling can increase)\n")
cat("All deterministic checks passed. No operating characteristics were rerun.\n")
