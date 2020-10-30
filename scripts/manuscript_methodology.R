pep_tmp<- rep(0,21)
# 8-mers
for(j in 4:11){
  pep_tmp[j:(j+7)]<- pep_tmp[j:(j+7)]+1
}
# 9-mers
for(j in 3:11){
  pep_tmp[j:(j+8)]<- pep_tmp[j:(j+8)]+1
}
# 10-mers
for(j in 2:11){
  pep_tmp[j:(j+9)]<- pep_tmp[j:(j+9)]+1
}
# 11-mers
for(j in 1:11){
  pep_tmp[j:(j+10)]<- pep_tmp[j:(j+10)]+1
}

# plot
svglite::svglite("results/figs/manuscript_methodology.svg")
# par(mfrow=c(2,1))
plot(pep_tmp,type="b",axes=F,xlab=NA,ylab="# peptides containing aa",ylim=c(0,80))
axis(1,at=1:21,1:21)
axis(2)

# Add peptides
n<- 80
points(1:21,y=rep(n,21))

# 8-mers
for(j in 4:11){
  n<- n-1
  points(j:(j+7),y=rep(n,8))
}
# 9-mers
for(j in 3:11){
  n<- n-1
  points(j:(j+8),y=rep(n,9))
}
# 10-mers
for(j in 2:11){
  n<- n-1
  points(j:(j+9),y=rep(n,10))
}
# 11-mers
for(j in 1:11){
  n<- n-1
  points(j:(j+10),y=rep(n,11))
}

dev.off()
