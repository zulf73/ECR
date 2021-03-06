library(LaplacesDemon)
library(dplyr)
library(plyr)
# Romantic factors
factor1 <- c(-0.2176663 ,0.0000000,0.2357452,
             0.0000000, -0.2175266, 0.0000000,
             -0.2372089, 0.0000000, -0.2411613,
             0.0000000, -0.1949819, 0.0000000,
             -0.2291032,  0.0000000,  0.2213134,  
             0.0000000, -0.2311514, 0.0000000,
             0.2102744, 0.0000000, -0.1902973,
             0.0000000, -0.2435696,  0.0000000,
             0.2391180, 0.0000000, 0.2315982,
             0.0000000, 0.2077698, 0.1310756,
             0.2239868,  0.1013154, 0.2230355,
             0.0000000,  0.2291235,  0.0000000 )
factor2<-c(
  0.0000000, -0.2617577,  0.0000000, -0.2488826,
 -0.1169330, -0.2632161,  0.0000000, -0.2621093,
 -0.1061264, -0.2467630, -0.1542195, -0.1965699,
 -0.1140653, -0.2332962,  0.0000000, -0.1950029,
 -0.1101189, -0.2506772,  0.0000000, -0.2051451,
  0.0000000,  0.2233301,  0.0000000, -0.2160262,
  0.0000000, -0.2127988,  0.0000000, -0.1826593,
  0.0000000, -0.1952212,  0.0000000, -0.1812358,
  0.0000000, -0.1910284,  0.0000000, -0.2034946)

desc<-c(
  'Q1. prefer not to show how I feel',
  'Q2. worry feeling abandoned',
  'Q3. comfortable with lovers',
  'Q4. worry re rels',
  'Q5. pull away when getting close',
  'Q6. worry lack of reciprocation',
  'Q7. uncomfortable too close',
  'Q8. worry losing lover',
  'Q9. uncomfortable opening up',
  'Q10. want reciprocation',
  'Q11. i want cannot get close',
  'Q12. i want to merge scares away lover',
  'Q13. nervous when partner comes too close',
  'Q14. worry about being alone',
  'Q15. comfortable sharing private thoughts feelings',
  'Q16. my desire to be close scares lovers',
  'Q17. avoid getting too close',
  'Q18. need assurance i am loved',
  'Q19. find it easy to get close',
  'Q20. sometimes force my partners to show more feelings',
  'Q21. difficult to depend on partner',
  'Q22. do not worry about abandonment',
  'Q23. prefer distance from partner',
  'Q24. upset or angry when not interest',
  'Q25. tell just about everything',
  'Q26. partners do not want to be as close as i do',
  'Q27. discuss problems and concerns',
  'Q28. anxious and insecure out of rel',
  'Q29.feel comfortable depending on partner',
  'Q30. partner not around as much as i would like',
  'Q31. comfort advice help',
  'Q32. frustrating when partner unavailable',
  'Q33. in need helps to turn to partner',
  'Q34. on disapproval i feel bad about myself',
  'Q35. comfort and reassurance',
  'Q36.resent time spent away'
  )

probs<-function(x){
  z<-hist(x, probability=T)$density
  z<-z[z>0]
  z/sum(z)
}

sq_density<-function(x,y){
  M <- matrix(data=0,ncol=5,nrow=5)
  n <- length(x)
  for ( j in 1:n) {
    M[x[j],y[j]] <- M[x[j],y[j]] + 1
  }
  M/sum(M)
}

udiff<-function( mtx ){
KLD( mtx, matrix(data=1,nrow=5,ncol=5)/25.0 )$mean.sum.KLD
}


matrix_udiff<-function( i, j ){
  mtx <-sq_density( br[,i], br[,j])
  udiff(mtx)
}

if (FALSE) {
kld_unif_mtx <- matrix(data=0, nrow=36, ncol=36)
for (p in 1:36){
  for (q in 1:36){
    kld_unif_mtx[p,q]<-matrix_udiff( p, q )
  }
}
}

###########################################################
# We will use estimated Markov Probability Density
# as the uninformative background rather than uniform 
# this will give us valuable insight

markovProbMatrixData <- c(
  0.09350204, 0.1638430, 0.1094789, 0.1573022, 0.06916148,
  0.27154551, 0.1742725, 0.1377139, 0.2241656, 0.18989142,
  0.10610550, 0.2898345, 0.1658414, 0.3124668, 0.12283513,
  0.10984837, 0.2494882, 0.1987590, 0.3105267, 0.12822478,
  0.10078791, 0.2543676, 0.1732982, 0.3311191, 0.13703650,
  0.15960731, 0.1947235, 0.1546190, 0.2396075, 0.24789529
)

markovProbMatrix <- matrix(data=markovProbMatrixData, nrow=5, ncol=5)
# need to transform this to a joint probability matrix
# P( A | B ) = P( A, B )/P(B)
# markovProbMatrix[ i, j ] =: P( i | j ) = P( i, j )/P(j)
# We better record P(j)

absolute.probs<-rep(0,5)
absolute.probs[1] <- sum(br[,1:36]==1)
absolute.probs[2] <- sum(br[,1:36]==2)
absolute.probs[3] <- sum(br[,1:36]==3)
absolute.probs[4] <- sum(br[,1:36]==4)
absolute.probs[5] <- sum(br[,1:36]==5)
absolute.probs <- absolute.probs/sum(absolute.probs)


markovJointProbs <- markovProbMatrix
for (j in 1:5){
markovJointProbs[,j] < markovJointProbs[,j] * absolute.probs[j]
}

mdiff<-function( mtx ){
  KLD( mtx, markovJointProbs)$mean.sum.KLD
}


matrix_mdiff<-function( i, j ){
  mtx <-sq_density( br[,i], br[,j])
  mdiff(mtx)
}

if (FALSE) {
  kld_markov_mtx <- matrix(data=0, nrow=36, ncol=36)
  for (p in 1:36){
    for (q in 1:36){
      kld_markov_mtx[p,q]<-matrix_mdiff( p, q )
    }
  }
}

  
soft.thresh<-function(x,t){
  y<-x-matrix(data=1,nrow=nrow(x),ncol=ncol(x))*t
  y[y<0]<-0
  y
}


#########################################
# qpairs_desc should include
# x,y, discrepancy, question_x, question_y
# Let's create it

qp <- which(kld_markov_mtx>0.3,arr.ind=T)
n <- dim(qp)[1]

d<-rep(0,n)
qx<-rep(0,n)
qy<-rep(0,n)
qa<-rep(0,n)
qb<-rep(0,n)
for ( jj in 1:n){
  a <- qp[jj,1]
  b <- qp[jj,2]
  qx[jj]<-a
  qy[jj]<-b
  d[jj] <- kld_markov_mtx[a,b]
  qa[jj] <- desc[a] 
  qb[jj] <- desc[b]
}
qpairs_desc<-data.frame(d=d,a=qx,b=qy,qa=qa,qb=qb)

informative_5x5_matrices <-list()
J <- list() # hate typing
for (q in 1:500){
  print(q);
  print(qpairs_desc[q,'qa']);
  print(qpairs_desc[q,'qb']);
  qx <- qpairs_desc[q,'a']
  qy <- qpairs_desc[q,'b']
  mtx <- sq_density(br[,qx],br[,qy])
  J[[q]]<-mtx
}

#################################################
# Now we have all the interesting square matrices
# this is real data without any processing of
# joint distributions

get_joint_distribution <- function( a, b ){
  index <- which( qpairs_desc[,'a']==a & qpairs_desc[,'b']==b)
  print(index)
  
  cmtx <-matrix(0,nrow=2,ncol=2)
  bmtx <-J[[index]]*100
  cmtx[1,1]<-sum(bmtx[1:2,1:2])
  cmtx[1,2]<-sum(bmtx[1:2,4:5])
  cmtx[2,1]<-sum(bmtx[4:5,1:2])
  cmtx[2,2]<-sum(bmtx[4:5,4:5])
  list( metadata = qpairs_desc[index,], mtx = J[[index]],cmtx=cmtx )
}

qpd<-arrange(qpairs_desc,desc(d))
# Calculate the correlations of
# marginals 
n<-min(dim(qpd)[1],500)
cor_mgnl<-rep(0,n)
b11<-rep(0,n)
b12<-rep(0,n)
b21<-rep(0,n)
b22<-rep(0,n)
dom_quad<-rep(0,n)
maxb<-rep(0,n)
bc<-rep(0,n)
quad_conc<-rep(0,n)
for (r in 1:n){
  try({
  a <- qpd[r,2]
  b <- qpd[r,3]
  lX <- get_joint_distribution( a,b)
  v <- cor(colSums(lX$mtx),rowSums(lX$mtx))
  cor_mgnl[r]<- v
  b11[r] <- lX$cmtx[1,1]
  b12[r] <- lX$cmtx[1,2]
  b21[r] <- lX$cmtx[2,1]
  b22[r] <- lX$cmtx[2,2]
  bij <-c(b11[r],b12[r],b21[r],b22[r])
  dom_quad[r]<-which.max(bij)
  maxb[r] <- max(bij)
  bc[r] <- mean(bij[bij != maxb[r] ])
  quad_conc[r]<-(maxb[r]-bc[r])/maxb[r]
  })
}

ecr_pairs <- qpd[1:500,]
ecr_pairs[,'cmgnl']<-cor_mgnl
ecr_pairs[,'b11']<-b11
ecr_pairs[,'b12']<-b12
ecr_pairs[,'b21']<-b21
ecr_pairs[,'b22']<-b22
ecr_pairs[,'maxb']<-maxb
ecr_pairs[,'bc']<-bc
ecr_pairs[,'quad_conc']<-quad_conc

library(ggplot2)
gghistogram<-function( df, var) {
# Create a Histogram
ggplot(data = df, aes(x = var)) + 
  geom_histogram(binwidth = 10, .color = "midnightblue") + 
  theme(legend.position = "top")
}

################################################
# we want functionality where we filter
# by conditions on the answers and then 
# produce probability distributions for the
# remaining questions

ecr_probabilities <- function(dataset){
  cn <- colnames(dataset)
  out_data <- matrix(0,ncol=length(cn),nrow=5)
  for (j in 1:length(cn)){
    y <- plyr::count( dataset, cn[j])
    z <- y[1:5,2]
    out_data[,j]<-z/sum(z)
  }
  ans <- data.frame(out_data)
  colnames(ans)<-cn
  ans
}

ecr_cond_probs<-function( df, conditions ){
  dataset <- df %>% filter(conditions)
  ans <- ecr_probabilities( dataset )
  ans
}

############################################
# How fast do the (A_1,...,A_{36}) converge
# to a delta function, i.e. concentration
# is almost all at one point?
# We need a measure of concentration at a
# point for 5x5 matrix and then we need to
# take the average over all the variables

mtx_concentration <- function( mtx ){
  maxm <-max(mtx)
  dm <- dim(mtx)
  sz <- dm[1]*dm[2]
  tot <- sum(mtx)
  ans <- (2*maxm-tot)/(sz-1)
  ans
}
vec_concentration <- function( mtx ){
  maxm <-max(mtx)
  dm <- dim(mtx)
  sz <- length(mtx)
  tot <- sum(mtx)
  ans <- (maxm-(tot-maxm)/(sz-1))/maxm
  ans
}

###############################
# Monte Carlo simulate
# Answering questions

mean_steps_to_delta<-function( nMonte ){
  all_conc<-matrix(0,ncol=10,nrow=nMonte)
  all_qs <- matrix(0,ncol=10,nrow=nMonte)
  all_as <- matrix(0,ncol=10,nrow=nMonte)
  for (k in 1:nMonte){
    qs <-sample(seq(1,36),size=10,replace=F)
    cond = T;
    as <-sample(seq(1,5),size=10,replace=T)
    aconc <- rep(0,10)
    for (j in 1:10 ){
      
      nloop<-1
      while(T) {
        new_a = sample(seq(1,5))
        cond_test <- cond & (br[,qs[j]]==new_a)
        X <-ecr_cond_probs( br[,1:36], cond_test)
        if (!all(is.na(X))){
          as[j] <- new_a
          cond <- cond_test
          break
        }
        nloop <- nloop + 1
        if (nloop>20){
          break
        }
      }
      X <-ecr_cond_probs( br[,1:36], cond)
      
      nX <- ncol(X)
      tot <- 0
      ct <- 0
      for (r in 1:nX){
        if ( r %in% qs ){
          next
        }
        v <- X[,r]
        if (!is.na(sum(v))){
          ct <- ct + 1
          tot <- tot + vec_concentration(v)
        }
      }
      if ( ct > 0 ){
      avg_conc <- tot/ct
      aconc[j]<-avg_conc
      print(length(aconc))
      all_conc[k,]<-aconc
      all_qs[k,]<-qs
      all_as[k,]<-as
      print(paste('step',ct,' conc=',avg_conc))
      }
    }
  }
  list(c=all_conc,q=all_qs,a=all_as)
}

triple_dist<-function(a,b,c){
  xtbl<-matrix(0, nrow=5,ncol=25)
  dataset<-br[,1:36]
  for (r in 1:nrow(dataset)){
    ia <- dataset[r,a]
    ib <- dataset[r,b]
    ic <- max(dataset[r,c],1)
    #print(paste(ia,ib,ic))
    xtbl[ia,(ic-1)*5 + ib]<-xtbl[ia,(ic-1)*5+ib]+1.
  }
  xtbl<-xtbl*100/sum(xtbl)
  sparsity <- sum(xtbl<0.3)/125.
  list(X=xtbl,sp=sparsity)
}

Jt <- list()

list_triple_dists<-function(K){
  as<-sample(seq(1,36),K,replace=T)
  bs<-sample(seq(1,36),K,replace=T)
  cs<-sample(seq(1,36),K,replace=T)
  sps<-rep(0,K)
  for (k in 1:K){
    TD<-triple_dist(as[k],bs[k],cs[k])
    sps[k]<-TD$sp
    Jt[[k]]<-TD$X
  }
  ans<-data.frame(a=as,b=bs,c=cs,sp=sps)
  ans
}

if (FALSE){
triple_df <- list_triple_dists(2000)
}

enumerate_pairs<-function( v ){
  n<-length(v)
  pairs <- list()
  for (k in 2:n){
    for (j in 1:k){
      pairs <- append(pairs, c(v[j],v[k]))
    }
  }
  pairs
}

get_trp_dist<-function( ia, ib, ic ){
  t<-triple_df
  idx <- which( (t$a==ia) & (t$b==ib) & (t$c==ic))
  if (length(idx)==0){
    return(NULL)
  }
  print(idx)
  Jt[[idx]] 
}

feasible_path<- function( path, new_a ){
  if (length(path)>36){
    return(F)
  }
  valid <- T
  pp <- enumerate_pairs( path )
  for (pair in pp){
    
    fdist<-get_trp_dist(pair[1], pair[2], length(path)+1)
    # value at path[pair[1]], path[pair[2]], length(path)+1
    if (length(fdist)>0){
      ia <- path[pair[1]]
      ib <- path[pair[2]]
      ic <- length(path)+1
      val <- fdist[ia, (ic-1)*5+ib ]
      if (val < thresh){
        valid<-F
        break
      }
    }
  }
  valid
}
ecrp<-ecr_probabilities(br[,1:36])

sim_erc<-function(){
  start<-rep(1,3)
  start[1] <- sample(seq(1,5),1, replace=T, 
                     prob=ecrp[,1])
  start[2] <- sample(seq(1,5),1, replace=T, 
                     prob=ecrp[,2])
  start[3] <- sample(seq(1,5),1, replace=T, 
                     prob=ecrp[,3])
  path <- start
  for (k in 4:36){
    new_a <- sample(seq(1,5),1,prob=ecrp[,k])
    if (feasible_path( path, new_a)){
      path <- append(path,new_a)
      next
    }
  }
  path
}

#############################
# classes determined from examination
# of higher than 0.8 correlation
# of Markov-Density adjusted marginal
# probability correlations

ecr_classes <-c( 1,
                 2,
                 2,
                 2,
                 1,
                 2,
                 1,
                 2,
                 1,
                 2,
                 1,
                 1,
                 1,
                 2,
                 2,
                 1,
                 1,
                 2,
                 2,
                 1,
                 1,
                 1,
                 1,
                 2,
                 1,
                 1,
                 2,
                 1,
                 1,
                 2,
                 2,
                 2,
                 2,
                 2,
                 2,
                 1
)

########################################
# The two factors let us leave at 
# Zulf's Facor 1
# Zulf's factor 2
# Let us consider candidates for description
#
# 1. Anxiety of Loneliness
# 2. Demands of Togetherness
# these are not good names yet but provisional
