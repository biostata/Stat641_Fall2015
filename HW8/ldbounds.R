##  Package ldbounds.R (unreleased version)  ##


"bounds" <-
function(x,t2=x,iuse=1,asf=NULL,alpha=0.05,phi=rep(1,length(alpha)),sides=2,ztrun=rep(8,length(alpha))){
   if (!is.numeric(x)){
      stop("'x' must be a vector of analysis times or the number of analysis times")
   }
   if (length(x)==1){
      t <- (1:x)/x
      if (t2==x){
         t2 <- t
      }
   }
   else{
      t <- x
   }
   if (length(t) != length(t2)){
      stop("Original and second time scales must be vectors of the same length.")
    }
   if ({min(t) < 0.0000001}|{max(t) > 1.0000001}|{min(t2) < 0.0000001}){
      stop("Analysis times must be in (0,1].  Second time scale values must be positive.")
    }
   t3 <- t2
   t2 <- t2/max(t2)
   if ({sum({t-c(0,t[-length(t)]) < 0.0000001}) > 0}|{sum({t2-c(0,t2[-length(t)]) < 0.0000001}) > 0}){
      stop("Analysis times must be ordered from smallest to largest.")
    }
   if ({sum(alpha < 0.0000001) > 0}|{sum(alpha) > 1.0000001}){
      stop("Each component of alpha must be positive and their sum cannot exceed 1.")
    }
   if (length(iuse) != length(alpha)&(length(iuse) != length(asf))){
      stop("For two-sided bounds, the lengths of the iuse and alpha vectors must both be 2.")
    }
   if (length(asf)==2){
      if ({class(asf[[1]])!="function"}|{class(asf[[2]])!="function"}){
         stop("Alpha spending function must be of class 'function'.")
       }
      alpha[1] <- asf[[1]](1)
      alpha[2] <- asf[[2]](1)
    }
   if (length(asf)==1){
      if (class(asf)!="function"){
         stop("Alpha spending function must be of class 'function'.")
       }
      alpha <- asf(1)
    }
   if (sum(iuse==5)<length(asf)){
      stop("Can't specify 2 spending functions unless iuse=c(5,5).")
    }
   if ({sum(iuse==5)>0}&{length(asf)==0}){
      stop("If iuse=5, must specify spending function.")
    }
   if (sum({iuse==3}|{iuse==4}) > length(phi)){
      stop("Phi must be specified for each boundary that uses spending function 3 or 4.")
    }
   if (sum({iuse==3}&{phi <= 0}) > 0){
      stop("For power family (iuse=3), phi must be positive.")
    }
   if (sum({iuse==4}&{phi==0}) > 0){
      stop("For Hwang-Shih-DeCani family (iuse=4), phi cannot be 0.")
    }
   if (length(phi)==1) phi <- rep(phi,2)
   if (!sides%in%c(1,2)){
      stop("Sides must be 1 or 2.")
   }
   if ({length(alpha)==1}|{{length(alpha)==2}&{alpha[1]==alpha[2]}&{iuse[1]==iuse[2]}&{length(asf)!=2}&{ztrun[1]==ztrun[2]}&{{length(phi)==1}|{phi[1]==phi[2]}}}){
      if ({length(alpha)==1}&{sides==1}){
         type <- 1
         alph <- alpha
       }
      if ({length(alpha)==1}&{sides==2}){
         type <- 2
         alph <- alpha
       }
      if (length(alpha)==2){
         type <- 2
         alph <- 2*alpha[1]
       }
      ld <- landem(t,t2,sides,iuse[1],asf,alph,phi[1],ztrun[1])
      ubnd <- ld$upper.bounds
      lbnd <- ld$lower.bounds
      epr <- ld$exit.pr
      dpr <- ld$diff.pr
      spend <- ld$spend
    }
   else{
      type <- 3
      ld1 <- landem(t,t2,1,iuse[1],asf[[1]],alpha[1],phi[1],ztrun[1])
      ld2 <- landem(t,t2,1,iuse[2],asf[[2]],alpha[2],phi[2],ztrun[2])
      lbnd <- -ld1$upper.bounds
      ubnd <- ld2$upper.bounds
      epr <- ld1$exit.pr+ld2$exit.pr
      dpr <- ld1$diff.pr+ld2$diff.pr
      spend <- c(ld1$spend,ld2$spend)
    }
   nom.alpha <- 1-pnorm(ubnd)+pnorm(lbnd)
   ans <- list(bounds.type=type,spending.type=spend,time=t,time2=t3,alpha=alpha,overall.alpha=sum(alpha),lower.bounds=lbnd,upper.bounds=ubnd,exit.pr=epr,diff.pr=dpr,nom.alpha=nom.alpha)
   class(ans) <- "bounds"
   return(ans)
 }

"alphas" <-
function(iuse,asf,alpha,phi,side,t){
   tol <- 10^(-13)
   if (iuse==1){
      pe <- 2*(1-pnorm(qnorm(1-(alpha/side)/2)/sqrt(t)))
      spend <- "O'Brien-Fleming"
    }
   else if (iuse==2){
      pe <- (alpha/side)*log(1+(exp(1)-1)*t)
      spend <- "Pocock"
    }
   else if (iuse==3){
      pe <- (alpha/side)*t^phi
      spend <- "Power Family: alpha * t^phi"
    }
   else if (iuse==4){
      pe <- (alpha/side)*(1-exp(phi*t))/(1-exp(-phi))
      spend <- "Hwang-Shih-DeCani Family"
    }
   else if (iuse==5){
     if(missing(alpha)) alpha <- asf(1)
     if(any(diff(asf(t))<=0.0000001))
       stop("Alpha Spending function must an increasing function.")
     if(asf(1)>1 ) stop("Alpha Spending function must be less than or equal to 1.")
     spend <- "User-specified spending function"
     pe <- (1/side)*asf(t)
    }
   else stop("Must choose 1, 2, 3, 4, or 5 as spending function.")
   pe <- side*pe
   pd <- pe-c(0,pe[-length(pe)])
   if (sum(as.integer({pd<0.0000001*(-1)}|{pd>1.0000001})) >= 1){
      warning("Spending function error")
      pd <- min(1,pd)
      pd <- max(0,pd)
    }
   for (j in 1:length(pd)){
      if (pd[j] < tol){
         warning("Type I error spent too small for analysis #",j,"\n",
                 "Zero used as approximation for ",pd[j])
         pd[j] <- 0
       }
    }
   ans <- list(pe=pe,pd=pd,spend=spend)
   return(ans)
 }

"drift" <-
function(x,za=NULL,zb=NULL,t2=x,pow=NULL,drft=NULL,conf=NULL,pval=NULL,pvaltime=NULL,zval=zb[length(zb)]){
   if (inherits(x, "bounds")){
      za <- x$lower.bounds
      zb <- x$upper.bounds
      t <- x$time
      t2 <- x$time2
   }
   else if (is.numeric(x)){
      t <- x
   }
   else{
      stop("'x' must be a vector of analysis times or inherit from class \"bounds\"")
   }
   if ({{!is.null(za)}&{length(za) != length(zb)}}|{length(t) != length(zb)}|{length(t) != length(t2)}){
      stop("Analysis times, upper bounds, and lower bounds must all be vectors of the same length.")
    }
   if ({min(t) <= 0}|{max(t) > 1}|{min(t2) <= 0}){
      stop("Analysis times must be in (0,1].  Second time scale values must be positive.")
    }
   if (sum({t-c(0,t[-length(t)]) <= 0}|{t2-c(0,t[-length(t2)]) <= 0}) > 0){
      stop("Analysis times must be ordered from smallest to largest.")
    }
   if ({is.null(za)}&{!is.null(zb)})
       za <- -zb
   t3 <- t2
   t2 <- t2/max(t2)
   if (!is.null(pow)+!is.null(drft)+!is.null(conf)+!is.null(pval)!=1){
      stop("Only one of power, drift, confidence level, or p-value ordering can be given.")
    }
   else if (is.null(pow)&is.null(drft)&is.null(conf)&is.null(pval)){
     drft=0
   }
   drift1 <- NULL
   if (!is.null(pow)){
      if ({pow <= 0}|{pow > 1}){
         stop("Power must be in (0,1].")
       }
      type <- 1
      drift1 <- adrift(t2,za,zb,pow)
    }
   if (!is.null(drft)){
      type <- 2
      drift1 <- drft
    }
   if (!is.null(drift1)){
      gl <- glan(t2,za,zb,drift1)
      if (!is.null(drft)) pow <- gl$pr
      ans <- list(type=type,time=t,time2=t3,lower.bounds=za,upper.bounds=zb,power=pow,
                  drift=drift1,lower.probs=gl$qneg,upper.probs=gl$qpos,
                  exit.probs=gl$qneg+gl$qpos,cum.exit=cumsum(gl$qneg+gl$qpos))
    }
   if (!is.null(conf)){
      if (zval < 0){
         stop("Confidence interval is only for nonnegative final Z value.")
       }
      conf.limit <- ci(conf,zval,t2,za,zb)
      ans <- list(type=3,time=t,time2=t3,lower.bounds=za,upper.bounds=zb,
                  conf.level=conf,final.zvalue=zval,conf.interval=conf.limit)
    }
   if (!is.null(pval)){
      if (zval < 0){
         stop("P-value is only for nonnegative Z value.")
       }
      p.value <- adj.p(pval,pvaltime,zval,t2,zb)
      ans <- list(type=4,time=t,time2=t3,lower.bounds=za,upper.bounds=zb,
                  conf.level=conf,analysis.time=pvaltime,final.zvalue=zval,p.ordering=pval,p.value=p.value)
    }
   class(ans) <- "drift"
   return(ans)
 }

"adj.p" <-
function(pval,pvaltime,zval,t,up.bound){
    if (!pval%in%c("SW","LR")){
        stop("Possible p-value orderings are stagewise (SW) and likelihood ratio (LR).")
    }
    if (is.null(pvaltime)){
        stop("P-value time must correspond to one of the analysis times.")
    }
    if (!is.null(pvaltime)){
        if (pvaltime>length(up.bound)){
           stop("P-value time must correspond to one of the analysis times.")
       }
    }    
    if (pval=="SW"){
        p.drift <- drift(zb=c(up.bound[1:(pvaltime-1)],zval),za=rep(-10,3),t=t[1:pvaltime],drft=0)
        p.value <- summary(p.drift)$bounds1[,'Cum exit pr.'][pvaltime]
    }
    else{
        lr.exit <- rep(0,length(up.bound))
        maxval1 <- max(up.bound[1],zval)
        lr1 <- drift(zb=maxval1,za=-10,t=t[1],drft=0)
        lr.exit[1] <- lr1$exit[1]
        for (j in 1:(length(up.bound)-1)){
            maxval <- max(up.bound[j+1],zval)
            lr <- drift(zb=c(up.bound[1:j],maxval),za=rep(-10,j+1),t=t[1:(j+1)],drft=0)
            lr.exit[j+1] <- lr$exit[j+1]
        }
        p.value <- sum(lr.exit)
    }
    return(p.value)
}
"adrift" <-
function(t,za,zb,pow){
   dr <- (zb[length(t)]+qnorm(pow))/sqrt(t[length(t)])
   drft <- bisect(t,za,zb,pow,dr)
   return(drft)
 }

"bisect" <-
function(t,za,zb,target,drft=0,upper=FALSE){
   tol <- 0.000001
   dl <- 0.25
   gotlo <- 0
   gothi <- 0
   prev <- 0
   pr <- 0
   while ({abs(pr-target) > tol}&{abs(drft-prev) > tol/10}){
      glan.out <- glan(t,za,zb,drft)
      if (upper){
         pr <- sum(glan.out$qpos)
       }
      if (!upper){
         pr <- glan.out$pr
       }
      if (pr > target+tol){
         hi <- drft
         drft <- drft-dl
         gothi <- 1
       }
      if (pr < target-tol){
         lo <- drft
         drft <- drft+dl
         gotlo <- 1
       }
      if ({gothi==1}&{gotlo==1}){
         prev <- drft
         drft <- (lo+hi)/2
       }
    }
   if ({abs(drft-prev) <= tol/10}&{abs(pr-target) > tol}){
      warning("Convergence problem")
    }
   return(drft)
 }

"bsearch" <-
function(last,nints,i,pd,stdv,ya,yb){
   tol <- 10^(-7)
   del <- 10
   uppr <- yb[i-1]
   q <- qp(uppr,last,nints[i-1],ya[i-1],yb[i-1],stdv)
   while (abs(q-pd) > tol){
      del <- del/10
      incr <- 2*as.integer(q > pd+tol)-1
      j <- 1
      while (j <= 50){
         uppr <- uppr+incr*del
         q <- qp(uppr,last,nints[i-1],ya[i-1],yb[i-1],stdv)
         if ({abs(q-pd) > tol}&{j==50}){
            stop("Error in search: not converging")
          }
         else if ({{incr==1}&{q <= pd+tol}}|{{incr==-1}&{q >= pd-tol}}){
            j <- 50
          }
         j <- j+1
       }
    }
   ybval <- uppr
   return(ybval)
 }

"ci" <-
function(conf,value,t,za,zb){
   zb[length(t)] <- value
   zcrit <- qnorm(1-(1-conf)/2)
   limit <- (value+c(-1,1)*zcrit)/sqrt(t[length(t)])
   target <- c(0,1)*conf+(1-conf)/2
   lim1 <- bisect(t,za,zb,target[1],limit[1],upper=TRUE)
   lim2 <- bisect(t,za,zb,target[2],limit[2],upper=TRUE)
   lim <- list(lower.limit=lim1,upper.limit=lim2)
   return(lim)
 }

"commonbounds" <-
function(looks,t=(1:looks)/looks,t2=t,iuse="OF",alpha=0.05,sides=2){
   if ({!is.null(looks)}&{!is.numeric(looks)}){
      stop("'looks' must be an integer.")
   }
   if (sum(t==(1:length(t))/length(t))<length(t)){
      warning("Time points are not equally spaced.")
   }
   if (length(t) != length(t2)){
      stop("Original and second time scales must be vectors of the same length.")
    }
   if ({min(t) < 0.0000001}|{max(t) > 1.0000001}|{min(t2) < 0.0000001}){
      stop("Analysis times must be in (0,1].  Second time scale values must be positive.")
    }
   t3 <- t2
   t2 <- t2/max(t2)
   if ({sum({t-c(0,t[-length(t)]) < 0.0000001}) > 0}|{sum({t2-c(0,t2[-length(t)]) < 0.0000001}) > 0}){
      stop("Analysis times must be ordered from smallest to largest.")
    }
   if (sum(!iuse%in%c("PK","OF","HP"))>0){
      stop("Boundary type (iuse) must be \"PK\" or \"OF\".")
   }
   if ({sum(alpha < 0.0000001) > 0}|{sum(alpha) > 1.0000001}){
      stop("Each component of alpha must be positive and their sum cannot exceed 1.")
    }
   if (length(iuse) != length(alpha)){
      stop("For two-sided bounds, the lengths of the iuse and alpha vectors must both be 2.")
    }
   if (!sides%in%c(1,2)){
      stop("Sides must be 1 or 2.")
   }
   if ({length(alpha)==1}|{{length(alpha)==2}&{alpha[1]==alpha[2]}&{iuse[1]==iuse[2]}}){
      if ({length(alpha)==1}&{sides==2}){
          alph <- alpha/2
      }
      else{
          alph <- alpha
      }
      if (iuse[1]=="PK"){
          root <- uniroot(search.glan.pocock,c(1.5,2.3+0.05*looks),k=looks,alpha=alph)$root
          ubnd <- rep(root,looks)
          spend <- "Pocock"
      }
      if (iuse[1]=="OF"){
          root <- uniroot(search.glan.obrien,c(1.5,2+0.05*looks),k=looks,alpha=alph)$root
          ubnd <- root/sqrt((1:looks)/looks)
          spend <- "O'Brien-Fleming"
      }
      if ({length(alpha)==1}&{sides==1}){
         type <- 4
         lbnd <- rep(-8,length(ubnd))
       }
      if ({length(alpha)==2}|{{length(alpha)==1}&{sides==2}}){
         type <- 5
         lbnd <- -1*ubnd
       }
      drift.for.probs <- drift(za=lbnd,zb=ubnd,t=t2,drft=0)
      dpr <- drift.for.probs$upper.probs
      epr <- cumsum(dpr)
    }
   else{
      type <- 6
      spend <- c("","")
      if (iuse[1]=="PK"){
          root <- uniroot(search.glan.pocock,c(1.5,2.3+0.05*looks),k=looks,alpha=alpha[1])$root
          lbnd <- -1*rep(root,looks)
          spend[1] <- "Pocock"
      }
      if (iuse[1]=="OF"){
          root <- uniroot(search.glan.obrien,c(1.5,2+0.05*looks),k=looks,alpha=alpha[1])$root
          lbnd <- -1*root/sqrt((1:looks)/looks)
          spend[1] <- "O'Brien-Fleming"
      }
      if (iuse[2]=="PK"){
          root <- uniroot(search.glan.pocock,c(1.5,2.3+0.05*looks),k=looks,alpha=alpha[2])$root
          ubnd <- rep(root,looks)
          spend[2] <- "Pocock"
      }
      if (iuse[2]=="OF"){
          root <- uniroot(search.glan.obrien,c(1.5,2+0.05*looks),k=looks,alpha=alpha[2])$root
          ubnd <- root/sqrt((1:looks)/looks)
          spend[2] <- "O'Brien-Fleming"
      }
      drift.for.probs <- drift(za=lbnd,zb=ubnd,t=t2,drft=0)
      dpr <- drift.for.probs$upper.probs+drift.for.probs$lower.probs
      epr <- cumsum(dpr)
    }
   nom.alpha <- 1-pnorm(ubnd)+pnorm(lbnd)
   ans <- list(bounds.type=type,spending.type=spend,time=t,time2=t3,alpha=alpha,overall.alpha=sum(alpha),lower.bounds=lbnd,upper.bounds=ubnd,exit.pr=epr,diff.pr=dpr,nom.alpha=nom.alpha)
   class(ans) <- "bounds"
   return(ans)
 }

"cprob" <-
function(last,nints,ya,yb,i,stdv){
   hlast <- (yb[i-1]-ya[i-1])/nints[i-1]
   grid <- seq(ya[i-1],yb[i-1],length=nints[i-1]+1)
   pupr <- (1-pnorm(yb[i],mean=grid,sd=stdv))*last
   plow <- pnorm(ya[i],mean=grid,sd=stdv)*last
   tqpos <- 0.5*hlast*(2*sum(pupr)-pupr[1]-pupr[length(pupr)]) # This is "trap"
   tqneg <- 0.5*hlast*(2*sum(plow)-plow[1]-plow[length(plow)]) # This is "trap"
   ans <- list(qpos=tqpos,qneg=tqneg)
   return(ans)
 }

"fcab" <-
function(last,nints,yam1,h,x,stdv){
  f <- last*dnorm(h*c(0:nints)+yam1,mean=matrix(rep(x,nints+1),nints+1,length(x),byrow=TRUE),sd=stdv)
  area <- 0.5*h*(2*colSums(f)-f[1,]-f[nrow(f),]) # This is "trap"
  return(area)
}

"glan" <-
function(t,za,zb,drft){
   h <- 0.05
   stdv <- sqrt(t-c(0,t[-length(t)])) # These are subroutine "sd"
   sdproc <- sqrt(t)                  # These are subroutine "sd"
   yb <- zb*sdproc-drft*t
   ya <- za*sdproc-drft*t
   nints <- ceiling((yb-ya)/(h*stdv))
   qneg1 <- pnorm(za[1],mean=drft*t[1]/stdv[1])
   qpos1 <- 1-pnorm(zb[1],mean=drft*t[1]/stdv[1])
   cp <- matrix(0,length(t),2)
   cp[1,] <- c(qpos1,qneg1)
   if (length(t) >= 2){
      grid <- seq(ya[1],yb[1],length=nints[1]+1) # These are "first"
      last <- dnorm(grid,mean=0,sd=stdv[1])      # These are "first"
      for (i in 2:length(t)){
         cpr <- cprob(last,nints,ya,yb,i,stdv[i])
         cp[i,] <- c(cpr[[1]],cpr[[2]])
         if (i < length(t)){
            hlast <- (yb[i-1]-ya[i-1])/nints[i-1]                 # These are "other"
            x <- seq(ya[i],yb[i],length=nints[i]+1)               # These are "other"
            last <- fcab(last,nints[i-1],ya[i-1],hlast,x,stdv[i]) # These are "other"
          }
       }
    }
   pr <- sum(cp)
   ans <- list(pr=pr,qpos=cp[,1],qneg=cp[,2])
   return(ans)
 }

"landem" <-
function(t,t2,side,iuse,asf,alpha,phi,ztrun){
   h <- 0.05
   zninf <- -8
   tol <- 0.0000001
   stdv <- sqrt(t2-c(0,t2[-length(t2)])) # These are subroutine "sd"
   sdproc <- sqrt(t2)                    # These are subroutine "sd"
   alph <- alphas(iuse,asf,alpha,phi,side,t)
   za <- zb <- ya <- yb <- nints <- rep(0,length(t))
   pd <- alph$pd
   pe <- alph$pe
   if (pd[1]==0){
      zb[1] <- -zninf
      if (zb[1] > ztrun){
         zb[1] <- ztrun
         pd[1] <- side*(1-pnorm(zb[1]))
         pe[1] <- pd[1]
         if (length(t) > 1) pd[2] <- pe[2]-pe[1]
       }
      yb[1] <- zb[1]*stdv[1]
    }
   else if (pd[1] < 1){
      zb[1] <- qnorm(1-pd[1]/side)
      if (zb[1] > ztrun){
         zb[1] <- ztrun
         pd[1] <- side*(1-pnorm(zb[1]))
         pe[1] <- pd[1]
         if (length(t) > 1) pd[2] <- pe[2]-pe[1]
       }
      yb[1] <- zb[1]*stdv[1]
    }
   if (side==1){
      za[1] <- zninf
      ya[1] <- za[1]*stdv[1]
    }
   else if (side != 1){
      za[1] <- -zb[1]
      ya[1] <- -yb[1]
    }
   nints[1] <- ceiling((yb[1]-ya[1])/(h*stdv[1]))
   if (length(t) >= 2){
      grid <- seq(ya[1],yb[1],length=nints[1]+1) # These are "first"
      last <- dnorm(grid,mean=0,sd=stdv[1])      # These are "first"
      for (i in 2:length(t)){
         if ({pd[i] < 0}|{pd[i] > 1}){
            warning("Possible error in spending function.  May be due to truncation.")
            pd[i] <- min(1,pd[i])
            pd[i] <- max(0,pd[i])
          }
         if (pd[i] < tol){
            zb[i] <- -zninf
            if (zb[i] > ztrun){
               zb[i] <- ztrun
               pd[i] <- side*qp(zb[i]*sdproc[i],last,nints[i-1],ya[i-1],yb[i-1],stdv[i])
               pe[i] <- pd[i]+pe[i-1]
               if (i < length(t)) pd[i+1] <- pe[i+1]-pe[i]
             }
            yb[i] <- zb[i]*sdproc[i]
          }
         else if (pd[i]==1) zb[i] <- yb[i] <- 0
         else if ({pd[i] >= tol}&{pd[i] < 1}){
            yb[i] <- bsearch(last,nints,i,pd[i]/side,stdv[i],ya,yb)
            zb[i] <- yb[i]/sdproc[i]
            if (zb[i] > ztrun){
               zb[i] <- ztrun
               pd[i] <- side*qp(zb[i]*sdproc[i],last,nints[i-1],ya[i-1],yb[i-1],stdv[i])
               pe[i] <- pd[i]+pe[i-1]
               if (i < length(t)){
                  pd[i+1] <- pe[i+1]-pe[i]
                }
             }
            yb[i] <- zb[i]*sdproc[i]
          }
         if (side==1){
            ya[i] <- zninf*sdproc[i]
            za[i] <- zninf
          }
         else if (side==2){
            ya[i] <- -yb[i]
            za[i] <- -zb[i]
          }
         nints[i] <- ceiling((yb[i]-ya[i])/(h*stdv[i]))
         if (i < length(t)){
            hlast <- (yb[i-1]-ya[i-1])/nints[i-1]                 # These are "other"
            x <- seq(ya[i],yb[i],length=nints[i]+1)               # These are "other"
            last <- fcab(last,nints[i-1],ya[i-1],hlast,x,stdv[i]) # These are "other"
          }
       }
    }
   ans <- list(lower.bounds=za,upper.bounds=zb,exit.pr=pe,diff.pr=pd,spend=alph$spend)
   return(ans)
 }

"plot.bounds" <-
function(x, scale = "z", main = NULL, xlab = NULL, ylab = NULL,
         xlim, ylim, las=1, pch=19, type="o",add=F,...){
    if (!((inherits(x, "bounds"))|(inherits(x, "drift"))))
     stop("'x' must inherit from class \"bounds\" or \"drift\"")
    if (!scale%in%c("z","b"))
     stop("Scale must be either \"z\" (z-value) or \"b\" (b-value)")
    if (is.null(main))
       main <- "Sequential boundaries using the Lan-DeMets method"
    if (is.null(xlab))
       xlab <- "Time"
    if (is.null(ylab)){
       if (scale=="z"){
          ylab <- "Z"
       }
       else{
          ylab <- "B"
       }
    }
    z <- c(0,x$time)
    r <- rep(0,length(z))
    if(missing(xlim)) xlim <- c(0,z[length(z)])
    if ({inherits(x, "bounds")}&{x$bounds.type==1}){
       u <- c(NA,x$upper.bounds)
       if (scale=="b"){
          u <- u*sqrt(z)
       }
       if(missing(ylim)) ylim <- c(0,max(u,na.rm=T))
       if(add) lines(z,u, pch=pch, type=type,...)
       else plot(z,u, main = main, xlab = xlab, ylab = ylab, xlim=xlim,
                 ylim=ylim, las=las, pch=pch, type=type,...)
       points(z,r, ...)
       lines(z,r,...)
    }                                 
    else{
       u <- c(NA,x$upper.bounds)
       l <- c(NA,x$lower.bounds)
       if (scale=="b"){
          u <- u*sqrt(z)
          l <- l*sqrt(z)
       }
       if(missing(ylim)) ylim <- c(min(l,na.rm=T),max(u,na.rm=T))
       if(add) lines(z,u, pch=pch, type=type,...)
       else plot(z,u, main = main, xlab = xlab, ylab = ylab, xlim=xlim,
            ylim=ylim, las=las, pch=pch, type=type,...)
       points(z,l,pch=pch, ...)
       lines(z,l,...)
       points(z,r, ...)
       lines(z,r,...)
    }
}


"plot.drift" <-
function(x, scale = "z", main = NULL, xlab = NULL, ylab = NULL,
         xlim, ylim, las=1, pch=19, type="o",add=F, ...){
    if (!((inherits(x, "bounds"))|(inherits(x, "drift"))))
     stop("'x' must inherit from class \"bounds\" or \"drift\"")
    if (!scale%in%c("z","b"))
     stop("Scale must be either \"z\" (z-value) or \"b\" (b-value)")
    if (is.null(main))
       main <- "Sequential boundaries using the Lan-DeMets method"
    if (is.null(xlab))
       xlab <- "Time"
    if (is.null(ylab)){
       if (scale=="z"){
          ylab <- "Z"
       }
       else{
          ylab <- "B"
       }
    }
    z <- c(0,x$time)
    r <- rep(0,length(z))
    if(missing(xlim)) xlim <- c(0,z[length(z)])
    if ({inherits(x, "bounds")}&{x$bounds.type==1}){
       u <- c(NA,x$upper.bounds)
       if (scale=="b"){
          u <- u*sqrt(z)
       }
       if(missing(ylim)) ylim <- c(0,max(u,na.rm=T))
       if(add) lines(z,u, pch=pch, type=type,...)
       else plot(z,u, main = main, xlab = xlab, ylab = ylab, xlim=xlim, ylim=ylim,
            las=las, pch=pch, type=type,...)
       points(z,r, ...)
       lines(z,r,...)
    }                                 
    else{
       u <- c(NA,x$upper.bounds)
       l <- c(NA,x$lower.bounds)
       if (scale=="b"){
          u <- u*sqrt(z)
          l <- l*sqrt(z)
       }
       if(missing(ylim)) ylim <- c(min(l,na.rm=T),max(u,na.rm=T))
       if(add) lines(z,u, pch=pch, type=type,...)
       else plot(z,u, main = main, xlab = xlab, ylab = ylab, xlim=xlim, ylim=ylim, las=las, pch=pch,
            type=type,...)
       points(z,l,pch=19, ...)
       lines(z,l,...)
       points(z,r, ...)
       lines(z,r,...)
    }
}
"print.bounds" <-
function(object, ...)
{
  z <- object
  if (!inherits(z, "bounds"))
    stop("'object' must inherit from class \"bounds\"")
  p <- length(z$time)
  if (identical(z$time,z$time2)){
    b <- matrix(NA, p, 3)
    b[,1:3] <- c(z$time, z$lower.bounds, z$upper.bounds) 
    colnames(b) <- c("Time", "Lower", "Upper")
  }
  else{
    b <- matrix(NA, p, 4)
    b[,1:4] <- c(z$time, z$time2, z$lower.bounds, z$upper.bounds) 
    colnames(b) <- c("Time", "Time 2", "Lower", "Upper")
  }
  ans <- list()
  ans$type <- z$bounds.type
  ans$spending <- z$spending.type
  ans$n <- p
  ans$alpha <- z$alpha    
  ans$oalpha <- z$overall.alpha
  ans$bounds <- b
  rownames(ans$bounds) <- rownames(ans$bounds, do.NULL = FALSE, prefix = "")
  if (ans$type%in%(1:3)){
      cat("\nLan-DeMets bounds for a given spending function \n", "\nn = ", ans$n, "\nOverall alpha: ", ans$oalpha, "\n")
  }
  if (ans$type%in%(4:6)){
      cat("\nGroup sequential boundaries  \n", "\nn = ", ans$n, "\nOverall alpha: ", ans$oalpha, "\n")
  }
  if (ans$type%in%c(1,4)){
    if (ans$type==1){
        cat("\nType: One-Sided Bounds", "\nalpha: ", ans$alpha, "\nSpending function:", ans$spending, "\n", "\nBoundaries:\n")
    }
    if (ans$type==4){
        cat("\nType: One-Sided Bounds", "\nalpha: ", ans$alpha, "\nBoundary type (non-alpha-spending):", ans$spending, "\n", "\nBoundaries:\n")
    }
    if (ncol(ans$bounds)==3) 
      print.default(ans$bounds[,-2], digits = 5, quote = FALSE, print.gap = 2, ...)
    else
      print.default(ans$bounds[,-3], digits = 5, quote = FALSE, print.gap = 2, ...)
  cat("\n")
  }
  else{
    if (ans$type==2){
        if (length(ans$alpha)==2){
            cat("\nType: Two-Sided Symmetric Bounds", "\nLower alpha: ", ans$alpha[1], "\nUpper alpha: ", ans$alpha[2], "\nSpending function: ", ans$spending, "\n")
        }
        else{
            cat("\nType: Two-Sided Symmetric Bounds", "\nLower alpha: ", ans$alpha/2, "\nUpper alpha: ", ans$alpha/2, "\nSpending function: ", ans$spending, "\n")
        }            
    }                                               
    if (ans$type==5){
        if (length(ans$alpha)==2){
          cat("\nType: Two-Sided Symmetric Bounds", "\nLower alpha: ", ans$alpha[1], "\nUpper alpha: ", ans$alpha[1], "\nBoundary type (non-alpha-spending): ", ans$spending, "\n")
        }
        else{
          cat("\nType: Two-Sided Symmetric Bounds", "\nLower alpha: ", ans$alpha/2, "\nUpper alpha: ", ans$alpha/2, "\nBoundary type (non-alpha-spending): ", ans$spending, "\n")
        }
    }                                               
    if (ans$type==3){
      cat("\nType: Two-Sided Asymmetric Bounds", "\nLower alpha: ", ans$alpha[1], "\nSpending function for the lower boundary: ", ans$spending[1], "\nUpper alpha: ", ans$alpha[2], "\nSpending function for the upper boundary: ", ans$spending[2], "\n")
    }
    if (ans$type==6){
      cat("\nType: Two-Sided Asymmetric Bounds", "\nLower alpha: ", ans$alpha[1], "\nType of (non-alpha-spending) lower boundary: ", ans$spending[1], "\nUpper alpha: ", ans$alpha[2], "\nType of (non-alpha-spending) upper boundary: ", ans$spending[2], "\n")
    }
    cat("\nBoundaries:\n")
    print.default(ans$bounds, quote = FALSE, print.gap = 2, ...)
    cat("\n")
  }
}
    
"print.drift" <-
function(x, digit = 5, ...)
{
   z <- x
   if (!inherits(z, "drift")) 
     stop("'x' must inherit from class \"drift\"")
   ans <- list()
   ans$type <- z$type
   ans$n <- length(z$time)
   if ((ans$type==1)|(ans$type==2)){
     ans$power <- z$power
     ans$drift <- z$drift
     if (identical(z$time,z$time2)){           
       b <- matrix(NA, ans$n, 3)
       b[,1:3] <- c(z$time, z$lower.probs, z$upper.probs)
       colnames(b) <- c("Time", "Lower probs", "Upper probs")  
       ans$bounds1 <- b
     }
     else{           
       b <- matrix(NA, ans$n, 4)
       b[,1:4] <- c(z$time, z$time2, z$lower.probs, z$upper.probs)
       colnames(b) <- c("Time", "Time 2", "Lower probs", "Upper probs")
       ans$bounds1 <- b
     }   
   }     
   if (ans$type==3){
     ans$level <- z$conf.level
     ans$fzvalue <- z$final.zvalue
     ans$interval <- z$conf.interval        
   }
   if (ans$type==3){
     ans$level <- z$conf.level
     ans$fzvalue <- z$final.zvalue
     ans$interval <- z$conf.interval        
   }
   if (ans$type==4){
     if (z$p.ordering=="SW"){
         ans$p.ordering <- "Stage-wise"
     }
     if (z$p.ordering=="LR"){
         ans$p.ordering <- "Likelihood ratio "
     }
     ans$fzvalue <- z$final.zvalue
     ans$analysis.time <- z$analysis.time
     ans$p.value <- z$p.value        
   }
   if (identical(z$time,z$time2)){
     ans$bounds <- matrix(c(z$time, z$lower.bounds, z$upper.bounds), ncol=3, dimnames = list(NULL,c("Time", "Lower", "Upper")))
   }   
   else{
     ans$bounds <- matrix(c(z$time, z$time2, z$lower.bounds, z$upper.bounds), ncol=4, dimnames = list(NULL,c("Time", "Time 2", "Lower", "Upper"))) 
   }
   rownames(ans$bounds) <- rownames(ans$bounds, do.NULL = FALSE, prefix = "")
   cat("\nLan-DeMets method for group sequential boundaries \n", "\nn = ", ans$n, "\n")
   cat("\nBoundaries: \n") 
   if ((ans$type==1)|(ans$type==2)){
     rownames(ans$bounds1) <- rownames(ans$bounds1, do.NULL = FALSE, prefix = "")
     print.default(cbind(ans$bounds,ans$bounds1[,-1]), quote = FALSE, print.gap = 2, ...)
     cat("\nPower : ", ans$power, "\n","\nDrift: ", ans$drift, "\n\n")
   }
   if (ans$type==3){
     low <- ans$interval$lower.limit
     up <- ans$interval$upper.limit
     cat("\nConfidence interval at the end of the trial: \n", "\nConfidence level: ", ans$level, "\nLast Z value: ", ans$fzvalue, "\n", 100*ans$level, "% confidence interval: (", low, ",", up, ")\n") 
   }
   if (ans$type==4){
     cat("\nAdjusted p-value: \n", "\nOrdering method: ", ans$p.ordering, "\nLook: ", ans$analysis.time, "\nZ value observed at that time: ", ans$fzvalue, "\n", "P-value: ", ans$p.value, "\n") 
   }
}

"print.summary.bounds" <-
function(x, digit = 5, ...)
{
   z <- x
   if (!inherits(z, "summary.bounds")) 
     stop("'x' must inherit from class \"summary.bounds\"")
   rownames(z$bounds) <- rownames(z$bounds, do.NULL = FALSE, prefix = "")
   if (z$type%in%(1:3)){
     cat("\nLan-DeMets bounds for a given spending function \n", "\nn = ", z$n, "\nOverall alpha: ", z$oalpha, "\n")
   }
   if (z$type%in%(4:6)){
      cat("\nGroup sequential boundaries  \n", "\nn = ", z$n, "\nOverall alpha: ", z$oalpha, "\n")
   }
   if (z$type%in%c(1,4)){
   if (z$type==1){
      cat("\nType: One-Sided Bounds", "\nalpha: ", z$alpha, "\nSpending function:", z$spending, "\n", "\nBoundaries:\n")
   }
   if (z$type==4){
      cat("\nType: One-Sided Bounds", "\nalpha: ", z$alpha, "\nBoundary type (non-alpha-spending):", z$spending, "\n", "\nBoundaries:\n")
    }
   if (ncol(z$bounds)==6) 
      print.default(z$bounds[,-2], digits = 5, quote = FALSE, print.gap = 2)
   else
      print.default(z$bounds[,-3], digits = 5, quote = FALSE, print.gap = 2)
   }
   else{
      if (z$type==2){
        if (length(z$alpha)==2){
          cat("\nType: Two-Sided Symmetric Bounds", "\nLower alpha: ", z$alpha[1], "\nUpper alpha: ", z$alpha[1], "\nSpending function: ", z$spending, "\n")
        }
        else{
          cat("\nType: Two-Sided Symmetric Bounds", "\nLower alpha: ", z$alpha/2, "\nUpper alpha: ", z$alpha/2, "\nSpending function: ", z$spending, "\n")
        }
      }                                               
      if (z$type==5){
        if (length(z$alpha)==2){
          cat("\nType: Two-Sided Symmetric Bounds", "\nLower alpha: ", z$alpha[1], "\nUpper alpha: ", z$alpha[1], "\nBoundary type (non-alpha-spending): ", z$spending, "\n")
        }
        else{
          cat("\nType: Two-Sided Symmetric Bounds", "\nLower alpha: ", z$alpha/2, "\nUpper alpha: ", z$alpha/2, "\nBoundary type (non-alpha-spending): ", z$spending, "\n")
        }
      }                                               
      if (z$type==3){
 cat("\nType: Two-Sided Asymmetric Bounds", "\nLower alpha: ", z$alpha[1], "\nSpending function for the lower boundary: ", z$spending[1], "\nUpper alpha: ", z$alpha[2], "\nSpending function for the upper boundary: ", z$spending[2], "\n")
      }
      if (z$type==6){
        cat("\nType: Two-Sided Asymmetric Bounds", "\nLower alpha: ", z$alpha[1], "\nType of (non-alpha-spending) lower boundary: ", z$spending[1], "\nUpper alpha: ", z$alpha[2], "\nType of (non-alpha-spending) upper boundary: ", z$spending[2], "\n")
      }
      cat("\nBoundaries:\n")
      print.default(z$bounds, digits = digit, quote = FALSE, print.gap = 2)
   }
}

"print.summary.drift" <-
function(x, digit = 5, ...)
{
   z <- x
   if (!inherits(z, "summary.drift")) 
     stop("'x' must inherit from class \"summary.drift\"")
   rownames(z$bounds) <- rownames(z$bounds, do.NULL = FALSE, prefix = "")
   cat("\nLan-DeMets method for group sequential boundaries \n", "\nn = ", z$n, "\n")
   cat("\nBoundaries: \n") 
   print.default(z$bounds, digits = digit, quote = FALSE, print.gap = 2)
   if ((z$type==1)|(z$type==2)){
      cat("\nPower : ", z$power, "\n","\nDrift: ", z$drift, "\n", "\n")
      rownames(z$bounds1) <- rownames(z$bounds1, do.NULL = FALSE, prefix = "")
      print.default(z$bounds1, digits = digit, quote = FALSE, print.gap = 2)
   }
   if (z$type==3){
     low <- z$interval$lower.limit
     up <- z$interval$upper.limit
     cat("\nConfidence interval at the end of the trial: \n", "\nConfidence level: ", z$level, "\nLast Z value: ", z$fzvalue, "\n", 100*z$level, "% confidence interval: (", low, ",", up, ")\n") 
   }
   if (z$type==4){
     cat("\nAdjusted p-value: \n", "\nOrdering method: ", z$p.ordering, "\nLook: ", z$analysis.time, "\nZ value observed at that time: ", z$fzvalue, "\n", "P-value: ", z$p.value, "\n") 
   }
}

"qp" <-
function(xq,last,nints,yam1,ybm1,stdv){
   hlast <- (ybm1-yam1)/nints
   grid <- seq(yam1,ybm1,length=nints+1)
   fun <- last*pnorm(grid,mean=xq,sd=stdv)
   qp <- 0.5*hlast*(2*sum(fun)-fun[1]-fun[length(fun)]) # This is "trap"
   return(qp)
 }

"search.glan.obrien" <-
function(k,c,alpha){
   return(glan((1:k)/k,rep(-8,k),c/sqrt((1:k)/k),0)$pr-alpha)
}
"search.glan.pocock" <-
function(k,c,alpha){
   return(glan((1:k)/k,rep(-8,k),rep(c,k),0)$pr-alpha)
}
"summary.bounds" <-
function (object, digits=5, ...) 
{
    z <- object
    if (!inherits(z, "bounds")) 
     stop("'object' must inherit from class \"bounds\"")
    p <- length(z$time)
    if (identical(z$time,z$time2)){
       b <- matrix(NA, p, 6)
       b[,1:6] <- c(z$time, z$lower.bounds, z$upper.bounds, z$exit.pr, z$diff.pr, z$nom.alpha) 
       colnames(b) <- c("Time", "Lower", "Upper", "Exit pr.", "Diff. pr.", "Nominal Alpha")
    }
    else{
       b <- matrix(NA, p, 6)
       b[,1:7] <- c(z$time, z$time2, z$lower.bounds, z$upper.bounds, z$exit.pr, z$diff.pr, z$nom.alpha) 
       colnames(b) <- c("Time", "Time 2", "Lower", "Upper", "Exit pr.", "Diff. pr.", "Nominal Alpha")
    }
    ans <- list()
    ans$type <- z$bounds.type
    ans$spending <- z$spending.type
    ans$n <- p
    ans$alpha <- z$alpha    
    ans$oalpha <- z$overall.alpha
    ans$bounds <- b
    class(ans) <- "summary.bounds"
    return(ans)
}

"summary.drift" <-
function (object, ...) 
{
    z <- object
    if (!inherits(z, "drift")) 
     stop("'object' must inherit from class \"drift\"")
    ans <- list()
    ans$type <- z$type
    ans$n <- length(z$time)
    if ((ans$type==1)|(ans$type==2)){
        ans$power <- z$power
        ans$drift <- z$drift
        if (identical(z$time,z$time2)){           
           b <- matrix(NA, ans$n, 5)
           b[,1:5] <- c(z$time, z$lower.probs, z$upper.probs, z$exit.probs, z$cum.exit)
           colnames(b) <- c("Time", "Lower probs", "Upper probs", "Exit pr.", "Cum exit pr.")  
           ans$bounds1 <- b
        }
        else{           
           b <- matrix(NA, ans$n, 6)
           b[,1:6] <- c(z$time, z$time2, z$lower.probs, z$upper.probs, z$exit.probs, z$cum.exit)
           colnames(b) <- c("Time", "Time 2", "Lower probs", "Upper probs", "Exit pr.", "Cum exit pr.")
           ans$bounds1 <- b
        }   
    }     
    if (ans$type==3){
 ans$level <- z$conf.level
        ans$fzvalue <- z$final.zvalue
        ans$interval <- z$conf.interval        
    }
   if (ans$type==4){
     if (z$p.ordering=="SW"){
         ans$p.ordering <- "Stage-wise"
     }
     if (z$p.ordering=="LR"){
         ans$p.ordering <- "Likelihood ratio "
     }
     ans$fzvalue <- z$final.zvalue
     ans$analysis.time <- z$analysis.time
     ans$p.value <- z$p.value        
   }
    if (identical(z$time,z$time2)){
        ans$bounds <- matrix(c(z$time, z$lower.bounds, z$upper.bounds), ncol=3, dimnames = list(NULL,c("Time", "Lower", "Upper")))
    }   
    else{
        ans$bounds <- matrix(c(z$time, z$time2, z$lower.bounds, z$upper.bounds), ncol=4, dimnames = list(NULL,c("Time", "Time 2", "Lower", "Upper"))) 
    }
    class(ans) <- "summary.drift"
    return(ans)
}

