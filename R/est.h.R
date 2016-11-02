est.h <-
function(x.sort,x,y,z,h,par,ui,ci){
    
    n.z = ncol(z);
    n.x = length(x.sort);
    alpha = par[,1];
    beta = as.matrix(par[,2:(n.z+1)]);
    gamma = par[,n.z+2];
    phi = par[,n.z+3];
    n = length(x);
    
    # estimate h
    f = function(h){
      h1 = c(0,h,1);
      hhat = H.est(x.sort,h1,x);
      mu = inv.logit(alpha+beta%*%t(z)+outer(gamma,hhat));
      p = phi * mu;
      q = phi * (1-mu);
      f = sum(dbeta(t(y),p,q,log=TRUE));
      return(-f);
    }
    gr = function(h){
      h1 = c(0,h,1);
      hhat = H.est(x.sort,h1,x);
      mu = inv.logit(alpha+beta%*%t(z)+outer(gamma,hhat));
      p = phi * mu;
      q = phi * (1-mu);
      y.tmp = t(log(y/(1-y)));
      tmp = phi * gamma * mu * (1-mu) * (y.tmp - digamma(p) + digamma(q));
      dev = rep(NA,length(x.sort)-2);
      for(i in 1:length(dev)){
        dev[i] = sum(t(tmp)*((x-x.sort[i])/(x.sort[i+1]-x.sort[i])*(x>x.sort[i])*(x<=x.sort[i+1])
                             +(1-(x-x.sort[i+1])/(x.sort[i+2]-x.sort[i+1]))*(x>x.sort[i+1])*(x<=x.sort[i+2])));
      }
      return(-dev);
    }
    theta0=h[-c(1,n.x)];
    res = constrOptim(theta=theta0,f=f,grad=gr,ui=ui,ci=ci,method='BFGS',control=(list(maxit=1e4)));
    return(c(0,res$par,1));
  }
