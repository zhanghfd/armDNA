est.par <-
function(y,z,h.est,par,cl){
    
    n.z = ncol(z);
    estimation = function(yj,z,h.est)
    {
    # estimate alpha, beta, gamma, phi
      fn = function(par){
          alpha = par[1];
          beta = par[2:(n.z+1)];
          gamma = par[n.z+2];
          phi = par[n.z+3];
          mu = exp(alpha+as.vector(z%*%beta)+gamma*h.est)/(1+exp(alpha+as.vector(z%*%beta)+gamma*h.est));
          p = phi * mu;
          q = phi * (1-mu);
          res = - sum(dbeta(yj,p,q, log = TRUE));
          return(res);
      }
      gr = function(par){
          alpha = par[1];
          beta = par[2:(n.z+1)];
          gamma = par[n.z+2];
          phi = par[n.z+3];
          mu = exp(alpha+as.vector(z%*%beta)+gamma*h.est)/(1+exp(alpha+as.vector(z%*%beta)+gamma*h.est));
          p = phi * mu;
          q = phi * (1-mu);
          digamma.p = digamma(p);
          digamma.q = digamma(q);
          digamma.phi = digamma(phi);
          y.tmp = log(yj/(1-yj));
          tmp = phi * mu * (1-mu) * (y.tmp - digamma.p + digamma.q);
          di.alpha = sum(tmp);
          di.beta = as.vector(tmp%*%z);
          di.gamma = sum(tmp*h.est);
          di.phi = sum(mu*(y.tmp-digamma.p+digamma.q)+log(1-yj)+digamma.phi-digamma.q);
          return(-c(di.alpha,di.beta,di.gamma,di.phi));
      }
      res = optim(par=c(rep(0,n.z+2),1),fn=fn,gr=gr,method='BFGS',control=(list(maxit=1e4)));
      est = res$par;
      return(est);
    }

    y.list = unclass(data.frame(y));
    estpar = parSapply(cl,y.list,estimation,z=z,h.est=h.est)

    return(t(estpar));
}
