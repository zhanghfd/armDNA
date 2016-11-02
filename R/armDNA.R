armDNA <-
  function(dat,pow=NULL,age.num=5,cl.cores=1){
    
    p = ncol(dat$z);
    
    cl.cores = as.integer(cl.cores);
    age.num = as.integer(age.num);
    if(!is.null(pow)){
      pow = as.numeric(pow);
    }
    if(cl.cores < 1){
      stop("The number of cores, cl.cores, should be a positive integer.")      
    }

    m = ncol(dat$y); # m: number of mDNAs
      
    alpha = gamma = rep(0,m);
    beta =  matrix(0,m,p);   # covariate effect
    phi = rep(1,m);
    par = cbind(alpha,beta,gamma,phi); # par: initial values of parameters alpha, beta, gamma, phi
    
    epsilon = 1e-6;
    x = dat$x;
    y = dat$y;
    z = dat$z;

    row.names(y) = NULL;
    x = round(x,2);
    y = epsilon/2 + y * (1-epsilon);
    n = nrow(y); # sample size
    
    xmin = round(min(x),2);
    xmax = round(max(x),2);
    
    par1 = par;
    
    if(is.null(pow)){
      if(age.num<2){
        stop("The number of age nodes, age.num, should be an integer greater than 1.")
      }
      
      x.sort = round(as.vector(quantile(dat$x,0:age.num/age.num)),2);
      n.x = length(x.sort); # number of distinct ages
      
      ui = matrix(0,n.x-1,n.x-2);
      diag(ui) = 1;
      for(i in 1:(n.x-2)){
        ui[i+1,i] = -1;
      }
      ci = rep(0,n.x-1);
      ci[n.x-1] = -1;
      
      ini = c(0.2,0.5,1,2,5); # candidate initial seeds
      y.list = unclass(data.frame(y));
      x.norm = (x-xmin)/(xmax-xmin);
      cl = makeCluster(cl.cores);

      naive.res = matrix(NA,10,m);
      for(i in 1:5){
        naive.res[c(i,i+5),] = parSapply(cl,y.list,lmfit,x.norm^ini[i],z)
      }
      rss = rowSums(naive.res[6:10,]);  
      index = which(rss==max(rss)); # initial seed selected
      used = (naive.res[index,]<0.05); # preliminary screening of marks
      
      h0 = ((x.sort-xmin)/(xmax-xmin))^ini[index];
      h0.est = H.est(x.sort,h0,x);
      par2 = est.par(y[,used],z,h0.est,par1[used,],cl);
      h1 = h0;
      h2 = h1 + 1;
      
      # iteration algorithm to estimate h_hat
      while(max(abs(h2-h1))>1e-3){
        par3 = par2;
        h2 = h1;
        h1 = est.h(x.sort,x,y[,used],z,h2,par3,ui,ci);
        h1.est = H.est(x.sort,h1,x);
        par2 = est.par(y[,used],z,h1.est,par3,cl);
      }

      par1 = est.par(y,z,h1.est,par1,cl);
      sd = hess(y,z,h1.est,par1,cl);
      xh = data.frame(cbind(x.sort,h1));
      names(xh) = c('age','h.est');
      stopCluster(cl);
      
    }else{
      if(!is.null(pow) | pow>0){
        cl = makeCluster(cl.cores);    
        h1.est = ((x-xmin)/(xmax-xmin))^pow; # initial h
        par1 = est.par(y,z,h1.est,par1,cl);
        sd = hess(y,z,h1.est,par1,cl);
        xh = NULL;
        stopCluster(cl);
      }else{
        stop("The power parameter used in the transformation function, pow, should be positive.")
      }
    }

    id = c(p+2,2:(p+1));
    par = par1[,id];
    se = sd[,id];
    p.value = 1-pchisq((par/se)^2,1);

    res = list(par=par,se=se,p.value=p.value,age.h=xh);
    
    return(res);
    
  }
    