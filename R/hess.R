hess <-
function(y,z,h.est,par,cl){

    grr = function(pary,z,h.est){
        n.z = ncol(z);
        alpha = pary[1];
        beta = pary[2:(n.z+1)];
        gamma = pary[n.z+2];
        phi = pary[n.z+3];
        yj = pary[-(1:(n.z+3))];
        mu = exp(alpha+as.vector(z%*%beta)+gamma*h.est)/(1+exp(alpha+as.vector(z%*%beta)+gamma*h.est));
        p = phi * mu;
        q = phi * (1-mu);
        digamma.p = digamma(p);
        digamma.q = digamma(q);
        digamma.phi = digamma(phi);
        trigamma.p = trigamma(p);
        trigamma.q = trigamma(q);
        trigamma.phi = trigamma(phi);
        y.tmp = log(yj/(1-yj));
        tmp1 = phi * (mu * (1-mu) * (1-2*mu) * (y.tmp - digamma.p + digamma.q) + phi*mu^2*(1-mu)^2*(-trigamma.p-trigamma.q));
        tmp2 = mu*(1-mu)*(y.tmp-digamma.p+digamma.q-phi*mu*trigamma.p+phi*(1-mu)*trigamma.q);
    
        alpha.alpha = sum(tmp1);
        alpha.beta = as.vector(tmp1%*%z);
        alpha.gamma = sum(tmp1*h.est);
        alpha.phi = sum(tmp2);
        beta.beta = (t(tmp1*z)%*%z);
        beta.gamma = as.vector((tmp1*h.est) %*% z);
        beta.phi = as.vector(tmp2%*%z);
        gamma.gamma = sum(tmp1*h.est^2);
        gamma.phi = sum(tmp2*h.est);
        phi.phi = sum(-mu^2*trigamma.p-(1-mu)^2*trigamma.q+trigamma.phi);
        mat1 = c(alpha.alpha,alpha.beta,alpha.gamma,alpha.phi);
        mat2 = cbind(alpha.beta,beta.beta,beta.gamma,beta.phi);
        mat3 = c(alpha.gamma,beta.gamma,gamma.gamma,gamma.phi);
        mat4 = c(alpha.phi,beta.phi,gamma.phi,phi.phi);
        mat = rbind(mat1,mat2,mat3,mat4);
        
        res = sqrt(abs(diag(solve(-mat))));
        return(res)
    }
    
    pary = rbind(t(par),y);
    pary.list = unclass(data.frame(pary));
    hes = parSapply(cl,pary.list,grr,z=z,h.est=h.est);
    return(t(hes));
    
}
