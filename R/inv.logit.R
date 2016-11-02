inv.logit <-
function(x){
    tmp = exp(x);
    return(tmp/(1+tmp));
}
