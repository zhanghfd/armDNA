H.est <-
function(x.sort,h,x){
    res = rep(NA,length(x));
    tmp = (x <= x.sort[2]);
    res[tmp] = (x[tmp] - x.sort[1]) / (x.sort[2] - x.sort[1]) * h[2];
    for(j in 2:(length(h)-1)){
        tmp = (x > x.sort[j] & x <= x.sort[j+1]);
        res[tmp] = (x[tmp] - x.sort[j]) / (x.sort[j+1] - x.sort[j]) * (h[j+1] - h[j]) + h[j];
    }
    return(res);
}
