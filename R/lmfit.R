lmfit = function(y,x,z){
  fit = lm(y ~ x + z);
  pv = summary(fit)$coefficients[2,4];
  rs = summary(fit)$r.squared;
  return(c(pv,rs));
}
