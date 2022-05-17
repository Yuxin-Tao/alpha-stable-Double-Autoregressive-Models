### The derivative of f(x) on alpha, with improvements from Matsui(2006)

### The following formulas follow from the expression in Matsui(2006)
formula26 <- function(x,alpha,k){
  s = 0
  for (i in 1:k){
    s = s+digamma((2*i+1)/alpha+1)*gamma((2*i+1)/alpha+1)/
      factorial(2*i)*(-x^2)^i
  }
  return(-s/(pi*alpha^2))
}

formula27 <- function(x,alpha,k){
  s = 0
  for (i in 1:k){
    s = s+digamma(i*alpha+1)*gamma(i*alpha+1)/factorial(i-1)*(-1)^(i-1)*
      sin(pi*alpha*i/2)*x^(-i*alpha-1)+
      gamma(i*alpha+1)/factorial(i-1)*(-1)^(i-1)*(pi/2*cos(pi*alpha*i/2)-
      log(x)*sin(pi*alpha*i/2))*x^(-i*alpha-1)
  }
  return(s/pi)
}

formula28 <- function(x,alpha){
  logx2 = log(1+x^2)
  x2 = 1+x^2
  f1_a = 1/pi*((x^2-1)/(x2^2)*(1-C-logx2/2)+2*x/(x2^2)*atan(x))
  f1_aa = (x^4-6*x^2+1)/(pi*x2^3)*(pi^2/6+(1-C-logx2/2)^2-1-(atan(x))^2)+
    8*x*(x^2-1)/(pi*x2^3)*atan(x)*(1.5-C-logx2/2)+
    2/(pi*x2^3)*((1-3*x^2)*(1-C-logx2/2)-x*x2*atan(x))
  s = f1_a + f1_aa*(alpha-1)
  return(s)
}

#### Revised numerical calculation of df(x)/dalpha of symmetric stable dist
#### For formula(23), directly integrate without dividing into three parts
#### Numerically more stable
dstable_da_new <- function(xx, alpha, tol=64*.Machine$double.eps, subdivisions=1000, log=FALSE)
{
  # if (alpha==2) {
  #   return(1/(2*sqrt(pi))*exp(-xx^2/4)*(-xx/2))
  # } else 
  if (alpha==1){
    return(((1-C)*(xx^2-1)-(xx^2-1)/2*log(xx^2+1)+2*xx*atan(xx))/(pi*(xx^2+1)^2))
  }
  
  out = NULL
  for (k in 1:length(xx)){
    x = abs(xx[k])       # df(x)/dalpha is even
    ### when alpha is close to 1
    if ((alpha>0.989)&(alpha<=1.021)){ ## a bit different from Matsui(2006)
      out[k] = formula28(x, alpha); next
    }
    
    ### when x is close to 0
    if ((alpha>=0.1)&&(alpha<=0.2) && (x<=1e-16)){
      out[k] = formula26(x,alpha,1); next
    } 
    if ((alpha>0.2)&&(alpha<=0.3)&&(x<1e-7) || 
        (alpha>0.3)&&(alpha<=0.99)&&(x<1e-5)){
      out[k] = formula26(x,alpha,5); next
    }
    if ((alpha>1.01)&&(alpha<=1.9999)&&(x<1e-5)){
      out[k] = formula26(x,alpha,10); next
    }
    if ((alpha>1.9999)&&(alpha<=2)&&(x<=8)){
      out[k] = formula26(x,alpha,85); next
    }    
    
    ### when x tends to infinity
    if (((alpha>=0.1)&&(alpha<=0.99) || 
         (alpha>1.01)&&(alpha<=1.9999)) && (x>10^(3/(1+alpha)))){
      out[k] = formula27(x,alpha,10); next
    }
    if ((alpha>1.9999)&&(alpha<=2)&&(x>8)){
      out[k] = formula27(x,alpha,20); next
    } 
    
    ### Other cases, following formula (23)
    zeta.tol = .4e-15
    f.zeta = -1/(pi*alpha^2)*digamma(1+1/alpha)*gamma(1+1/alpha)
    if(is.finite(x) && x <= zeta.tol * (zeta.tol+ x)) { ## if x is too small
      out[k] = f.zeta
      next
    }

    ## Function to integrate: f(..) = c2 * \int_{0}^{\pi/2} g*h*(1-g)exp(-g) du
    ## Without splitting the integral
    
    a_1 <- alpha - 1
    g <- function(th) {
      r <- th
      i.bnd <- abs(pi/2 -sign(a_1)*th) < 64*.Machine$double.eps
      r[i.bnd] <- 0
      th <- th[io <- !i.bnd]
      att <- alpha*th ## = alpha*(theta)
      r[io] <- (cos(th) * (x/sin(att))^alpha)^(1/a_1) * cos(att-th)
      r
    }
    g_deriv <- function(th) {
      att <- alpha*th
      r = (x*cos(th)/sin(att))^(alpha/a_1) * ((-1/a_1^2*log(x*cos(th)/sin(att))-
                                                 alpha/a_1*cos(att)*th/sin(att))*cos(att-th)/cos(th)-sin(att-th)*th/cos(th))
      if(any(nax <- is.na(th)))
        r[nax] <- NA_real_
      # if(any(lrg <- !nax & th > .large.exp.arg))# e.g. th == Inf
      #   r[lrg] <- 0
      r
    }
    
    if((alpha >= 1 &&
        ((!is.na(g. <- g( pi2	)) && g. > .large.exp.arg) || identical(g(0), 0))) ||
       (alpha  < 1 &&
        ((!is.na(g. <- g(0)) && g. > .large.exp.arg) || identical(g(pi2), 0)))) {
      ## g() is numerically too large *or* 0 even where it should be inf
      ## ===>	 g() * exp(-g()) is 0 everywhere
      out[k] = 0
      next
    }
    
    ## Function to integrate: f(..) = c2 * \int_{-\theta_0}^{\pi/2} g1(u) du
    g1 <- function(th) {
      r <- g_deriv(th)*exp(-g(th))*(1-g(th))
      if(any(nax <- is.na(th)))
        r[nax] <- NA_real_
      if(any(lrg <- !nax & g(th) > .large.exp.arg))# e.g. g(th) == Inf
        r[lrg] <- 0
      r
    }
    c2 = alpha / (pi*abs(a_1)*x)
    
    I = .integrate2(g1,0,pi/2,subdivisions=subdivisions, rel.tol= tol, abs.tol= tol)
    ans = c2 * I
    ans = -1/(a_1*alpha)*dstable(x,alpha=alpha,beta=0) + ans
    out[k] = ans
  }
  
  ## a bit different from Matsui(2006)
  i0 <- (out == 0) & (abs(xx)>1e10) # --> we can do better using asymptotic:
  if (any(i0)) {   ### (xx[i0]>0)*2-1 is equal to sign(xx[i0])
    out[i0] <- (digamma(alpha)*gamma(alpha)/pi*sin(alpha*pi2)+
                  gamma(alpha)/2 * cos(alpha*pi2))* alpha*abs(xx[i0])^(-(1+alpha))+
      gamma(alpha)/pi * sin(alpha*pi2)*(1-alpha*log(abs(xx[i0])))*abs(xx[i0])^(-(1+alpha))
  }
  
  return(out)
}
