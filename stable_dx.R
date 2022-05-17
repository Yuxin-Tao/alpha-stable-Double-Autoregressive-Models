### The derivative of f(x) on x, with improvements from Matsui(2006)

C = 0.57721566490153286060     # Euler constant
.large.exp.arg <- -(.Machine$double.min.exp * log(2)) ## == 708.396...
pi2 = pi/2
##  x*exp(-x)  numerically stably, with correct limit 0 for x --> Inf
x.exp.m.x <- function(x) {
  r <- x*exp(-x)
  if(any(nax <- is.na(x)))
    r[nax] <- NA_real_
  if(any(lrg <- !nax & x > .large.exp.arg))# e.g. x == Inf
    r[lrg] <- 0
  r
}
x2.exp.m.x <- function(x) {
  r <- x^2 * exp(-x)
  if(any(nax <- is.na(x)))
    r[nax] <- NA_real_
  if(any(lrg <- !nax & x > .large.exp.arg))# e.g. x == Inf
    r[lrg] <- 0
  r
}
.e.plus <- function(x, eps) x + eps* abs(x)
.e.minus<- function(x, eps) x - eps* abs(x)
pi2.. <- function(eps) pi2 * (1 - eps) ## == .e.minus(pi/2, eps), slight more efficiently

.integrate2 <- function(f, lower, upper, ..., subdivisions, rel.tol, abs.tol,
                        stop.on.error = FALSE)
{### Numerically Integrate, no errors, but warnings
  ri <- integrate(f, lower, upper, ..., subdivisions=subdivisions,
                  rel.tol=rel.tol, abs.tol=abs.tol, stop.on.error=stop.on.error)
  if((msg <- ri[["message"]]) != "OK")
    warning(msg) ## NB: "roundoff error ..." happens many times
  ri[["value"]]
}

C.stable.tail <- function(alpha, log = FALSE) {
  stopifnot(0 <= alpha, alpha <= 2)
  r <- alpha
  i0 <- alpha == 0
  r[i0] <- if(log) -log(2) else 0.5
  r[!i0] <- gamma(alpha)/pi * sin(alpha*pi2)
  if(any(a2 <- alpha == 2)) r[a2] <- if(log) -Inf else 0
  r
}

dPareto <- function(x, alpha, beta, log = FALSE) {
  if(any(neg <- x < 0)) { ## left tail
    x[neg] <- -x[neg]
    beta <- rep(beta, length.out=length(x))
    beta[neg] <- -beta[neg]
  }
  alpha*(1+beta)* C.stable.tail(alpha)* x^(-(1+alpha))
}

### The following formulas follow from the expression in Matsui(2006)
formula16 <- function(x,alpha,k){
  s = 0
  for (i in 1:k){
    s = s+gamma((2*i+1)/alpha)/factorial(2*i-1) * (-1)^i * x^(2*i-1)
  }
  return(s/(pi*alpha))
}

formula17 <- function(x,alpha,k){
  s = 0
  for (i in 1:k){
    s = s+gamma(i*alpha+2)/factorial(i)*(-1)^(i)*sin(pi*alpha*i/2)*x^(-i*alpha-2)
  }
  return(s/pi)
}

formula18 <- function(x,alpha){
  logx2 = log(1+x^2)
  x2 = 1+x^2
  f1 = -2*x/(pi*x2^2)
  f1_a = 1/pi*((-2*x^3+6*x)/(x2^3)*(1.5-C-logx2/2)+(2-6*x^2)/(x2^3)*atan(x))
  f1_aa = -2*x*(x^4-14*x^2+9)/(pi*x2^4)*(pi^2/6+(1-C-logx2/2)^2-1-(atan(x))^2)-
    8*(3*x^4-8*x^2+1)/(pi*x2^4)*atan(x)*(1.5-C-logx2/2)-
    2*x*(x^4-22*x^2+17)/(pi*x2^4)*(1-C-logx2/2)-
    4*(x^4-6*x^2+1)/(pi*x2^4)*atan(x)+8*x*(x^2-1)/(pi*x2^4)
  s = f1+f1_a*(alpha-1)+1/2*f1_aa*(alpha-1)^2
  return(s)
}


#### Revised numerical calculation of df(x)/dx of density of symmetric stable dist
dstable_dx_new <- function(xx, alpha, tol=64*.Machine$double.eps, subdivisions=1000, log=FALSE)
{
  if (alpha==2) {
    return(1/(2*sqrt(pi))*exp(-xx^2/4)*(-xx/2))
  } else if (alpha==1){
    return(-2*xx/(pi*(1+xx^2)^2))
  }
  
  out = NULL
  for (k in 1:length(xx)){
    x = abs(xx[k])       # f(x) is even
    sign_x = sign(xx[k]) # df/dx is odd
    ### when alpha is close to 1
    if ((alpha>0.989)&(alpha<=1.021)){ ## a bit different from Matsui(2006)
      out[k] = sign_x*formula18(x, alpha); next
    }
    
    ### when x is close to 0
    if ((alpha>=0.1)&&(alpha<=0.25)&&(x<1e-8) ||  ## a bit different from Matsui(2006)
        (alpha>0.25)&&(alpha<=0.3)&&(x<1e-6) ||
        (alpha>0.3)&&(alpha<=0.99)&&(x<1e-5)){
      out[k] = sign_x*formula16(x,alpha,5); next
    }
    if ((alpha>1.01) && (x<1e-3)){
      out[k] = sign_x*formula16(x,alpha,10); next
    }    

    ### when x tends to infinity
    if (((alpha>=0.1)&&(alpha<=0.99) ||    ## a bit different from Matsui(2006)
         (alpha>1.01)&&(alpha<=1.99999)) && (x>10^(3/(1+alpha)))){
      out[k] = sign_x*formula17(x,alpha,10); next
    }
    
    ### Other cases, following formula (15)
    zeta.tol = .4e-15
    theta0 = 0
    f.zeta = 0
    if(is.finite(x) && x <= zeta.tol * (zeta.tol+ x)) { ## if x is too small
      out[k] = f.zeta
      next
    }

    a_1 <- alpha - 1
    ##' g() is strictly monotone -- Nolan(1997) ["3. Numerical Considerations"]
    ##'     alpha >= 1  <==>  g() is falling, ie. from Inf --> 0;  otherwise growing from 0 to +Inf
    g <- function(th) {
      r <- th
      i.bnd <- abs(pi/2 -sign(a_1)*th) < 64*.Machine$double.eps
      r[i.bnd] <- 0
      th <- th[io <- !i.bnd]
      att <- alpha*th ## = alpha*(theta)
      r[io] <- (cos(th) * (x/sin(att))^alpha)^(1/a_1) * cos(att-th)
      r
    }
    ## Function to integrate: f(..) = c2 * \int_{-\theta_0}^{\pi/2} g1(u) du
    g1 <- function(th) {
      ## g1 :=  g(.)^2 * exp(-g(.))
      x2.exp.m.x( g(th) )
    }
    c2 <- alpha^2 / (pi*abs(a_1)*(a_1)*x^2)
    
    ## Now, result = c2 * \int_{-t0}^{pi/2}  g1(u) du  ,  we "only" need the integral
    ## where however, g1(.) may look to be (almost) zero almost everywhere and just have a small peak
    ## ==> Find the peak, split the integral into two parts of for intervals  (t0, t_max) + (t_max, pi/2)
    
    ## NB: g() is monotone, see above
    if((alpha >= 1 &&
        ((!is.na(g. <- g( pi2	)) && g. > .large.exp.arg) || identical(g(0), 0))) ||
       (alpha  < 1 &&
        ((!is.na(g. <- g(0)) && g. > .large.exp.arg) || identical(g(pi2), 0)))) {
      ## g() is numerically too large *or* 0 even where it should be inf
      ## ===>	 g() * exp(-g()) is 0 everywhere
      out[k] = 0
      next
    }
    
    g. <- if (alpha >= 1) g(.e.plus(0, 1e-6)) else g(pi2..(1e-6))
    if (is.na(g.)) # g() is not usable --- FIXME rather use *asymptotic dPareto()?
      if (x < .01){
        out[k] = f.zeta
        next
      }
    
    Int <- function(a,b){.integrate2(g1, lower = a, upper = b,
                                     subdivisions=subdivisions, rel.tol= tol, abs.tol= tol)}
    
    ## We know that the maximum of g1(.) is = 4*exp(-2) = 0.5413  "at" g(.) == 2
    ## find that by uniroot :
    ## g(.) == 2  <==>  log(g(.)-1) == -log(2)
    ## However, it can be that the maximum is at the boundary,  and
    ## g(.) > 2 everywhere or  g(.) < 2  everywhere  {in that case we could revert to optimize..}
    
    if ((alpha >= 1 && !is.na(g. <- g(pi2)) && g. > 2) ||
        (alpha <	1 && !is.na(g. <- g(pi2)) && g. < 2)){
      g1.th2 <- g1( theta2 <- pi2..(1e-6) )
    } else if((alpha <  1 && g(0) > 2) ||
              (alpha >= 1 && g(0) < 2)){
      g1.th2 <- g1( theta2 <- .e.plus(0, 1e-6) )
    } else {
      ## when alpha ~=< 1 (0.998 e.g.),  g(x) is == 0 (numerically) on a wide range;
      ## uniroot is not good enough, and we should *increase* -theta0
      ## or decrease pi2 such that it can find the root:
      l.th <- 0
      u.th <- pi2
      if(alpha < 1) { ## g() is *increasing from 0 ..
        while ((g.t <- g(.th <- (l.th + pi2)/2)) == 0) l.th <- .th
        if(g.t == 2)# decrease upper limit {needed, e.g. for alpha = 1e-20}
          while ((g.t <- g(.th <- (l.th + u.th)/2)) == 1) u.th <- .th
        if(abs(u.th - l.th) < 1e-13) {# do not trust g()
          out[k] = if(log)-Inf else 0
          next
        }
      }
      
      ur1 <- uniroot(function(th) g(th) - 2,
                     lower = l.th, upper = u.th, tol = .Machine$double.eps)
      ur2 <- tryCatch(uniroot(function(th) log(g(th))-log(2),
                              lower = l.th, upper = u.th, tol = .Machine$double.eps),
                      error=function(e)e)
      g.1 <- x2.exp.m.x(ur1$f.root+2)
      g.2 <- if(inherits(ur2, "error")) -Inf else x2.exp.m.x(2*exp(ur2$f.root))
      if(g.1 >= g.2) {
        theta2 <- ur1$root
        g1.th2 <- g.1 ## == g1(theta2)
      } else {
        theta2 <- ur2$root
        g1.th2 <- g.2
      }
    }
    
    ## now, because g1()'s peak (at th = theta2) may be extreme, we find two more intermediate values
    ## NB: Theoretically: Max = 0.5413 = g1(theta2)  ==> 1e-4 is a very small fraction of that
    eps <- 1e-4
    if((do1 <- g1.th2 > eps && g1(0) < eps))    ## to the left:
      th1 <- uniroot(function(th) g1(th) - eps, lower = 0, upper = theta2,
                     tol = tol)$root
    if((do4 <- g1.th2 > eps && g1(pi2) < eps))  ## to the right:
      th3 <- uniroot(function(th) g1(th) - eps, lower = theta2, upper = pi2,
                     tol = tol)$root
    
    if(do1) {
      r1 <- Int(0, th1)
      r2 <- Int(   th1, theta2)
    } else {
      r1 <- 0
      r2 <- Int(0,      theta2)
    }
    if(do4) {
      r3 <- Int(        theta2, th3)
      r4 <- Int(                th3, pi2)
    } else {
      r3 <- Int(        theta2,      pi2)
      r4 <- 0
    }
    ans = c2*(r1+r2+r3+r4)
    ans = 1/(a_1*x)*dstable(x,alpha=alpha,beta=0)-ans
    ans = sign_x*ans
    out[k] = ans
  }
  
  ## a bit different from Matsui(2006)
  i0 <- (out == 0) & (abs(xx)>1e10) # --> we can do better using asymptotic:
  if (any(i0)) {   ### (xx[i0]>0)*2-1 is equal to sign(xx[i0])
    out[i0] <- -((xx[i0]>0)*2-1)*(1+alpha)*alpha * C.stable.tail(alpha)*abs(xx[i0])^(-(2+alpha))
  }
  
  return(out)
}
