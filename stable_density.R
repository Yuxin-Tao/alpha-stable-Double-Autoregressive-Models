### rebuild the package code "dstable", with improvements from Matsui(2006)
### give virtually the same results as "dstable" for most values of x and alpha

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
formula8 <- function(x,alpha,k){
  s = 0
  for (i in 0:k){
    s = s+gamma((2*i+1)/alpha)/factorial(2*i)*(-1)^i*x^(2*i)
  }
  return(s/(pi*alpha))
}

formula9 <- function(x,alpha,k){
  s = 0
  for (i in 1:k){
    s = s+gamma(i*alpha+1)/factorial(i)*(-1)^(i-1)*sin(pi*alpha*i/2)*x^(-i*alpha-1)
  }
  return(s/pi)
}

G <- function(v,y){
  logy2 = log(1+y^2); y2 = 1+y^2
  z = atan(y)
  sin = sin(v*z); cos = cos(v*z)
  phi = digamma(v); phi_d = trigamma(v); phi_dd = psigamma(v,2)
  ga = gamma(v); ga_d = phi*ga; ga_dd = ga*(phi_d+phi^2);
  ga_ddd = ga*(phi_dd+3*phi_d*phi+phi^3)
  g = y2^(-v/2)*(ga/4*logy2^2-ga_d*logy2+ga_dd)*(cos*(phi-logy2/2)-z*sin(v))+
    y2^(-v/2)*(-ga*logy2+2*ga_d)*(-z*sin*(phi-logy2/2)-z^2*cos+cos*phi_d)+
    y2^(-v/2)*ga*(-z^2*cos*(phi-logy2/2)-2*z*sin*phi_d+z^3*sin+cos*phi_dd)
  return(g)
}

formula10 <- function(x,alpha){
  logx2 = log(1+x^2)
  x2 = 1+x^2
  f1 = 1/(pi*x2)
  f1_a = 1/pi*((x^2-1)/(x2^2)*(1-C-logx2/2)+2*x/(x2^2)*atan(x))
  f1_aa = (x^4-6*x^2+1)/(pi*x2^3)*(pi^2/6+(1-C-logx2/2)^2-1-(atan(x))^2)+
    8*x*(x^2-1)/(pi*x2^3)*atan(x)*(1.5-C-logx2/2)+
    2/(pi*x2^3)*((1-3*x^2)*(1-C-logx2/2)-x*x2*atan(x))
  f1_aaa = 1/pi*(-G(4,x)+3*G(3,x)-G(2,x))
  s = f1+f1_a*(alpha-1)+1/2*f1_aa*(alpha-1)^2+1/6*f1_aaa*(alpha-1)^3
  return(s)
}

formula14 <- function(x,alpha){
  return(gamma(alpha+1)/2*x^(-alpha-1)*(2-alpha))
}

#### Revised numerical calculation of density of symmetric stable dist
dstable_new <- function(xx, alpha, tol=64*.Machine$double.eps, subdivisions=1000, log=FALSE)
{
  if (alpha==2) {
    return(dnorm(xx, mean = 0, sd = sqrt(2)))
  } else if (alpha==1){
    return(dcauchy(xx))
  }
  
  out = NULL
  for (k in 1:length(xx)){
    x = abs(xx[k])  # f(x) is even
    ### when alpha is close to 1
    if ((alpha>0.99)&(alpha<=1.01)){
      out[k] = formula10(x, alpha); next
    }
    
    ### when x is close to 0
    if ((alpha>=0.1) && (alpha<=0.2) && (x<1e-16)){
      out[k] = formula8(x,alpha,1); next
    }
    if ((alpha>0.2)&&(alpha<=0.5)&&(x<1e-8) || 
        (alpha>0.5)&&(alpha<=0.99)&&(x<1e-5)){
      out[k] = formula8(x,alpha,5); next
    }
    if ((alpha>1.01) && (alpha<=1.99999) && (x<1e-5)){
      out[k] = formula8(x,alpha,10); next
    }    
    if ((alpha>1.99999) && (alpha<=2) && (x<7)){
      out[k] = formula8(x,alpha,85); next
    }
    ### when x tends to infinity
    if (((alpha>=0.1)&&(alpha<=0.99) || 
         (alpha>1.01)&&(alpha<=1.99999)) && (x>10^(3/(1+alpha)))){
      out[k] = formula9(x,alpha,10); next
    }
    
    ### other cases
    zeta.tol = .4e-15
    theta0 = 0
    f.zeta = gamma(1+1/alpha)/pi
    if(is.finite(x) && x <= zeta.tol * (zeta.tol+ x)) { ## if x is too small
      out[k] = f.zeta
      next
    }
    ## the real check should be about the feasibility of g() below, or its integration
    
    a_1 <- alpha - 1
    ##' g() is strictly monotone -- Nolan(1997) ["3. Numerical Considerations"]
    ##'     alpha >= 1  <==>  g() is falling, ie. from Inf --> 0;  otherwise growing from 0 to +Inf
    g <- function(th) {
      r <- th
      # r[which(r==-0)] = +0
      ## g(-pi/2) or g(pi/2) could become  NaN --> work around
      i.bnd <- abs(pi/2 -sign(a_1)*th) < 64*.Machine$double.eps
      r[i.bnd] <- 0
      th <- th[io <- !i.bnd]
      att <- alpha*th ## = alpha*(theta)
      r[io] <- (cos(th) * (x/sin(att))^alpha)^(1/a_1) * cos(att-th)
      r
    }
    ## Function to integrate: dstable(..)= f(..) = c2 * \int_{-\theta_0}^{\pi/2} g1(u) du
    g1 <- function(th) {
      ## g1 :=  g(.) exp(-g(.))
      x.exp.m.x( g(th) )
    }
    c2 <- ( alpha / (pi*abs(a_1)*x) )
    
    ## Now, result = c2 * \int_{-t0}^{pi/2}  g1(u) du  ,  we "only" need the integral
    ## where however, g1(.) may look to be (almost) zero almost everywhere and just have a small peak
    ## ==> Find the peak, split the integral into two parts of for intervals  (t0, t_max) + (t_max, pi/2)
    
    ## However, this may still be bad, e.g., for dstable(71.61531, alpha=1.001, beta=0.6),
    ## or  dstable(1.205, 0.75, -0.5)
    ##   the 2nd integral was "completely wrong" (basically zero, instead of ..e-5)
    
    ## NB: g() is monotone, see above
    if((alpha >= 1 &&
        ((!is.na(g. <- g( pi2	)) && g. > .large.exp.arg) || identical(g(0), 0))) ||
       (alpha  < 1 &&
        ((!is.na(g. <- g(0)) && g. > .large.exp.arg) || identical(g(pi2), 0)))) {
      ## g() is numerically too large *or* 0 even where it should be inf
      ## ===>	 g() * exp(-g()) is 0 everywhere
      out[k] = if(log)-Inf else 0
      next
    }
    
    g. <- if(alpha >= 1) g(.e.plus(0, 1e-6)) else g(pi2..(1e-6))
    if (is.na(g.))# g() is not usable --- FIXME rather use *asymptotic dPareto()?
      if (x < .01) {
        out[k] = f.zeta
        next
      }

    Int <- function(a,b){.integrate2(g1, lower = a, upper = b,
                                     subdivisions=subdivisions, rel.tol= tol, abs.tol= tol)}
    
    ## We know that the maximum of g1(.) is = exp(-1) = 0.3679  "at" g(.) == 1
    ## find that by uniroot :
    ## g(.) == 1  <==>  log(g(.)) == 0   --- the latter is better conditioned,
    ##                                       e.g., for (x = -1, alpha = 0.95, beta = 0.6)
    ## the former is better for  dstable(-122717558, alpha = 1.8, beta = 0.3, pm = 1)
    ## However, it can be that the maximum is at the boundary,  and
    ## g(.) > 1 everywhere or  g(.) < 1  everywhere  {in that case we could revert to optimize..}
    
    if((alpha >= 1 && !is.na(g. <- g(pi2)) && g. > 1) ||
       (alpha <	 1 && !is.na(g. <- g(pi2)) && g. < 1))
      g1.th2 <- g1( theta2 <- pi2..(1e-6) )
    else if((alpha <  1 && g(0) > 1) ||
            (alpha >= 1 && g(0) < 1))
      g1.th2 <- g1( theta2 <- .e.plus(0, 1e-6) )
    else {
      ## when alpha ~=< 1 (0.998 e.g.),  g(x) is == 0 (numerically) on a wide range;
      ## uniroot is not good enough, and we should *increase* -theta0
      ## or decrease pi2 such that it can find the root:
      l.th <- 0
      u.th <- pi2
      if(alpha < 1) { ## g() is *in*creasing from 0 ..
        while ((g.t <- g(.th <- (l.th + pi2)/2)) == 0) l.th <- .th
        if(g.t == 1)# decrease upper limit {needed, e.g. for alpha = 1e-20}
          while ((g.t <- g(.th <- (l.th + u.th)/2)) == 1) u.th <- .th
        if(abs(u.th - l.th) < 1e-13) {# do not trust g()
          out[k] = if(log)-Inf else 0
          next
        }
      }
      
      ur1 <- uniroot(function(th) g(th) - 1,
                     lower = l.th, upper = u.th, tol = .Machine$double.eps)
      ## consider using safeUroot() [ ~/R/Pkgs/copula/R/safeUroot.R ] !!
      ur2 <- tryCatch(uniroot(function(th) log(g(th)),
                              lower = l.th, upper = u.th, tol = .Machine$double.eps),
                      error=function(e)e)
      g.1 <- x.exp.m.x(ur1$f.root+1)
      g.2 <- if(inherits(ur2, "error")) -Inf else x.exp.m.x(exp(ur2$f.root))
      if(g.1 >= g.2) {
        theta2 <- ur1$root
        g1.th2 <- g.1 ## == g1(theta2)
      } else {
        theta2 <- ur2$root
        g1.th2 <- g.2
      }
    }
    ## now, because g1()'s peak (at th = theta2) may be extreme, we find two more intermediate values
    ## NB: Theoretically: Max = 0.3679 = g1(theta2)  ==> 1e-4 is a very small fraction of that
    ## to the left:
    eps <- 1e-4
    if((do1 <- g1.th2 > eps && g1(0) < eps))
      th1 <- uniroot(function(th) g1(th) - eps, lower = 0, upper = theta2,
                     tol = tol)$root
    if((do4 <- g1.th2 > eps && g1(pi2) < eps))
      ## to the right:
      th3 <- uniroot(function(th) g1(th) - eps, lower = theta2, upper = pi2,
                     tol = tol)$root
    
    if(do1) {
      r1 <- Int(0, th1)
      r2 <- Int(         th1, theta2)
    } else {
      r1 <- 0
      r2 <- Int(0,      theta2)
    }
    if(do4) {
      r3 <- Int(              theta2, th3)
      r4 <- Int(                      th3, pi2)
    } else {
      r3 <- Int(              theta2,      pi2)
      r4 <- 0
    }
    ans = c2*(r1+r2+r3+r4)
    out[k] = ans
    
    ### when alpha is close to 2
    if (alpha>1.99999) {
      out[k] = max(ans, formula14(x,alpha))
    }
  }

  return(out)
}
