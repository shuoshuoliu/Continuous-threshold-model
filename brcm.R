library(survival)
library(MASS)
###  generic function for Biomarker Rectifier Continuous Models(brcm) using the
###  Rectifier Linear Unit (ReLU).
brcm = function(x, ...) {
  UseMethod("brcm")
}

### the follow function take input like
###
###       brcm(Surv(time, status) ~ trt + z1 + z2 + z3 + w)
###
### and will fit a Cox PH model with
###
###       lambda_0(t) exp(b1*trt + b2*z1 + b3*z2 + b4*z3 +
###                        b5*max(w-c, 0) + b6*trt*max(w-c, 0))
###
brcm.formula = function(formula, data, subset, robust = FALSE,
                        method = c("gradient", "profile"), epsilon =  NULL, ...) {
  method = match.arg(method)

  mf <- match.call()
  m <- match(c("formula", "data", "subset"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())

  #mf = model.frame(formula = formula, data = data)
  x = model.matrix(attr(mf, "terms"), data = mf)
  y = model.response(mf)

  if (class(y) == "Surv") {
    family = "surv"
    st = sort(y[, 1], decreasing = TRUE, index.return = TRUE)
    idx = st$ix
    y = y[idx, ]
    x = x[idx, ]
    ## remove the intercept
    x = x[, -1]
  }

  ### Do box cox transformation here if w > 0.
  p = length(x[1, ])
  w = x[, p]

  #if(min(w)<=0) stop("Please adjust biomarker to positive for box-cox transformation")

  if(min(w)>0) {
    bx=boxcox(w~1, plotit=FALSE)
    lambda = bx$x[which.max(bx$y)]
    x[, p] = ((w)^lambda-1)/lambda
  } else
    lambda = 2

  w = x[, p]

  ### Use 20 points for first pass search
  if(is.null(epsilon)) epsilon = (max(w) - min(w))/20;

  ### seq_c for search grid points
  seq_c=seq(min(w), max(w), epsilon)

  fit = brcm.default(x, y, method, robust = FALSE, seq_c, ...)
  fit$call = match.call()
  fit$formula = formula
  #fit$lambda = lambda
  fit$method = method

  return(fit)
}

brcm.default = function (x, y, method, robust = FALSE, seq_c, ...) {
  x = as.matrix(x)

  method = method

  ### first pass search for initial value
  fit0 = reluProFit(x, y, seq_c)
  theta0 = c(fit0$beta, fit0$c.max)

  if (method == "gradient")
    fit = reluGradFit(x, y, theta0)
  else {
    lgc = as.matrix(fit0$log_c)

    ### second pass: refined search around possible MLE
    c2 = lgc[lgc[, 2] > (fit0$logLik - 3.84), 1]

    epsilon = seq_c[2]-seq_c[1]
    epsilon2 = (max(c2) - min(c2) + epsilon)/60
    seq_c2 = seq(min(c2)- 0.5*epsilon, max(c2) + 0.5*epsilon, epsilon2)
    fit = reluProFit(x, y, seq_c2)
    tmp = rbind(lgc, fit$log_c)
    fit$log_c = tmp[order(tmp[, 1]), ]
  }

  ### Find s.e
  grdInfo = reluGradient(fit$theta, x, y, info = TRUE)
  pi_theta = grdInfo$pi_theta

  #score in specified model
  sV = ginv(pi_theta)

  #numerical method in specified
  #I_numerical = .reluInfo(fit$theta, x, y)
  #I_numerical_inv = ginv(I_numerical)

  #robust var by second derivative in misspecified
  second = grdInfo$second
  robust_second = ginv(second)%*%pi_theta%*%ginv(second)

  if (isFALSE(robust)) fit$sd.score = sqrt(diag(sV))         ######
  else fit$sd.score = sqrt(diag(robust_second))

  fit$call = match.call()

  class(fit) = "brcm"
  return(fit)
}

########### help functions ###########################
.evalLP = function(theta, x) {
  ### this code is re-written to handle model of
  ###
  ###            Surv(time, status) ~ trt + x2 + x3 + ... + biomarker
  ###
  ### input 'x' is a nxp matrix with 2 or more columns:
  ### x1 = treatment z,
  ### x2, ..., x_{p-1} are other covariates,
  ### xp = biomarker variable w.
  p = length(x[1, ])
  p1 = p-1
  x1 = x[, 1]
  z = x[, 1:p1]
  w = x[, p]

  ### beta = theta1, ... theta[p+1], cut point = theta[p+2]
  beta = theta[1:(p+1)]
  cx   = theta[p+2]

  phi = ifelse(w >= cx, w-cx, 0)
  eta = ifelse(w >= cx, -(beta[p]+beta[p+1]*x1), 0)
  eta_i = -(beta[p]+beta[p+1]*x1)

  ### output X is a n x (p+1) matrix: [x1, ... x_{p-1}, phi, phi*x1]
  X = cbind(z, phi, phi*x1)

  ### Linear predict vector lp
  lp = (X%*%beta)[, 1]
  return(list(X = X, lp = lp, eta = eta, eta_i = eta_i, w = w, cx = cx, phi = phi))
}


### Numerical method for information matrix
.reluInfo = function(theta, x, y){
  p = length(theta)
  h = 0.0001
  Info = matrix(0, p, p)
  for(i in 1:p) {
    theta1 = theta
    theta2 = theta
    theta1[i] = theta[i] - h
    theta2[i] = theta[i] + h
    Info[i, ] = .5*(reluGradient(theta2, x, y) - reluGradient(theta1, x, y))/h  #should be theta1-theta2
  }
  return(Info)
}

###################################################
# this program was validated on September 12, 2018
###################################################
reluLglik=function(theta, x, y){
  ev = .evalLP(theta, x)
  lp  = ev$lp
  elp = exp(lp)

  status = y[, 2]
  ### find cost as negative of log likelihood function
  J = -sum(status*(lp-log(cumsum(elp))))
  return(J)
}

############ To be Validated for V ########################
reluGradient = function(theta, x, y, info = FALSE){
  ev = .evalLP(theta, x)
  lp = ev$lp
  elp = exp(lp)
  status = y[, 2]

  elp = exp(lp)

  Xa = cbind(ev$X, ev$eta)
  s0 = cumsum(elp)
  s1 = apply(Xa*elp, 2, cumsum)
  U = Xa - s1/s0  ### U is a n*p matrix, each row ui
  ### grd was check to be correct
  grd = -colSums(status * U)

  n = nrow(Xa)
  p = ncol(Xa)
  z2 = cbind(rep(1,length(ev$X[, 1])), ev$X[, 1])
  I = ifelse(ev$w >= ev$cx, 1, 0)
  z2_I = z2*I

  if(info) {
    mtm = function(m) {return(m%*%t(m))}
    ### add code for information matrix here
    ### these code shall be check to make sure they are correct
    UU = apply(U, 1, mtm) ### p*p*n each with u_i*t(u_i)
    U2 = aperm(array(UU, c(p, p, n)), c(3, 1, 2))  ### now n*p*p
    pi_theta = colSums(status * U2, dims = 1)

    ### Find 2nd order derivative,
    ### code is validated on Sept 18, 2018 with numerical method
    ## Sm is a lower triangular matrix of 1 because time is sorted by decending
    Sm = matrix(1, n, n)
    Sm[upper.tri(Sm)] = 0

    XX = array(apply(Xa, 1, mtm), c(p, p, n))

    # multiply each XX(p, p, i) with elp[i], by change elp to a p*p*n array
    # with each of i th pxp matrix = elp[i]
    X2 = XX * array(rep(elp, each = p*p), c(p, p, n))

    ### calculate s2, a n*p*p array, X2%*%Sm like cumsum over riskset
    s2 = apply(X2, c(1, 2), function(x, y){return(y%*%x)}, Sm)

    ### calculate s1_i*t(s1_i)
    # aperm changes a p1*p2*p3 array to a p3*p1*p2 array with c(3, 1, 2)
    s1_2 = aperm(array(apply(s1, 1, mtm), c(p, p, n)),c(3, 1, 2))

    I_theta = -colSums(status*(s2/c(s0)-s1_2/c(s0)^2), dims = 1)

    #for v_2rc
    s1rc = apply(z2_I*elp, 2, cumsum)
    V_2rc = status%*%(s1rc/(s0)-z2_I)

    #for v_2cc
    den=density(ev$w)

    pt1=which(den$x>=theta[p])[1]
    f.w.c = 0.5*(den$y[pt1]+den$y[pt1-1])

    s1cc = apply(z2%*%theta[(p-2):(p-1)]*elp, 2, cumsum)
    V_2cc = f.w.c*(status%*%(z2%*%theta[(p-2):(p-1)]-s1cc/(s0)))

    #file V.2
    V.2=matrix(0,nrow=p,ncol=p)
    V.2[p, p]= V_2cc
    V.2[p,p-2] = V_2rc[1]
    V.2[p,p-1] = V_2rc[2]
    V.2[p-2,p] = V_2rc[1]
    V.2[p-1,p] = V_2rc[2]

    second = I_theta + V.2

    upi = list(U = grd, pi_theta = pi_theta, I_theta = I_theta, second = second)
    return(upi)
  } else {
    return(grd)
  }
}

########## main fit functions ###########
reluProFit = function(x, y, seq_c) {
  p = length(x[1, ])
  ### p main effect, 1 interaction, 1 threshold
  theta0 = rep(0, p+2)
  c1 = seq_c

  nc = length(c1)
  log_c = matrix(NaN, nc, 2)
  log_c[, 1] = c1

  lgLik = -1e10
  for (j in 1:nc) {
    theta0[p+2] = c1[j]
    X = .evalLP(theta0, x)$X

    ### Make sure matrix X is not singular with too many no effective biomarker
    if(mean(X[, p+1]>0) < 0.05) next
    fit = coxph(y~X)
    log_c[j, 2]= fit$loglik[2]

    if (log_c[j, 2] > lgLik) {
      beta = fit$coefficients
      lgLik = log_c[j, 2]
      c.max = c1[j]
    }
  }
  theta = c(beta, c.max)
  log_c = as.matrix(log_c[!is.nan(log_c[, 2]), ])

  if (ncol(log_c) == 1)
    log_c=t(log_c)
  else log_c = log_c

  x = x
  y = y
  status = y[, 2]

  int_names = paste(colnames(x)[1], ":", colnames(x)[p], sep="")
  var_names = c(colnames(x)[1:p], int_names, "c.max")

  return(list(beta = beta, c.max = c.max, theta = theta, logLik = lgLik,
              x = x, y = y, status = status, log_c = log_c, var_names = var_names))
}

reluGradFit = function(x, y, theta0) {
  lmax = nlminb(theta0, reluLglik, reluGradient, x = x, y = y,
                lower = c(rep(-10, 4)), upper = c(rep(10, 4)))
  theta = lmax$par
  #print(as.vector(theta))
  p = length(theta)
  beta = theta[1:(p-1)]
  #print(beta)
  c.max = theta[p]
  logLik = -lmax$objective
  x = x
  y = y

  m = length(x[1, ])
  status = y[, 2]
  int_names = paste(colnames(x)[1], ":", colnames(x)[m], sep="")
  var_names = c(colnames(x)[1:m], int_names, "c.max")
  return(list(beta = beta, c.max = c.max, theta = theta, logLik = logLik,
              x = x, y = y, status = status, var_names = var_names))
}

print.brcm = function(x, ...) {
  cat("brcm: Cox PH model with ReLU\n")
  cat("Call:\n")
  print(x$call)
  cat("\n")

  cat('Coefficients and threshold parameter\n')
  bc = cbind(unname(t(x$beta)), x$c.max)
  colnames(bc) = t(x$var_names)
  print(bc, digits=4)
}

summary.brcm = function(object,...){
  zscore = qnorm(1-.05/2)

  sd = object$sd.score
  theta = t(cbind(unname(t(object$beta)), object$c.max))

  length = length(theta)
  p = length(object$object[1, ])
  w = object$object[, p]
  qtl = cbind(theta-zscore*sd, theta+zscore*sd)

  TAB1 = cbind(theta, sd, theta/sd, qtl[,1], qtl[,2])
  p_value = 2*pnorm(-abs(theta[1:length-1]/sd[1:length-1]))
  P_c = 2*pnorm(-abs((theta[length]-min(w))/sd[length]))
  p_values = c(p_value,P_c)
  TAB2 = cbind(TAB1, p_values)
  colnames(TAB2) = c("Estimate", "Std.Err", "Z value", "95% CI(lower)", "95% CI(upper)", "Pr(>z)")
  rownames(TAB2) = object$var_names

  results = list(call=object$call,TAB2=TAB2)
  class(results) = "summary.brcm"
  return(results)
}

print.summary.brcm = function(x,...){
  cat("Call:\n")
  print(x$call)
  cat("\n")

  cat('Summary table of coefficients and threshold parameter\n')
  printCoefmat(x$TAB2, digits=4)
}

#plot.brcm = function(x, ...){
  #method = x$method
  #if (method == "profile"){
    #c1 = x$log_c[, 1]
    #log_c = x$log_c[, 2]
    #plot(c1, log_c, type="o", main="Loglikelihood vs c by Profile", xlab="Threshold Parameter c",
         #ylab="Loglikelihood")
    #abline(v=c1[which.max(log_c)], h=max(log_c), lty=2, lwd=2, col="blue")
    #points(c1[which.max(log_c)],max(log_c), col='red', pch=19)
  #}
  #else {
    #length = length(x$logLik)
    #plot(1:length, x$logLik, main="Loglikelihood vs iterations by Gradient",
         #xlab="Number of iterations", ylab="Loglikelihood")
  #}
#}

print.residuals.brcm = function(x, type, ...) {
  ### code for martingle residuals
  p = length(x$x[1, ])
  p1 = p-1
  x1 = x$x[, 1]
  z = x$x[, 1:p1]
  w = x$x[, p]
  cx = x$c.max
  beta = x$beta

  phi = ifelse(w >= cx, w-cx, 0)

  X = cbind(z, phi, phi*x1)
  s0 = cumsum(exp(X%*%beta))
  H0 = rev(cumsum(rev(x$status/s0)))                  #cumulative baseline hazard
  #H0 = sum(x$status/s0)
  H = H0*exp(X%*%beta)    #cumulative          hazard
  r = x$status-H          #martingale residuals

  if (type == "Martingale"){
    results = list(r=r, w=w)
    class(results) = "residuals.brcm"
    return(results)
  }
  else {
    sv = survfit(Surv(H, x$status)~1)
    #plot(sv$time, -log(sv$surv), ylab="H.hat(ri)", xlab="ri")
    t=sv$time
    H_e = -log(sv$surv)
    results = list(t=t, H_e=H_e)
    class(results) = "residuals.brcm"
    return(results)
  }
}

plot.brcm = function(x, type, ...) {
  if (type == "Martingale"){
    plot(x$w, x$r, xlab="biomarker", ylab="Martingale residuals")
    lw = lowess(x$r ~ x$w)
    lines(lw, lwd=3, col="red")
  }
  else {
    plot(x$t, x$H_e, xlab="Cox-Snell Residuals", ylab="Estimated cumulative hazards")
    lw = lowess(x$H_e ~ x$t)
    lines(lw, lwd=3, col="red")
  }
}


predict.brcm = function(x, ...){
  lp = x$lp
  return(lp)
}

###### Generate survival data ################
surv.gendat =  function(n, c0, beta){
  ### biomarkers
  x = rnorm(n, 0.2, 2)
  ### code below to generate data with ReLU
  x1 = ifelse(x >= c0, x-c0, 0)
  z = rbinom(n,1,0.5)  ####used binomial to determine 0 or 1####
  zx = z*x1
  X = cbind(z, x1, zx)
  h0 = .5
  h = h0*exp(X%*%beta)

  stime = rexp(n, h)               #Failure time.
  endstudy = runif(n, 0, 5)
  cstatus = ifelse(stime>endstudy, 0, 1) ##Censored at end of study time.
  #cat('\nCensoring: ', 1-mean(cstatus), '\n')
  stime = ifelse(stime>endstudy, endstudy, stime)
  dat = cbind(time=stime, status=cstatus, z=z, x=x)
  dat = data.frame(dat)
  return(dat)
}
