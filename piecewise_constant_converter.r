### Make grid of parameters (in log space for scale and probably also for shape)
require(dplyr); require(FAdist); require(compiler)


### Compute LLs on grid with one pass through data (batch for now, or batch updates "quasi-stream")
### Select scale and shape
######

### Convert log logistic multiplier into piecewise-constant model
# question, do you ever need to query conditional CDF?

### define unit time
### get quantile times
### select nearest unit times
### approximate hazard (piecewise constant (poisson), piecewise linear? (cox?))

# TODO ret_subtimeline_ is not using offset appropriately from event times data
# returns a piecewise constant approximated disribution of p (or equiv q) with parameters lb, ub, unit, offset, qtile
piecewise_constant_subtimeline = (function(p, q, ...,
                                          t0=0, lb=0, ub=Inf,
                                          unit=0.01, offset = 0.000, qtile = 0.01, 
                                          givenpoints=NULL, as.hazard=F, inf.default = 0, notinf.default = 0) {
  dots = list(...)
  if(is.null(givenpoints)) {
    # if(is.infinite(t0)) t0 = 0
    if(is.infinite(t0)) return(data.frame(lb=lb,ub=ub,logRateFactor=inf.default))
    if(lb < t0) stop("lb must be >= t0")
    if(lb >= ub) stop("lb must be < ub")
    if(is.na(offset)) {
      if(lb + qtile - (lb %% qtile) < ub) offset = qtile - (lb %% qtile)
      else offset = 0
    }
    # get exact qtiles from rounded (grid) positions
    qstarttile = do.call("p", list(lb-t0) %>% c(dots))
    qendtile = do.call("p", list(ub-t0) %>% c(dots))
    if(qtile > qendtile) qtile = qendtile
    distub = do.call("q", list(seq(qstarttile,qendtile,qtile)) %>% c(dots)) %>%
      (function(x) { ( ( (( x - offset ) / unit) %>% round() ) * unit + offset ) } ) %>%
      (function(x) { if (all.equal(x[length(x)],ub-t0) %>% isTRUE()) x else c(x, ub-t0) } )
    distub = distub[distub>lb-t0] # = max(distub[distub<=lb],lb)
  } else {
    distub = givenpoints
    qstarttile = 0
  } # TODO there is an unchecked case where distub is empty with positive ub-lb, lb==distub (which is silently passing through and giving null vector (incorrect I believe))
  
  # qtdf = table(distub) %>% data.frame() %>% tbl_df() %>%
  #   mutate(distub=(as.numeric(levels(distub))[distub])) %>%
  #   arrange(distub) %>%
  #   select(-Freq) %>%
  qtdf = data.frame(distub=sort(unique(distub))) %>% tbl_df() %>%
    mutate(qtiles = do.call("p", list(distub) %>% c(dots)),
           ub = distub+t0) %>%
    mutate(lb = lag(get("ub"), default = lb),
           qtileprev = lag(get("qtiles"), default = qstarttile),
           logRate = log( -log(1-(qtiles-qtileprev))/(ub-lb) ) + notinf.default
           ) %>%
    (function(x) if(as.hazard) x %>% mutate(logRate = logRate - log(1-qtileprev)) else x)(.) %>%
     select(lb, ub, logRate) %>%
    rename(logRateFactor=logRate) %>% (function(x) {if(nrow(x)>0 && x[1,"ub"] == x[1,"lb"]) x[-1,] else x})
   #qtdf$rate[dim(qtdf)[1]] = qtdf$rate[dim(qtdf)[1]-1]
  return(qtdf)
}) %>% cmpfun()

# test
plot_pcs = function(pcs) {
  require(ggplot2)
  ggplot(pcs, aes(x=ub, y=logRateFactor)) +
    geom_segment(aes(y=logRateFactor, yend=logRateFactor,x=lb, xend=ub)) + xlab("Time") + theme_get()
}

demo = F
if(demo) {
piecewise_constant_subtimeline(p=pllog, q=qllog, scale=log(10), shape=0.3) %>% plot_pcs()


# written out
unit = 0.01
offset = 0.001 # will change to align t at evidence onto unit time grid
qtile = 0.01
pars = list()
pars$scale = log(10)
pars$shape = 0.1
qtimes = qllog(seq(qtile,1,qtile), scale = pars$scale, shape=pars$shape)
qtimes = round((qtimes-offset)/unit)*unit+offset
qtimes[qtimes<=0] = max(qtimes[qtimes<=0],offset) # in case binding is to 
# get exact qtiles from rounded positions
qtdf = table(qtimes) %>% data.frame() %>% tbl_df()
qtdf$qtimes = with(qtdf, as.numeric(levels(qtimes))[qtimes])
qtdf = qtdf %>%
  arrange(qtimes) %>%
  select(-Freq) %>%
  mutate(qtiles = pllog(qtimes, scale= pars$scale, shape=pars$shape)) %>%
  mutate(qprev = lag(qtimes, default = 0)) %>%
  mutate(qtileprev = lag(qtiles, default = 0)) %>%
  mutate(rate = -log(1-(qtiles-qtileprev))/(qtimes-qprev))
#qtdf$rate[dim(qtdf)[1]] = -log(qtdf$qtiles[dim(qtdf)[1]]-qtdf$qtileprev[dim(qtdf)[1]])/qtdf$qprev[dim(qtdf)[1]]
qtdf$rate[dim(qtdf)[1]] = qtdf$rate[dim(qtdf)[1]-1]
qtdf = qtdf %>% 
  select(qprev,rate) %>%
  rename(qtimes=qprev)
#qtdf$rate[1] = qtdf$rate[1] + 0 #weird dexp behavior NaN from 0
}

### Piecewise-exponential density and distribution functions
# while it looks like msm dpexp could work, that density is based on the rates accounting for survival,
# whereas the rates computed here are unadjusted (they form a valid prob dist by concatenation of CDFs)
dpexp = function(x, rates, times) {
  rates[rowSums(outer(x,times,">="))]
}
ppexp = function(q, rates, times) {
  if(times[1]!=0) stop("first time must be 0")
  ind = rowSums(outer(q,times,">="))
  cdfs = c(cumsum(1-exp(-rates[-length(rates)]*(times[-1]-times[-length(times)]))),1)
  qcdfs = 1-exp(-rates[ind]*(q - times[ind]))
  qcdfs[cdfs[ind]==1] = qcdfs[cdfs[ind]==1]*(1-cdfs[length(cdfs)-1])
  qcdfs + c(0,cdfs)[ind]
}


if(demo) {
### demo
xs = seq(0,50,0.1)
ds = dllog(xs,scale=pars$scale, shape=pars$shape) # density
hs = ds / (1-pllog(xs,scale=pars$scale, shape=pars$shape)) # hazard
ps = pllog(xs,scale=pars$scale, shape=pars$shape) # hazard
plot(x=xs,y=ds, ylim=c(0, max(ds,na.rm=T)))
segments(x0=qtdf$qtimes,x1=c(qtdf$qtimes[-1],max(xs)),y0=qtdf$rate, y1=qtdf$rate, col=2, pch=18, lwd=4)

plot(x=xs,y=hs, ylim=c(0, max(hs,na.rm=T)))
points(x=xs,
       y=dpexp(x=xs,rates=qtdf$rate,times=qtdf$qtimes)/(1-ppexp(q=xs,rates=qtdf$rate,times = qtdf$qtimes)),
       pch=18,cex=0.5,
       col=2, type='l')

plot(x=xs,y=ps, ylim=c(0,max(ps,na.rm=T)))
points(x=xs,
       y=ppexp(q=xs,rates=qtdf$rate,times = qtdf$qtimes),
       pch=18,cex=0.5,
       col=2, type='l')
}