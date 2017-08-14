### Partitioned MLE for piece-wise constant log logistic addition

### This code illustrates the identification of the log logistic MLE estimate for censored data 
### It functions for log-logistic right-censored. Not left nor interval censored.
### All we need for the algorithm is conditional log logistic with right censoring.

source("piecewise_constant_converter.r")

require(FAdist)
require(dplyr)
require(compiler)

# likelihood function (parameters, data) {}
#   data = times to events
#   parameters = partition scale multipliers, scale and shape
#
#   input representation:
#     [partition multiplier], [(tte, event?, partition index)]  ***create fn to transform into latter***
#         or
#     [(tte, event?, partition multiplier)]  ***this one***

#   given tbl_df of (tte, event?, partition multiplier), do:
#   
pcllll = (function(df) {
  df %>% mutate(ll=llll(c(1,1),df))
}) %>% cmpfun()

# note dllog uses scale ~ log(median), so add multiplier in log space #
# parameterization: 
#   shape: mu = log(alpha)
#   scale: s = 1/beta
llll = (function(pars, df) { # ?? don't you want to sum over unique likelihood-per-row?
  -sum(with(df, (1-event) * llsurvivell(scale=pars[1]+logmultiplier, shape=pars[2], data=tte) +
           event * lleventll(scale=pars[1]+logmultiplier, shape=pars[2], data=tte) ))
}) %>% cmpfun()

llsurvivell = (function(scale, shape, data) {  # negative ll of survival #
  return(log( 1/(1+(data/exp(scale))^(1/shape)) ) )
}) %>% cmpfun()

lleventll = (function(scale, shape, data) {  # ll of density #
  return(log( dllog(data, shape = shape, scale = scale) ))
}) %>% cmpfun()

bound_range = (function(x, pcmax, pcmin, check=F) {
  if(check) {
    if(any(x!=pmin(pmax(x, pcmin),pcmax)))
      browser()
  }
  pmin(pmax(x, pcmin),pcmax)
}) %>% cmpfun()

### Now develop conditional survival:
### - (lb, ub, event, logmultiplier, lastTime)
### - pc must be ordered (within groups) and non-overlapping; df can be; need function to map df to pc (if no pc specified)
### - df and pc must specify groupId; df must specify subgroupId -- else assumes all different
llcll = (function(pars, df, pc, groupId=c(), pcmax = 1e7, pcmin = -1e3) {
  if(identical(df,pc))
    dfpc = T
  if(length(pars)==3)
    pc$logmultiplier = pars[3] + pc$logmultiplier # OLD min(pars[3], min) limit multiplicative increase ability NEW use constrOptim
  if(pars[2]<1e-4) pars[2] = 1e-4

  #pc$logmultiplier = rep(0,nrow(pc)) #probably can clean up a lot if we don't need logmultiplier in here # but we do!
  # groups: df durations can overlap but pc durations cannot
  
  #empiricaldist = cumsum(exp(c(pc$logmultiplier,-10) + pllog(c(pc$lb,Inf), scale=pars[1], shape=pars[2], log.p = T))) %>%
  #  (function(x) x/max(x))(.) # Z vector
  # browser()
  # empiricaldist = pc %>% 
  #   group_by_(.dots=groupId %>% lapply(as.symbol)) %>% 
  #   mutate(empiricaldist= 
  #     cumsum(exp(c(logmultiplier,-100))*
  #              (pllog(c(ub,Inf), scale=pars[1],shape=pars[2]) -
  #                 pllog(c(lb,ub[length(logmultiplier)]), scale=pars[1],shape=pars[2])
  #              )
  #            ) %>% 
  #       #(function(x) x/max(x))(.) %>% 
  #       .[-(length(logmultiplier)+1)],
  #     lbemp = lag(empiricaldist, default=0)) # %>% ungroup() %>% .[["lbemp"]]

  if(dfpc) {
    # empiricaldist = empiricaldist %>% ungroup() %>%
    ints = 1:nrow(pc) #1:nrow(.)
    ints2 =1:nrow(pc) #1:nrow(.)
  } else if("lastTime" %in% groupId) { #TODO generalize to some given "time" variable
    pcints = df %>% group_by_(.dots=groupId %>% lapply(as.symbol)) %>% #TODO this is wrong for groups beyond "lastTime"
      do(intsub = findInterval(.$ub, pc$lb[pc$lastTime==.$lastTime[1]],left.open = T, rightmost.closed = T) +
           min(which(pc$lastTime==.$lastTime[1]))-1, 
         intslb = findInterval(.$lb, pc$lb[pc$lastTime==.$lastTime[1]],left.open = T, rightmost.closed = T) +
           min(which(pc$lastTime==.$lastTime[1]))-1) %>% ungroup() %>%
      do(ub = .$intsub %>% unlist(),
         lb = .$intslb %>% unlist())
    ints = pcints$ub[[1]]
    ints2 = pcints$lb[[1]]
  } else { # for optimization that does not have a lastTime and thus no empiricaldist shift
    ints = findInterval(df$ub, pc$lb, left.open = T, rightmost.closed = T)
    ints2 = findInterval(df$lb, pc$lb,left.open = T, rightmost.closed = T)
  }
  
  # if(any(empiricaldist$lbemp>1)) {
  #   browser()
  #   return(Inf)
  # }
  
  csurv = (1-pllog(df$ub, scale=pars[1],shape=pars[2]))/
                (1-pllog(pc$lb[ints], scale=pars[1],shape=pars[2]))
  if(any(is.nan(csurv)) || any(csurv<=0)) {
    return(pcmax)
    browser()
  }
  
  eventll = 
    # df$event * (pc$logmultiplier[ints] +
    # # df$event * (pc$logmultiplier[findInterval(df$ub,pc$lb)] +
    #               lleventll(scale=pars[1], shape=pars[2], data=df$ub))  # ll corresponding to unnormalized density
    # hazard components:
    #   1. pc$logmultiplier[ints], fixed, not needed for optimization (but remember to in likelihood calc.)
    #   2. lleventll(scale=pars[1], shape=pars[2], data=df$ub) * Survival(t)
      (pc$logmultiplier[ints] + lleventll(scale=pars[1], shape=pars[2], data=df$ub) - log(csurv)
            )[as.logical(df$event)] %>% bound_range(pcmax,pcmin)
    
  # TODO these "conditional survivals" are not technically correct because they might not integrate to 1; they integrate to 1-Surv(t_0);
  # where t_0 is the first active time and might not be the 0 time (lastTime) of the distribution.
  # also, they are not normalized--works often because f(t)*g(t) is generally << 1 (not for coinciding pointy distributions though)  
  # trouble on the horizon here?
  
  #which pc$ub am I above?
  #F(nearest lower ub) + F(residual above nearest lower ub)
  #ints = findInterval(df$ub, pc$lb, left.open = T, rightmost.closed = T) # should be within group
  # ints = pcints$ub[[1]]
  # survs = 1 - (empiricaldist$lbemp[ints] +
  #                exp(pc$logmultiplier[ints])*
  #                (pllog(df$ub, scale=pars[1],shape=pars[2])-pllog(pc$lb[ints], scale=pars[1],shape=pars[2]))
  #              )
  
  # if(any((1-df$event)*survs<0)) {
  # # if(any(survs<0)) {  # f(t) = h(t)*S(t) modification, count the survival no matter what
  #   # browser() # TODO what to do about par[3] ??
  #   # Proposal: bound at some large positive number (to discourage choosing parameter there).
  #   #           other data can override this constant bound (multiple negative likelihood components)
  #   #           which means you also should provide a bound on negative likelihood components
  #   #           (bounding shape is a start but insufficient: log(shape f)+log(shape f)+... is still unbounded)
  #   #           this will have the effect of not choosing to further optimize already high likelihood regions, which in most
  #   #           cases seem desirable.
  #   #           this procedurally should work, but on unseen test data there is always the possibility that the empiricaldist > 1
  #   #           suffice to throw a warning and attribute a constant bound low likelihood in this case.
  #   return(pcmax)
  #   return(Inf)
  # }
    #ints2 = findInterval(df$lb, pc$lb,left.open = T, rightmost.closed = T)
    # ints2 = pcints$lb[[1]]
  # survs2 = 1-empiricaldist$lbemp[ints2]
    # noeventll = (1-df$event)*(log(survs) - log(survs2)) # vectorized --> supports NaNs which throw errors
  
  # noeventll = pmax((log(survs[as.logical(1-df$event)]) - log(survs2[as.logical(1-df$event)])),pcmin)
  
  noeventll = (exp(pc$logmultiplier[ints])*log(csurv)) %>% bound_range(pcmax,pcmin)
  # count survival time regardless of event: interval, with y={0,1}: L(t) = f(t) = h(t)^y*S(t)
  
  # noeventll = pmax((log(survs) - log(survs2)),pcmin)
  
  
  # noeventll = 
  #   (1-df$event) * (
  #     (function(e,du,pl,mult)
  #     {ints = findInterval(du, pl,left.open = T, rightmost.closed = T); log(1- 
  #        (e[ints] + exp(mult[ints])*(pllog(du, scale=pars[1],shape=pars[2])-pllog(pl[ints], scale=pars[1],shape=pars[2])))
  #        )
  #     })(empiricaldist, df$ub, pc$lb, pc$logmultiplier) -
  #     
  #     (function(e,dl,pl,mult)
  #       {ints = findInterval(dl, pl,left.open = T, rightmost.closed = T); log(1- e[ints] + 0)}
  #       )(empiricaldist, df$lb, pc$lb, pc$logmultiplier)
  #   )  # ll corresponding to unnormalized survival (but survival is of empirical dist), approximate
  if(is.nan(sum(eventll) + sum(noeventll))) {
    browser()
  }
  # browser()
  
  -(sum(eventll) + sum(noeventll))
  
  #OLD 2 (1-df$event) * (log(1-empiricaldist[findInterval(df$ub, pc$lb)] + 0) -
                    #OLD 2 log(1-empiricaldist[findInterval(df$lb,pc$lb)] + 0 ))
  # OLD (1-df$event) * (weightedcomponents + llcsurvivell(scale=pars[1], shape=pars[2], data=df))
}) %>% cmpfun()

# will not work for noncontinuous durations?
make_pc_from_blank_df = (function(df) {
  pc = df %>% arrange(lb) %>% (function(x) fix_timelineCPP(x$lb, x$ub))(.)
  if(!("logmultiplier" %in% df))
    pc$logmultiplier = 0
  else
    pc$logmultiplier = df %>% mutate_active_duration_indices(pc) %>%
      value_apply_onto_durations(pc, "logmultiplier")
  pc
}) %>% cmpfun()

apply_constraints_to_mlps = (function(pars, maxIncrease=10, minShape=0.01) {
  if(length(pars)>=2 && pars[2] < 0.01) pars[2] = 0.01
  if(length(pars)>=3 && pars[3] > 10) pars[3] = 10
  pars
}) %>% cmpfun()

# survival at ub given survival at lb = survival at ub / survival at lb
llcsurvivell = (function(scale, shape, data) {
  return(
    log( 1/(1+(data$ub/exp(scale))^(1/shape)) ) -
      log( 1/(1+(data$lb/exp(scale))^(1/shape)) )
  )
}) %>% cmpfun()

llceventll = (function(scale, shape, data) {
  return(
    log( dllog(data$ub, shape = shape, scale = scale) ) - 
      log( 1/(1+(data$lb/exp(scale))^(1/shape)) )
    )
}) %>% cmpfun()

demo = F
if(demo) {
### demo
# illustration of singular point MLE 
# dat = data.frame(matrix(c(c(1,1,1),c(2,0,5)),nrow = 2, byrow = T)) %>% tbl_df()
# names(dat) = c("tte", "event", "logmultiplier")
# dat[3:100,] = dat[1,] + runif(300)*0.01
# dat[,2] = c(1,0,rep(1,98))
# dat[,3] = rep(0,100)
# dat[5:7,1] = 2
# dat[5:7,3] = log(2)

# opt = optim(c(1,1),llll, df=dat)
# plot(x=seq(0,2,0.002),y=dllog(seq(0,2,0.002), scale = opt$par[1], shape=opt$par[2]), pch=18)

# illustration of non-event survival balanced with singular point mle 
# dat2 = dat
# dat2[50:100,"tte"] = dat[50:100,"tte"]+1
# opt2 = optim(c(1,1),llll, df=dat2)
# plot(x=seq(0,2,0.002),y=dllog(seq(0,2,0.002), scale = opt2$par[1], shape=opt2$par[2]), pch=18)
# 
# # illustration of differing start points in survival function on MLE
# dat3 = dat %>% mutate(lb = 0, ub = tte) %>% select(-tte)
# dat3[50:nrow(dat3),c("lb","ub")] = dat3[50:nrow(dat3),c("lb","ub")]+1
# opt3 = optim(c(1,1),llcll, df=dat3)
# plot(x=seq(0,2,0.002),y=dllog(seq(0,2,0.002), scale = opt3$par[1], shape=opt3$par[2]), pch=18)
# points(x=seq(0,2,0.002),y=dllog(seq(0,2,0.002), scale = opt$par[1], shape=opt$par[2]), col=2, pch=18)
# points(x=seq(0,2,0.002),y=dllog(seq(0,2,0.002), scale = opt2$par[1], shape=opt2$par[2]), col=3, pch=18)

### debug: FIXED ensure you understand shape and scale parameterization so that multiplier integrates correctly

prior4 = piecewise_constant_subtimeline(pllog, qllog, scale = 0, shape = 0.6, ub=10, qtile = 0.001) %>%
  mutate(event = 0) %>% rename(logmultiplier = logRateFactor)
lack = 20
dat4 = data.frame(
  event = c(1,1,1,1,0,0,rep(0,lack)),
  lb = c(0,0,0,1,0,0,rep(0,lack)),
  ub = c(2,1,2.1,2,4,1,rep(3,lack))
) %>% tbl_df() %>%
  mutate(logmultiplier = prior4$logmultiplier[findInterval(ub,prior4$ub)])
opt4 = optim(c(1,1),llcll, df=dat4, pc=prior4)

plot(x=seq(0,10,0.005),y=dllog(seq(0,10,0.005), scale = opt4$par[1], shape=opt4$par[2]), pch=18, cex=0.2,
     )#log="y")
segments(x0 = prior4$lb, x1=prior4$ub,y0=exp(prior4$logmultiplier), col=2)

k= 0.5
dat4mixed = piecewise_constant_subtimeline(pllog, qllog, scale = opt4$par[1], shape = opt4$par[2], givenpoints = c(0,prior4$ub)) %>%
  rename(logmultiplier = logRateFactor) %>% mutate(combined = k*logmultiplier+prior4$logmultiplier)
with(dat4mixed, segments(x0=lb,x1=ub, y0=exp(combined)*exp(max(dat4mixed$logmultiplier)-max(dat4mixed$combined)), y1=exp(combined)*exp(max(dat4mixed$logmultiplier)-max(dat4mixed$combined)), col="blue"))

# log concavity problem for shape<1
# plot(y=dllog(seq(0,5,0.01),log = T, scale=0,shape=0.005)+dllog(seq(0,5,0.01),log = T, scale=log(1.5),shape=1.05)+dllog(seq(0,5,0.01),log = T, scale=log(3),shape=0.005), x=seq(0,5,0.01), pch=18)

# can optimize over arbitrary prior?
prior5 = piecewise_constant_subtimeline(pllog, qllog, scale = 0, shape = 0.6, ub=10, qtile = 0.001) %>%
  mutate(event = 0) %>% rename(logmultiplier = logRateFactor) %>% mutate(logmultiplier=logmultiplier*(1+sin((1:nrow(.))/40*pi))/2)
lack = 100
dat5 = data.frame(
  event = c(1,1,1,1,0,0,rep(0,lack)),
  lb = c(0,0,0,1,0,0,rep(0,lack)),
  ub = c(2,1,2.1,2,4,1,rep(3,lack))
) %>% tbl_df() %>%
  mutate(logmultiplier = prior5$logmultiplier[findInterval(ub,prior5$ub)])
opt5 = optim(c(1,1),llcll, df=dat5, pc=prior5) %>%
  (function(x) {x$par = apply_constraints_to_mlps(x$par); x})(.)

plot(x=seq(0,10,0.005),y=dllog(seq(0,10,0.005), scale = opt5$par[1], shape=opt5$par[2]), pch=18, cex=0.8,
     )#log="y")
points(x=seq(0,10,0.005),y=dllog(seq(0,10,0.005), scale = opt4$par[1], shape=opt4$par[2]), pch=18, cex=0.2, col="grey") #dat4 and dat5 same but opt different because of different logmultipliers
segments(x0 = prior5$lb, x1=prior5$ub,y0=exp(prior5$logmultiplier), col=2)

dat5mixed = piecewise_constant_subtimeline(pllog, qllog, scale = opt5$par[1], shape = opt5$par[2], givenpoints = c(0,prior5$ub)) %>%
  rename(logmultiplier = logRateFactor) %>% mutate(combined = logmultiplier+prior5$logmultiplier)
with(dat5mixed, segments(x0=lb,x1=ub, y0=exp(combined)*exp(max(dat5mixed$logmultiplier)-max(dat5mixed$combined)), y1=exp(combined)*exp(max(dat5mixed$logmultiplier)-max(dat5mixed$combined)), col="blue"))

### Demo moving from exponential prior for period data; cut or not cut
prior6 = piecewise_constant_subtimeline(pexp, qexp, rate=1/3, ub=10, qtile = 0.001) %>%
  mutate(event = 0) %>% rename(logmultiplier = logRateFactor)
dat6 = data.frame(
  event = c(1,1,1,0),
  lb = c(0,0,0,0),
  ub = c(3,3,3,100)
)
opt6 = constrOptim(theta=c(1,1,-1),
                   f=llcll,
                   grad=NULL,
                   ui=rbind( # shape must be > 0.01, maxIncrease must be less than 5
                     c(0,1,0),
                     c(0,0,-1)),
                   ci=c(0.01,-5),
                   df=dat6,
                   pc=prior6)

plot(x=seq(0,10,0.005),y=dllog(seq(0,10,0.005), scale = opt6$par[1], shape=opt6$par[2]), pch=18, cex=0.8,
)#log="y")
segments(x0 = prior6$lb, x1=prior6$ub,y0=exp(prior6$logmultiplier), col=2)
  
# TODO need a function than can take the 3 pars and produce a tte distribution
dat6mixed = piecewise_constant_subtimeline(pllog, qllog, scale = opt6$par[1], shape = opt6$par[2], givenpoints = c(0,prior6$ub)) %>%
  rename(logmultiplier = logRateFactor) %>% mutate(combined = logmultiplier+prior6$logmultiplier)
with(dat6mixed, segments(x0=lb,x1=ub, y0=exp(combined)*exp(max(logmultiplier)-max(combined)), y1=exp(combined)*exp(max(logmultiplier)-max(combined)), col="blue"))

opt6b =  constrOptim(theta=c(1,1,-1),
                   f=llcll,
                   grad=NULL,
                   ui=rbind( # shape must be > 0.01, maxIncrease must be less than 5
                     c(0,1,0),
                     c(0,0,-1)),
                   ci=c(0.01,-5),
                   df=dat6,
                   pc=make_pc_from_blank_df(dat6))
points(x=seq(0,10,0.005),y=dllog(seq(0,10,0.005), scale = opt6b$par[1], shape=opt6b$par[2]), pch=18, cex=0.4, col=3
)#log="y")

# DOTO need to calibrate (because multiplicative density is unnormalized) --> could put great weight into exp tail!
#   added par[3] which is the multiplier across the distribution

### DOTO won't dat be exponential size with division of trajectory into spacings with different logmultipliers? --> use grid
###   yes I believe you do
###   consider batch updates; what happens is that you'll get suboptimal ML parameters;
###   could get specify piecewise approximation based on log logistic parameters, get ML approximation
###     and use to LSS estimate log logistic parameters as update?

###   how about radix approximation, i.e. you have a fixed number of logmultipliers e.g. 100k of them?
###   this appears more promising
###   yes you can just find largest LL on grid of parameters (won't be optimal, but probably not too bad?)

### DOTO this parameterization ignores likelihood as product survival densities: p(t|t>T) with p log-logistic --> llcll instead of llll
### edge cases will exist near observations 
}

