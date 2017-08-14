### generate plots illustrating concept of log logistic entry into piecewise constant hazard

source("piecewise_constant_converter.r")
require(dplyr)
require(ggplot2)
require(FAdist)
df = data.frame(lb = c(0,1,3,8), ub=c(1,3,8,10), logRateFactor=c(log(0.05),log(0.25),log(0.5),log(0.05))) %>% tbl_df()

get_time_to_event = function(lbubhazard, sequence) {
  expandedLbubHazard = insert_events_into_expansion(lbubhazard, sequence) %>% select(-M) %>% rowwise() %>%
    mutate(logRateFactor = lbubhazard[[findInterval(lb, lbubhazard$lb),"logRateFactor"]]) %>% ungroup()
  expandedLbubHazard %>%
    mutate(survpart = exp(-exp(logRateFactor)*(ub-lb))) %>%
    mutate(surv = lag(cumprod(survpart),default=1)) %>%
    mutate(tte = logRateFactor + log(surv)) %>% select(-survpart)
  
}

sequence = seq(0,10,0.01)
llshape = 0.25; llscale = log(6)
heightGamma = 2
#llsequence = dllog(sequence, shape=llshape, scale = llscale, log=T)

graphdf = get_time_to_event(df, sequence) %>% 
  mutate(llRateFactor = dllog(ub, shape=llshape, scale = llscale, log=T)-log(1-pllog(ub, shape=llshape, scale = llscale))) %>%
  mutate(lltte = dllog(ub, shape=llshape, scale = llscale, log=T)) %>%
  mutate(combinedRateFactor = logRateFactor+llRateFactor+heightGamma)
combinedtte = graphdf %>% select(lb,ub, logRateFactor=combinedRateFactor) %>% get_time_to_event(sequence)
graphdf = graphdf %>% bind_cols(combinedtte = data.frame(combinedtte=combinedtte$tte))

cbbPalette <- c("#000000", "#56B4E9", "#009E73", "#E69F00", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

ggpc = ggplot(data = graphdf, aes(x=ub)) + 
  geom_line(aes(y=exp(tte), color="Time to event")) + 
  geom_line(aes(y=exp(logRateFactor), color="Hazard")) +
  theme_classic() + theme(axis.title.x = element_blank()) + ylab("Density; hazard") +
  scale_colour_manual(values=cbbPalette, name="Piecewise")  

ggll = ggplot(data = graphdf, aes(x=ub)) + 
  geom_line(aes(y=exp(lltte), color="Time to event")) +
  geom_line(aes(y=exp(llRateFactor), color="Hazard")) +
  theme_classic() + theme(axis.title.x = element_blank()) + ylab("Density; hazard") +
  scale_colour_manual(values=cbbPalette, name="Log logistic")  

ggco = ggplot(data = graphdf, aes(x=ub)) + 
  #geom_line(aes(y=exp(combinedRateFactor), color="Combined hazard")) +
  geom_line(aes(y=exp(combinedtte), color="Time to event")) +
  theme_classic() + ylab("Density") +
  xlab("Time") + scale_colour_manual(values=cbbPalette[-1], name="Combined model")  

llap = piecewise_constant_subtimeline(p=pllog,
                               q=qllog,
                               scale=llscale,
                               shape=llshape,
                               t0=0,
                               lb=0,
                               ub=10,
                               unit=0.01,
                               offset=0,
                               as.hazard = T)
graphdf = graphdf %>% bind_cols(data.frame(combinedApproximateRateFactor = product_of_timelines(llap, graphdf %>% select(lb,ub,logRateFactor))$logRateFactor + heightGamma))

ggap = ggplot(data = graphdf, aes(x=ub)) + 
  geom_line(aes(y=exp(combinedRateFactor), color="Hazard")) +
  geom_segment(aes(y=exp(combinedApproximateRateFactor), yend=exp(combinedApproximateRateFactor), x=lb,xend=ub, color="Approximate hazard")) +
  ylab("Hazard") + theme_classic() + 
  scale_colour_manual(values=cbbPalette[c(4,1)], name="Combined model") + theme(axis.title.x = element_blank())
ggplot_build(ggap)

library(grid)
pdf(file = "simulations/conceptGraph.pdf",width = 5,height=7)
grid.newpage()
grid.draw(rbind(ggplotGrob(ggpc), ggplotGrob(ggll), ggplotGrob(ggap), ggplotGrob(ggco), size = "last"))
dev.off()
