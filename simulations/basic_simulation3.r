source("forest_modifications.r")
source("learn_one_random_proposal.r")

set.seed(1298653)
pt_exp_timeline = function(lb=0, ub=110, rate = 0.1, overshoot=10) {
  data.frame(
    time = rexp(ceiling((ub-lb)*rate*overshoot), rate = 0.1) %>% cumsum() %>%
      (function(.) .[.+lb<ub]+lb)
  ) %>% tbl_df()
}

npts3 = 50;
dats3 = data.frame(id = 1:npts3 %>% as.character()) %>%
  tbl_df() %>%
  group_by(id) %>% 
  do(pts = data.frame(time = pt_exp_timeline()) %>% tbl_df() %>% mutate(event=rep("A",nrow(.))))
lbs = data.frame( lb = rep(0,npts3), ub = rep(110,npts3), id = 1:npts3 %>% as.character()) %>%
  group_by(id) %>%
  do(lbubs=data.frame(.) %>% select(-id))

dats3 = dats3 %>% left_join(by="id", lbs)

tree = data.frame(
  condition=character(),
  timeFrame=numeric(),
  trueNode=numeric(),
  falseNode=numeric(),
  distribution=numeric(),
  par0=character(),
  par1=numeric(),
  par2=numeric(),
  modifiable=numeric(),
  stringsAsFactors = F
) %>% tbl_df()

tree[dim(tree)[1]+1,] = c("" ,NA,NA,NA,1 ,NA,0 ,NA, 1)
finish_tree_initialization = function(tree) {
  tree$timeFrame = as.numeric(tree$timeFrame)
  tree$distribution = as.numeric(tree$distribution)
  tree$par0 = as.character(tree$par0)
  tree$par1 = as.numeric(tree$par1)
  tree$par2 = as.numeric(tree$par2)
  tree$trueNode = as.numeric(tree$trueNode)
  tree$falseNode = as.numeric(tree$falseNode)
  tree$modifiable = as.numeric(tree$modifiable)
  tree
}
tree = finish_tree_initialization(tree)
forest = list(tree=tree, roots = 1)

set.seed(44); output = learn_pc(forest, dats3, "A", iter=1,
                                pf = generate_one_proposal_po, 
                                maintainBlankTree = F,
                                pars=list(
                                  type="default",
                                  localityFactor = 1,
                                  maxLength=100,
                                  probSplit=0
                                )
)

set.seed(44); outputAForest = learn_pc(forest, dats3, "A", iter=10,
                                       pf = generate_one_proposal_po, 
                                       pars=list(
                                         type="default",
                                         localityFactor = 1,
                                         maxLength=100,
                                         probSplit=0.5
                                       )
) # works

set.seed(44); outputAForest = learn_pc(forest, dats3, "A", iter=10,
                                       pf = generate_one_proposal_po, 
                                       pars=list(
                                         type="default",
                                         localityFactor = 1,
                                         maxLength=100,
                                         probSplit=0
                                       )
) # DOTO crashes, now runs but there is (1) an offset issue and (2) replacedIndex is erroneously NA
# forest with both exp splits and ll mods

set.seed(44); outputBForest = learn_pc(forest, dats3, "A", iter=30,
                                       pf = generate_one_proposal_po, 
                                       pars=list(
                                         type="default",
                                         localityFactor = 1,
                                         maxLength=100,
                                         probSplit=0.5
                                       )
)

# Rprof(tmp <- tempfile(), line.profiling=T)
# output = learn_pc(forest,dats,"A",iter=2)
# Rprof()
# summaryRprof(tmp,lines = "show")

plot_pt_rate_graph = function(dats, model, index=1) {
  require(ggplot2)
  df = dats %>% 
    mutate(sts = list(subtimeline_forest_preexpansion(model[["forest"]][["tree"]],
                                                      pts,
                                                      model[["forest"]][["roots"]],
                                                      lbubs,
                                                      1))) %>% # DOTO shows an improper expansion!
    .[["sts"]] %>% 
    .[[index]] %>% 
    data.frame
  ggplot(df, aes(x=ub, y=logRateFactor)) +
    geom_point(data=dats$pts[[index]], aes(x=time,y=0, color=event)) +
    geom_segment(aes(y=logRateFactor, yend=logRateFactor,x=lb, xend=ub)) +
    xlab("Time")
}
pt_rate_table = function(dats, model, index = 1) {
  subtimeline_forest_preexpansion(model$forest$tree,
                                  dats$pts[[index]],
                                  model$forest$roots,
                                  dats$lbubs[[index]],
                                  1)
}

plot_pt_rate_graph(dats3,output,7)

plot_pt_rate_graph(dats3,outputAForest,2)
pt_rate_table(dats3, outputAForest,7)

ggout2 = plot_pt_rate_graph(dats3,outputBForest,2)
pdf(file="simulations/poisson.pdf", width=8, height=2)
ggplot_build(ggout2)
dev.off()