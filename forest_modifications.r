# remember to set the working directory with setwd(...)
source("subtimeline_maker.r")
source("likelihood_functions.r")
source("partitioned_pcll.r")

library(tidyr)
library(compiler)

## Forest modification proposals

## sketch
# proposals occur at:
#   1. leaf nodes
#   2. new trees (stump)
#   3. switch from leaf to distribution

# proposal types:
# a. splits ([event], time)
# b. nuisance modeling (parameters)
# c. 

# ...
# aa. ghost tree (jump across forest)
# ab. product distribution at external leaf

# replace tree[active,] with
# [1] split
# [2] trueNode
# [3] falseNode
replace_leaf = (function(tree, active, internalTrueFalseTree, markReplace=F) {
  if(!exists("replaceLeafCols")) replaceLeafCols = c("condition", "timeFrame", "trueNode","falseNode","distribution",
                                                     "par0","par1","par2","modifiable")
  if(internalTrueFalseTree[[1,"trueNode"]] != 2 || internalTrueFalseTree[[1,"falseNode"]] != 3)
    stop("internalTrueFalseTree must point trueNode to second row and falseNode to third row")
  internalTrueFalseTree[1,c("trueNode","falseNode")] = (dim(tree)[1] + 1:2)
  tree[active,] = internalTrueFalseTree[1,replaceLeafCols] 
  tree = tree %>% bind_rows(internalTrueFalseTree[-1,replaceLeafCols])
  
  if(markReplace) {
    tree = tree %>%
      mutate(replaced=rep(F,nrow(.)))
    tree[c(active,nrow(tree)-1, nrow(tree)),"replaced"] = T
    tree
  } else {
    tree
  }
}) %>% cmpfun()

modify_leaf = (function(tree, active, leaf) {
  if(!exists("replaceLeafCols")) replaceLeafCols = c("condition", "timeFrame", "trueNode","falseNode","distribution",
                                                     "par0","par1","par2","modifiable")
  tree[active,] = leaf[replaceLeafCols]
  tree
}) %>% cmpfun()

# Adds a blank tree to the end of the forest;
#   modifiable: 0 --> fixed, 1 --> modifiable, NA --> raise node (modified through raise (par[3]) only)
add_blank_tree = (function(forest, modifiable = 1) {
  forest$tree = forest$tree %>% rbind(list("",NA,NA,NA,1,NA,0,NA,modifiable))
  forest$roots = c(forest$roots, nrow(forest$tree))
  forest
}) %>% cmpfun()

# Get a blank tree
get_blank_tree = (function(modifiable = 1) {
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
  
  tree[dim(tree)[1]+1,] = c("" ,NA,NA,NA,1 ,NA,0 ,NA, modifiable)
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
  tree
}) %>% cmpfun()

# Get a blank forest
get_blank_forest = (function(modifiable = 1) {
  forest = list(tree = get_blank_tree(modifiable), roots = 1)
  forest
}) %>% cmpfun()

# The raise node adjusts the logmultiplier by par1, which is set by par[3] in non-exp optimization
add_raise_node = (function(forest) {
  # if(sum(is.na(forest[["tree"]][,"modifiable"]) &
  #        forest[["tree"]][,"distribution"]==1, na.rm=T)==0) {
    add_blank_tree(forest, modifiable = NA)
  # } else {
  #   warning("forest already has raise node, no raise node added")
  #   forest
  # }
}) %>% cmpfun()

# generate proposals
# 1. random generation
# 2. identified-pattern generation
# 3. exhaustive

# generate structure
# identify parameter estimates
# select best
# repeat


# intermediate representation: compute unexpanded forest
# then conduct lesion as necessary

add_internaltruefalsetrees = (function(tree, events, timeFrames) {
  n = length(events)
  if(n==0) return(tree)
  tree %>% 
    rbind(list(c(NA           ,"","") %>% rep(n) %>% replace(seq(1,3*n,3),events %>% as.character()),
               c(NA           ,NA,NA) %>% rep(n) %>% replace(seq(1,3*n,3),timeFrames),
               (c(dim(tree)[1],NA,NA) %>% rep(n) ) + 0 + 1:(3*n),
               (c(dim(tree)[1],NA,NA) %>% rep(n) ) + 1 + 1:(3*n),
               c(NA           ,1,1)   %>% rep(n), 
               c(NA           ,NA,NA) %>% rep(n),
               c(NA           ,0,0)   %>% rep(n),
               c(NA           ,NA,NA)  %>% rep(n),
               c(NA           ,1, 1) %>% rep(n)
               ))
    #rbind(list(""   ,        NA,NA,NA,1 ,0 ,NA)) %>%
    #rbind(list(""   ,        NA,NA,NA,1 ,0 ,NA))
}) %>% cmpfun()
add_modificationleaves = (function(tree, events, distributions, modifiable=1) {
  n = length(distributions)
  chosenEvents = sample(events,size = n)
  tree %>% 
    rbind(list(rep("",n) %>% (function(.) {.[distributions==2]=as.character(chosenEvents[distributions==2]); .})(.),
               rep(NA,n),
               rep(NA,n),
               rep(NA,n),
               distributions,
               chosenEvents, # par0: for distribution == 2, this is the raise value; it is NA for distribution == 1 or NA
               rep(0,n), # par1
               rep(NA,n), # par2
               rep(modifiable, n) # modifiable
               ))
}) %>% cmpfun()

generate_parameters_from_data_distribution = (function(target, dats, probSplit=0.5, localityFactor = 2) {
  list(probSplit = probSplit,
       type = "distribution",
       dots = list(localityFactor = localityFactor, dats=dats, target=target)
  )
}) %>% cmpfun()
generate_parameters_from_gapped_data_distribution = (function(target, dats, probSplit=0.5, localityFactor = 2, gap) {
  list(probSplit = probSplit,
       type = "gap",
       dots = list(localityFactor = localityFactor, dats=dats, target=target, gap = gap)
  )
  #gap is a df with columns size and gap
}) %>% cmpfun()

reordered_sample = function(ns, dots) {
  do.call("sample", list(x=dots[["t"]], size=ns))
}
noisy_reordered_sample = function(ns, dots) {
  do.call("sample", list(x=dots[["t"]], size=ns)) * 2 * runif(ns)^dots$localityFactor
}

generate_random_proposals = (function(tree,
                                     events,
                                     proposals=1,
                                     pars = list(
                                       probSplit=0.5,
                                       type="default",
                                       localityFactor = 2,
                                       maxLength=100,
                                       d = NA, # distribution function to sample from
                                       dots = NA # parameters to pass to 'd', with nsplits preceding those arguments
                                       )
                                     ) {
  # browser()
  # tis = sample(which(tree[["condition"]]==""), size = proposals, replace=T)
  if(is.null(pars))
    pars = list(
      probSplit=0.5,
      type="default",
      localityFactor = 2,
      maxLength=100,
      d = NA, # distribution function to sample from
      dots = NA # parameters to pass to 'd', with nsplits preceding those arguments
    )
  nsplits = floor(pars$probSplit*proposals)
  nmods = floor((1-pars$probSplit)*proposals)
  if(nsplits + nmods != proposals ) {
    if(proposals*(pars$probSplit-nsplits/proposals) > runif(1))
      nsplits = nsplits + 1
    else
      nmods = nmods + 1
  }

  if(pars$type == "default") {
    proposalTrees = tree[1,] %>%
      add_internaltruefalsetrees(sample(events,nsplits, replace=T),
                                 pars$maxLength^(runif(nsplits)^pars$localityFactor)-1
                                 ) %>%
      .[-1,] %>%
      mutate(proposalId = rep(seq(1,length.out = nsplits), each=3),
             proposalType = 1)
  } else if(pars$type == "distribution" | pars$type =="gap") {
    proposedEvents = sample(events,nsplits, replace=T)
    times = pars$dots$dats  %>%
      mutate(ets = list(pts %>% filter(event==proposedEvents[1] | event==pars$dots$target))) %>%  sample_n(5)%>%
      select(id,ets) %>% unnest() %>% 
      group_by(id) %>% mutate(nt = lag(cumsum(event==pars$dots$target),default = 0)) %>%
      group_by(id, nt) %>% mutate(wantedTimes = max(time)-time) %>%
      ungroup() %>% .[["wantedTimes"]] %>% .[.!=0] %>%
      (function(.) if(length(.)==0) 1e4 else .)(.)
    sampledTime = sample(times,1) *2*(runif(1)^pars$dots$localityFactor)
    if(pars$type=="gap") {
      minTime = pars$dots$gap %>% mutate(first = cumsum(nrow(tree)<size)) %>% filter(first == 1) %>% (function(.) if(nrow(.)==0) 0 else .[[1,"gap"]])(.)
      tries = 0
      while(sampledTime < minTime & tries < 1000) {
        times = pars$dots$dats  %>%
          mutate(ets = list(pts %>% filter(event==proposedEvents[1] | event==pars$dots$target))) %>%  sample_n(5)%>%
          select(id,ets) %>% unnest() %>% 
          group_by(id) %>% mutate(nt = lag(cumsum(event==pars$dots$target),default = 0)) %>%
          group_by(id, nt) %>% mutate(wantedTimes = max(time)-time) %>%
          ungroup() %>% .[["wantedTimes"]] %>% .[.!=0] %>%
          (function(.) if(length(.)==0) 1e4 else .)(.)
        sampledTime = sample(times,1)*2*(runif(1)^pars$dots$localityFactor)
        tries = tries + 1
      }
      if(tries == 1000)
        warning("tried 1000 times, no available sampling")
    }
    proposalTrees = tree[1,] %>%
      # add_internaltruefalsetrees(proposedEvents,
      #                            do.call(pars$d,c(nsplits,pars["dots"])) #only works with functions who take lists as inputs
      # ) %>%
      add_internaltruefalsetrees(proposedEvents,
                                 sampledTime   
      ) %>%
      .[-1,] %>%
      mutate(proposalId = rep(seq(1,length.out = nsplits), each=3),
             proposalType = 1)      
  }
  
  
  proposalMods = tree[1,] %>%
    add_modificationleaves(events, rep(2,nmods))%>%#, replace = T)) %>%
    .[-1,] %>% 
    mutate(proposalId = seq(3*nsplits+1,length.out = nmods),
           proposalType = 2)
    
  proposalTrees %>% bind_rows(proposalMods)
}) %>% cmpfun()

generate_proposal_cuts = (function(pars, n=1) {
  pars$maxLength^(runif(n)^pars$localityFactor)-1
}) %>% cmpfun()
# given a list of proposals, produce a map object for the reducer to use.
#   the map object will contain all the precompute statistics necessary for the reduce function
#   - for exponential proposal, this is mlss
#   - for loglogistic proposal, this is relativeTimeLLInput
#   proposal should be a list (or df) containing:
#   - $proposalType
#   - $existingTree
#   - $roots
#   - $proposalLocation (in tree)
#   - $proposedTree
evaluate_proposals_map = (function(proposal, od, target) {
  od = od %>%
    group_by(id) %>%
    mutate(withProposal = list(
         get_withproposal_object(tree = proposal[["existingForest"]][["tree"]],
                                 proposedTree = proposal[["proposedForest"]][["tree"]],
                                 roots = proposal[["existingForest"]][["roots"]],
                                 proposalLocation = proposal$proposalLocation,
                                 data = pts[[1]],
                                 lbubs = lbubs[[1]])
         )) %>%
    ungroup()
                                        
  # TODO ensure the max and min logRates enter into the MLE (MAP of exp, but REDUCE of ll)
  if(proposal$proposalType == 1 || 
     (proposal$proposalType == 2 && 
      proposal[["proposedForest"]][["tree"]][proposal$proposedForest$tree$replaced==T,"distribution"]==1)
     ) { #exp
    #return mlss
    od %>% 
      group_by(id) %>%
      mutate(mapobject = 
           list(get_mlss(proposal[["proposedForest"]][["tree"]], proposal[["proposedForest"]][["roots"]],
                         withProposal[[1]], pts[[1]], lbubs[[1]], target))) %>%
      #TODO caution only takes first of ids in group_by (so problem if nonunique ids)
      ungroup()
  } else if(proposal$proposalType == 2) {
    od %>% 
      group_by(id) %>% 
      mutate(mapobject = 
           list(get_ll_mapobject(proposedTree = proposal[["proposedForest"]][["tree"]],
                 proposalLocation = proposal$proposalLocation, 
                 withProposal = withProposal[[1]],
                 data = pts[[1]],
                 target = target) 
                )
           ) %>% 
      ungroup()
  }
}) %>% cmpfun()

get_withproposal_object = (function(tree, proposedTree, roots, proposalLocation, data, lbubs) {
  # 1. make additional window splits for proposals affecting durations
  # 2. if split only affects one tree, combine factors for other trees
  #    b. compute this for all held out trees
  # 3. per proposal:
  #    filter out unaffected durations
  #    apply proposal change
  
  
  #withProposal[1,"treeIndex"] = 5 # to test having distribution type 2 (log logistic) overwriting
  withProposal = proposedTree %>%
    subtimeline_forest_preexpansion(data,roots,lbubs, 4) %>% 
    do((function(x) if(!("replacedIndex" %in% names(x))) 
      mutate(x, replacedIndex=
               ifelse(treeIndex==proposalLocation,
                      yes = proposalLocation, 
                      no=NA)
             ) else x)(.)) %>%
    mutate(oldIndexInNewTree = (function(ri,ti) { ninsa = !is.na(ri); ti[ninsa] = ri[ninsa]; ti })(replacedIndex,treeIndex)) %>% 
    group_by(treeNumber,treeIndex,replacedIndex,oldIndexInNewTree) %>%
    do((function(x) {ret_subtimeline_tree(tree = tree,
                                         data = data,
                                         active = x$oldIndexInNewTree,
                                         lbubs = x,
                                         leafType = 1)})(.)) %>%
    ungroup() %>%
    select(-oldIndexInNewTree) %>%
    arrange(treeNumber,lb,ub)
    
    #withProposal %>% do((function(x) data.frame(la=mean(x$treeIndex),lo=T))(.)) #combines into tibble if same type/size of tbl
  withProposal
}) %>% cmpfun()

get_mlss = (function(tree, roots, withProposal, data, lbubs, target) {
  # versus
  #tree %>% subtimeline_forest_preexpansion(data,roots,lbubs, 4) %>% arrange(treeNumber,lb, lb-ub)

  # expanded computation of sum(logRateFactor)
  # collect events: Mp_hat and qp*Tp_hat
  expansion = withProposal %>% arrange(lb) %>% (function(.) fix_timelineCPP(.$lb, .$ub))
  factors = withProposal %>% select(lb,ub,logRateFactor) %>% 
    mutate_active_duration_indices(expansion) %>%
    value_sum_onto_durations(expansion, value="logRateFactor")
  expansionCounts = findInterval((data %>% filter(event==target) %>% .[["time"]] - 1e-14), expansion$lb, left.open=T, rightmost.closed = T) %>%
    vector_of_counts(nrow(expansion))
  # note fudge factor to deal with double precision problems of R's findInterval(...); this limits the precision to 1e-14 (finer will cause uncaught errors!)
  
  ### Proposal type 1: branch exponential calculation
  liCalculation = expansion %>%
    mutate(factors=factors) %>%
    mutate(qptp_hat = exp(factors)*(ub-lb)) %>% # may need to conduct proposal expansion prior to here
    mutate(Mp_hat = expansionCounts) %>% tbl_df() # may need to conduct proposal expansion prior to here

  # convert expansion into summaries in collapsed form (i.e. collapse)
  collapsedCalculation = withProposal %>% mutate_active_duration_indices(expansion) %>%
    mutate(qptp_hat = value_sum_from_durations(.,liCalculation,"qptp_hat"),
           Mp_hat = value_sum_from_durations(.,liCalculation,"Mp_hat")
    )
  
  # filter the replacedIndex intervals, summarise, and get qp', and thus dLL.  
  mlss = collapsedCalculation %>% filter(!is.na(replacedIndex)) %>% # TODO can we move up this filter command to make faster?
    group_by(treeIndex, replacedIndex) %>%
    summarise(M = sum(Mp_hat), QT = sum(qptp_hat)) %>% ungroup() # possible double counting here if overlapping mods, watch out
  
  # maximum likelihood is just those by log qp' and (qp'-1)
  mlss = mlss %>% mutate(log_qp_prime = (function(.) { x = log(.$M)-log(.$QT); x[is.infinite(x)]=-20; x})(.), dll = M*log_qp_prime - expm1(log_qp_prime)*QT)
  mlss 
}) %>% cmpfun()

get_ll_mapobject = (function(proposedTree, proposalLocation, withProposal, data, target) {
  ### Proposal type 2: log logistic calculation
  # get logmultipliers across active areas in withProposal expansions + data expansions
  # these go into llcll(df = active areas, pc = logmultipliers in active areas)
  # these need to be aligned with respect to lastTime -- couldn't the logmultipliers be different for
  # different active areas? yes ... so then calculation should be over logmultiplier-specific event,
  # i.e. no overlapping intervals. How to compute from this (single df, not {df, pc})?
  activeLogLogistic = withProposal %>% 
    filter(!is.na(replacedIndex)) # i.e. the durations to covers
  if(nrow(activeLogLogistic)==0)
    return(withProposal %>%
             mutate(M=0, lastTime=-Inf) %>%
             select(-replacedIndex,-treeIndex,-treeNumber) %>%
             .[0,])
    
  # expand for the model (multiple trees preexpanded and overlapping) and
  #   sum the logRateFactors *excluding* replaceIndex?
  expansion = withProposal %>% arrange(lb) %>% (function(.) fix_timelineCPP(.$lb, .$ub))
  factors = withProposal %>% 
    filter(is.na(replacedIndex)) %>% # exclude logRateFactors that will be replaced
    select(lb,ub,logRateFactor) %>% 
    mutate_active_duration_indices(expansion) %>%
    value_sum_onto_durations(expansion, value="logRateFactor")
  expansion = expansion %>% mutate(summedLogRateFactor = factors)

  # expand further to account for data events    
  
  expansionOnData = expansion %>%
    insert_events_into_expansion(data %>% filter(event==target) %>% .[["time"]])
  
  expansionLogLogistic = expansionOnData %>%
    mutate(logRateFactor=value_apply_onto_durations(
      expansion %>% mutate_active_duration_indices(expansionOnData), ., "summedLogRateFactor"))
  
  activeLogLogistic = activeLogLogistic %>% 
    mutate_active_duration_indices(expansionLogLogistic)
  
  relativeTimeLLInput = expansionLogLogistic %>%
    # DOTO you can't just throw away overlapping intervals in other trees
    # now filter down to active expanded intervals (active is in withProposal)
    mutate(replacedIndex=value_apply_onto_durations(activeLogLogistic, ., "replacedIndex")) %>%
    filter(replacedIndex > 0) %>%
    mutate(lastTime = (function(lb,eventTimes)
      c(-Inf, eventTimes, Inf)[1+findInterval(lb, eventTimes)]
      )(.$lb,filter(data, data$event==proposedTree[[proposalLocation,"condition"]]) %>% .[["time"]]) ) %>% #parameter to manipulate here # old par0->condition
    mutate(lb = ifelse(is.finite(lb-lastTime), lb-lastTime, lb), 
           ub = ifelse(is.finite(ub-lastTime), ub-lastTime, ub)) %>% 
    #      if "never", do not overwrite times with -Inf, instead let the next function deal with the non-shift
    select(-replacedIndex)

  # browser()
  relativeTimeLLInput
}) %>% cmpfun()

# od should contain column of mapobjects
evaluate_proposals_reduce = (function(od, proposal) {
  if(proposal$proposalType == 1 ||
     (proposal$proposalType == 2 &&
      proposal[["proposedForest"]][["tree"]][proposal[["proposedForest"]][["tree"]]$replaced==T,"distribution"]==1)
     ) {
    reduce_exp_mapobjects(od)
  } else if(proposal$proposalType == 2) {
    reduce_ll_mapobjects(od)
  }
}) %>% cmpfun()

reduce_exp_mapobjects = (function(od) {
  od %>% select(mapobject) %>% 
    unnest() %>% 
    group_by(treeIndex, replacedIndex) %>% 
    summarise(M = sum(M), QT = sum(QT)) %>%
    ungroup() %>%
    mutate(log_qp_prime = (function(.) { x = log(.$M)-log(.$QT); x[is.infinite(x)]=-100; x})(.),
           dll = M*log_qp_prime - expm1(log_qp_prime)*QT)
}) %>% cmpfun()

reduce_ll_mapobjects = (function(od) {
  rtlli = od %>% select(id,mapobject) %>% unnest()
  if(nrow(rtlli) == 0)
    return(data.frame(scale = 0,
                      shape = 0,
                      raise = 0,
                      dll = 0) %>% tbl_df() %>% .[0,])
  # opt = constrOptim(theta=c(1,1,-1), # Linear constraints Nelder-Mead (slow)
  #                   f=llcll,
  #                   grad=NULL,
  #                   ui=rbind( # shape must be > 0.01, maxIncrease must be less than 5
  #                     c(0,1,0),
  #                     c(0,0,-1)),
  #                   ci=c(0.01,-5),
  #                   df=rtlli %>% rename(event=M, logmultiplier=logRateFactor),
  #                   pc=rtlli %>% rename(event=M, logmultiplier=logRateFactor))
  
  dfpc = rtlli %>% rename(event=M, logmultiplier=logRateFactor)
  
  opt = optim(par=c(1,1,-1), # Box-constraint optimization (box on shape and maxIncrease)
              f=llcll,
              method = "L-BFGS-B",
              gr=NULL, # shape must be > 0.01, maxIncrease must be less than 5
              lower = c(-5, 0.01, -100), # TODO parameterize constraints, require li finite, use "safely"?
              upper = c(23,Inf,0),
              # control = list(trace=4),
              df=dfpc,
              pc=dfpc,
              groupId = c("id","lastTime"))
  # browser()
  
  dll = -(opt$value -
    (rtlli %>%
       mutate(ll = -M*logRateFactor + exp(logRateFactor)*(ub-lb)) %>% select(ll) %>% sum()
    ))
  # dll, larger is better, TODO is this responsible to compare likelihoods across models?
  # compare opt$value to likelihood without (which is just the exp methods of calculation)
  
  
  data.frame(scale = opt$par[1],
             shape = opt$par[2],
             # raise = 0,
             raise = opt$par[3],
             dll = dll) %>% tbl_df()
}) %>% cmpfun()

  # DOTO check that logRateFactor does not contain the proposal or the old one? I think it includes one?
  
  # TODO formulate how what it means to add a multiplicative factor during optimization
  #      because it is not changing scale, which is what adding to logmultiplier does
  
  
  # get attachment of lastTime appropriate (fix precedingtime(...)?)
  
  # now convert lbubs to "time-since-last-event";
  #   but you did the expansion based on 
  #   specify defaults:
  #     "never seen" default: t=0
  #     "just saw" default: for periodics, reset t=0
  
# auto-adjust with fully-applied exp rate?
    # include in ll calc? include in optim?
    
  # liCalculation has everything to calculate pc_llll,
  #   but not really:
  #   (1) we have to filter to applicable subtimeline
  #   (2) we need to expand on when the applicable event occurs (Mp_hat is inexact)
  #   so push up above expansionCounts
  #   and define what to do when the event happens (default to unobserved? reset to t=0?)
  #     as unobserved: will make events likely, i.e. kickstarter, switch on
  #     as rest to t=0: will make events periodic, i.e. clock
  #       this requires expansions (or just overwriting of lb): (because either parent or child event resets)
  #     can do both and select/record as a parameter
  #     semantically it doesn't really matter; we are just using the method to get a instaneous rate with
  #       with effective rate parameterizations
  #     we will need to account for this in the splitting
  #     can bypass this and just choose the easier one for now
  

  # make_learning_data_structure = function() {}
  # TODO
  ### If exponential proposal, then:
  # Objective: count duration and events in each interval for treeIndex and replace index
  # y. have a row for all-but-treeNumber-old-value-multiplicative-effect-of-q.
  #    (sum all trees, subtract out unused tree number)
  #      --> this is an expansion of size O(kn), k expansion size, n number of data durations
  #    (to calculate q*T you need the multiplicative effect)
  # z. keep affected [duration x replacedIndex] only
  # 0. if replaceIndex is not expanded, expand it
  # 1. use indicator column and real value (lb-ub) column
  # 2. sum by groups of treeIndex and replaceIndex
  # 3. compute MLE and dLL.
  # 4. select proposal or stop.
  
  
  
  ### If log logistic proposal, then:
  # y. have a row for all but treeNumber old value multiplicative effect of q.
  #    (sum all trees, subtract out unused tree number)
  #    (to calculate q*T you need the multiplicative effect)
  # z. keep affected [duration x replacedIndex] only
  # 0. if replaceIndex is not expanded, expand it
  # 1. use indicator column and real value (lb-ub) column
  # 2. sum by groups of treeIndex and replaceIndex
  #    (use this as your data set to compute log logistic parameters)
  # 3. compute MLE and dLL.
  # 4. select proposal or stop.
  
  # is the the duration expansion some sort of a database operation on durations?

runTests = F
if(runTests)
  source("tests/forest_modifications_tests.r")
