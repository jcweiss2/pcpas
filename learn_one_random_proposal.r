### learn_one_random_proposal.r
require(dplyr)
require(compiler)

proposal_location_generator = (function(tree) {
  # browser()
  data.frame(index = 1:nrow(tree),
             modifiableValue = exp(tree$modifiable)) %>% 
    as_tibble() %>%
    filter(!is.na(modifiableValue)) %>% mutate(modifiableValue = cumsum(modifiableValue)/sum(modifiableValue)) %>%
    mutate(belowIndex = (runif(1) < modifiableValue)*index) %>% .[["belowIndex"]] %>% (function(.) min(.[.>0]))(.)
}) %>% cmpfun()

# generates a prop(osal) object
generate_one_proposal_po = (function(forest, dats, ...) {
  # proposalLocation = which(forest$tree$modifiable > 0) %>% sample(1)
  proposalLocation = proposal_location_generator(forest$tree)
  dots = list(...)[names(list(...)) %in% names(formals(generate_random_proposals))]
  proposal = do.call("generate_random_proposals",
                     list(forest[["tree"]],
                       dats[[1,"pts"]][["event"]]%>%unique()%>%as.character(),
                       proposals = 1,
                       dots$pars)
  )

  
  # Proposal type 1: branch
  if(proposal$proposalType[1] == 1)
    proposedTree = forest$tree %>% replace_leaf(active = proposalLocation[1],
                                       internalTrueFalseTree = proposal[1:3,], T)
  # TODO reconcile # Proposal type 2: distribution override
  else if(proposal$proposalType[1] == 2)
    proposedTree = forest$tree %>%
      modify_leaf(active = proposalLocation, proposal[1,]) %>%
      mutate(replaced=(1:nrow(.)==proposalLocation))
  
  proposedForest = list(tree = proposedTree, roots = forest$roots)

  # TODO: a default mlp for having no data enter a "true/false" node --> currently unrepresented in mlss 
  
  #   proposal should be a list (or df) containing:
  #   - $proposalType
  #   - $existingTree
  #   - $roots
  #   - $proposalLocation (in tree)
  #   - $proposalTree
  prop = list()
  prop$proposalType = proposal$proposalType[1]
  prop$existingForest = forest
  prop$proposalLocation = proposalLocation
  prop$proposedForest = proposedForest
  # TODO sometimes trueNode or falseNode is not represented because no hits; need to give an MLP
  
  prop
}) %>% cmpfun()

prop_log_logistic = (function(forest, dats, ...) {
  dots = list(...)
  pars = dots$pars
  proposedForest = list(
    tree = add_modificationleaves(forest$tree,pars$proposalEvent,2) %>% mutate(replaced=F) %>% (function(.) {.$replaced[nrow(.)]=T;.})(.),
    roots = c(forest$roots,nrow(forest$tree)+1)
  )
  #   proposal should be a list (or df) containing:
  #   - $proposalType
  #   - $existingTree
  #   - $roots
  #   - $proposalLocation (in tree)
  #   - $proposalTree
  prop = list()
  prop$proposalType = 2
  prop$existingForest = forest %>% add_blank_tree()
  prop$proposalLocation = nrow(forest$tree)+1
  prop$proposedForest = proposedForest
  prop
}) %>% cmpfun()

learn_pc = (function(forest, od, target, iter=10, pf = generate_one_proposal_po, maintainBlankTree=T,...) {
  model = list()
  model$forest = forest
  learn_pc_model(model, od, target, iter, pf, maintainBlankTree, ...)
}) %>% cmpfun()
learn_pc_model = (function(model, od, target, iter, pf = generate_one_proposal_po, maintainBlankTree,...) {
  if(iter==0) return(model)
  cat(paste0(iter," "))
  
  prop = pf(model$forest, od, ...)
  
  proposalResults = evaluate_proposals_map(prop, od %>% ungroup(), target) %>%
    evaluate_proposals_reduce(prop)
  
  dll = proposalResults$dll %>% sum()
  # browser()
  
  if(is.nan(dll)) {
    # browser()
    warning("proposal caused dll to be NaN")
  } else if(dll > 1) {# AIC criterion (1 par --> 2 pars, net 1)
    # if(dll > 900)
    #   browser()
    # accept proposal
    if(prop$proposalType == 1 ||
       (prop$proposalType==2 && prop$proposedForest$tree[prop$proposedForest$tree$replaced==T,"distribution"]==1)) { #TODO debug this case
      log_qp = model$forest$tree[proposalResults$replacedIndex,"par1"]
      
      prop$proposedForest$tree = prop$proposedForest$tree %>% select(-replaced)
      prop$proposedForest$tree[proposalResults$treeIndex,"par1"] = proposalResults$log_qp_prime + log_qp
      model$forest = prop$proposedForest
      # model$forest$tree = prop$proposedTree %>% select(-replaced)
      # model$forest$tree[proposalResults$treeIndex, "par1"] = proposalResults$log_qp_prime + log_qp
    } else if(prop$proposalType == 2) {
      prop$proposedForest$tree[prop$proposedForest$tree$replaced==T, c("par0","par1","par2","modifiable")] =
        c(proposalResults$raise, proposalResults$scale, proposalResults$shape, 1)
      #TODO do something with shape
      prop$proposedForest$tree = prop$proposedForest$tree %>% select(-replaced)
      
      
      # OBSOLETE, now use par0 as raise when hits ll node # Raise node: to adjust for ll distribution multiplicative hazard
      # # TODO But you need multiple raise nodes potentially (for non-root log-logistic distribution adjustments, so let's turn off par3)
      # if(sum(is.na(prop[["proposedForest"]][["tree"]][,"modifiable"]) &
      #        prop[["proposedForest"]][["tree"]][,"distribution"]==1, na.rm = T
      #        ) == 0)
      #   prop$proposedForest = prop$proposedForest %>% add_raise_node()
      # if(sum(is.na(prop[["proposedForest"]][["tree"]][,"modifiable"]) &
      #        prop[["proposedForest"]][["tree"]][,"distribution"]==1, na.rm = T
      #        ) == 1) { # raise node
      #   raiseNode = which(is.na(prop[["proposedForest"]][["tree"]][,"modifiable"]) &
      #                     prop[["proposedForest"]][["tree"]][,"distribution"]==1)
      #   prop[["proposedForest"]][["tree"]][raiseNode,"par1"] =
      #     # prop[["proposedForest"]][["tree"]][raiseNode,"par1"] +  # TODO even with this commented out raise is still way off
      #     proposalResults$raise
      # }
      model$forest = prop$proposedForest
      
    }
    
    if(maintainBlankTree) {
      # browser()
      if(sum(
            model[["forest"]][["tree"]][ model[["forest"]][["roots"]], "modifiable" ] > 0 &
            model[["forest"]][["tree"]][ model[["forest"]][["roots"]], "distribution" ]==1,
            na.rm = T
            ) == 0
      ) {
        dots = list(...)[names(list(...)) %in% names(formals(add_blank_tree))]
        if(length(dots)==0)
          model[["forest"]] = model[["forest"]] %>% add_blank_tree()
        else
          model[["forest"]] = model[["forest"]] %>% add_blank_tree(dots$modifiable)
      }
    }
    
    if(with(model,!exists("stats")))
      model$stats = data.frame(iter = iter, dll = dll) %>% tbl_df() %>% mutate(forest = list(model$forest))
    else
      model$stats = model$stats %>% bind_rows(data.frame(iter=iter, dll=dll) %>% tbl_df() %>% mutate(forest=list(model$forest)))
  }
  learn_pc_model(model, od, target, iter-1, pf, maintainBlankTree, ...)

}) %>% cmpfun()
