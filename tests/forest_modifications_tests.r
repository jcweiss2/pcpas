#forest_modifications_test.r
source("forest_modifications.r")
source("tests/subtimeline_maker_population_tests.r")


#tests
set.seed(42)
generate_random_proposals(tree, dat[[1]]%>%unique(),proposals = 1)

replace_leaf(tree,
             active = 3,
             internalTrueFalseTree = generate_random_proposals(tree, dat[[1]]%>%unique(),proposals = 1))

tree2 = replace_leaf(tree,
                     active = 3,
                     internalTrueFalseTree = generate_random_proposals(tree, dat[[1]]%>%unique(),proposals = 1))
tree2 %>%
  subtimeline_forest_preexpansion(od$pts[[1]],rep(c(1,2,3),1),od$lbubs[[1]], 3) %>% arrange(lb, lb-ub) %>% data.frame

tree %>% subtimeline_forest_preexpansion(od$pts[[1]],3,od$lbubs[[2]], 2) %>% arrange(lb, lb-ub)
tree2 %>% subtimeline_forest_preexpansion(od$pts[[1]],3,od$lbubs[[2]], 2) %>% arrange(lb, lb-ub)

tree %>% subtimeline_forest_preexpansion(od$pts[[1]],3,od$lbubs[[1]], 3) %>% arrange(lb, lb-ub)
tree2 %>% subtimeline_forest_preexpansion(od$pts[[1]],3,od$lbubs[[1]], 3) %>% arrange(lb, lb-ub)

tree %>% subtimeline_forest_preexpansion(od$pts[[1]],3,od$lbubs[[1]], 4) %>% arrange(lb, lb-ub)
tree2 %>% subtimeline_forest_preexpansion(od$pts[[1]],3,od$lbubs[[1]], 4) %>% arrange(lb, lb-ub)

tree %>% subtimeline_forest_preexpansion(od$pts[[1]],3,od$lbubs[[1]], 1) %>% arrange(lb, lb-ub)
tree2 %>% subtimeline_forest_preexpansion(od$pts[[1]],3,od$lbubs[[1]], 1) %>% arrange(lb, lb-ub) %>% data.frame()

tree %>% subtimeline_forest_preexpansion(od$pts[[1]],3,od$lbubs[[2]], 1) %>% arrange(lb, lb-ub) %>% ll(od$pts[[1]],"A")
tree2 %>% subtimeline_forest_preexpansion(od$pts[[1]],3,od$lbubs[[2]], 1) %>% arrange(lb, lb-ub) %>% ll(od$pts[[1]],"A")

###

tree %>% subtimeline_forest_preexpansion(od$pts[[1]],rep(c(1,2,3),1),od$lbubs[[1]], 1) %>% arrange(lb, lb-ub) %>% ll(od$pts[[1]], "A")
tree2 %>% subtimeline_forest_preexpansion(od$pts[[1]],rep(c(1,2,3),1),od$lbubs[[1]], 1) %>% arrange(lb, lb-ub) %>% ll(od$pts[[1]],"A")



### test: generate proposals, do map-reduce through forest_modification functions

# proposal_location = tree %>%
#   select(condition) %>%
#   mutate(index=1:(dim(tree)[1])) %>%
#   filter(condition == "") %>% select(index) %>% .[[1]] %>% sample(1)
proposalLocation = which(!is.na(tree$modifiable) & tree$modifiable>0) %>% sample(2)
proposal = generate_random_proposals(tree, od$pts[[1]][[1]]%>%unique(),proposals = 2)
tree %>% subtimeline_forest_preexpansion(od4pts[[1]],3,od$lbubs[[1]], 4) %>%
  arrange(lb, lb-ub) %>% mutate(replace=treeIndex %in% proposalLocation)

# Proposal type 1: branch
proposedTree = tree %>% replace_leaf(active = proposalLocation[1],
                                     internalTrueFalseTree = proposal[1:3,], T)
withProposal = get_withproposal_object(tree,
                        proposedTree,
                        roots = c(1,3),
                        proposalLocation = proposalLocation[1],
                        data = od$pts[[1]],
                        lbubs = od$lbubs[[1]])
# DOTO fix bug here -- expansion for withProposal loses columns with data expansion;
#      workaround using more group variables in group_by

get_mlss(proposedTree, c(1,3), withProposal, od$pts[[1]], od$lbubs[[1]], target="A")
# DOTO: a default mlp for having no data enter a "true/false" node --> currently unrepresented in mlss;
#       doesn't appear necessary for updates

#   proposal should be a list (or df) containing:
#   - $proposalType
#   - $existingTree
#   - $roots
#   - $proposalLocation (in tree)
#   - $proposalTree
prop = list()
prop$proposalType = 1
prop$existingForest = list(tree=tree, roots = 1)
prop$proposalLocation = proposalLocation[1]
prop$proposedForest = list(tree=proposedTree, roots = 1)
evaluate_proposals_map(prop, od %>% ungroup(), target="A") %>% select(mapobject) %>% unnest()
# TODO sometimes trueNode or falseNode is not represented because no hits; need to give an MLP

evaluate_proposals_map(prop, od %>% ungroup(), target="A") %>%
  evaluate_proposals_reduce(prop)


# Proposal type 2: distribution override
proposedTree = tree %>%
  modify_leaf(active = proposalLocation[2], proposal[4,]) %>%
  # [4,] because generate 2 proposals of length 3 and 1
  mutate(replaced=(1:nrow(.)==proposalLocation[2]))
withProposal = get_withproposal_object(tree,
                        proposedTree,
                        roots = c(1,3),
                        proposalLocation = proposalLocation[2],
                        data = od$pts[[1]],
                        lbubs = od$lbubs[[1]])

get_ll_mapobject(proposedTree = proposedTree,
                 proposalLocation = proposalLocation[2], 
                 withProposal = withProposal,
                 data = od$pts[[1]],
                 target = "A") 

#   proposal should be a list (or df) containing:
#   - $proposalType
#   - $existingTree
#   - $roots
#   - $proposalLocation (in tree)
#   - $proposalTree
prop = list()
prop$proposalType = 2
prop$existingForest = list(tree=tree,roots=1)
prop$proposalLocation = proposalLocation[2]
# prop$proposalLocation = 3
prop$proposedForest = list(tree=proposedTree, roots = 1)

evaluate_proposals_map(prop, od %>% ungroup(), "A") # %>% select(mapobject) %>% unnest()

evaluate_proposals_map(prop, od %>% ungroup(), "A") %>%
  evaluate_proposals_reduce(prop)

#simpler:
# llcll(pars = c(1,1,-1), df = rtlli %>% rename(event=M, logmultiplier=logRateFactor),
#       pc = rtlli %>% rename(event=M, logmultiplier=logRateFactor),
#       groupId = c("id","lastTime"))

# TODO do appropriate grouping by lastTime
# plot(x=seq(0,20,0.1),y=dllog(seq(0,20,0.1), scale = opt$par[1], shape=opt$par[2]), pch=18, log="y", ylim=c(1e-10,1))
# 
# testmixed = with(rtlli,
#                  data.frame(ub = ub, lb = lb,
#                             qlb = pllog(rtlli$lb, scale = opt$par[1], shape = opt$par[2]),
#                             qub = pllog(rtlli$ub, scale = opt$par[1], shape = opt$par[2])) %>%
#                    mutate(logRateFactor = log( -log(1-(qub-qlb))/(ub-lb)))) %>%
#   rename(logmultiplier = logRateFactor) %>%
#   mutate(combined = logmultiplier+rtlli$logRateFactor)
# segments(x0 = rtlli$lb, x1=rtlli$ub,
#          y0=exp(rtlli$logRateFactor), col=2)
# with(testmixed, segments(x0=lb,x1=ub,
#                          y0=exp(combined)*exp(max(testmixed$logmultiplier)-max(testmixed$combined)),
#                          y1=exp(combined)*exp(max(testmixed$logmultiplier)-max(testmixed$combined)),
#                          col="blue")) #TODO fix (doesn't use par[3])
# legend(x = "bottomright", legend = c("loglogistic", "prior","posterior"), col=c("black","red","blue"), lwd=2)
# abline(v=1)
# abline(v=8)

doOld = F
if(doOld) {
rt = relativeTimeLLInput %>% rename(event=M, logmultiplier=logRateFactor)
opt = constrOptim(theta=c(1,1,-1),
                  f=llcll,
                  grad=NULL,
                  ui=rbind( # shape must be > 0.01, maxIncrease must be less than 5
                    c(0,1,0),
                    c(0,0,-1)),
                  ci=c(0.01,-5),
                  df=rt,
                  pc=rt)

dll = -opt$value -
  (relativeTimeLLInput %>%
     mutate(ll = M*logRateFactor - exp(logRateFactor)*(ub-lb)) %>% select(ll) %>% sum()
  ) # larger is better

#simpler:
llcll(pars = c(1,1,-1), df = relativeTimeLLInput %>% rename(event=M, logmultiplier=logRateFactor),
      pc = relativeTimeLLInput %>% rename(event=M, logmultiplier=logRateFactor),
      groupId = "lastTime")

# TODO do appropriate grouping by lastTime
plot(x=seq(0,20,0.1),y=dllog(seq(0,20,0.1), scale = opt$par[1], shape=opt$par[2]), pch=18, log="y", ylim=c(1e-20,1))

testmixed = with(relativeTimeLLInput,
                 data.frame(ub = ub, lb = lb,
                            qlb = pllog(relativeTimeLLInput$lb, scale = opt$par[1], shape = opt$par[2]),
                            qub = pllog(relativeTimeLLInput$ub, scale = opt$par[1], shape = opt$par[2])) %>%
                   mutate(logRateFactor = log( -log(1-(qub-qlb))/(ub-lb)))) %>%
  rename(logmultiplier = logRateFactor) %>%
  mutate(combined = logmultiplier+relativeTimeLLInput$logRateFactor)
segments(x0 = relativeTimeLLInput$lb, x1=relativeTimeLLInput$ub,
         y0=exp(relativeTimeLLInput$logRateFactor), col=2)
with(testmixed, segments(x0=lb,x1=ub,
                         y0=exp(combined)*exp(max(testmixed$logmultiplier)-max(testmixed$combined)),
                         y1=exp(combined)*exp(max(testmixed$logmultiplier)-max(testmixed$combined)),
                         col="blue")) #TODO fix (doesn't use par[3])
legend(x = "bottomright", legend = c("loglogistic", "prior","posterior"), col=c("black","red","blue"), lwd=2)
abline(v=1)
abline(v=8)
}
