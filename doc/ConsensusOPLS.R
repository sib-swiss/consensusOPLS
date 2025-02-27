## ----setup, class.source='fold-hide'------------------------------------------
#install.packages("knitr")
library(knitr)
opts_chunk$set(echo = TRUE)

# To ensure reproducibility
set.seed(12)

## ----packages_installation, eval=FALSE, class.source='fold-hide'--------------
#  pkgs <- c("ggplot2", "ggrepel", "DT", "psych")
#  sapply(pkgs, function(x) {
#      if (!requireNamespace(x, quietly = TRUE)) {
#          install.packages(x)
#      }
#  })
#  if (!require("BiocManager", quietly = TRUE))
#      install.packages("BiocManager")
#  if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
#      BiocManager::install("ComplexHeatmap")
#  }

## ----packages_load, warning=FALSE, message=FALSE, class.source='fold-hide'----
library(ggplot2) # to make beautiful graphs
library(ggrepel) # to annotate ggplot2 graph
library(DT) # to make interactive data tables
library(psych) # to make specific quantitative summaries
library(ComplexHeatmap) # to make heatmap with density plot

## ----theme_ggplot2, class.source='fold-hide'----------------------------------
theme_graphs <- theme_bw() + theme(strip.text = element_text(size=14),
                                   axis.title = element_text(size=16),
                                   axis.text = element_text(size=14),
                                   plot.title = element_text(size=16),
                                   legend.title = element_text(size=14))

## ----load_ConsensusOPLS_package-----------------------------------------------
# Detaching and uninstalling a package
# detach("package:ConsensusOPLS", unload=T)
# remove.packages("ConsensusOPLS")

# Installing a package
#devtools::install_github("sib-swiss/consensusOPLS/codes/ConsensusOPLS")
library(ConsensusOPLS)

## ----import_demo_3_Omics------------------------------------------------------
data(demo_3_Omics, package = "ConsensusOPLS")

## ----check_dims, class.source='fold-hide'-------------------------------------
# Check dimension
BlockNames <- c("MetaboData", "MicroData", "ProteoData")
nbrBlocs <- length(BlockNames)
dims <- lapply(X=demo_3_Omics[BlockNames], FUN=dim)
names(dims) <- BlockNames
dims

# Remove unuseful object for the next steps
rm(dims)

## ----check_orders, class.source='fold-hide'-----------------------------------
rns <- do.call(cbind, lapply(X=demo_3_Omics[BlockNames], rownames))
rns.unique <- apply(rns, 1, function(x) length(unique(x)))
stopifnot(max(rns.unique) == 1)

## ----describe_data_by_Y_function, class.source='fold-hide'--------------------
require(psych)
require(DT)
describe_data_by_Y <- function(data, group) {
    bloc_by_Y <- describeBy(x = data, group = group,
                            mat = TRUE)[, c("group1", "n", "mean", "sd",
                                            "median", "min", "max", "range",
                                            "se")]
    bloc_by_Y[3:ncol(bloc_by_Y)] <- round(bloc_by_Y[3:ncol(bloc_by_Y)], 
                                          digits = 2)
  return (datatable(bloc_by_Y))
}

## ----describe_data_by_Y_bloc1-------------------------------------------------
describe_data_by_Y(data = demo_3_Omics[[BlockNames[1]]],
                   group = demo_3_Omics$ObsNames[,1])

## ----describe_data_by_Y_bloc2-------------------------------------------------
describe_data_by_Y(data = demo_3_Omics[[BlockNames[2]]],
                   group = demo_3_Omics$ObsNames[,1])

## ----describe_data_by_Y_bloc3-------------------------------------------------
describe_data_by_Y(data = demo_3_Omics[[BlockNames[3]]],
                   group = demo_3_Omics$ObsNames[,1])

## ----scale_data, class.source='fold-hide'-------------------------------------
# Save not scaled data
demo_3_Omics_not_scaled <- demo_3_Omics

# Scaling data
demo_3_Omics[BlockNames] <- lapply(X=demo_3_Omics[BlockNames], FUN=function(x)
    scale(x, center = TRUE, scale = TRUE)
)

## ----heatmap_function, message = FALSE, class.source='fold-hide'--------------
heatmap_data <- function(data, bloc_name, factor = NULL){
  if(!is.null(factor)){
    ht <- Heatmap(
      matrix = data, name = "Values",
      row_dend_width = unit(3, "cm"),
      column_dend_height = unit(3, "cm"),
      column_title = paste0("Heatmap of ", bloc_name),
      row_split = factor,
      row_title = "Y = %s",
      row_title_rot = 0
    )
  } else{
    ht <- Heatmap(
      matrix = data, name = "Values",
      row_dend_width = unit(3, "cm"),
      column_dend_height = unit(3, "cm"),
      column_title = paste0("Heatmap of ", bloc_name)
    )
  }
  return(ht)
}

## ----heatmap_no_scale, message = FALSE, class.source='fold-hide'--------------
# Heat map for each data block
lapply(X = 1:nbrBlocs,
       FUN = function(X){
         bloc <- BlockNames[X]
         heatmap_data(data = demo_3_Omics_not_scaled[[bloc]],
                      bloc_name = bloc,
                      factor = demo_3_Omics_not_scaled$Y[,1])})

## ----heatmap_scale, message = FALSE, class.source='fold-hide'-----------------
# Heat map for each data block
lapply(X = 1:nbrBlocs,
       FUN = function(X){
         bloc <- BlockNames[X]
         heatmap_data(data = demo_3_Omics[[bloc]],
                      bloc_name = bloc,
                      factor = demo_3_Omics$Y[,1])})

## ----heatmap_density, class.source='fold-hide'--------------------------------
# Heatmap with density for each data bloc
lapply(X = 1:nbrBlocs,
       FUN = function(X){
         bloc <- BlockNames[X]
         factor <- demo_3_Omics$Y[, 1]
         densityHeatmap(t(demo_3_Omics[[bloc]]),
                        ylab = bloc,
                        column_split  = factor,
                        column_title = "Y = %s")})

## ----rm_unscale_data, class.source='fold-hide'--------------------------------
# Remove unscaled data
rm(demo_3_Omics_not_scaled)

## ----define_cv_parameters-----------------------------------------------------
# Number of predictive component(s)
LVsPred <- 1

# Maximum number of orthogonal components
LVsOrtho <- 3

# Number of cross-validation folds
CVfolds <- nrow(demo_3_Omics[[BlockNames[[1]]]])
CVfolds

## ----run_consensusOPLSmodel---------------------------------------------------
copls.da <- ConsensusOPLS(data = demo_3_Omics[BlockNames],
                          Y = demo_3_Omics$Y,
                          maxPcomp = LVsPred,
                          maxOcomp  = LVsOrtho,
                          modelType = "da",
                          nperm = 100,
                          cvType = "nfold",
                          nfold = 14,
                          nMC = 100,
                          cvFrac = 4/5,
                          kernelParams = list(type = "p", 
                                              params = c(order = 1)),
                          mc.cores = 1)

## ----outputs_model------------------------------------------------------------
print("The main results of the ConsensusOPLS analysis:")
copls.da

print("List of available outputs in the ConsensusOPLS object:")
summary(attributes(copls.da))

## ----print_main_results_R2, class.source='fold-hide'--------------------------
position <- copls.da@nOcomp

paste0('R2: ', round(copls.da@R2Y[position], 4))

## ----print_main_results_Q2, class.source='fold-hide'--------------------------
paste0('Q2: ', round(copls.da@Q2[position], 4))

## ----print_main_results_DQ2, class.source='fold-hide'-------------------------
paste0('DQ2: ', round(copls.da@DQ2[position], 4))

## ----extract_VIP, class.source='fold-hide'------------------------------------
# Compute the VIP
VIP <- copls.da@VIP

# Multiply VIP * sign(loadings for predictive component)
VIP_plot <- lapply(X = 1:nbrBlocs,
                   FUN = function(X){
                       signe_loadings <- sign(copls.da@loadings[[X]][, "p_1"])
                       result <- VIP[[X]][, "p_1"]*signe_loadings
                       return(sort(result, decreasing = TRUE))})
names(VIP_plot) <- BlockNames

## ----plot_VIP, class.source='fold-hide'---------------------------------------
# Metabo data
ggplot(data = data.frame(
  "variables" = factor(names(VIP_plot[[1]]),
                       levels=names(VIP_plot[[1]])[order(abs(VIP_plot[[1]]), 
                                                         decreasing=T)]), 
  "valeur" = VIP_plot[[1]]), 
  aes(x = variables, y = valeur)) +
  geom_bar(stat = "identity") +
  labs(title = paste0("Barplot of ", names(VIP_plot)[1])) +
  xlab("Predictive variables") +
  ylab("VIP x loading sign") +
  theme_graphs +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) 

# Microarray data
ggplot(data = data.frame(
  "variables" = factor(names(VIP_plot[[2]]),
                       levels=names(VIP_plot[[2]])[order(abs(VIP_plot[[2]]), 
                                                         decreasing=T)]), 
  "valeur" = VIP_plot[[2]]), 
  aes(x = variables, y = valeur)) +
  geom_bar(stat = "identity") +
  labs(title = paste0("Barplot of ", names(VIP_plot)[2])) +
  xlab("Predictive variables") +
  ylab("VIP x loading sign") +
  theme_graphs +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) 

# Proteo data
ggplot(data = data.frame(
  "variables" = factor(names(VIP_plot[[3]]),
                       levels=names(VIP_plot[[3]])[order(abs(VIP_plot[[3]]), 
                                                         decreasing=T)]), 
  "valeur" = VIP_plot[[3]]), 
  aes(x = variables, y = valeur)) +
  geom_bar(stat = "identity") +
  labs(title = paste0("Barplot of ", names(VIP_plot)[3])) +
  xlab("Predictive variables") +
  ylab("VIP x loading sign") +
  theme_graphs +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) 

## ----ggplot_score_data, class.source='fold-hide', warning=FALSE---------------
ggplot(data = data.frame("p_1" = copls.da@scores[, "p_1"],
                                  "o_1" = copls.da@scores[, "o_1"],
                                  "Labs" = as.matrix(unlist(demo_3_Omics$ObsNames[, 1]))),
                aes(x = p_1, y = o_1, label = Labs, 
                    shape = Labs, colour = Labs)) +
  xlab("Predictive component") +
  ylab("Orthogonal component") +
  ggtitle("ConsensusOPLS Score plot")+
  geom_point(size = 2.5) + 
  geom_text_repel(size = 4, show.legend = FALSE) + 
  theme_graphs+
  scale_color_manual(values = c("#7F3C8D", "#11A579"))

## ----ggplot_data_pred_compo, class.source='fold-hide'-------------------------
ggplot(
  data = data.frame("Values" = copls.da@blockContribution[, "p_1"],
                    "Blocks" = as.factor(labels(demo_3_Omics[1:nbrBlocs]))),
  aes(x = Blocks, y = Values,
      fill = Blocks, labels = Values)) +
  geom_bar(stat = 'identity') + 
  geom_text(aes(label=round(Values, 2)), vjust=-0.3) +
  ggtitle("Block contributions to the predictive component")+
  xlab("Data blocks") +
  ylab("Weight") +
  theme_graphs +
  scale_color_manual(values = c("#1B9E77", "#D95F02", "#7570B3"))

## ----ggplot_data_1st_ortho_compo, class.source='fold-hide'--------------------
ggplot(
  data = 
    data.frame("Values" = copls.da@blockContribution[, "o_1"],
               "Blocks" = as.factor(labels(demo_3_Omics[1:nbrBlocs]))),
  aes(x = Blocks, y = Values,
      fill = Blocks, labels = Values)) +
  geom_bar(stat = 'identity') + 
  geom_text(aes(label=round(Values, 2)), vjust=-0.3) +
  ggtitle("Block contributions to the first orthogonal component") +
  xlab("Data blocks") +
  ylab("Weight") +
  theme_graphs +
  scale_color_manual(values = c("#1B9E77", "#D95F02", "#7570B3"))

## ----ggplot_data_pred_ortho, message = FALSE, class.source='fold-hide'--------
data_two_plots <- data.frame("Values" = copls.da@blockContribution[, "p_1"],
                             "Type" = "Pred",
                             "Blocks" = labels(demo_3_Omics[1:nbrBlocs]))
data_two_plots <- data.frame("Values" = c(data_two_plots$Values,
                                          copls.da@blockContribution[, "o_1"]),
                             "Type" = c(data_two_plots$Type,
                                        rep("Ortho", times = length(copls.da@blockContribution[, "o_1"]))),
                             "Blocks" = c(data_two_plots$Blocks,
                                          labels(demo_3_Omics[1:nbrBlocs])))

ggplot(data = data_two_plots,
       aes(x = factor(Type), 
           y = Values, 
           fill = factor(Type))) +
  geom_bar(stat = 'identity') + 
  ggtitle("Block contributions to each component")+
  geom_text(aes(label=round(Values, 2)), vjust=-0.3) +
  xlab("Data blocks") +
  ylab("Weight") +
  facet_wrap(. ~ Blocks)+
  theme_graphs+
  scale_fill_discrete(name = "Component")+
  scale_fill_manual(values = c("#7F3C8D", "#11A579"))

## ----plot_bloc_PredVSOrtho_bis, message = FALSE, class.source='fold-hide'-----
ggplot(data = data_two_plots,
                           aes(x = Blocks, 
                               y = Values, 
                               fill = Blocks)) +
    geom_bar(stat = 'identity') +
    geom_text(aes(label=round(Values, 2)), vjust=-0.3) +
    ggtitle("Block contributions to each component") +
    xlab("Components") +
    ylab("Weight") +
    facet_wrap(. ~ factor(Type)) +
    theme_graphs +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          plot.title = element_text(hjust = 0.5, 
                                    margin = margin(t = 5, r = 0, b = 0, l = 100))) +
    scale_fill_manual(values = c("#1B9E77", "#D95F02", "#7570B3"))

## ----ggplot_data_pred_vs_ortho, message = FALSE, warning = FALSE, class.source='fold-hide'----
ggplot(data = data.frame("Pred" = copls.da@blockContribution[, "p_1"],
                         "Ortho" = copls.da@blockContribution[, "o_1"],
                         "Labels" = labels(demo_3_Omics[1:nbrBlocs])),
       aes(x = Pred, y = Ortho, label = Labels, 
           shape = Labels, colour = Labels)) +
  xlab("Predictive component") +
  ylab("Orthogonal component") +
  ggtitle("Block contributions predictive vs. orthogonal") +
  geom_point(size = 2.5) + 
  geom_text_repel(size = 4, show.legend = FALSE) + 
  theme_graphs +
  scale_color_manual(values = c("#1B9E77", "#D95F02", "#7570B3"))

## ----create_data_loadings-----------------------------------------------------
loadings <- copls.da@loadings
data_loads <- sapply(X = 1:nbrBlocs,
                     FUN = function(X){
                       data.frame("Pred" = 
                                    loadings[[X]][, grep(pattern = "p_",
                                                               x = colnames(loadings[[X]]),
                                                               fixed = TRUE)],
                                  "Ortho" = 
                                    loadings[[X]][, grep(pattern = "o_",
                                                               x = colnames(loadings[[X]]),
                                                               fixed = TRUE)],
                                  "Labels" = labels(demo_3_Omics[1:nbrBlocs])[[X]])
                     })
data_loads <- as.data.frame(data_loads)

## ----ggplot_data_loadings, class.source='fold-hide'---------------------------
ggplot() +
  geom_point(data = as.data.frame(data_loads$V1),
             aes(x = Pred, y = Ortho, colour = Labels), 
             size = 2.5, alpha = 0.5) + 
  geom_point(data = as.data.frame(data_loads$V2),
             aes(x = Pred, y = Ortho, colour = Labels),
             size = 2.5, alpha = 0.5) +
  geom_point(data = as.data.frame(data_loads$V3),
             aes(x = Pred, y = Ortho, colour = Labels),
             size = 2.5, alpha = 0.5) +
  xlab("Predictive component") +
  ylab("Orthogonal component") +
  ggtitle("Loadings plot on first orthogonal and predictive component")+
  theme_graphs+
  scale_color_manual(values = c("#1B9E77", "#D95F02", "#7570B3"))

## ----create_data_loadings_VIP-------------------------------------------------
loadings <- do.call(rbind.data.frame, copls.da@loadings)
loadings$block <- do.call(c, lapply(names(copls.da@loadings), function(x) 
    rep(x, nrow(copls.da@loadings[[x]]))))
loadings$variable <- gsub(paste(paste0(names(copls.da@loadings), '.'), 
                                collapse='|'), '', 
                          rownames(loadings))

VIP <- do.call(rbind.data.frame, copls.da@VIP)
VIP$block <- do.call(c, lapply(names(copls.da@VIP), function(x) 
    rep(x, nrow(copls.da@VIP[[x]]))))
VIP$variable <- gsub(paste(paste0(names(copls.da@VIP), '.'), 
                                collapse='|'), '', 
                          rownames(VIP))
        
loadings_VIP <- merge(x = loadings[, c("p_1", "variable")], 
                      y = VIP[, c("p_1", "variable")], 
                      by = "variable", all = TRUE)
colnames(loadings_VIP) <- c("variable", "loadings", "VIP")
loadings_VIP <- merge(x = loadings_VIP, 
                      y = loadings[, c("block", "variable")], 
                      by = "variable", all = TRUE)
loadings_VIP$label <- ifelse(loadings_VIP$VIP > 1, loadings_VIP$variable, NA)

## ----ggplot_data_loadings_VIP, class.source='fold-hide'-----------------------
ggplot(data = loadings_VIP,
       aes(x=loadings, y=VIP, col=block, label = label)) +
  geom_point(size = 2.5, alpha = 0.5) + 
  xlab("Predictive component") +
  ylab("Variable Importance in Projection") +
  ggtitle("VIP versus loadings on predictive components")+
  theme_graphs+
  scale_color_manual(values = c("#1B9E77", "#D95F02", "#7570B3"))

## ----run_permutations, warning=FALSE------------------------------------------
PermRes <- copls.da@permStats

## ----plot_R2_perm, class.source='fold-hide'-----------------------------------
ggplot(data = data.frame("R2Yperm" = PermRes$R2Y),
       aes(x = R2Yperm)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30,
                 color="grey", fill="grey") +
  geom_density(color = "blue", linewidth = 0.5) +
  geom_vline(aes(xintercept=PermRes$R2Y[1]), 
             color="blue", linetype="dashed", size=1) +
  xlab("R2 values") +
  ylab("Frequency") +
  ggtitle("R2 Permutation test")+
  theme_graphs

## ----plot_Q2_perm, message = FALSE, class.source='fold-hide'------------------
ggplot(data = data.frame("Q2Yperm" = PermRes$Q2Y),
       aes(x = Q2Yperm)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30,
                 color="grey", fill="grey") +
  geom_density(color = "blue", size = 0.5) +
  geom_vline(aes(xintercept=PermRes$Q2Y[1]), 
             color="blue", linetype="dashed", size=1) +
  xlab("Q2 values") +
  ylab("Frequency") +
  ggtitle("Q2 Permutation test")+
  theme_graphs

## ----plot_DQ2_perm, message = FALSE, class.source='fold-hide'-----------------
ggplot(data = data.frame("DQ2Yperm" = PermRes$DQ2Y),
       aes(x = DQ2Yperm)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30,
                 color="grey", fill="grey") +
  geom_density(color = "blue", size = 0.5) +
  geom_vline(aes(xintercept=PermRes$DQ2Y[1]), 
             color="blue", linetype="dashed", size=1) +
  xlab("DQ2 values") +
  ylab("Frequency") +
  ggtitle("DQ2 Permutation test")+
  theme_graphs

## ----reproducibility----------------------------------------------------------
sessionInfo()

