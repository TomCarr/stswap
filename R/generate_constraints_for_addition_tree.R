#' @export

##########################################
###define main function for constraints###
##########################################
generate_constraints_for_addition_tree <- function(species_trees, 
backbone_trees, 
clean_backbone_tree,
addition_trees, 
addition_trees_best, 
n_sites, 
refer_to_backbone,
do_dating,
clade_name_pos,
continual_search,
existing_species_position,
backbone_order_position,
backbone_species_position,
addition_species_position,
frequent_root,
outgroup_is_pruned){

clean_backbone <- read.tree(clean_backbone_tree)

for (z in 1:length(list.files(species_trees))){ ### generate a new dated species tree (for the specific clade) for each overall dated species tree  
done <- 0
while (done == 0){
if (continual_search == FALSE){
done <- 1
}
###trees###
###########
existing_species_tree <- read.tree(paste(species_trees, list.files(species_trees)[[z]], sep=""))
run <- unlist(strsplit(list.files(species_trees)[[z]], split=".tre_dated"))[[1]] ### get the run number for that species tree
###
backbone_tree <-  read.tree(paste(backbone_trees, list.files(backbone_trees, pattern = paste("run_", run, "_", sep="")), sep="")) ### get the corresponding backbone tree
backbone_tree <- drop.tip(backbone_tree, which(backbone_tree$tip.label %in% clean_backbone$tip.label == FALSE))
backbone_tree_bootstrap_number <- unlist(strsplit(list.files(backbone_trees, pattern = paste("run_", run, "_", sep=""))[[1]], split="_"))[[2]]
###
if (length(list.files(addition_trees)) == 1){ ### if its one file
addition_species_tree <- read.tree(paste(addition_trees, list.files(addition_trees)[[1]], sep=""))
} else { ### or if its lost of files
addition_species_tree <- lapply(paste(addition_trees, list.files(addition_trees), sep=""), read.tree)
}
if (class(addition_species_tree) != "list"){ ### better in a list for the next steps
addition_species_tree <- list(addition_species_tree)
}
species_sample <- sample(seq(1,length(addition_species_tree),1),1)
addition_species_tree <- addition_species_tree[[species_sample]]

###get clade name###
####################
clade <- paste(unlist(strsplit(addition_trees,split=""))[clade_name_pos],collapse="")

#############################
###deal with the outgroups###
#############################
best_tree <- read.tree(paste(addition_trees_best, list.files(addition_trees_best), sep="")[[1]])
if (outgroup_is_pruned == FALSE){
good_to_go <- "no"    
###this first step gets rid of the outgroups###
###############################################
side_one <- extract.clade(addition_species_tree, length(addition_species_tree$tip.label) + 2)$tip.label
side_two <- addition_species_tree$tip.label[-which(addition_species_tree$tip.label %in% side_one == TRUE)]
if (length(side_one) > length(side_two)){
addition_species_tree <- drop.tip(addition_species_tree, side_two)
} else {
addition_species_tree <- drop.tip(addition_species_tree, side_one)
}
###now check if the extracted ingroup is the correct size###
############################################################
best_side_one <- extract.clade(best_tree, length(best_tree$tip.label) + 2)$tip.label
best_side_two <- best_tree$tip.label[-which(best_tree$tip.label %in% best_side_one == TRUE)]
if (length(best_side_one) > length(best_side_two)){
best_extraction <- drop.tip(best_tree, best_side_two)
} else {
best_extraction <- drop.tip(best_tree, best_side_one)
}
if (length(which(best_extraction$tip.label %in% addition_species_tree$tip.label == FALSE)) == 0){
if (length(which(addition_species_tree$tip.label %in% best_extraction$tip.label == FALSE)) == 0){
good_to_go <- "yes"
}
}
} else {
good_to_go <- "yes"
}
#####################
#####################
#####################

if (good_to_go == "yes"){
###get species information###
#############################
backbone_tree_species <- paste(unlist(lapply(lapply(lapply(backbone_tree$tip.label, strsplit, "_"), "[[", 1),"[[",backbone_species_position[[1]])), unlist(lapply(lapply(lapply(backbone_tree$tip.label, strsplit, "_"), "[[", 1),"[[",backbone_species_position[[2]])), sep="_")
addition_species_tree_species <-  paste(unlist(lapply(lapply(lapply(addition_species_tree$tip.label, strsplit, "_"), "[[", 1),"[[",addition_species_position[[1]])), unlist(lapply(lapply(lapply(addition_species_tree$tip.label, strsplit, "_"), "[[", 1),"[[",addition_species_position[[2]])), sep="_")
existing_species_tree_species <-  paste(unlist(lapply(lapply(lapply(existing_species_tree$tip.label, strsplit, "_"), "[[", 1),"[[",existing_species_position[[1]])), unlist(lapply(lapply(lapply(existing_species_tree$tip.label, strsplit, "_"), "[[", 1),"[[",existing_species_position[[2]])), sep="_")

###if the species tree to add is represented in the backbone tree###
####################################################################
if ((length(which(addition_species_tree_species %in% backbone_tree_species)) > 0) & (refer_to_backbone == TRUE)){
generate_tree_pl_inputs_complex(backbone_tree, addition_species_tree, existing_species_tree, run, backbone_tree_bootstrap_number, backbone_species_position, addition_species_position, existing_species_position, do_dating, clade, species_sample,frequent_root,backbone_order_position)
} else {
### or if it is not ###
generate_tree_pl_input_simple(existing_species_tree,addition_species_tree, run,backbone_tree_bootstrap_number, existing_species_tree_species, addition_species_position,do_dating,clade, species_sample)
}
if (length(list.files(paste(clade, "_to_date/", sep=""), pattern=paste("run_", run, "_bootstrap_tree_to_date.tre_dated", sep=""))) > 0){
done <- 1
}
}
}
}
}
