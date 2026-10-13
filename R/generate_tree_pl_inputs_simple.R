#' @export

generate_tree_pl_input_simple <- function(existing_species_tree,
addition_species_tree,
run, 
backbone_tree_bootstrap_number, 
existing_species_tree_species, 
addition_tree_species_poss, 
do_dating, 
clade, 
species_sample){

if (length(addition_species_tree) > 0){
###work through on the species tree###
######################################
message(paste("working on config file for ", clade, sep=""))
if (dir.exists(paste(clade, "_to_date/", sep="")) == FALSE){
dir.create(paste(clade, "_to_date/", sep=""))
}

species_to_date <- addition_species_tree
ages <- vector(mode="numeric", length=0)
constraints <- vector("list", length=0)
constraints_mini <- vector("list", length=0)
age_addition <- vector(mode="numeric", length=0)
constraint_addition <- vector(mode="numeric", length=0)
species_to_date <- addition_species_tree

###now generate the constraint file###
######################################
side_one <- extract.clade(species_to_date, (length(species_to_date$tip.label)+2))
side_one_species <- paste(unlist(lapply(lapply(lapply(side_one$tip.label, strsplit, "_"), "[[", 1), "[", addition_tree_species_poss[[1]])), unlist(lapply(lapply(lapply(side_one$tip.label, strsplit, "_"), "[[", 1), "[", addition_tree_species_poss[[2]])), sep="_")
side_two <- species_to_date$tip.label[which(species_to_date$tip.label %in% side_one$tip.label == FALSE)]
side_two_species <- paste(unlist(lapply(lapply(lapply(side_two, strsplit, "_"), "[[", 1), "[", addition_tree_species_poss[[1]])), unlist(lapply(lapply(lapply(side_two, strsplit, "_"), "[[", 1), "[", addition_tree_species_poss[[2]])), sep="_")

if (dir.exists(paste(clade, "_to_date/", sep="")) == FALSE){
dir.create(paste(clade, "_to_date/", sep=""))
}
write(paste("treefile = ", clade, "_to_date/backbone_", backbone_tree_bootstrap_number, "_bs_", species_sample, "_run_", run, "_bootstrap_tree_to_date.tre", sep=""), paste(clade, "_to_date/", run, ".treePL", sep=""), append=TRUE)
write("numsites=10000", paste(clade, "_to_date/", run, ".treePL", sep=""), append=TRUE)
if ((length(which(side_one_species %in% existing_species_tree_species)) > 0) & (length(which(side_two_species %in% existing_species_tree_species)) > 0)){
######
write(paste("mrca = crown", " ", species_to_date$tip.label[[1]], " ", species_to_date$tip.label[length(species_to_date$tip.label)], sep=""), paste(clade, "_to_date/", run, ".treePL", sep=""), append=TRUE)
write(paste("max = crown", " ", max(node.depth.edgelength(extract.clade(existing_species_tree, findMRCA(existing_species_tree, which(existing_species_tree_species %in% c(side_one_species, side_two_species)), "node")))), sep=""), paste(clade, "_to_date/", run, ".treePL", sep=""), append=TRUE)
write(paste("min = crown", " ", max(node.depth.edgelength(extract.clade(existing_species_tree, findMRCA(existing_species_tree, which(existing_species_tree_species %in% c(side_one_species, side_two_species)), "node")))), sep=""), paste(clade, "_to_date/", run, ".treePL", sep=""), append=TRUE)
} else {
if (length(which(c(side_one_species, side_two_species) %in% existing_species_tree_species)) > 1){
write(paste("mrca = crown", " ", species_to_date$tip.label[[1]], " ", species_to_date$tip.label[length(species_to_date$tip.label)], sep=""), paste(clade, "_to_date/", run, ".treePL", sep=""), append=TRUE)
write(paste("max = crown", " ", max(node.depth.edgelength(extract.clade(existing_species_tree, existing_species_tree[[1]][,1][[which(existing_species_tree[[1]][,2] == findMRCA(existing_species_tree, which(existing_species_tree_species %in% c(side_one_species, side_two_species)), "node"))]]))), sep=""), paste(clade, "_to_date/", run, ".treePL", sep=""), append=TRUE)
}
if (length(which(c(side_one_species, side_two_species) %in% existing_species_tree_species)) == 1){
write(paste("mrca = crown", " ", species_to_date$tip.label[[1]], " ", species_to_date$tip.label[length(species_to_date$tip.label)], sep=""), paste(clade, "_to_date/", run, ".treePL", sep=""), append=TRUE)
write(paste("max = crown", " ", max(node.depth.edgelength(extract.clade(existing_species_tree, existing_species_tree[[1]][,1][[which(existing_species_tree[[1]][,2] == which(existing_species_tree_species %in% c(side_one_species, side_two_species)), "node")]]))), sep=""), paste(clade, "_to_date/", run, ".treePL", sep=""), append=TRUE)
}
}
write("pliter=2", paste(clade, "_to_date/", run, ".treePL", sep=""), append=TRUE)
write("smooth=10", paste(clade, "_to_date/", run, ".treePL", sep=""), append=TRUE)
write("plsimaniter=10000", paste(clade, "_to_date/", run, ".treePL", sep=""), append=TRUE)
write("nthreads=8", paste(clade, "_to_date/", run, ".treePL", sep=""), append=TRUE)
write("thorough", paste(clade, "_to_date/", run, ".treePL", sep=""), append=TRUE)
write(paste("outfile = ", paste(clade, "_to_date/backbone_", backbone_tree_bootstrap_number, "_bs_", species_sample, "_run_", run, "_bootstrap_tree_to_date.tre", "_dated", sep="")), paste(clade, "_to_date/", run, ".treePL", sep=""), append=TRUE)
write("prime", paste(clade, "_to_date/", run, ".treePL", sep=""), append=TRUE)
write.tree(species_to_date, paste(clade, "_to_date/backbone_", backbone_tree_bootstrap_number, "_bs_", species_sample, "_run_", run, "_bootstrap_tree_to_date.tre", sep=""))    
if (do_dating == TRUE){
run_tree_pl(paste(clade, "_to_date/", sep=""), run)
}
}
}
