#' @export

graft_into_new_tree <- function(species_trees, addition_trees_dated, existing_species_position, addition_species_position, type){

if (type == "dated"){
for (z in 1:length(list.files(species_trees))){
existing_species_tree <- read.tree(paste(species_trees, list.files(species_trees)[[z]], sep=""))
run <- unlist(strsplit(list.files(species_trees)[[z]], split=".tre_dated"))[[1]]
if (list.files(addition_trees_dated, pattern=paste("_run_", run, "_bootstrap_tree_to_date.tre_dated", sep=""))[[1]] %in% list.files(addition_trees_dated)){
addition_species_tree <- read.tree(paste(addition_trees_dated, list.files(addition_trees_dated, pattern=paste("_run_", run, "_bootstrap_tree_to_date.tre_dated", sep=""))[[1]], sep=""))
graft_dated_multiple_complex(existing_species_tree, addition_species_tree, run, existing_species_position, addition_species_position)
message(paste("done tree ", z, sep=""))
}
}
}

if (type == "undated"){
for (z in 1:length(list.files(species_trees))){
existing_species_tree <- read.tree(paste(species_trees, list.files(species_trees)[[z]], sep=""))
run <- unlist(strsplit(list.files(species_trees)[[z]], split=".tre"))[[1]]
if (list.files(addition_trees_dated, pattern=paste("_run_", run, "_bootstrap_tree_to_date.tre", sep=""))[[1]] %in% list.files(addition_trees_dated)){
addition_species_tree <- read.tree(paste(addition_trees_dated, list.files(addition_trees_dated, pattern=paste("_run_", run, "_bootstrap_tree_to_date.tre", sep=""))[[1]], sep=""))
graft_undated_multiple_complex(existing_species_tree, addition_species_tree, run, existing_species_position, addition_species_position)
message(paste("done tree ", z, sep=""))
}
}
}
}

