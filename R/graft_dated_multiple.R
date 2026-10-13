#' @export

graft_dated_multiple_complex <- function(existing_species_tree, addition_species_tree, run, existing_species_position, addition_species_position){ 

###get information###
#####################

###get the species###
existing_species_tree_species_only <- paste(unlist(lapply(lapply(lapply(existing_species_tree$tip.label, strsplit, "_"), "[[", 1),"[",existing_species_position[[1]])), unlist(lapply(lapply(lapply(existing_species_tree$tip.label, strsplit, "_"), "[[", 1),"[",existing_species_position[[2]])), sep="_")
addition_species_tree_species_only <- paste(unlist(lapply(lapply(lapply(addition_species_tree$tip.label, strsplit, "_"), "[[", 1),"[",addition_species_position[[1]])), unlist(lapply(lapply(lapply(addition_species_tree$tip.label, strsplit, "_"), "[[", 1),"[",addition_species_position[[2]])), sep="_")

###figure out the graft positions###
####################################
if (length(which(addition_species_tree_species_only %in% existing_species_tree_species_only)) > 1){ ### if the tree already has more than one of the species of what is being added
if (round(max(node.depth.edgelength(extract.clade(existing_species_tree, findMRCA(existing_species_tree, which(existing_species_tree_species_only %in% addition_species_tree_species_only), "node")))),4) == round(max(node.depth.edgelength(addition_species_tree)),4)){ ### if its calibrated to the crown of what it overlaps with
mrcas <- extract.clade(existing_species_tree, findMRCA(existing_species_tree, which(existing_species_tree_species_only %in% addition_species_tree_species_only), "node"))$tip.label
position <- 0
} else { ### if its not grafted onto a crown node - it could be grafting within a genus for example
if (round(max(node.depth.edgelength(extract.clade(existing_species_tree, existing_species_tree[[1]][,1][which(existing_species_tree[[1]][,2] == findMRCA(existing_species_tree, which(existing_species_tree_species_only %in% addition_species_tree_species_only), "node"))]))),4) == round(max(node.depth.edgelength(addition_species_tree)),4)){
mrcas <- extract.clade(existing_species_tree, findMRCA(existing_species_tree, which(existing_species_tree_species_only %in% addition_species_tree_species_only), "node"))$tip.label
position <- (existing_species_tree$edge.length[which(existing_species_tree[[1]][,2] == findMRCA(existing_species_tree, which(existing_species_tree_species_only %in% addition_species_tree_species_only), "node"))] - 0.0001)   
} else {   
mrcas_init <- extract.clade(existing_species_tree, findMRCA(existing_species_tree, which(existing_species_tree_species_only %in% addition_species_tree_species_only), "node"))$tip.label  ### removing the entiure genus because relationship in relation to additional things in genus is unknown - but could make same argument for initial step
mrcas <- mrcas_init[[1]]
existing_species_tree <- drop.tip(existing_species_tree, mrcas_init[-1])
position <- max(node.depth.edgelength(addition_species_tree))-0.0001
keep_searching <- 0
while (keep_searching == 0){
if (round(existing_species_tree$edge.length[which(existing_species_tree[[1]][,2] == which(existing_species_tree$tip.label == mrcas))], 4) >= position){
keep_searching <- 1
} else {
removal <- existing_species_tree[[1]][,1][which(existing_species_tree[[1]][,2] == which(existing_species_tree$tip.label == mrcas))]
existing_species_tree <- drop.tip(existing_species_tree, extract.clade(existing_species_tree, removal)$tip.label[-which(extract.clade(existing_species_tree, removal)$tip.label == mrcas)])
}
}
}
}
} else {
if (round(existing_species_tree$edge.length[which(existing_species_tree[[1]][,2] == which(existing_species_tree_species_only %in% addition_species_tree_species_only))], 4) == round(max(node.depth.edgelength(addition_species_tree)),4)){ ### why would this ever not be satisfied
mrcas <- existing_species_tree$tip.label[which(existing_species_tree_species_only %in% addition_species_tree_species_only)]
position <- max(node.depth.edgelength(addition_species_tree))-0.0001
} else {
mrcas <- existing_species_tree$tip.label[which(existing_species_tree_species_only %in% addition_species_tree_species_only)]
position <- max(node.depth.edgelength(addition_species_tree))-0.0001
keep_searching <- 0
while (keep_searching == 0){
if (round(existing_species_tree$edge.length[which(existing_species_tree[[1]][,2] == which(existing_species_tree$tip.label == mrcas))], 4) >= position){
keep_searching <- 1
} else {
stem_clade <- extract.clade(existing_species_tree, existing_species_tree[[1]][,1][which(existing_species_tree[[1]][,2] == which(existing_species_tree$tip.label == mrcas))])
removal <- stem_clade$tip.label[-which(stem_clade$tip.label == mrcas)]
existing_species_tree <- drop.tip(existing_species_tree, removal)
}
}
}
}
###do the actual grafting###
############################
big_final <- existing_species_tree
if (length(mrcas) > 1){
if (position <= big_final$edge.length[which(big_final[[1]][,2] == findMRCA(big_final, mrcas, "node"))]){
big_final <- bind.tree(big_final, addition_species_tree, findMRCA(big_final, mrcas, "node"), position=position)
} else {
big_final <- bind.tree(big_final, addition_species_tree, findMRCA(big_final, mrcas, "node"), position=position)
}
} else {
if (position <= big_final$edge.length[which(big_final[[1]][,2] == which(big_final$tip.label == mrcas))]){
big_final <- bind.tree(big_final, addition_species_tree, which(big_final$tip.label == mrcas), position=position)
} else {
big_final <- bind.tree(big_final, addition_species_tree, which(big_final$tip.label == mrcas), position=position)
}
}
big_final <- drop.tip(big_final, mrcas)
write.tree(big_final, paste("NEWTREE_run_", run, "_bootstrap_tree_to_date.tre_dated", sep=""))
}



