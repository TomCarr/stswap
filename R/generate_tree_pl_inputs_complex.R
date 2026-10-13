#' @export

generate_tree_pl_inputs_complex <- function(backbone_tree,
addition_species_tree, 
existing_species_tree,
run, 
backbone_tree_bootstrap_number, 
backbone_species_position, 
addition_species_position,
existing_species_species_position, 
do_dating,
clade, 
species_sample,
frequent_root,
backbone_order_position){


if (length(addition_species_tree) > 0){

###get backbone variables###
############################
backbone_species_only <- paste(unlist(lapply(lapply(lapply(backbone_tree$tip.label, strsplit, "_"), "[[", 1),"[",backbone_species_position[[1]])), unlist(lapply(lapply(lapply(backbone_tree$tip.label, strsplit, "_"), "[[", 1),"[",backbone_species_position[[2]])), sep="_") 
backbone_genus_only <- unlist(lapply(lapply(lapply(backbone_tree$tip.label, strsplit, "_"), "[[", 1),"[",backbone_species_position[[1]]))
if (backbone_order_position != FALSE){
backbone_order_only <- unlist(lapply(lapply(lapply(backbone_tree$tip.label, strsplit, "_"), "[[", 1),"[",backbone_order_position))
}

###work through on the species tree###
######################################
message(paste("working on config file for ", clade, sep=""))
if (dir.exists(paste(clade, "_to_date/", sep="")) == FALSE){
dir.create(paste(clade, "_to_date/", sep=""))
}
already_estimated <- "no"

###do check on whether we already have the subclade constrained to that backbone###
###################################################################################
if ((length(list.files(paste(clade, "_to_date/", sep=""),pattern=paste("backbone_", backbone_tree_bootstrap_number, "_bs_", species_sample, "_run_", sep=""))) > 0) & (round(max(node.depth.edgelength(backbone_tree)), 4) == round(frequent_root, 4))){
run_numbers <- vector(mode="numeric", length=0)
for (b in 1:length(list.files(paste(clade, "_to_date/", sep=""),pattern=paste("backbone_", backbone_tree_bootstrap_number, "_bs_", species_sample, "_run_", sep="")))){
run_numbers <- append(run_numbers, unlist(strsplit(list.files(paste(clade, "_to_date/", sep=""),pattern=paste("backbone_", backbone_tree_bootstrap_number, "_bs_", species_sample, "_run", sep="")), split="_"))[[6]])
}
run_numbers <- unique(run_numbers)
good_run_numbers <- vector(mode="numeric", length=0)
for (b in 1:length(run_numbers)){
other_bs_backbone_trees_constraining_this_clade_runs_ages <- max(node.depth.edgelength(read.tree(paste(backbone_trees, list.files(backbone_trees, pattern = paste("run_", run_numbers[[b]], "_bootstrap_tree_to_date.tre_dated", sep=""))[[1]], sep=""))))
if (round(other_bs_backbone_trees_constraining_this_clade_runs_ages, 4) == round(frequent_root, 4)){
good_run_numbers <- append(good_run_numbers, run_numbers[[b]])
}
}
if (length(good_run_numbers) > 0){
good_other_one <- good_run_numbers[[1]]
files <- list.files(paste(clade, "_to_date/", sep=""),pattern=paste("backbone_", backbone_tree_bootstrap_number, "_bs_", species_sample, "_run_", good_other_one, sep=""))
file_select <- grepl("_dated", files)
if (length(which(file_select == TRUE)) > 0){
file <- files[which(file_select == TRUE)[[1]]]
file_undated <- files[[1]]
file.copy(paste(clade, "_to_date/", file, sep=""), paste(clade, "_to_date/", "backbone_", backbone_tree_bootstrap_number, "_bs_", species_sample, "_run_", run, "_bootstrap_tree_to_date.tre_dated", sep=""))
file.copy(paste(clade, "_to_date/", file_undated, sep=""), paste(clade, "_to_date/", "backbone_", backbone_tree_bootstrap_number, "_bs_", species_sample, "_run_", run, "_bootstrap_tree_to_date.tre", sep=""))
keep_searching <- "no"
already_estimated <- "yes" 
do_dating <<- "no"
}
}
}

species_to_date <- addition_species_tree
ages <- vector(mode="numeric", length=0)
constraints <- vector("list", length=0)
constraints_mini <- vector("list", length=0)
age_addition <- vector(mode="numeric", length=0)
constraint_addition <- vector(mode="numeric", length=0)

###checking number of orders in species level tree###
#####################################################
if (already_estimated == "no"){

message("building new config file")
orders <- vector(mode="numeric", length=0)
species_genus <- unlist(lapply(lapply(lapply(species_to_date$tip.label, strsplit, "_"), "[[", 1),"[",addition_species_position[[1]])) ### threshold for matching now elevated to genus

################################
###get the actual constraints###
################################

###define key species tree variables for deriving appropriate constraints###
############################################################################
species_species <- paste(unlist(lapply(lapply(lapply(species_to_date$tip.label, strsplit, "_"), "[[", 1),"[",addition_species_position[[1]])), unlist(lapply(lapply(lapply(species_to_date$tip.label, strsplit, "_"), "[[", 1),"[",addition_species_position[[2]])), sep="_")
species_genus_only <- unlist(lapply(lapply(lapply(species_to_date$tip.label, strsplit, "_"), "[[", 1),"[",addition_species_position[[1]]))
message("got species and genera in addition")

###work out which clade you are doing and get overlapping genera###
###################################################################
species_equiv <- drop.tip(species_to_date, which(species_genus_only %in% backbone_genus_only == FALSE)) 
species_equiv_species <- paste(unlist(lapply(lapply(lapply(species_equiv$tip.label, strsplit, "_"), "[[", 1),"[",addition_species_position[[1]])), unlist(lapply(lapply(lapply(species_equiv$tip.label, strsplit, "_"), "[[", 1),"[",addition_species_position[[2]])), sep="_")
message("got species equivalent to backbone")

###################################################################################################################################
if (length(which(backbone_genus_only %in% species_genus_only)) > 1){ ###have to add a provision if theres a single tip for a genus###
###################################################################################################################################
clade_backbone_tree <- extract.clade(backbone_tree, findMRCA(backbone_tree, which(backbone_genus_only %in% species_genus_only), "node"))
clade_backbone_tree <- drop.tip(clade_backbone_tree, which(unlist(lapply(lapply(lapply(clade_backbone_tree$tip.label, strsplit, "_"), "[[", 1),"[",backbone_species_position[[1]])) %in% species_genus_only == FALSE))
###the aim is to remove extra things from the clade backbone tree, unless they are the only thing that is sampled for that genus
if (length(which(paste(unlist(lapply(lapply(lapply(clade_backbone_tree$tip.label, strsplit, "_"), "[[", 1),"[",backbone_species_position[[1]])), unlist(lapply(lapply(lapply(clade_backbone_tree$tip.label, strsplit, "_"), "[[", 1),"[",backbone_species_position[[2]])), sep="_") %in% species_species == FALSE))>0){
message("clade backbone tree cleaning 1 start")
missing <- clade_backbone_tree$tip.label[which(paste(unlist(lapply(lapply(lapply(clade_backbone_tree$tip.label, strsplit, "_"), "[[", 1),"[",backbone_species_position[[1]])), unlist(lapply(lapply(lapply(clade_backbone_tree$tip.label, strsplit, "_"), "[[", 1),"[",backbone_species_position[[2]])), sep="_") %in% species_species == FALSE)]
message("got missing")
missing_genera <- unlist(lapply(lapply(lapply(missing, strsplit, "_"), "[[", 1),"[",backbone_species_position[[1]]))
message("got missing genera")
clade_backbone_tree_genera <- unlist(lapply(lapply(lapply(clade_backbone_tree$tip.label, strsplit, "_"), "[[", 1),"[",backbone_species_position[[1]]))
removal <- vector(mode="numeric", length=0)
removal_just_genera <- vector(mode="numeric", length=0)
for (i in 1:length(missing_genera)){
if (length(which(clade_backbone_tree_genera == missing_genera[[i]])) > 1){ # do nothing if the missing thing is the sole representative of a genus
removal <- append(removal, missing[[i]])
removal_just_genera <- append(removal_just_genera, missing_genera[[i]])
}
}
######
removal_from_removal <- vector(mode="numeric", length=0)
if (length(removal) > 0){
for (i in 1:length(unique(removal_just_genera))){
if (length(which(clade_backbone_tree_genera %in% unique(removal_just_genera)[[i]])) > 1){ ### fine just sends it down to restore it if its just one
if (length(which(clade_backbone_tree$tip.label %in% removal[which(removal_just_genera == unique(removal_just_genera)[[i]])])) == length(which(clade_backbone_tree_genera %in% unique(removal_just_genera)[[i]]))){ ### if all the representatives for the genus are in removal
if (length(extract.clade(clade_backbone_tree, findMRCA(clade_backbone_tree, which(clade_backbone_tree_genera == unique(removal_just_genera)[[i]]), "node"))$tip.label) == length(which(clade_backbone_tree_genera == unique(removal_just_genera)[[i]]))){ ###if they form a clade
removal_from_removal <- append(removal_from_removal, removal[which(removal_just_genera == unique(removal_just_genera)[[i]])[[1]]])
} else { ### removing things from removal if outside addition and in paftol - specifically the oldest, and will set as a max onstraint on the root node of the whole thing
removal_stem_ages <- vector(mode="numeric", length=0)
for (j in which(removal_just_genera %in% unique(removal_just_genera[[i]]))){ 
ancestor_run <- Ancestors(clade_backbone_tree, which(clade_backbone_tree$tip.label == removal[[j]]))
keep_rising <- "yes"
k=1
while(keep_rising == "yes"){
stem <- extract.clade(clade_backbone_tree, ancestor_run[[k]])$tip.label
stem_genera <- unlist(lapply(lapply(lapply(stem, strsplit, "_"), "[[", 1),"[",backbone_species_position[[1]]))
if (length(which(stem_genera %in% unique(removal_just_genera[[i]]))) != length(stem_genera)){
keep_rising <- "no"
removal_stem_ages <- append(removal_stem_ages, max(node.depth.edgelength(extract.clade(clade_backbone_tree, ancestor_run[[k]]))))
}
k=k+1
}
}
removal_from_removal <- append(removal_from_removal, removal[which(removal_just_genera %in% unique(removal_just_genera)[[i]])][which(removal_stem_ages == max(removal_stem_ages))]) # correct
}
}
} else {
removal_from_removal <- append(removal_from_removal, removal[which(removal_just_genera %in% unique(removal_just_genera)[[i]])])
}
}
}
if (length(removal_from_removal) > 0){
removal <- removal[-which(removal %in% removal_from_removal)]
}
######
if (length(removal) > 0){
clade_backbone_tree <- drop.tip(clade_backbone_tree, removal)
}
}

if (length(clade_backbone_tree$tip.label) > 1){

###get the constraints###
#########################
message("now figuring constraints")
for (i in (length(clade_backbone_tree$tip.label)+1):((length(clade_backbone_tree$tip.label)*2)-1)){ ### going through each node in the clade backbone tree
backbone_node_genera <- unlist(lapply(lapply(lapply(extract.clade(clade_backbone_tree, i)$tip.label, strsplit, "_"), "[[", 1), "[",backbone_species_position[[1]])) ### get genera that are descended from the node

###testing and constraining between genera###
#############################################
if (length(unique(backbone_node_genera)) > 1){ ### only do constraints between genera if there is more than one genus descended from the backbone node
if (length(which(unique(species_genus_only) %in% unique(backbone_node_genera))) > 1){ ###  and if the species tree has more than one of those genera
species_node_genera <- unlist(lapply(lapply(lapply(extract.clade(species_to_date, findMRCA(species_to_date, which(species_genus_only %in% unique(backbone_node_genera)), "node"))$tip.label, strsplit, "_"), "[[", 1), "[", 1)) ### get the descendant genera from the equivalent node in the species tree
backbone_matchers <- species_node_genera[which(species_node_genera %in% backbone_genus_only)] ### this is accross all backbone genera 
clade_matchers <- species_node_genera[which(species_node_genera %in% backbone_node_genera)] ### this is accross just the backbone node genera
if (length(backbone_matchers[-which(backbone_matchers %in% clade_matchers)]) > 0){ ### if its matching other genera within the mrce of the backbone node genera there is an incongruence
message("incongruence1")
} else {
constraints[[length(constraints)+1]] <- extract.clade(species_to_date, findMRCA(species_to_date, which(species_genus_only %in% unique(backbone_node_genera)), "node"))$tip.label 
ages <- append(ages, max(node.depth.edgelength(extract.clade(clade_backbone_tree, i))))
}
###testing and constraining within genera###
############################################
} 
} else { ### if there is only one genus
backbone_node_species <- paste(unlist(lapply(lapply(lapply(extract.clade(clade_backbone_tree, i)$tip.label, strsplit, "_"), "[[", 1), "[",backbone_species_position[[1]])), unlist(lapply(lapply(lapply(extract.clade(clade_backbone_tree, i)$tip.label, strsplit, "_"), "[[", 1), "[",backbone_species_position[[2]])), sep="_") ### get the species in the backbone tree
if (length(which(species_species %in% unique(backbone_node_species))) > 1){ ### if the species tree has more than one of the backbone node species
species_node_species <- extract.clade(species_to_date, findMRCA(species_to_date, which(species_species %in% unique(backbone_node_species)), "node"))$tip.label ### get the species descended from the mrca of the overlapping species
species_node_species <- paste(unlist(lapply(lapply(lapply(species_node_species, strsplit, "_"), "[[", 1),"[",1)), unlist(lapply(lapply(lapply(species_node_species, strsplit, "_"), "[[", 1),"[",addition_species_position[[2]])), sep="_")
backbone_matchers <- species_node_species[which(species_node_species %in% backbone_species_only)]
clade_matchers <- species_node_species[which(species_node_species %in% backbone_node_species)]
if (length(backbone_matchers[-which(backbone_matchers %in% clade_matchers)]) > 0){ ### if its matching other species in backbone, there is an incongruence
message("incongruence2")
} else { 
constraints[[length(constraints)+1]] <- extract.clade(species_to_date, findMRCA(species_to_date, which(species_species %in% unique(backbone_node_species)), "node"))$tip.label 
ages <- append(ages, max(node.depth.edgelength(extract.clade(clade_backbone_tree, i))))
}
species_node_genera <- unlist(lapply(lapply(lapply(extract.clade(species_to_date, findMRCA(species_to_date, which(species_genus_only %in% unique(backbone_node_genera)), "node"))$tip.label, strsplit, "_"), "[[", 1), "[", 1)) ### get the descendant genera from the equivalent node in the species tree
}
}
}

###add in another layer of constraints here### ### sets max age constraints where possible
############################################## ### seems to be monophyletic things in the backbone tree that aren't necessarily monophyletyic in the species level tree

backbone_clade_genera <- unlist(lapply(lapply(lapply(clade_backbone_tree$tip.label, strsplit, "_"), "[[", 1), "[",backbone_species_position[[1]])) ### get genera that are descended from the node
table_genera <- table(species_genus_only)
if (length(table_genera) > 1){
for (i in 1:length(table_genera)){
if (table_genera[[i]] > 1){ ### if there is more than one tip for the genus
if (names(table_genera)[[i]] %in% backbone_clade_genera){ ### if its in the backbone tree
if (findMRCA(species_to_date, which(species_genus_only == names(table_genera)[[i]]), "node") != length(species_to_date$tip.label)+1){ ###if its not the root node of the species tree
if (length(which(backbone_clade_genera == names(table_genera)[[i]])) > 1){ ###if there's more than one tip for the genus in the backbone tree
if (length(extract.clade(clade_backbone_tree, findMRCA(clade_backbone_tree, which(backbone_clade_genera %in% names(table_genera)[[i]]), "node"))$tip.label) == length(which(backbone_clade_genera %in% names(table_genera)[[i]]))){ ### if its monophyletic in the backbone clade tree 
if (findMRCA(clade_backbone_tree, which(backbone_clade_genera %in% names(table_genera)[[i]]), "node") != length(clade_backbone_tree$tip.label) + 1){ ### if we are not at the root node of the backbone clade tree
constraints[[length(constraints)+1]] <- c(extract.clade(species_to_date, findMRCA(species_to_date, which(species_genus_only == names(table_genera)[[i]]), "node"))$tip.label, "max") 
ages <- append(ages, max(node.depth.edgelength(extract.clade(clade_backbone_tree, clade_backbone_tree[[1]][,1][[which(clade_backbone_tree[[1]][,2] == findMRCA(clade_backbone_tree, which(backbone_clade_genera %in% names(table_genera)[[i]]), "node"))]]))))
}
}
} else {
if (clade_backbone_tree[[1]][,1][[which(clade_backbone_tree[[1]][,2] == which(backbone_clade_genera %in% names(table_genera)[[i]]))]] != length(clade_backbone_tree$tip.label) + 1){
constraints[[length(constraints)+1]] <- c(extract.clade(species_to_date, findMRCA(species_to_date, which(species_genus_only == names(table_genera)[[i]]), "node"))$tip.label, "max") 
ages <- append(ages, max(node.depth.edgelength(extract.clade(clade_backbone_tree, clade_backbone_tree[[1]][,1][[which(clade_backbone_tree[[1]][,2] == which(backbone_clade_genera %in% names(table_genera)[[i]]))]]))))
}
}
}
}
}
}
}

###additional max constraint search in comparison to existing species tree###
#############################################################################

existing_species_tree_genera <- unlist(lapply(lapply(lapply(existing_species_tree$tip.label, strsplit, "_"), "[[", 1),"[[",existing_species_species_position[[1]]))
addition_species_tree_genera <- unlist(lapply(lapply(lapply(addition_species_tree$tip.label, strsplit, "_"), "[[", 1),"[[",addition_species_position[[1]]))
additional_age <- max(node.depth.edgelength(extract.clade(existing_species_tree, existing_species_tree[[1]][,1][[which(existing_species_tree[[1]][,2] == findMRCA(existing_species_tree, which(existing_species_tree_genera %in% addition_species_tree_genera), "node"))]])))
if (length(which(ages > additional_age)) > 0){
ages <- append(ages, additional_age)
constraints[[length(constraints)+1]] <- c(species_to_date$tip.label, "max")
constraints_trimmer <- vector(mode="numeric", length=0)
if (length(which(ages > ages[length(ages)])) > 0){
constraints_trimmer <- which(ages > ages[length(ages)])
constraints <- constraints[-constraints_trimmer]
ages <- ages[-constraints_trimmer]
} 
}

###remove duplicated constraints - we have been through every node on the backbone equivalent clade###
####################################################################################################

###first filter removes any same constraints 
if (length(constraints) > 0){
constraints_trimmer <- vector(mode="numeric", length=0)
for (x in 1:length(constraints)){
if ("max" %in% constraints[[x]] == FALSE){
for (y in (1:length(constraints))[-x]){
if (setequal(constraints[[x]],constraints[[y]]) == TRUE){
ages_one <- ages[[x]]
ages_two <- ages[[y]]
if (ages_one > ages_two){
constraints_trimmer <- append(constraints_trimmer, y)
} else {
constraints_trimmer <- append(constraints_trimmer, x)
}
}
if ("max" %in% constraints[[y]]){
tester <- constraints[[y]][-length(constraints[[y]])]
if (setequal(constraints[[x]],tester) == TRUE){
constraints_trimmer <- append(constraints_trimmer, y)
}
}
}
}
###second filter removes any nested equal or older constraints
if ("max" %in% constraints[[x]]){
tester <- constraints[[x]][-length(constraints[[x]])]
decs <- getDescendants(species_to_date, findMRCA(species_to_date, tester, "node"), curr=NULL)
}
if ("max" %in% constraints[[x]]==FALSE){
tester <- constraints[[x]]
decs <- getDescendants(species_to_date, findMRCA(species_to_date, tester, "node"), curr=NULL)
}
for (y in (1:length(constraints))[-x]){
if ("max" %in% constraints[[y]] == FALSE){
constraints_node <- findMRCA(species_to_date, constraints[[y]], "node")
}
if ("max" %in% constraints[[y]] == TRUE){
constraints_node <- findMRCA(species_to_date, constraints[[y]][-length(constraints[[y]])], "node")
}
if (constraints_node %in% decs){
if (ages[[x]] <= ages[[y]]){
constraints_trimmer <- append(constraints_trimmer, x)
}
}
}
}
constraints_trimmer <- unique(constraints_trimmer)
if (length(constraints_trimmer) > 0){
message(paste(length(constraints_trimmer), " weird constraints got made and are being removed", sep=""))
constraints <- constraints[-constraints_trimmer]
ages <- ages[-constraints_trimmer]
}

###add in a root constraint if necessary###
##############################################this check is for cases where the species_to_date tree has extra taxa outside the backbone clade of interest
if (findMRCA(species_to_date, constraints[[1]], "node") != length(species_to_date$tip.label)+1){ ### if the first constraint isn't at the root node of the species tree
if (length(which(species_genus_only[-which(species_genus_only %in% unlist(lapply(lapply(lapply(constraints[[1]], strsplit, "_"), "[[", 1),"[",addition_species_position[[1]])))] %in% backbone_genus_only == FALSE)) == length(species_genus_only[-which(species_genus_only %in% unlist(lapply(lapply(lapply(constraints[[1]], strsplit, "_"), "[[", 1),"[",addition_species_position[[1]])))])){ ###if the species genus only genera outside the constraint and not in backbone are the same length as those outside the constraint - this does make sense
constraint_addition <- species_to_date$tip.label
age_addition <- max(node.depth.edgelength(extract.clade(backbone_tree, backbone_tree[[1]][,1][[which(backbone_tree[[1]][,2] == findMRCA(backbone_tree, clade_backbone_tree$tip.label, "node"))]])))-0.001 ### basically constraint max onto stem node
} else { 
constraint_addition <- species_to_date$tip.label
age_addition <- c(max(node.depth.edgelength(clade_backbone_tree)), "and min") ### constraint to clade backbone tree if somehow there are other overlapping backbone genera outside - seems impossible
message("weird root node correction")
}
}

###do check on whether both sides of constraint addition have backbone genus - even if some of ones outside main clade are backbone, the first branching event may be non backbone###
###############################################################################################################################################################################
if (length(constraint_addition) > 0){ ### basically a check on constraint additions
constraint_tree <- extract.clade(species_to_date, findMRCA(species_to_date, constraint_addition, "node"))
side_one <- extract.clade(constraint_tree, length(constraint_tree$tip.label)+2)$tip.label
side_one_genera <- unlist(lapply(lapply(lapply(side_one, strsplit, "_"), "[[", 1),"[",addition_species_position[[1]]))

side_two <- constraint_tree$tip.label[-which(constraint_tree$tip.label %in% side_one == FALSE)]
side_two_genera <- unlist(lapply(lapply(lapply(side_two, strsplit, "_"), "[[", 1),"[",addition_species_position[[1]]))

if (length(which(side_one_genera %in% backbone_genus_only)) == 0 || length(which(side_two_genera %in% backbone_genus_only)) == 0){
age_addition <- max(node.depth.edgelength(extract.clade(backbone_tree, backbone_tree[[1]][,1][[which(backbone_tree[[1]][,2] == findMRCA(backbone_tree, clade_backbone_tree$tip.label, "node"))]])))-0.001
}
}

###now an additional check for when the constraint is at the root of the species tree, but has the wrong date due to topological issue###
#########################################################################################################################################

if (findMRCA(species_to_date, constraints[[1]], "node") == length(species_to_date$tip.label)+1){
if (ages[[1]] != max(node.depth.edgelength(clade_backbone_tree))){
ages[[1]] <- max(node.depth.edgelength(clade_backbone_tree))
}
}
} else {
constraint_addition <- species_to_date$tip.label
age_addition <- max(node.depth.edgelength(extract.clade(backbone_tree, backbone_tree[[1]][,1][[which(backbone_tree[[1]][,2] == findMRCA(backbone_tree, clade_backbone_tree$tip.label, "node"))]])))-0.001
}
}
}

####################################################################################################################################
if (length(ages) == 0){ ###have to add a provision if theres a single tip for a genus###
####################################################################################################################################
if (length(which(backbone_genus_only %in% species_genus_only)) == 1){
ages <- max(node.depth.edgelength(backbone_tree)) - node.depth.edgelength(backbone_tree)[backbone_tree[[1]][,1][which(backbone_tree[[1]][,2] == which(backbone_genus_only %in% species_genus_only))]] - 0.001
constraints_mini[[length(constraints_mini)+1]] <- species_to_date$tip.label
}
if (length(clade_backbone_tree$tip.label) <= 1){
backbone_clade_genera <- unlist(strsplit(clade_backbone_tree$tip.label[[1]], split = "_"))[[3]]
if (length(which(backbone_genus_only %in% backbone_clade_genera)) >1){
backbone_clade_genera <- unlist(lapply(lapply(lapply(clade_backbone_tree$tip.label, strsplit, "_"), "[[", 1), "[",backbone_species_position[[1]])) ### get genera that are descended from the node
ages <- max(node.depth.edgelength(backbone_tree)) - node.depth.edgelength(backbone_tree)[backbone_tree[[1]][,1][which(backbone_tree[[1]][,2] == findMRCA(backbone_tree, which(backbone_genus_only %in% backbone_clade_genera), "node"))]] - 0.001
constraints_mini[[length(constraints_mini)+1]] <- species_to_date$tip.label
} else {
backbone_clade_genera <- unlist(lapply(lapply(lapply(clade_backbone_tree$tip.label, strsplit, "_"), "[[", 1), "[",backbone_species_position[[1]])) ### get genera that are descended from the node
ages <- max(node.depth.edgelength(backbone_tree)) - node.depth.edgelength(backbone_tree)[backbone_tree[[1]][,1][which(backbone_tree[[1]][,2] == which(backbone_genus_only %in% backbone_clade_genera))]] - 0.001
constraints_mini[[length(constraints_mini)+1]] <- species_to_date$tip.label
}
}
}

###write it to a file###
########################
write(paste("treefile = ", clade, "_to_date/backbone_", backbone_tree_bootstrap_number, "_bs_", species_sample, "_run_", run, "_bootstrap_tree_to_date.tre", sep=""), paste(clade, "_to_date/", run, ".treePL", sep=""), append=TRUE)
write("numsites=10000", paste(clade, "_to_date/", run, ".treePL", sep=""), append=TRUE)
if (length(constraints) > 0){
for (c in 1:length(constraints)){
if (constraints[[c]][[length(constraints[[c]])]] != "max"){
write(paste("mrca = crown_", c, " ", constraints[[c]][[1]], " ", constraints[[c]][[length(constraints[[c]])]], sep=""), paste(clade, "_to_date/", run, ".treePL", sep=""), append=TRUE)
write(paste("max = crown_", c, " ", ages[[c]], sep=""), paste(clade, "_to_date/", run, ".treePL", sep=""), append=TRUE)
write(paste("min = crown_", c, " ", ages[[c]], sep=""), paste(clade, "_to_date/", run, ".treePL", sep=""), append=TRUE)
} else {
write(paste("mrca = crown_", c, " ", constraints[[c]][[1]], " ", constraints[[c]][[length(constraints[[c]])-1]], sep=""), paste(clade, "_to_date/", run, ".treePL", sep=""), append=TRUE)
write(paste("max = crown_", c, " ", ages[[c]], sep=""), paste(clade, "_to_date/", run, ".treePL", sep=""), append=TRUE)
}
}
}
if (length(constraint_addition) > 0 & length(age_addition) == 1){
write(paste("mrca = crown_extra", " ", constraint_addition[[1]], " ", constraint_addition[[length(constraint_addition)]], sep=""), paste(clade, "_to_date/", run, ".treePL", sep=""), append=TRUE)
write(paste("max = crown_extra", " ", age_addition, sep=""), paste(clade, "_to_date/", run, ".treePL", sep=""), append=TRUE)
}
if (length(constraint_addition) > 0 & length(age_addition) == 2){
write(paste("mrca = crown_extra", " ", constraint_addition[[1]], " ", constraint_addition[[length(constraint_addition)]], sep=""), paste(clade, "_to_date/", run, ".treePL", sep=""), append=TRUE)
write(paste("max = crown_extra", " ", age_addition[[1]], sep=""), paste(clade, "_to_date/", run, ".treePL", sep=""), append=TRUE)
write(paste("min = crown_extra", " ", age_addition[[1]], sep=""), paste(clade, "_to_date/", run, ".treePL", sep=""), append=TRUE)
}
if (length(constraints_mini) > 0){
write(paste("mrca = stem", " ", constraints_mini[[1]][[1]], " ", constraints_mini[[1]][[length(constraints_mini[[1]])]], sep=""), paste(clade, "_to_date/", run, ".treePL", sep=""), append=TRUE)
write(paste("max = stem", " ", ages[[1]], sep=""), paste(clade, "_to_date/", run, ".treePL", sep=""), append=TRUE)
}
write("pliter=2", paste(clade, "_to_date/", run, ".treePL", sep=""), append=TRUE)
write("smooth=10", paste(clade, "_to_date/", run, ".treePL", sep=""), append=TRUE)
write("plsimaniter=10000", paste(clade, "_to_date/", run, ".treePL", sep=""), append=TRUE)
write("nthreads=8", paste(clade, "_to_date/", run, ".treePL", sep=""), append=TRUE)
write("thorough", paste(clade, "_to_date/", run, ".treePL", sep=""), append=TRUE)
write(paste("outfile = ", clade, "_to_date/backbone_", backbone_tree_bootstrap_number, "_bs_", species_sample, "_run_", run, "_bootstrap_tree_to_date.tre_dated", sep=""), paste(clade, "_to_date/", run, ".treePL", sep=""), append=TRUE)
write("prime", paste(clade, "_to_date/", run, ".treePL", sep=""), append=TRUE)
write.tree(species_to_date, paste(clade, "_to_date/backbone_", backbone_tree_bootstrap_number, "_bs_", species_sample, "_run_", run, "_bootstrap_tree_to_date.tre", sep=""))    
if (do_dating == TRUE){
run_tree_pl(paste(clade, "_to_date/", sep=""), run)
}
}
}
}
