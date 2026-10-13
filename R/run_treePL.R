#' @export

run_tree_pl <- function(input_direc, run){ 
system(paste('treePL ', input_direc, run, '.treePL > ', input_direc, 'primer_', run, '.txt', sep=""))
###
primer <- readLines(paste(input_direc, "primer_", run, ".txt", sep=""))
if ("PLACE THE LINES BELOW IN THE CONFIGURATION FILE" %in% primer){
starter <- which(primer == "PLACE THE LINES BELOW IN THE CONFIGURATION FILE")
optimisations <- vector(mode="numeric", length=0)
for (i in (starter+1):length(primer)){
optimisations <- append(optimisations, primer[[i]])
}
###
treePL_input <- readLines(paste(input_direc, run, ".treePL", sep=""))
treePL_input <- treePL_input[-which(treePL_input == "prime")]
for (i in 1:length(optimisations)){
treePL_input <- append(treePL_input, optimisations[[i]])
}
write(treePL_input, paste(input_direc, run, ".treePL", sep=""))
system(paste('treePL ', input_direc, run, '.treePL > ', input_direc, 'actual_', run, '.txt', sep=""))
}
}

