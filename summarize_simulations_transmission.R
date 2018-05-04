# This script summarizes the average transmission probabilities for simulations. 
# Input: (command line arguments)
#	- directories: strings giving the folders where the simulation results 
#		are. These should be the parent folders where the simulation 
#		script saved its output, because the subfolders' names are used
#		to determine the parameters of the simulations.
#
# Output:
#	- files summarizing the average horizontal and vertical transmission
#		probabilities in patches where the symbiont is parasitic and
#		mutualistic. Simulations in the same parent folder with the same
#		parameters but different initial values will be grouped as columns
#		of a single file. Rows of the file indicate the state of
#		simulations at different time points.

library("plyr")

# Read and check input

directories <- commandArgs(trailingOnly=T)
if (length(directories) < 1) {
	writeLines(paste("This script summarizes the average transmission probabilites for simulations.",
		"\tInput: (command line arguments)",
		"\t - directories: strings giving the folders where the simulation results",
		"\t   are. These should be the parent folders where the simulation",
		"\t   script saved its output, because the subfolders' names are used",
		"\t   to determine the parameters of the simulations.\n",
		"\tOutput:",
		"\t - files summarizing the average horizontal and vertical transmission",
		"\t   probabilities in patches where the symbiont is parasitic and mutualistic.",
		"\t   Simulations in the same parent folder with the same parameters",
		"\t   but different initial values will be grouped as columns of a",
		"\t   single file. Rows of the file indicate the state of simulations",
		"\t   at different time points.\n", sep="\n"))

	stop("No files given for analysis.")
}

if (class(directories) != "character") {
	stop("File name(s) must be strings.")
}
if (!all(file.exists(directories))) {
	stop(paste0("The following directories cannot be found:\n",
		paste(directories[!file.exists(directories)], collapse="\n")))
}

# Summarize 

parentdir <- getwd()

for (mydir in directories) {
	setwd(parentdir)
	setwd(mydir)
	for (fs in list.dirs(full.names=F, recursive=F)) {
		for (d in list.dirs(fs, full.names=F, recursive=F)) {
			for (np in list.dirs(paste(fs, d, sep="/"), full.names=F, recursive=F)) {

				# Get the parameters
				patches <- as.numeric(strsplit(np, " = |, ")[[1]][4])
				patchesM <- (patches/2)+1:patches # Patches where the symbiont is mutualistic (M-patches)
				sim_files <- list.files(paste(fs, d, np, sep="/"))
				num_hv_sims <- length(sim_files) # number of simulations w/ different starting transmission probabilities
				time_points <- sort(unique(read.table(paste(fs, d, np, list.files(paste(fs, d, np, sep="/")), sep="/")[[1]], header=T, sep=",", strip.white=T)$time))
				patchP_h <- patchM_h <- patchP_v <- patchM_v <- matrix(nrow=length(time_points)+1, ncol=num_hv_sims+1)
				patchP_h[,1] <- patchM_h[,1] <- patchP_v[,1] <- patchM_v[,1] <- c(0, time_points)

				for (hv in 1:num_hv_sims) {
					sim_data <- read.table(paste(fs, d, np, sim_files[[hv]], sep="/"),
						header=T, sep=",", strip.white=T)
					sim_data$symb <- factor(sim_data$symb, levels=1:2)
					sim_data$patch <- as.factor(sim_data$patch)
					sim_data$patchM <- sim_data$patch %in% patchesM # gives 1 if the patch is an M-patch, 0 otherwise

					avgs <- ddply(sim_data[, c("time", "h", "v", "patchM")], c(.(time), .(patchM)), colwise(mean, na.rm=T), .drop=F)

					patchP_h[,hv+1] <- c(as.numeric(strsplit(sim_files[[hv]], " = |, |.csv")[[1]][2]), avgs[!(avgs$patchM), "h"])
					patchM_h[,hv+1] <- c(as.numeric(strsplit(sim_files[[hv]], " = |, |.csv")[[1]][2]), avgs[avgs$patchM, "h"])
					patchP_v[,hv+1] <- c(as.numeric(strsplit(sim_files[[hv]], " = |, |.csv")[[1]][4]), avgs[!(avgs$patchM), "v"])
					patchM_v[,hv+1] <- c(as.numeric(strsplit(sim_files[[hv]], " = |, |.csv")[[1]][4]), avgs[avgs$patchM, "v"])
				}
				# Record average transmission probabilities
				write.csv(patchP_h, file=paste0("Horiz trans evo patch P, ", fs, ", ", np, ", ", d, ".csv"))
				write.csv(patchM_h, file=paste0("Horiz trans evo patch M, ", fs, ", ", np, ", ", d, ".csv"))
				write.csv(patchP_v, file=paste0("Vert trans evo patch P, ", fs, ", ", np, ", ", d, ".csv"))
				write.csv(patchM_v, file=paste0("Vert trans evo patch M, ", fs, ", ", np, ", ", d, ".csv"))
			}
		}
	}
}
