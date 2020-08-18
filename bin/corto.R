#!/labs/genut/bin/Rscript

library(corto)
library(argparse)


parser <- ArgumentParser(description="Generate correlation-based DPI networks using corto")
parser$add_argument("--RDS", required=TRUE, help="RDS object containing expression matrix[genes,cells]")
parser$add_argument("--transpose", default=FALSE, help="if matrix[cells,genes], the transpose to matrix[genes,cells]")
parser$add_argument("--TFs", default="", type = "character", help="list of transcription factors")
parser$add_argument("--nbootstraps", default=100, help="Number of bootstraps to be performed. Default is 100")
parser$add_argument("--pvalue", default=1e-30, help="The p-value threshold for correlation significance (default 1e-30)")
parser$add_argument("--nthreads", default=1, help="number of threads to use for bootstrapping. Default is 1")
parser$add_argument("--sif", default=FALSE, help="Save network in SIF format")
parser$add_argument("--output", default="./", help="name of the Rdata output folder")
args <- parser$parse_args()


# Read matrix
inmat <- readRDS(args$RDS)
# transpose matrix
ifelse(test = args$transpose == TRUE, yes = inmat <- t(inmat), no = NULL)
# read centroids
ifelse(test = args$TFs == "", 
	yes = centroids <- rownames(inmat), 
	no = centroids <- row.names(inmat)[grep(gsub('.{1}$', '', paste0("^",read.table(file = args$TFs, stringsAsFactors = FALSE,)[[1]], "$", "|", collapse = "")), rownames(inmat))])

# Get regulon

regulon <- corto(inmat = inmat, 
              centroids = centroids, 
              nbootstraps = args$nbootstraps, 
              p = args$pvalue, 
              nthreads=args$nthreads)


# sabe network

name <- stringr::str_remove(string = basename(args$RDS), pattern = ".RDS")

if (args$sif == FALSE) {
	ifelse(test = args$output == "./", yes = outpath <- paste0(args$output, name, ".RData"), no = outpath <- args$output)
	save(regulon, file = outpath)
} else {
	ifelse(test = args$output == "./", yes = outpath <- paste0(args$output, name, ".sif"), no = outpath <- args$output)
	regulon2sif <- function(regulon){
				sif <- data.frame(TF = c(), target = c(), interaction = c())
				for (name in names(regulon)) {
						source = rep(name, length(regulon[[name]]$tfmod));
						psif    = cbind(source, names(regulon[[name]]$tfmod), regulon[[name]]$tfmod);
						sif <- rbind(sif, psif)}
				colnames(sif) <- c("TF", "target", "int")
				rownames(sif) <- NULL
				return(sif)
	}
	sif <- regulon2sif(regulon)
	write.table(sif, file = outpath, quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
}



