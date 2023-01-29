
library(metacell)
library("Matrix")
library(viridis)
library(zoo)
library(tgstat)

kidney_norns_stable_genes_cells_F = "./score/kidney_norns_stable_genes_cells_F.txt"
mouse2humanfile <- "./score/mouse2human.txt"


read_large_umis_n <- function(mat_id, bs = 1e4, cells = NULL, norm_to= 2000, do_log = F) {
	mat = scdb_mat(mat_id)
	if (is.null(cells)) { cells = mat@cells}
	ncells = length(cells)
	umis = NULL
	for (i in seq_len(ncells %/% bs + 1)) {
		from = (i - 1) * bs + 1
		to = min(i * bs, ncells)
		mt <- as.matrix(mat@mat[, cells[from:to]])
		mtn <- sweep(mt, 2, colSums(mt), "/") * norm_to
		if (do_log == T){
			umis = cbind(umis, log(mtn+1, 2))
		} else {
			umis = cbind(umis, mtn)
		}
		mt <- NULL
		mtn <- NULL
	}
	umis
}


image.2 = function(X, col = colorRampPalette(c("blue", "white", "red"))(1000), balance=F, annotate="both", zlim = NULL, hct = NULL, vct = NULL, lwd=1, lty=1, cex = 1, text=F, text_mat = X) 
{
	if (is.null(zlim)) {
		if (balance) {
			zlim = max(abs(X), na.rm=T) * c(-1,1)
		} else {
			zlim = quantile(X, c(0,1), na.rm=T)
		}
	}
	hcls = NULL; vcls = NULL; hticks = NULL; vticks = NULL
	nrow = nrow(X); ncol = ncol(X)
	if (!is.null(hct)) {
		X = X[order(hct),]; text_mat = text_mat[order(hct),]
		hticks = seq(-1 / (2 * (nrow-1)),1 + 1 / (2 * (nrow-1)), length.out = nrow + 1)
		hcls = cumsum(table(hct)) + 1
	}
	if (!is.null(vct)) {
		X = X[,order(vct)]; text_mat = text_mat[,order(vct)]
		vticks = seq(-1 / (2 * (ncol-1)),1 + 1 / (2 * (ncol-1)), length.out = ncol + 1)
		vcls = cumsum(table(vct)) + 1
	}
	message("zlim: ", zlim[1], "<>", zlim[2])
	image(t(X), axes = F, col = col, zlim = zlim)
	abline(h = hticks[hcls], v = vticks[vcls], lwd=lwd, lty=lty)
	if (annotate %in% c("rows", "both")) {
		mtext(rownames(X), las = 2, side=2, at = (1 - seq_len(nrow(X))) / (1 - nrow(X)), cex = cex)
		if (!is.null(hct)) {
			mtext(names(hcls), side = 4, las = 2, at = rowMeans(cbind(hticks[c(1,hcls[-length(hcls)])], hticks[hcls])), cex = cex)
		}
	}
	if (annotate %in% c("columns", "both")) {
		mtext(colnames(X), las = 2, side=1, at = (1 - seq_len(ncol(X))) / (1 - ncol(X)), cex = cex)
		if (!is.null(vct)) {
			mtext(names(vcls), side = 3, las = 2, at = rowMeans(cbind(vticks[c(1,vcls[-length(vcls)])], vticks[vcls])), cex = cex)
		}
	}
	if (text) {
	 hmed = seq(0,1,length.out=nrow); vmed = seq(0,1,length.out=ncol)
	 text(rep(vmed, each = nrow), rep(hmed, ncol), text_mat)
	}
}

score_on_mat <- function(id_use, mypath) {
	sc_mat = scdb_mat(id_use)
	metadata <- sc_mat@cell_metadata
	umis_n = read_large_umis_n(id_use)

	scorefilelist <- list.files(path=mypath, pattern="_F.txt")
	scorefilelist <- grep(pattern="score", scorefilelist, value=T)
	scorenames <- gsub(pattern="_score_genes_cells_F.txt", replacement="",gsub(pattern="kidney_norns_", replacement="", scorefilelist))
	scores = list()
	mouse2human <- read.table(file = mouse2humanfile, sep="\t", header=T)

	for(i in 1:length(scorenames)){
		add_genes <- read.table(file = paste0(mypath, scorefilelist[i]), sep="\t", header=T)[, "add_genes",drop=F]
		scores[scorenames[i]] <- data.frame(add_genes = unique(mouse2human[mouse2human[, "Gene.name"] %in% add_genes[,1], "Human.gene.name"]))
		scores[paste0(scorenames[i], "_V")] <- data.frame(rep(1, length(unique(mouse2human[mouse2human[, "Gene.name"] %in% add_genes[,1], "Human.gene.name"]))))
	}

	high_genes_stable_mouse <- read.table(file = kidney_norns_stable_genes_cells_F, sep="\t", header=T)[,]
	high_genes_stable <- unique(mouse2human[mouse2human[, "Gene.name"] %in% high_genes_stable_mouse, "Human.gene.name"])

	score_cutpoint <- as.data.frame(matrix(nrow=0, ncol=2))
	for(i in 1:length(scorenames)){
		score_cutpoint[scorenames[i],] <- read.table(file = paste0(mypath, "optimal_cutpoint_",scorenames[i],".txt"), sep="\t", header=T)[,1]
		gn <- as.character(unlist(scores[scorenames[i]]))
		gs <- as.numeric(unlist(scores[paste0(scorenames[i], "_V")]))
		if (length(setdiff(gn, rownames(umis_n))) > 0 ){score_cutpoint[scorenames[i],"V1"] = score_cutpoint[scorenames[i],"V1"]*(1 - (sum(gs[!gn %in% rownames(umis_n)]) / sum(gs)))  }
	}

	newdata = as.data.frame(colSums(umis_n[intersect(rownames(umis_n),high_genes_stable), ]))
	for(i in 1:length(scorenames)){
		newdata <- cbind(newdata, colSums(umis_n[intersect(rownames(umis_n),as.character(unlist(scores[scorenames[i]]))), ]))
	}
	colnames(newdata) <- c("high_genes_stable", scorenames)

	factortonorm <- median(newdata[,"high_genes_stable"])
	best_of_type <- as.data.frame(colnames(umis_n))
	best_of_type[, "type"] <- ""
	best_of_type[, "percent"] <- 0
	rownames(best_of_type) <- best_of_type[,1]
	for(test in rownames(score_cutpoint)){
		tempres <- round((newdata[, test])/(score_cutpoint[test, 1]/score_cutpoint[test, 2]*factortonorm+0.00001) * 100, 1)
		best_of_type <- apply(cbind(best_of_type, tempres), 1, function(X) {if(as.numeric(X[4]) > 100 & as.numeric(X[4]) > as.numeric(X[3])){c(X[1],test,X[4])} else {c(X[1],X[2],X[3])}})
		best_of_type <- t(best_of_type)
	}

	write.table(file=paste0(id_use, "_type_counts.txt"), sep="\t", table(best_of_type[,2]), quote=F)
	write.table(file=paste0(id_use, "_best_of_type.txt"), sep="\t", best_of_type, quote=F)

	res = sapply(rownames(score_cutpoint), function(test) round((newdata[, test])/(score_cutpoint[test, 1]/score_cutpoint[test, 2]*factortonorm+0.00001) * 100, 1))
	ceil_res = pmin(res, 200)

	anno = c(colnames(ceil_res), "none")[ max.col(cbind(ceil_res, 100))]

	png(file = paste0(id_use, "_type_cells.png"), width=2400, height=1500)
	par(mar = c(2, 15, 15, 1))
	image.2(t(ceil_res), col = colorRampPalette(c("white", "pink", "purple", "navyblue"))(1000), vct = factor(anno, levels = c(colnames(ceil_res), "none")))
	box()
	dev.off()
}
