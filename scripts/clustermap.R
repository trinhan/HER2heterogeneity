# R package for heatmap and hierarchical clustering visualization
# Author: Ole Christian Lingjærde
# Date: 05 Sep 2016 v.1.2.3

# plot.init		Call this first to set up correct margins etc
# plot.hmap 	Plot heatmap using your preferred colors
# plot.cbar 	Add color bar to one of the sides (to add extra info)
# plot.tree		Add column or row dendrogram to a heatmap
# hcluster     	Perform hierarchical clustering of rows or columns
# subclust		Manually or automatically identify subclusters
# get.subclust  Get subtype labels (for original row/col ordering of X)
# set.color     Set color coding for accompanying data
# plot.hmap.key Plot color key for heatmap
# plot.cbar.key	Plot color key for color bar (annotation bars)
# plot.text		Add text to one of the sides (for row/column labels)
# blocks		Generate an example data set with 110 rows and 40 columns

library(gclus)
library(foreach)
library(doParallel)
registerDoParallel(cores=7)

# Function to make a small example data set
blocks = function() {
	# Example data set with 110 rows and 40 columns
	A = cbind(matrix(rnorm(70*20),nrow=70), matrix(3+rnorm(70*20),nrow=70))
	B = cbind(-1 + matrix(rnorm(40*30),nrow=40), matrix(1+rnorm(40*10),nrow=40))
	X = rbind(A,B)

	# Additional information tracks for the columns
	y1 = c(rep("Lum A",14), rep("Lum B",6), rep("Her2",8), rep("Normal",4), rep("Basal",8))
	y2 = c(rep(1,10), rep(2,20), rep(3,10))
	y3 = runif(40, 1, 16)

	# Example text descriptor for rows and cols
	rownames(X) = paste("Gene", 1:110, sep="-")
	colnames(X) = paste("Patient", 1:40, sep="-")

	# Randomly permute rows and columns
	k.row = sample(110); k.col = sample(40)
	X = X[k.row, k.col]
	y1 = y1[k.col]; y2 = y2[k.col]; y3 = y3[k.col]

	# Collect it all in a list
	res = list(X=X, y1=y1, y2=y2, y3=y3)
	return(res)
}

pca.approx = function(X, ncomp) {
	if (ncomp > min(ncol(X),nrow(X))) {
		ncomp = min(ncol(X),nrow(X))
		message(paste("Note: only the first", ncomp, "principal components are used"))
	}
	z = svd(X)
	if (ncomp == 1) {
		d = matrix(z$d[1],1,1)
	} else {
		d = diag(z$d[1:ncomp])
	}
	X = z$u[,1:ncomp,drop=F] %*% d %*% t(z$v[,1:ncomp,drop=F])
	X
}

# Use the function below to perform clustering of a data matrix X. 
# To cluster both rows and columns, call it twice with clust="col" 
# and clust="row" respectively. The option reorder=T leads to a nice 
# arrangement of the clusters from left to right to ensure neighbors 
# are similar to each other

hcluster = function(X, clust="col", distance="euclidean", linkage="complete", reorder=T, ncomp=NA) {
	# clust: one of "col", "row"
	# dist: one of "euclidean", "pearson", "spearman"
	# linkage: one of "complete", "average", "single", "ward"
	# reorder: arrange the clustered items in optimized order
	# ncomp: if set to k, only the first k principal components are used in clustering
	if (!is.na(ncomp)) {
		X = pca.approx(X, ncomp)	
	}
	if (distance=="euclidean") {
		if (clust=="col") {
			d = dist(t(X), method="euclidean")
		} else {
			d = dist(X, method="euclidean")
		}
	} else if (distance=="pearson") {
		if (clust=="col") {
			d = as.dist((1-cor(X, method="pearson"))/2)
		} else {
			d = as.dist((1-cor(t(X), method="pearson"))/2)
		}
	} else if (distance=="spearman") {
		if (clust=="col") {
			d = as.dist((1-cor(X, method="spearman"))/2)
		} else {
			d = as.dist((1-cor(t(X), method="spearman"))/2)
		}
	}
	hc = hclust(d, method=linkage)
	if (reorder) {
		hc = reorder.hclust(hc, d)
	}
	tmp = .CLUSTERMAP
	if (clust=="row") {
		tmp$row.X = X
		tmp$row.ncomp = ncomp
		tmp$rowclust = hc
		tmp$row.linkage = linkage
		tmp$row.distance = distance
	} else if (clust=="col") {
		tmp$col.X = X
		tmp$col.ncomp = ncomp
		tmp$colclust = hc
		tmp$col.linkage = linkage
		tmp$col.distance = distance
	}
	assign(".CLUSTERMAP", tmp, envir=.GlobalEnv)
}

# Use the function below to identify subclusters in a dendrogram.
# The effect is that later calls to plot.tree will show the subclusters
# in different colors. You can either provide the desired number of 
# clusters using the argument k, or you can set k=NA to automatically
# estimate the number of cluster with one of two algorithms: gap or part. 
# Using gap (Tibshirani et al (2001), J Roy Statist Soc B, 63: 411-423)
# the result is an estimate of the optimal 'flat cut' of the dendrogram, 
# while using part (Nilsen et al (2012), Stat Appl Genet Mol Biol, 
# 12: 637-652) the result is an estimate of the optimal 'nonflat cut' of 
# the dendrogram; is essence the part algorithm is an extension of gap
# that applies gap recursively on subclusters previously identified. 
# B and min.size are used by PART (only if k=NA) and determines the
# number of permutations (B=25 is low, but computation time increases
# with increasing B; consider B=200 or higher for high-quality estimates)
# and the least number of elements in a single cluster (min.size=5 is
# quite arbitrary and you may want to change this). Note that PART may
# call some individuals as outliers and these will be shown in a separate
# color but may not be close to each other and will stick out as
# "clusters" of size potentially smaller than min.size.

subclust = function(k=NA, clust="col", method="gap", B=50, min.size=5, max.level=3) {
	tmp = .CLUSTERMAP
	if (clust=="col" && !is.na(k)) {
		tmp$colgroup = cutree(.CLUSTERMAP$colclust, k)
	} else if (clust=="col" && is.na(k) && method=="gap") {
		tmp$colgroup = gap.new(.CLUSTERMAP$col.X, clust, .CLUSTERMAP$col.distance, .CLUSTERMAP$col.linkage, B)
	} else if (clust=="col" && is.na(k) && method=="part") {
		tmp$colgroup = part.new(.CLUSTERMAP$col.X, clust, .CLUSTERMAP$col.distance, .CLUSTERMAP$col.linkage, B, min.size, max.level)
	} else if (clust=="row" && !is.na(k)) {
		tmp$rowgroup = cutree(.CLUSTERMAP$rowclust, k)
	} else if (clust=="row" && is.na(k) && method=="gap") {
		tmp$rowgroup = gap.new(.CLUSTERMAP$row.X, clust, .CLUSTERMAP$row.distance, .CLUSTERMAP$row.linkage, B)
	} else if (clust=="row" && is.na(k) && method=="part") {
		tmp$rowgroup = part.new(.CLUSTERMAP$row.X, clust, .CLUSTERMAP$row.distance, .CLUSTERMAP$row.linkage, B, min.size, max.level)
	} 
	assign(".CLUSTERMAP", tmp, envir=.GlobalEnv)
}

# Use the function below to retrieve cluster labels after identification
# of subclusters. By default, labels are returned in the same order as the 
# original ordering of columns/rows in the input data matrix X. To obtain
# labels in the order shown in the cluster tree, use order="tree". By 
# default, the function returns a vector of cluster labels. If you 
# supply row/column labels using the argument 'labels' the function returns
# a data frame with two columns, the first being the cluster labels and
# the second being the supplied row/column labels. 

get.subclust = function(clust="col", labels=c(), order="tree") {
	if (clust=="col") {
		group = .CLUSTERMAP$colgroup
		if (order=="tree") {
			group = group[.CLUSTERMAP$colclust$order]
			if (length(labels)>0) labels = labels[.CLUSTERMAP$colclust$order]
		}
	} else if (clust=="row") {
		group = .CLUSTERMAP$rowgroup
		if (order=="tree") {
			group = group[.CLUSTERMAP$rowclust$order]
			if (length(labels)>0) labels = labels[.CLUSTERMAP$rowclust$order]
		}
	}
	if (length(labels)>0) {
		res = data.frame(cluster=group, label=labels)
	} else {
		res = group
	}
	res
}

get.order = function(clust="col") {
	if (clust=="col") {
		return(.CLUSTERMAP$colclust$order)
	} else if (clust=="row") {
		return(.CLUSTERMAP$rowclust$order)
	}
}

# Use the function below to initiate a plot. The purpose is to make
# room on the four sides of the heatmap for additional elements such
# as cluster dendrograms, color bars with additional info, and labels.
# Supply as arguments the sides where you want such elements, with the
# coding side=1 (below), side=2 (left), side=3 (top), side=4 (right).
# If you want a tree (=dendrogram) on the left and top, let tree=c(2,3).
# If you want text labels on the right, let text=4. If you want color
# bars below the plot, let cbar=1. Examples:
#     plot.init()                Dendrogram on the left and on the top
#     plot.init(tree=c())        Make no room for additional elements
#     plot.init(tree=3, text=4)  Dendrogram on the top, labels on the right
# The last two arguments (inner and outer) do not overrule the other
# parameters, but allows the user to add (or even subtract) a specific
# amount of space in the inner and outer margins on a specified side of
# the heatmap. They may be used alone, or in combination with other
# parameters. They should always be vectors of length 4 giving the 
# desired margin extension for each of the four sides. An example of use
# would be to extend an inner margin by a certain amount because we
# plan to show many bar plots of accompanying data there that we want
# to give some space. Note that inner=c(0.4,0.4,0.4,0.4) means that each inner
# margin should be extended by a length corresponding to 0.4 times the
# total length (in that direction) of the heatmap itself. 

plot.init = function(tree=c(), text=c(), cbar=c(), ckey=T, inner=c(), outer=c()) {
	if (ckey) {
		layout(matrix(c(2,3,3,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1),6,3,byrow=TRUE))
		par(mar=c(2,2,2,2))
	}
	marg.outer = rep(0, 4)
	marg.inner = rep(0.01, 4)
	marg.outer[tree] = 0.3
	marg.outer[text] = 0.3
	marg.inner[cbar] = 0.1
	if (length(inner)==4) marg.inner = marg.inner + inner
	if (length(outer)==4) marg.outer = marg.outer + outer
	marg = list(inner=marg.inner, outer=marg.outer)
	tmp = list(margin=marg, rowclust=c(), colclust=c(), rowgroup=c(), 
			colgroup=c(), row.ncomp=NA, col.ncomp=NA, col.X=c(), row.X=c())
	assign(".CLUSTERMAP", tmp, envir = .GlobalEnv)
}

panel.init = function(mfrow=c(2,2)) {
	A = matrix(c(2,3,3,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1),6,3,byrow=TRUE)
	B = matrix(NA, 6*mfrow[1], 3*mfrow[2])
	for (i in 1:mfrow[1]) {
		for (j in 1:mfrow[2]) {
			B[(1+(i-1)*6):(i*6),(1+(j-1)*3):(j*3)] = A + 3*(j-1) + mfrow[2]*3*(i-1)
		}
	}
	layout(B)
	par(mar=c(0,0,0,0))
}

set.colgroup = function(group) {
  tmp = .CLUSTERMAP
  tmp$colgroup = group
  assign(".CLUSTERMAP", tmp, envir = .GlobalEnv)
}

set.colclust = function(clust) {
  tmp = .CLUSTERMAP
  tmp$colclust = list()
  tmp$colclust$order = clust
  assign(".CLUSTERMAP", tmp, envir = .GlobalEnv)
}

# For additional information that you want to put into color bars
# around the heatmap, you may either specify directly a color for
# each element, or if the information is numerical you may use one
# of the two functions below to convert to colors. See top section
# of Examplecode_clustermap.R for three examples of use.

set.color = function(x, type="discrete", label, color, na.color) {
	if (type=="continuous") {
		if (missing(label)) label = ""
		if (missing(color)) color = "green-white-red" # color=palette()
		if (missing(na.color)) na.color="grey"
		res = color.cont(x, label=label, color=color, na.color=na.color)
	} else if (type=="discrete") {
		if (missing(label)) label = ""
		if (missing(color)) color = palette() # color="green-white-red"
		if (missing(na.color)) na.color="grey"
		res = color.disc(x, label=label, color=color, na.color=na.color)
	}
	return(res)
}

color.cont = function(x, label, color, na.color) {
	M = max(abs(x),na.rm=T)
    z = rbind(as.numeric(x)/M)
    z.mis = is.na(z)
    z[z.mis] = mean(z[!z.mis])
    x.grid = seq(min(x,na.rm=T), max(x,na.rm=T), length=100)
    z.grid = rbind(as.numeric(x.grid)/M)
    c0 = col2rgb(unlist(strsplit(color,"-")))/255
    if (ncol(c0)==1) {
    	c0 = cbind(c(1,1,1), c0)
    }
   	if (ncol(c0)==2) {
    	z2 = (z-min(z))/(max(z)-min(z))
    	c1 = cbind(c0[,1]) %*% (1-z2) + cbind(c0[,2]) %*% z2
    	c2 = rgb(c1[1,], c1[2,], c1[3,])
    	z2.grid = as.numeric((z.grid-min(z.grid))/(max(z.grid)-min(z.grid)))
    	c1.grid = cbind(c0[,1]) %*% (1-z2.grid) + cbind(c0[,2]) %*% z2.grid
    	c2.grid = rgb(c1.grid[1,], c1.grid[2,], c1.grid[3,])
    } else if (ncol(c0)==3) {
    	c1 = matrix(c0[,2], 3, length(z), byrow=F) + 
    	cbind(c0[,1]-c0[,2]) %*% ifelse(z<0,-z,0) + 
    	cbind(c0[,3]-c0[,2]) %*% ifelse(z>0,z,0)
    	c2 = rgb(c1[1,], c1[2,], c1[3,])
    	z2.grid = as.numeric((z.grid-min(z.grid))/(max(z.grid)-min(z.grid)))    	
    	c1.grid = matrix(c0[,2], 3, length(z.grid), byrow=F) + 
    	cbind(c0[,1]-c0[,2]) %*% ifelse(z.grid<0,-z.grid,0) + 
    	cbind(c0[,3]-c0[,2]) %*% ifelse(z.grid>0,z.grid,0)
    	c2.grid = rgb(c1.grid[1,], c1.grid[2,], c1.grid[3,])
    }
    c2[z.mis] = "grey"
    key.x = x.grid #x[!z.mis]
    key.z = z2.grid #z2[!z.mis]
    key.c = c2.grid #c2[!z.mis]
    list(label=label, type="cont", x=x, color=c2, na.color=na.color, 
    	key=data.frame(x=key.x, z=key.z, color=key.c))
}

color.disc = function(x, label, color, na.color) {
    val = sort(unique(x))
    z = rep(NA, length(x))
    ncol = length(color)
    for (i in 1:length(val)) {
    	  z[x==val[i]] = color[1+(i-1) %% ncol]
    }
    z[is.na(x)] = na.color
	list(label=label, type="disc", x=as.factor(x), color=z, na.color=na.color, 
		key=data.frame(x=val, color=color[1+(0:(length(val)-1)) %% ncol]))
}

# The function below creates a color key, i.e. a graphical overview
# of what numerical values the different colors found in the heatmap 
# correspond to. 

plot.hmap.key = function() {
	colrange = COLOR_SCALE
	xmin = -max(abs(colrange$valuerange))
	xmax = max(abs(colrange$valuerange))
	plot(0, 0, type = "n", xlim = c(xmin, xmax), ylim = c(0,1), xlab = "", 
			ylab = "", xaxt = "n", yaxt = "n")
	axis(1, c(xmin,xmax), round(c(xmin,xmax) * 100)/100)
	axis(1, 0, 0)
	z1.plot = rbind(seq(-max(abs(colrange$valuerange)), 
		max(abs(colrange$valuerange)), length=100))
	z1 = rbind(seq(-1, 1, length=100))
	c0 = col2rgb(unlist(strsplit(colrange$colorscale,"-")))/255
	c1 = matrix(c0[,2],3,length(z1),byrow=F) + 
		  cbind(c0[,1]-c0[,2])%*%ifelse(z1<0,-z1,0) + 
			  cbind(c0[,3]-c0[,2])%*%ifelse(z1>0,z1,0)
	c2 = rgb(c1[1,],c1[2,],c1[3,])
	for (i in 1:(length(c2)-1)) {
		polygon(c(z1.plot[i], z1.plot[i], z1.plot[i+1], z1.plot[i+1]), c(0,1,1,0), 
			col = c2[i], border = FALSE)
	}
	if (sum(c0)>2) {abline(v=0,col="black")} else {abline(v=0,col="white")}
	#frame()
}

plot.cbar.key = function(mfrow=c(2,5), border=F, vsize=0.1, hsize=0.2) {
	frame()
	keys = .CLUSTERMAP$cbar.key
	type = .CLUSTERMAP$cbar.type
	nbar = length(keys)
	nrow = mfrow[1]
	ncol = mfrow[2] # ceiling(nbar/nrow)
	bty = ifelse(border, "o", "n")
	par(mfrow=c(nrow,ncol))
	par(mar=c(2,2,2,2))
	for (i in 1:nbar) {
		plot(0,0,xlim=c(0,1),ylim=c(0,1),type="n",xaxt="n",yaxt="n",bty=bty)
		if (type[i]=="disc") {
			nlevels = nrow(keys[[i]])
			text(0.1, 1, names(keys)[i], adj=c(0,1), cex=1.2, font=2)
			for (j in 1:nlevels) {
				polygon(c(0.1,0.1,0.1+hsize,0.1+hsize), 
					c(1-vsize-(j-1)*vsize, 1-vsize-j*vsize, 
						1-vsize-j*vsize, 1-vsize-(j-1)*vsize),
					col=as.character(keys[[i]][j,2]), border="black")
				text(0.1+2*hsize, 1-vsize-j*vsize+0.5*vsize, 
					as.character(keys[[i]][j,1]), adj=c(0,0.5))
			}
		} else if (type[i]=="cont") {
			xval = as.character(keys[[i]][,1])
			keep = !duplicated(xval)
			xval.unique = xval[keep]
			zval.unique = keys[[i]][keep,2]
			col.unique = as.character(keys[[i]][keep,3])
			nlevels = length(xval.unique)
			text(0.1, 1, names(keys)[i], adj=c(0,1), cex=1.2, font=2)
			for (j in 1:nlevels) {
				polygon(c(0.1,0.1,0.1+hsize,0.1+hsize), 
					c(1-vsize-0.8*(j-1)/nlevels, 1-vsize-0.8*j/nlevels, 
					  1-vsize-0.8*j/nlevels, 1-vsize-0.8*(j-1)/nlevels),
					col=col.unique[j], border=FALSE)
			}
			text(0.1+2*hsize, 1-vsize-0.5*0.8/nlevels, round(xval.unique[1],2), adj=c(0,1))
			text(0.1+2*hsize, 1-vsize-0.8+0.8*0.5/nlevels, round(xval.unique[length(xval.unique)],2), adj=c(0,0))	
		}
	}
}

# The function below is used to add labels to the sides of the heatmap.
# The first argument should be a vector of the same length as the number
# of columns in the heatmap (when side=1 or side=3), and a vector of the
# same length as the number of rows in the heatmap (when side=2 or side=4).
# You may change the size of the text with cex=... and you may truncate
# text labels to a certain maximal length (say 4) with maxchar=...

plot.text = function(txt, side=4, sep.outer=0, cex=0.5, maxchar=NA) {
	if (side==1 || side==3) {
	    if (length(.CLUSTERMAP$colclust)>0) {
		  	txt = txt[.CLUSTERMAP$colclust$order]
		}
	} else if (side==2 || side==4) {
		if (length(.CLUSTERMAP$rowclust)>0) {
			txt = txt[.CLUSTERMAP$rowclust$order]
		}
	}
	margin = .CLUSTERMAP$margin
    N = length(txt)
    if (!is.na(maxchar)) {
    		for (i in 1:N) txt[i] = substr(txt[i],1,maxchar)
    }
	if (side==1) {
		for (k in 1:N) {
			text((2*k-1)/(2*N), -margin$inner[1]-sep.outer, txt[k], 
				adj=c(0,0.5), srt=-90, cex=cex)
		}		
	} else if (side==2) {
		for (k in 1:N) {
			text(-margin$inner[2]-sep.outer, (2*k-1)/(2*N), txt[k], 
				adj=c(1,0.5), cex=cex)
		}		
	} else if (side==3) {
		for (k in 1:N) {
			text((2*k-1)/(2*N), 1+margin$inner[1]+sep.outer, txt[k], 
				adj=c(0,0.5), srt=90, cex=cex)
		}		
	} else if (side==4) {
		for (k in 1:N) {
			text(1+margin$inner[4]+sep.outer, (2*k-1)/(2*N), txt[k], 
				adj=c(0,0.5), cex=cex)
		}
	}
}
 
plot.cbar = function (..., side=3, labels=T, pvalue=F, pvalue.method="chisq", 
			cex=0.8, sep.outer=0.02, sep.inner=0.2, border=NA, lwd=0.5) {
    datatype = sapply(list(...), function(x) x$type)
    groups.x = sapply(list(...), function(x) x$x)
    groups = sapply(list(...), function(x) x$color)
    keys = lapply(list(...), function(x) x$key)
    labs = sapply(list(...), function(x) x$label)
    margin = .CLUSTERMAP$margin
    
	# For the jth color bar, we have a label (labels[j]), a set of values
	# (groups.x[,j]) and a set of colors (groups[,j]). The link between 
	# values and colors is given by keys[[j]]. 
    
  	tmp = .CLUSTERMAP
  	tmp$cbar.key = keys
  	tmp$cbar.type = datatype
  	if (labels) {
  		names(tmp$cbar.key) = labs
  	} else {
  		names(tmp$cbar.key) = rep("", length(keys))
  	}
 	assign(".CLUSTERMAP", tmp, envir = .GlobalEnv)
    
    if (side==1) {
    	  if (!is.null(.CLUSTERMAP$colclust)) {
    	  	groups = groups[.CLUSTERMAP$colclust$order,,drop=F]
    	  	groups.x = groups.x[.CLUSTERMAP$colclust$order,,drop=F]
    	  }
    	  for (j in 1:ncol(groups)) { # Number of bars
    		delta = margin$inner[1]
    		for (i in 1:nrow(groups)) {
      	  		coli = groups[i,j]
      	  		a0 = c(i-1, i-1, i, i)/nrow(groups)
      	  		b0 = c(-sep.outer-delta*(j-1)/ncol(groups),
       	  			-sep.outer-delta*(j-sep.inner)/ncol(groups),
       	  			-sep.outer-delta*(j-sep.inner)/ncol(groups),
       	  			-sep.outer-delta*(j-1)/ncol(groups))
       	  	    polygon(a0, b0, col = coli, border = border, lwd=lwd)
    		}
    		pval = NA
    	 	if (pvalue & length(.CLUSTERMAP$colgroup)>0 & datatype[j]=="disc") {
    	 		yj = .CLUSTERMAP$colgroup[.CLUSTERMAP$colclust$order]
    	 		xj = groups.x[,j]
    	 		if (length(unique(yj)) < 2 || length(unique(xj)) < 2) {
    	 			pval = 1
    	 		} else if (pvalue.method == "chisq") {
    	 			pval = chisq.test(xj,yj)$p.value
    	 		} else if (pvalue.method == "chisq.sim") {
    	 			pval = chisq.test(xj,yj,simulate.p.value=T,B=5000)$p.value
    	 		} else if (pvalue.method == "fisher") {
    	 			pval = fisher.test(xj,yj,simulate.p.value=T,B=5000)$p.value
    	 		} 
    	 	} else if (pvalue & length(.CLUSTERMAP$colgroup)>0 & datatype[j]=="cont") {
    	 		yj = .CLUSTERMAP$colgroup[.CLUSTERMAP$colclust$order]
    	 		xj = groups.x[,j]
    	 		if (length(unique(yj)) < 2 || length(unique(xj)) < 2) {
    	 			pval = 1
    	 		} else {
    	 			datj = data.frame(xj=xj, yj=yj)
    	 			pval = summary(aov(xj ~ yj, data=datj))[[1]][1,5]
    	 		}
    	 	}  	
    		info.txt = ""
    		if (labels & !is.na(pval)) {
    			info.txt = paste(labs[j], " (", paste("P =", 
    				format(pval, digits=4, nsmall=4)), ")", sep="")
    		} else if (labels & is.na(pval)) {
    			info.txt = labs[j]
    		} else if (!labels & !is.na(pval)) {
    			info.txt = paste("P =", format(pval, digits=4, nsmall=4))
    		}
    		text(1.03, -sep.outer-delta*(j-0.5)/ncol(groups), info.txt, 
    			adj=c(0,0.5), cex=cex)
  	  	}
    } else if (side==2) {
    	  if (!is.null(.CLUSTERMAP$rowclust)) {
    	     groups = groups[.CLUSTERMAP$rowclust$order,,drop=F]
    	  }
      for (j in 1:ncol(groups)) { # Number of bars
    		delta = margin$inner[2]
    		for (i in 1:nrow(groups)) {
      	  		coli = groups[i,j]
       	  		a0 = c(-sep.outer-delta*(j-1)/ncol(groups),
       	  			-sep.outer-delta*(j-sep.inner)/ncol(groups),
       	  			-sep.outer-delta*(j-sep.inner)/ncol(groups),
       	  			-sep.outer-delta*(j-1)/ncol(groups))
          		b0 = c(i-1, i-1, i, i)/nrow(groups)
          		polygon(a0, b0, col = coli, border = border, lwd=lwd)
        	}
        }
    } else if (side==3) {
    	  if (!is.null(.CLUSTERMAP$colclust)) {
    	    groups = groups[.CLUSTERMAP$colclust$order,,drop=F]
    	    groups.x = groups.x[.CLUSTERMAP$colclust$order,,drop=F]  
    	  }  	  
    	  for (j in 1:ncol(groups)) { # Number of bars
    		delta = margin$inner[3]
    		for (i in 1:nrow(groups)) {
      	  		coli = groups[i,j]
      	  		a0 = c(i-1, i-1, i, i)/nrow(groups)
       	  		b0 = c(1+sep.outer+delta*(j-1)/ncol(groups),
       	  			1+sep.outer+delta*(j-sep.inner)/ncol(groups),
       	  			1+sep.outer+delta*(j-sep.inner)/ncol(groups),
       	  			1+sep.outer+delta*(j-1)/ncol(groups))
       	  		polygon(a0, b0, col = coli, border = border, lwd=lwd)
    		}
    		pval = NA
    		if (pvalue & length(.CLUSTERMAP$colgroup)>0 & datatype[j]=="disc") {
    	 		yj = .CLUSTERMAP$colgroup[.CLUSTERMAP$colclust$order]
    	 		xj = groups.x[,j]
    	 		if (length(unique(yj)) < 2 || length(unique(xj)) < 2) {
    	 			pval = 1
    	 		} else if (pvalue.method == "chisq") {
    	 			pval = chisq.test(xj,yj)$p.value
    	 		} else if (pvalue.method == "chisq.sim") {
    	 			pval = chisq.test(xj,yj,simulate.p.value=T,B=5000)$p.value
    	 		} else if (pvalue.method == "fisher") {
    	 			pval = fisher.test(xj,yj,simulate.p.value=T,B=5000)$p.value
    	 		} 
    	 	} else if (pvalue & length(.CLUSTERMAP$colgroup)>0 & datatype[j]=="cont") {
    	 		yj = .CLUSTERMAP$colgroup[.CLUSTERMAP$colclust$order]
    	 		xj = groups.x[,j]
    	 		if (length(unique(yj)) < 2 || length(unique(xj)) < 2) {
    	 			pval = 1
    	 		} else {
    	 			datj = data.frame(xj=xj, yj=yj)
    	 			pval = summary(aov(xj ~ yj, data=datj))[[1]][1,5]
    	 		}
    	 	} 
    		info.txt = ""
    		if (labels & !is.na(pval)) {
    			info.txt = paste(labs[j], " (", paste("P =", 
    				format(pval, digits=4, nsmall=4)), ")", sep="")
    		} else if (labels & is.na(pval)) {
    			info.txt = labs[j]
    		} else if (!labels & !is.na(pval)) {
    			info.txt = paste("P =", format(pval, digits=4, nsmall=4))
    		}
    		text(1.03, 1+sep.outer+delta*(j-0.5)/ncol(groups), info.txt, 
    			adj=c(0,0.5), cex=cex) 		    		
    	  }
    } else if (side==4) {
    	  if (!is.null(.CLUSTERMAP$rowclust)) {
      	  groups = groups[.CLUSTERMAP$rowclust$order,,drop=F]
      }
    	  for (j in 1:ncol(groups)) { # Number of bars
    		delta = margin$inner[4]
    		for (i in 1:nrow(groups)) {
      	  		coli = groups[i,j]
          		a0 = c(1+sep.outer+delta*(j-1)/ncol(groups),
       	  			1+sep.outer+delta*(j-sep.inner)/ncol(groups),
       	  			1+sep.outer+delta*(j-sep.inner)/ncol(groups),
       	  			1+sep.outer+delta*(j-1)/ncol(groups))
          		b0 = c(i-1, i-1, i, i)/nrow(groups)
          		polygon(a0, b0, col = coli, border = border, lwd=lwd)
          	}
        }
    }
}

# Use this function to plot heatmaps of a numerical data matrix X. It is
# recommended that you print to a file (e.g. a pdf-file) instead to 
# increase speed (and use fast=F). Use the argument colorscale to control
# how numerical values are represented as colors in the heatmap. The
# argument should be three (valid R) colors separated by hyphens. 
# Some examples: "green-black-red", "yellow-black-blue", "black-white-black"
# (to ignore signs), and "white-grey-darkgrey".

plot.hmap = function (X, colorscale = "blue-white-red", pch = ".", cex = 10, col.na="grey") {
	if (length(.CLUSTERMAP$colclust)>0) {
		X = X[,.CLUSTERMAP$colclust$order]
	}
	if (length(.CLUSTERMAP$rowclust)>0) {
		X = X[.CLUSTERMAP$rowclust$order,]
	}
	margin = .CLUSTERMAP$margin
    plot(0, 0, xlim = c(-margin$inner[2]-margin$outer[2], 1+margin$inner[4]+margin$outer[4]), 
    	ylim = c(-margin$inner[1]-margin$outer[1], 1+margin$inner[3]+margin$outer[3]), 
    	type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
    Nrow = nrow(X)
    Ncol = ncol(X)
    x0 = seq(0, 1, length = Ncol)
    y0 = seq(0, 1, length = Nrow)
    x1 = rep(x0, Nrow)
    y1 = rep(y0, rep(Ncol, Nrow))
    z1 = rbind(as.numeric(t(X))/max(abs(X),na.rm=T))
    c0 = col2rgb(unlist(strsplit(colorscale,"-")))/255
    c1 = matrix(c0[,2],3,length(z1),byrow=F) + 
    	cbind(c0[,1]-c0[,2])%*%ifelse(is.na(z1),NA,ifelse(z1<0,-z1,0)) + 
    	cbind(c0[,3]-c0[,2])%*%ifelse(is.na(z1),NA,ifelse(z1>0,z1,0))
    	c2 = rep(col.na, length(z1))
    	ok = !is.na(c1[1,]) & !is.na(c1[2,]) & !is.na(c1[3,])
    c2[ok] = rgb(c1[1,ok],c1[2,ok],c1[3,ok])
    for (i in 1:Nrow) {
        for (j in 1:Ncol) {
           polygon(c(j - 1, j - 1, j, j)/Ncol, c(i - 1, 
              i, i, i - 1)/Nrow, col = c2[(i - 1) * Ncol + 
              j], border = FALSE)
        }
    }
    cscale = list(colorscale = colorscale, valuerange = range(X,na.rm=T))
    assign("COLOR_SCALE", cscale, envir=.GlobalEnv)
    invisible()
}

plot.tree = function (groups, side=3, sep=0.03, lwd=0.5, trim=T, colors=NA) {
	if (side==1 || side==3) {
		clust = .CLUSTERMAP$colclust
		if (missing(groups)) {
			if (length(.CLUSTERMAP$colgroup)>0) {
				groups = .CLUSTERMAP$colgroup[.CLUSTERMAP$colclust$order]
			} else {
				groups = rep(0, length(clust$order))
			}
		}
	} else if (side==2 || side==4) {
		clust = .CLUSTERMAP$rowclust
		if (missing(groups)) {
			if (length(.CLUSTERMAP$rowgroup)>0) {
				groups = .CLUSTERMAP$rowgroup[.CLUSTERMAP$rowclust$order]
			} else {
				groups = rep(0, length(clust$order))
			}
		}		
	}
	margin = .CLUSTERMAP$margin
    d = margin$inner[side]; d2 = margin$outer[side]; g = as.numeric(groups)
	if (side == 1) {
      	invisible(plot.fork.below(length(clust$height), clust, trim, g, d, d2, sep, lwd))
    } else if (side == 2) {
        invisible(plot.fork.left(length(clust$height), clust, trim, g, d, d2, sep, lwd))
    } else if (side == 3) {
    		invisible(plot.fork.above(length(clust$height), clust, trim, g, d, d2, sep, lwd))
    } else if (side == 4) {
        invisible(plot.fork.right(length(clust$height), clust, trim, g, d, d2, sep, lwd))
    }
}

plot.fork.left = function (i, clust, trim, groups, delta, delta2, sep, lwd) {
    N = length(clust$height)
    mrg = clust$merge
    if (trim) {
    		hgt = (clust$height[i]-min(clust$height))/(max(clust$height)-min(clust$height)) 
    } else {
    		hgt = clust$height[i]/max(clust$height)
    }
    a = rep(0, 2)
    group = rep(0, 2)
    b1 = -delta - sep - delta2 * hgt
    for (j in 1:2) {
        if (mrg[i, j] < 0) {
            k = which(clust$order == -mrg[i, j])
            a[j] = (2*k-1)/(2*(N+1))
            b0 = -delta - sep
            group[j] = groups[k]
        }
        else {
            tmp = plot.fork.left(mrg[i, j], clust, trim, groups, delta, delta2, sep, lwd)
            a[j] = tmp[1]
            b0 = -delta - sep - delta2 * tmp[2]
            group[j] = tmp[3]
        }
        if (group[j] != -1) {
            lines(c(b0, b1), c(a[j], a[j]), col = group[j] + 1, lend=1, lwd=lwd)
        }
        else {
            lines(c(b0, b1), c(a[j], a[j]), lend=1, lwd=lwd)
        }
    }
    if (diff(group) == 0 & group[1] != -1) {
        lines(c(b1, b1), c(a[1], a[2]), col = group[1] + 1, lend=2, lwd=lwd)
        invisible(return(c(mean(a), hgt, group[1])))
    }
    else {
        lines(c(b1, b1), c(a[1], a[2]), lend=2, lwd=lwd)
        invisible(return(c(mean(a), hgt, -1)))
    }
}

plot.fork.right = function (i, clust, trim, groups, delta, delta2, sep, lwd) {
    N = length(clust$height)
    mrg = clust$merge
    if (trim) {
    		hgt = (clust$height[i]-min(clust$height))/(max(clust$height)-min(clust$height)) 
    } else {
    		hgt = clust$height[i]/max(clust$height)
    }
    a = rep(0, 2)
    group = rep(0, 2)
    b1 = 1+delta+sep+delta2*hgt
    for (j in 1:2) {
        if (mrg[i, j] < 0) {
            k = which(clust$order == -mrg[i, j])
            a[j] = (2*k-1)/(2*(N+1))
            b0 = 1+delta+sep
            group[j] = groups[k]
        }
        else {
            tmp = plot.fork.right(mrg[i, j], clust, trim, groups, delta, delta2, sep, lwd)
            a[j] = tmp[1]
            b0 = 1+delta+sep+delta2*tmp[2]
            group[j] = tmp[3]
        }
        if (group[j] != -1) {
            lines(c(b0, b1), c(a[j], a[j]), col = group[j] + 1, lend=1, lwd=lwd)
        }
        else {
            lines(c(b0, b1), c(a[j], a[j]), lend=1, lwd=lwd)
        }
    }
    if (diff(group) == 0 & group[1] != -1) {
        lines(c(b1, b1), c(a[1], a[2]), col = group[1] + 1, lend=2, lwd=lwd)
        invisible(return(c(mean(a), hgt, group[1])))
    }
    else {
        lines(c(b1, b1), c(a[1], a[2]), lend=2, lwd=lwd)
        invisible(return(c(mean(a), hgt, -1)))
    }
}

plot.fork.above = function (i, clust, trim, groups, delta, delta2, sep, lwd) {
    N = length(clust$height)
    mrg = clust$merge
    if (trim) {
    		hgt = (clust$height[i]-min(clust$height))/(max(clust$height)-min(clust$height)) 
    } else {
    		hgt = clust$height[i]/max(clust$height)
    }
    a = rep(0, 2)
    group = rep(0, 2)
    b1 = 1 + delta + sep + delta2 * hgt
    for (j in 1:2) {
        if (mrg[i,j] < 0) { # Leaf node
            k = which(clust$order == -mrg[i,j])
            a[j] = (2*k-1)/(2*(N+1))
            b0 = 1 + delta + sep
            group[j] = groups[k]
        } else { # Inner node
            tmp = plot.fork.above(mrg[i,j], clust, trim, groups, delta, delta2, sep, lwd)
            a[j] = tmp[1]
            b0 = 1 + delta + sep + delta2 * tmp[2]
            group[j] = tmp[3]
        }
        if (group[j] != -1) {
            lines(c(a[j], a[j]), c(b0, b1), col = group[j] + 1, lend=1, lwd=lwd)
        } else {
            lines(c(a[j], a[j]), c(b0, b1), lend=1, lwd=lwd)
        }
    }
    if (diff(group) == 0 & group[1] != -1) {
        lines(c(a[1], a[2]), c(b1, b1), col = group[1] + 1, lend=2, lwd=lwd)
        invisible(return(c(mean(a), hgt, group[1])))
    } else {
        lines(c(a[1], a[2]), c(b1, b1), lend=2, lwd=lwd)
        invisible(return(c(mean(a), hgt, -1)))
    }
}

plot.fork.below = function (i, clust, trim, groups, delta, delta2, sep, lwd) {
    N = length(clust$height)
    mrg = clust$merge
    if (trim) {
    		hgt = (clust$height[i]-min(clust$height))/(max(clust$height)-min(clust$height)) 
    } else {
    		hgt = clust$height[i]/max(clust$height)
    }
    a = rep(0, 2)
    group = rep(0, 2)
    b1 = -delta-sep-delta2*hgt
    for (j in 1:2) {
        if (mrg[i, j] < 0) {
            k = which(clust$order == -mrg[i, j])
            a[j] = (2*k-1)/(2*(N+1))
            b0 = -delta-sep
            group[j] = groups[k]
        }
        else {
            tmp = plot.fork.below(mrg[i, j], clust, trim, groups, delta, delta2, sep, lwd)
            a[j] = tmp[1]
            b0 = -delta-sep-delta2*tmp[2]
            group[j] = tmp[3]
        }
        if (group[j] != -1) {
            lines(c(a[j], a[j]), c(b0, b1), col = group[j] + 1, lend=1, lwd=lwd)
        }
        else {
            lines(c(a[j], a[j]), c(b0, b1), lend=1, lwd=lwd)
        }
    }
    if (diff(group) == 0 & group[1] != -1) {
        lines(c(a[1], a[2]), c(b1, b1), col = group[1] + 1, lend=2, lwd=lwd)
        invisible(return(c(mean(a), hgt, group[1])))
    }
    else {
        lines(c(a[1], a[2]), c(b1, b1), lend=2, lwd=lwd)
        invisible(return(c(mean(a), hgt, -1)))
    }
}

gap.new = function(X, clust="col", distance="euclidean", linkage="complete", B=100) {
	if (clust=="col") {X = t(X)}
	if (distance=="euclidean") {
		res = gap.euclidean(X, linkage, B)$group
	} else if (distance=="pearson") {
		res = gap.pearson(X, linkage, B)$group
	} else if (distance=="spearman") {
		res = gap.spearman(X, linkage, B)$group
	}
	return(res)
}

part.new = function(X, clust="col", distance="euclidean", linkage="complete", B=100, minSize=10, maxLevel=3) {
	if (clust=="col") {X = t(X)}
	if (distance=="euclidean") {
     	res = part.euclidean(X, linkage, B, minSize, maxLevel)
  	} else if (distance=="pearson") {
		res = part.pearson(X, linkage, B, minSize, maxLevel)
	} else if (distance=="spearman") {
		res = part.spearman(X, linkage, B, minSize, maxLevel)
	}
    return(res)
}

part.euclidean = function(X, linkage, B, minSize, maxLevel, K=6) {
	if (nrow(X) < 2*minSize) {
		return(rep(1,nrow(X)))
	}
	if (maxLevel == 0) {
	    return(rep(1,nrow(X)))
	}
	K = min(nrow(X),K)
	k.opt = gap.euclidean(X, linkage, B, K)$k
	Dobs = dist(X, method="euclidean")
	Popt = cutree(hclust(Dobs, method=linkage), k=k.opt) # N x 1 matrix
	if (sum(table(Popt)>=minSize) < 2) {
		Popt = rep(1,nrow(X))
		k.opt = 1
	}
	group = rep(NA,nrow(X))
	if (k.opt > 1) {
		group.max = 0
		for (i in 1:k.opt) {
			group[Popt==i] = group.max + part.euclidean(X[Popt==i,,drop=F], linkage, B, minSize, maxLevel-1, K)
			group.max = max(group.max, group, na.rm=T)
		}
	} else {
		P2 = cutree(hclust(Dobs, method=linkage), k=2)
	    P2.1 = part.euclidean(X[P2==1,,drop=F], linkage, B, minSize, maxLevel-1, min(sum(P2==1), K))
		P2.2 = part.euclidean(X[P2==2,,drop=F], linkage, B, minSize, maxLevel-1, min(sum(P2==2), K))
		if (all(c(P2.1,P2.2) == 1)) {
			group[1:nrow(X)] = 1
		} else {
			group[P2==1] = P2.1
			group[P2==2] = max(P2.1,na.rm=T) + P2.2
		}
	}
	return(group)
}

part.pearson = function(X, linkage, B, minSize, maxLevel, K=6) {
	if (nrow(X) < 2*minSize) {
		return(rep(1,nrow(X)))
	}
	if (maxLevel == 0) {
	    return(rep(1,nrow(X)))
	}
	K = min(nrow(X),K)
	k.opt = gap.pearson(X, linkage, B, K)$k
	Dobs = as.dist((1-cor(t(X), method="pearson"))/2)
	Popt = cutree(hclust(Dobs, method=linkage), k=k.opt) # N x 1 matrix
	if (sum(table(Popt)>=minSize) < 2) {
		Popt = rep(1,nrow(X))
		k.opt = 1
	}
	group = rep(NA,nrow(X))
	if (k.opt > 1) {
		group.max = 0
		for (i in 1:k.opt) {
			group[Popt==i] = group.max + part.pearson(X[Popt==i,,drop=F], linkage, B, minSize, maxLevel-1, K)
			group.max = max(group.max, group, na.rm=T)
		}
	} else {
		P2 = cutree(hclust(Dobs, method=linkage), k=2)
	    P2.1 = part.pearson(X[P2==1,,drop=F], linkage, B, minSize, maxLevel-1, min(sum(P2==1), K))
		P2.2 = part.pearson(X[P2==2,,drop=F], linkage, B, minSize, maxLevel-1, min(sum(P2==1), K))
		if (all(c(P2.1,P2.2) == 1)) {
			group[1:nrow(X)] = 1
		} else {
			group[P2==1] = P2.1
			group[P2==2] = max(P2.1,na.rm=T) + P2.2
		}
	}
	return(group)
}

part.spearman = function(X, linkage, B, minSize, maxLevel, K=6) {
	if (nrow(X) < 2*minSize) {
		return(rep(1,nrow(X)))
	}
	if (maxLevel == 0) {
	    return(rep(1,nrow(X)))
	}
	K = min(nrow(X),K)
	k.opt = gap.spearman(X, linkage, B, K)$k
	Dobs = as.dist((1-cor(t(X), method="spearman"))/2)
	Popt = cutree(hclust(Dobs, method=linkage), k=k.opt) # N x 1 matrix
	if (sum(table(Popt)>=minSize) < 2) {
		Popt = rep(1,nrow(X))
		k.opt = 1
	}
	group = rep(NA,nrow(X))
	if (k.opt > 1) {
		group.max = 0
		for (i in 1:k.opt) {
			group[Popt==i] = group.max + part.spearman(X[Popt==i,,drop=F], linkage, B, minSize, maxLevel-1, K)
			group.max = max(group.max, group, na.rm=T)
		}
	} else {
		P2 = cutree(hclust(Dobs, method=linkage), k=2)
	    P2.1 = part.spearman(X[P2==1,,drop=F], linkage, B, minSize, maxLevel-1, min(sum(P2==1), K))
		P2.2 = part.spearman(X[P2==2,,drop=F], linkage, B, minSize, maxLevel-1, min(sum(P2==1), K))
		if (all(c(P2.1,P2.2) == 1)) {
			group[1:nrow(X)] = 1
		} else {
			group[P2==1] = P2.1
			group[P2==2] = max(P2.1,na.rm=T) + P2.2
		}
	}
	return(group)
}

gap.euclidean = function(X, linkage, B, K=6) {
	N = nrow(X)
	# Xobs = scale(X, scale=F)  ## REMOVED
	
	# Find W1,...,WK
	Wobs = rep(0, K)
	Dobs = dist(X, method="euclidean")
	Pobs = cutree(hclust(Dobs, method=linkage), 1:K) # N x K matrix
	for (k in 1:K) {
		Dk = as.matrix(Dobs)
		for (r in 1:k) {
			keep = (Pobs[,k]==r)
			Wobs[k] = Wobs[k] + sum(Dk[keep,keep])/(2*sum(keep))
		}
	}
		
	# Find W1*,...,WK* for B simulated data sets
	Xmean = outer(rep(1,nrow(X)), apply(X, 2, mean))
	Xobs = X - Xmean
	V = eigen(t(Xobs)%*%Xobs)$vectors  # svd(Xobs)$v
	R = apply(Xobs %*% V, 2, range)
	Wsim = matrix(0, B, K)
	# for (b in 1:B) {
	Wsim.list = foreach(b=1:B) %dopar% {
		print(paste("Permutation", b))
		# Simulate data and find partitions into 1,2,...,K clusters
		Xsim = Xmean + apply(R,2,function(x)runif(N,x[1],x[2])) %*% t(V)
		Dsim = dist(Xsim, method="euclidean") # Uses the most time
		Psim = cutree(hclust(Dsim, method=linkage), 1:K)
		wsim = rep(0,K)
		for (k in 1:K) {
			Dk = as.matrix(Dsim)
			for (r in 1:k) {
				keep = (Psim[,k]==r)
				# Wsim[b,k] = Wsim[b,k] + sum(Dk[keep,keep])/(2*sum(keep))
				wsim[k] = wsim[k] + sum(Dk[keep,keep])/(2*sum(keep))
			}
		}
		wsim
	}
	for (b in 1:B) Wsim[b,] = Wsim.list[[b]]
	
	# Calculate gap statistic
	res = apply(log(Wsim),2,mean) - log(Wobs)
	
	# Find optimal k
	sk = sqrt(1+1/B) * apply(log(Wsim),2,sd)
	k.opt = which(diff(res)-sk[-1] <= 0)[1]
	if (is.na(k.opt) || length(k.opt)==0) k.opt=1
	
	# Return result
	list(k=k.opt, gap=res, sk=sk, group=Pobs[,k.opt])
}

gap.pearson = function(X, linkage, B, K=6)  {
	N = nrow(X)
	# Xobs = scale(X, scale=F)
	
	# Find W1,...,WK
	Wobs = rep(0, K)
	Dobs = as.dist((1-cor(t(X), method="pearson"))/2)
	Pobs = cutree(hclust(Dobs, method=linkage), 1:K) # N x K matrix
	for (k in 1:K) {
		Dk = as.matrix(Dobs)
		for (r in 1:k) {
			keep = (Pobs[,k]==r)
			Wobs[k] = Wobs[k] + sum(Dk[keep,keep])/(2*sum(keep))
		}
	}
		
	# Find W1*,...,WK* for B simulated data sets
	Xmean = outer(rep(1,nrow(X)), apply(X, 2, mean))  ## NEW
	Xobs = X - Xmean  ## NEW
	V = eigen(t(Xobs)%*%Xobs)$vectors  # svd(Xobs)$v
	R = apply(Xobs %*% V, 2, range)
	Wsim = matrix(0, B, K)
	# for (b in 1:B) {
	Wsim.list = foreach(b=1:B) %dopar% {
	    print(paste("Permutation", b))
		# Simulate data and find partitions into 1,2,...,K clusters
		Xsim = Xmean + apply(R,2,function(x)runif(N,x[1],x[2])) %*% t(V)
		Dsim = as.dist((1-cor(t(Xsim), method="pearson"))/2)		
		Psim = cutree(hclust(Dsim, method=linkage), 1:K)
		wsim = rep(0,K)
		for (k in 1:K) {
			Dk = as.matrix(Dsim)
			for (r in 1:k) {
				keep = (Psim[,k]==r)
				# Wsim[b,k] = Wsim[b,k] + sum(Dk[keep,keep])/(2*sum(keep))
				wsim[k] = wsim[k] + sum(Dk[keep,keep])/(2*sum(keep))
			}
		}
		wsim
	}
	for (b in 1:B) Wsim[b,] = Wsim.list[[b]]
	
	# Calculate gap statistic
	res = apply(log(Wsim),2,mean) - log(Wobs)
	
	# Find optimal k
	sk = sqrt(1+1/B) * apply(log(Wsim),2,sd)
	k.opt = which(diff(res)-sk[-1] <= 0)[1]
	
	# Return result
	list(k=k.opt, gap=res, sk=sk, group=Pobs[,k.opt])
}

gap.spearman = function(X, linkage, B, K=6)  {
	N = nrow(X)
	# Xobs = scale(X, scale=F)
	
	# Find W1,...,WK
	Wobs = rep(0, K)
	Dobs = as.dist((1-cor(t(X), method="spearman"))/2)
	Pobs = cutree(hclust(Dobs, method=linkage), 1:K) # N x K matrix
	for (k in 1:K) {
		Dk = as.matrix(Dobs)
		for (r in 1:k) {
			keep = (Pobs[,k]==r)
			Wobs[k] = Wobs[k] + sum(Dk[keep,keep])/(2*sum(keep))
		}
	}
		
	# Find W1*,...,WK* for B simulated data sets
	Xmean = outer(rep(1,nrow(X)), apply(X, 2, mean))  ## NEW
	Xobs = X - Xmean  ## NEW
	V = eigen(t(Xobs)%*%Xobs)$vectors  # svd(Xobs)$v
	R = apply(Xobs %*% V, 2, range)
	Wsim = matrix(0, B, K)
	# for (b in 1:B) {
	Wsim.list = foreach(b=1:B) %dopar% {
		print(paste("Permutation", b))
		# Simulate data and find partitions into 1,2,...,K clusters
		Xsim = Xmean + apply(R,2,function(x)runif(N,x[1],x[2])) %*% t(V)
		Dsim = as.dist((1-cor(t(Xsim), method="spearman"))/2)		
		Psim = cutree(hclust(Dsim, method=linkage), 1:K)
		wsim = rep(0,K)
		for (k in 1:K) {
			Dk = as.matrix(Dsim)
			for (r in 1:k) {
				keep = (Psim[,k]==r)
				# Wsim[b,k] = Wsim[b,k] + sum(Dk[keep,keep])/(2*sum(keep))
				wsim[k] = wsim[k] + sum(Dk[keep,keep])/(2*sum(keep))
			}
		}
		wsim
	}
	for (b in 1:B) Wsim[b,] = Wsim.list[[b]]
	
	# Calculate gap statistic
	res = apply(log(Wsim),2,mean) - log(Wobs)
	
	# Find optimal k
	sk = sqrt(1+1/B) * apply(log(Wsim),2,sd)
	k.opt = which(diff(res)-sk[-1] <= 0)[1]
	
	# Return result
	list(k=k.opt, gap=res, sk=sk, group=Pobs[,k.opt])
}

