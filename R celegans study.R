## In this project my current R skills will be displayed. This is done by using data from an article;
## Jovic et al. 2017  (https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0189445)
## 
## The effect of heat stress on gene-expression in C. elegans will be investigated in this project.
## In short, this dataset contains expression levels of C. elegans before and during heat stress.
## Four replicates of the time points 0, 30, 60 and 120 minutes are used. 
##
## The following questions will be answered 
## 1) First the data will be checked and formatted if needed 
## 2) 


#################################### Start ##################################################################

######################## Checking and formatting the data ---------------------------------------------

### Load the expression data in file "Jovic_etal.txt"
expdat <- read.delim(file="Jovic_etal.txt",header = T,row.names = 1)
  ## Check the loaded data: structure, dimensions, row- and col-names and numbers.
  str(expdat)
  dim(expdat)
  nrow(expdat)
  ncol(expdat)
  rownames(expdat)
  colnames(expdat)

  ## To check whether the genIDs are ok
  head(expdat)

### Load the gene info in the file "Gene_info.txt"
  ginfo <- read.delim(file="Gene_info.txt",header = T,row.names = 1)
  ginfo
  ## To check the loaded data on structure, dimensions, row- and col-names and numbers.
  str(ginfo)
  dim(ginfo)
  nrow(ginfo)
  ncol(ginfo)
  rownames(ginfo)
  colnames(ginfo)
  
  ## To check whether geneIDs are ok
  head(ginfo)  
  
### To check if the order is the same for the expression dataa and gene info  
  rownames(ginfo) == rownames(expdat) 
  ## To check if sum of vector is equal to the number of rows
  check.row <- rownames(ginfo) == rownames(expdat)
  sum(check.row)
  nrow(ginfo)

  sum(rownames(ginfo) == rownames(expdat)) == nrow(ginfo) 


######################################################################################################
##################### Expression distribution ---------------------------------------------------------

### Find minimum, maximum, mean and median expression level per sample  
  min.per.samp <- apply(expdat,2,min)
  min.per.samp
  max.per.samp <- apply(expdat,2,max)
  max.per.samp
  mean.per.samp <- apply(expdat,2,mean)
  mean.per.samp
  med.per.samp <- apply(expdat,2,median)
  med.per.samp
  
  ## Display these values in a barplot() function 
  par(mfrow=c(2,2))
  barplot(min.per.samp,main = "Minimum")
  barplot(max.per.samp,main = "Maximum")
  barplot(mean.per.samp,main = "Mean")
  barplot(med.per.samp,main = "Median")
  

  tp <- rep(c(0,30,60,120),each=4)

### plot the expression values of the 513th gene against the timepoints first before doing the rest.
  x <- as.numeric(expdat[513,])  
  plot(tp,x,xlab = "Heat exposure (min)",ylab = "Gene expression (log2)",main = "Gene 513",cex=2,pch=19,col="orchid")

  boxplot(x~tp,col="grey",xlab = "Heat exposure (min)",ylab = "Gene expression (log2)",main = "Gene 513")


## make a 4 x4 plot with the histograms for the other samples
  par(mfrow=c(4,4),mar=c(1,2,3,1))
  for ( i in 1:16){
  hist(expdat[,i],breaks = 100,col="grey",main = tp[i],xlab = "Heat exposure (min)")
  }

  
  plot(density(expdat[,1]),xlab = "Heat exposure (min)",main = "Gene expression distribution")
  for ( i in 1:16){
    lines(density(expdat[,i]),col=rainbow(16)[i])
  }
  legend(15,0.35,legend = tp,fill = rainbow(16))
  

  m0 <- apply(expdat[,tp==0],1,mean)
  m30 <- apply(expdat[,tp==30],1,mean)
  m60 <- apply(expdat[,tp==60],1,mean)
  m120 <- apply(expdat[,tp==120],1,mean) 

  mexp <- data.frame(m0,m30,m60,m120)
  mexp

  rownames(mexp)
  
## to reduce the computational load, the genes that show variation over time will be used
  
  gvar <- apply(mexp,1,var)
  
  ## plot distribution of variances
  plot(density((gvar)))
  
  ## genes with a highe var of 0.1
  sum(gvar>0.1)

  var.selc <- gvar > 0.1
  
## plot the mean values of the t = 30 min timepoint against the t = 0 timepoint
  x <- mexp$m0[var.selc] ## timepoint 0, only those that have var > 0.1
  y <- mexp$m30[var.selc] ## timepoint 30, only those that have var > 0.1
  plot(x,y,xlab="Mean expression t=0",ylab = "Mean expression t=30",pch=19,cex=0.25)

## a 1 to 1 trendline 
  abline(coef = c(0,1),col="blue",lwd=3)
  
## Highlight genes that are up- or down-regulated and add trendlines for the +1 and -1 borders.
 
  abline(coef = c(1,1),col="red")
  abline(coef = c(-1,1),col="red")
  
## Get the names of the 10 most upregulated gene in this timepoint
  
  lograt <- log2(mexp$m30/mexp$m0)
  lograt

  top10.up<- ginfo[order(lograt,decreasing = T),][1:10,]
  top10.up
  
  
  x.text <- mexp[rownames(top10.up),]$m0
  y.text <- mexp[rownames(top10.up),]$m30
 
  use.text <- ginfo[rownames(top10.up),]$Public_name
  text(x.text,y.text,use.text)  
  
## most of the genes found in the top 10 are from the heat shock protein family, which makes sense in this case
  
## Restults of hsp-16.2 on t=0 and t= 30 minutes will be statistically tested 
  ## make two vectors one with the gene expression levels of hsp-16.2 on t=0 and one with the expression levels of t=30
  x0 <- expdat[ginfo$Public_name == "hsp-16.2",tp==0]
  x30 <- expdat[ginfo$Public_name == "hsp-16.2",tp==30]

  t.test(x0,x30)

  t.results <- t.test(x0,x30)
  t.results
 
  str(t.results)
 
  t.results$p.value
  t.results$estimate
  
  ## for() loop to calculate and record the p-value between t=0 and t=30
  
  t.res.mat <- matrix(data=NA,nrow = nrow(expdat),ncol = 3)
  rownames(t.res.mat) <- rownames(expdat)
  colnames(t.res.mat) <- c("p","est.x","est.y")
  for ( i in 1:nrow(expdat)){
    x0 <- expdat[i,tp==0]
    x30 <- expdat[i,tp==30]
    t.results <- t.test(x0,x30)
    t.res.mat[i,] <- c(t.results$p.value,t.results$estimate)
  }
  
  ## data.frame from the result matrix s
  
  t.res.mat <- data.frame(t.res.mat)
  
  ## to check the first few lines
  t.res.mat[1:5,]
  

  log.p <- -log10(t.res.mat$p)
  
  ## genes with a significantly different threshold at -log10(p) > 2
  sum(log.p>2)
  

  sig.genes <- ginfo[rownames(expdat)[log.p>2],]
  grep("hsp",sig.genes$Public_name)
  
  
  ## scatter plot
  var.selc <- gvar > 0.1
  x <- mexp$m0[var.selc] ## timepoint 0, only those that have var > 0.1
  y <- mexp$m30[var.selc] ## timepoint 30, only those that have var > 0.1
  plot(x,y,xlab="Mean expression t=0",ylab = "Mean expression t=30",pch=19,cex=0.25)
  abline(coef = c(0,1),col="blue",lwd=3)
  abline(coef = c(1,1),col="red")
  abline(coef = c(-1,1),col="red")
  sx <- mexp$m0[log.p>2]
  sy <- mexp$m30[log.p>2]
  points(sx,sy,col="red",pch=19,cex=1)
   
  legend(12,8,legend =c("significant","not-significant"),fill = c("red","black"))
  
  
  ## to explore the relation between effect and significance a volcano plot is suitable
  
  ## return ot orginial
 
  lograt <- log2(2^t.res.mat$est.y/2^t.res.mat$est.x)
  plot(lograt,log.p,xlab="Effect ratio (log2)",ylab = "Significance -log10(p)",pch=19,cex=0.5)
  
  abline(h=2,col="red",lty="dashed")
  
  
  do.t.test <- function(exp.val,grp1,grp2){
              x <- exp.val[grp1]
              y <- exp.val[grp2]
              res <- t.test(x,y)
              res <- c(res$p.value,res$estimate) 
              return(res)
  }
  
  ## test function with gene in row 513
  do.t.test(expdat[513,],tp==0,tp==30)
  do.t.test(expdat[513,],tp==0,tp==60)
  do.t.test(expdat[513,],tp==0,tp==120)
  do.t.test(expdat[513,],tp==30,tp==60)
  
  
  ## appy function to all gene tp=60
  res2 <- apply(expdat,1,do.t.test,tp==0,tp==60)
  
  dim(res2)
  ## flip the rows and columns

  res2 <- t(res2)
  
 
  colnames(res2) <- c("p","est.x","est.y")
  res2 <- data.frame(res2)
  res2[1:5,]
  
  ### Make a volcano plot 
  lograt <- log2(2^res2$est.y/2^res2$est.x)
  log.p <- -log10(res2$p)
  plot(lograt,log.p,xlab="Effect ratio (log2)",ylab = "Significance -log10(p)",pch=19,cex=0.5)
 abline(h=2,col="red",lty="dashed")
  

 ### generate heatmap
 
 ### select genes that are significantly different between t=0 and t=60 (-log10(p)>4)
 sig.genes <- rownames(res2)[log.p>4]
 ## for the signifianct genes get the expression data from all samples
 sig.data <- expdat[sig.genes,]
 

 sig.data <- as.matrix(sig.data)
 heatmap(sig.data)
 
 ### explore the settings of the heatmap() function, colors, layout etc.