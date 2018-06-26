#'Makes a sunbeam plot
#'Shows which SNP is expected to have the lowest p-value
#'
#'
#'@title sunbeam_whichtop
#'@export
#'@param CV1 The first of the true causal SNPs
#'@param CV2 The second of the true causal SNPs
#'@param Tag The Tag SNP
#'@param N0 The number of control samples
#'@param N1 The number of case samples
#'@param gammatrue The odds ratio in the true causal model
#'@param LD The LD reference dataset (a matrix)
#'@param freq The frequency reference dataset
#'@param ORmax The maximum OR we want to consider
#'@param ORgrid The spacing between ORs
#'@param z_score_sig The Z Score at which we conclude that no SNP is top
#'@param col_list Which colours to use in the plot
#'@return a sunbeam plot showing lowest expected p-value SNPs
#'@author Mary Fortune
sunbeam_whichtop<-function(CV1,CV2,Tag,N0,N1,gammatrue,LD,freq,ORmax=1,ORgrid=0.05,z_score_sig=-qnorm((0.1)/2),col_list=c('white', 'firebrick1', 'firebrick4', 'gold','blue')){
	#get the SNP names
	snps<-colnames(freq)[-ncol(freq)]
	#check that the same SNPs are present in the freq data, the LD matrix and the input SNPs
	if (!setequal(snps,colnames(LD))){stop("Reference data mismatch: the freq data and the LD data do not contain the same SNPs")}
	if (!is.element(CV1,snps)){stop("CV1 is not in the reference data")}
	if (!is.element(CV2,snps)){stop("CV2 is not in the reference data")}
	if (!is.element(Tag,snps)){stop("Tag is not in the reference data")}
	#set up the grid causal models which we will test
	gamma_list<-seq(-ORmax,ORmax,by=ORgrid)
	n_gamma<-length(gamma_list)
	#compute the list of probabilities of genotype given causal SNP genotype
	GenoProbList<-make_GenoProbList(snps,c(CV1,CV2),freq)
	#matrix giving which SNP will come top in expectation at each point on the grid
	top_SNP<-matrix(-1,n_gamma,n_gamma)
	#populate top_SNP
	for (ii in 1:n_gamma){
		for (jj in 1:n_gamma){
			#gamma at this point on the grid
			try_gamma<-c(gamma_list[ii],gamma_list[jj])
			#expected Z Score at this gamma
			est_z_score<-abs(est_statistic(N0,N1,snps,c(CV1,CV2),try_gamma,freq,GenoProbList))
			#which SNP comes top
			GWAS_hit<-which.max(est_z_score)
			if (max(est_z_score)<z_score_sig){
				top_SNP[ii,jj]<-0
			}else if(GWAS_hit==which(snps==CV1)){
				top_SNP[ii,jj]<-1
		      	}else if(GWAS_hit==which(snps==CV2)){
				top_SNP[ii,jj]<-2
		      	}else if(GWAS_hit==which(snps==Tag)){
				top_SNP[ii,jj]<-3
		      	}else{
				top_SNP[ii,jj]<-4
		      	}
		}
	}
	#which colour will correspond to each SNP being top
	col.l <- colorRampPalette(col_list[1:(max(top_SNP)+1)])
	#we only want to label the plot at a subset of points
	subset<-which(abs(round(gamma_list,1)-gamma_list)<0.001)
	#make the plot, including a point corresponding to gammatrue
	levelplot(top_SNP,scales=list(x=list(at=subset, labels=gamma_list[subset]),y=list(at=subset, labels=gamma_list[subset])),xlab=list("log OR for CV1",cex=1.5),ylab=list("log OR for CV2",cex=1.5),col.regions=col.l,colorkey=FALSE, 		panel = function(...){
		panel.levelplot(...)
		panel.abline(h = which(gamma_list==0))
            	panel.abline(v = which(gamma_list==0))
		panel.xyplot((gammatrue[1]+ORmax+ORgrid)/ORgrid,(gammatrue[2]+ORmax+ORgrid)/ORgrid, cex = 1, col = 1,pch=19)
        })
}
