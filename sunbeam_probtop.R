
#'Makes a sunbeam plot
#'Shows the probability of the input SNP having the lowest p-value
#'
#'
#'@title sunbeam_prob
#'@export
#'@param CV1 The first of the true causal SNPs
#'@param CV2 The second of the true causal SNPs
#'@param SNPint The SNP we are analysing 
#'@param N0 The number of control samples
#'@param N1 The number of case samples
#'@param gammatrue The odds ratio in the true causal model
#'@param LD The LD reference dataset (a matrix)
#'@param freq The frequency reference dataset
#'@param ORmax The maximum OR we want to consider
#'@param ORgrid The spacing between ORs
#'@param z_score_sig The Z Score at which we conclude that no SNP is top
#'@param col_list Which colours to use in the plot
#'@return a sunbeam plot showing the probability of SNP SNPint having the lowest p-value
#'@author Mary Fortune
sunbeam_probtop<-function(CV1,CV2,SNPint,N0,N1,gammatrue,LD,freq,ORmax=1,ORgrid=0.05,z_score_sig=-qnorm((0.1)/2),col_list=c('white', 'purple')){
	#get the SNP names
	snps<-colnames(freq)[-ncol(freq)]
	#check that the same SNPs are present in the freq data, the LD matrix and the input SNPs
	if (!setequal(snps,colnames(LD))){stop("Reference data mismatch: the freq data and the LD data do not contain the same SNPs")}
	if (!is.element(CV1,snps)){stop("CV1 is not in the reference data")}
	if (!is.element(CV2,snps)){stop("CV2 is not in the reference data")}
	if (!is.element(SNPint,snps)){stop("SNPint is not in the reference data")}
	#set up the grid causal models which we will test
	gamma_list<-seq(-ORmax,ORmax,by=ORgrid)
	n_gamma<-length(gamma_list)
	#compute the list of probabilities of genotype given causal SNP genotype
	GenoProbList<-make_GenoProbList(snps,c(CV1,CV2),freq)
	#matrix giving the probability that SNPint will be the top SNP at each point on the grid
	prob_SNPint_top<-matrix(-1,n_gamma,n_gamma)
	#populate top_SNP
	for (ii in 1:n_gamma){
		for (jj in 1:n_gamma){
			#gamma at this point on the grid
			try_gamma<-c(gamma_list[ii],gamma_list[jj])
			#expected Z Score at this gamma
			est_z_score<-abs(est_statistic(N0,N1,snps,c(CV1,CV2),try_gamma,freq,GenoProbList))
			#prob SNPint comes top
			if (max(est_z_score)<z_score_sig){
				prob_SNPint_top[ii,jj]<-0
			}else{
				prob_SNPint_top[ii,jj]<-prob_snp_top(which(snps==SNPint),est_z_score,LD)
		      	}
		}
	}
	#which colour will correspond to SNPint being top
	col.l <- colorRampPalette(col_list)
	#we only want to label the plot at a subset of points
	subset<-which(abs(round(gamma_list,1)-gamma_list)<0.001)
	#make the plot, including a point corresponding to gammatrue
	levelplot(prob_SNPint_top,scales=list(x=list(at=subset, labels=gamma_list[subset]),y=list(at=subset, labels=gamma_list[subset])),xlab=list("log OR for CV1",cex=1.5),ylab=list("log OR for CV2",cex=1.5),col.regions=col.l,panel = function(...){
		panel.levelplot(...)
		panel.abline(h = which(gamma_list==0))
            	panel.abline(v = which(gamma_list==0))
		panel.xyplot((gammatrue[1]+ORmax+ORgrid)/ORgrid,(gammatrue[2]+ORmax+ORgrid)/ORgrid, cex = 1, col = 1,pch=19)
        })
}


##' Internal function, prob_snp_top
##'
##' This function takes in a SNP and the model for the Z Score
##' and computes the probability of this SNP having the lowest p-value
##' @title prob_snp_top
##' @param whichsnp the SNP we are interested in
##' @param exp_z_score a vector giving the expected Z Score
##' @param LD a matrix giving the LD between SNPs in this region 
##' @return a vector giving the probability of this SNP coming top 
##' @author Mary Fortune
prob_snp_top<-function(whichsnp,exp_z_score,LD){
	num_top<-0
	num_itt<-1000
	sim_z_score<-abs(rmvnorm(n=num_itt,mean=exp_z_score,sigma=LD))
	num_top<-length(which(apply(sim_z_score,1,which.max)==whichsnp))
	return(num_top/num_itt)
}

