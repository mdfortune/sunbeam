#plotting th specific cases Jenn asked for on 15/1/18

#Sunbeam plot with SNPs 
#A = rs61839660, beta = -0.18649 (imm_10_6134703)
#D = rs62626317, beta = -0.22831 (imm_10_6150931)

library(lattice)
library(snpStats)
library(data.table)
library(snpStatsWriter)
library(plyr)
library(dplyr)
library(mlogitBMA)
library(nnet)
library(combinat)
library(boot)
library(Rcpp)
library(corpcor)
library(mvtnorm)
library(weights)
library(RColorBrewer)

library(devtools)
install_github("mdfortune/simGWAS")
library(simGWAS)

#load the freq file
load("~/freqdatasets/T1D_regions/D_chr10_6068495_6237542_freq.RData")

#and get the SNP names
snps<-colnames(freq)[-ncol(freq)]

#What significance level will be colour as white
#z_score_sig<--qnorm((5*10^-8)/2)
z_score_sig<--qnorm((0.1)/2)


#How large do we want the grid to be, and how fine?
ORmax<-1
ORgrid<-0.005


#MS analysis
#A SNP: rs61839660 == imm_10_6134703
#D SNP: rs62626317 == imm_10_6150931
#B (Tag) SNP: rs2104286 ==imm_10_6139051

CV1<-"imm_10_6134703"
CV2<-"imm_10_6150931"
Tag<-"imm_10_6139051"

#what is the "real" value of gamma?
gamma.true<-c(-0.18649,-0.22831)

#How many cases and how many controls? 
N0<-12747
N1<-4461

	#the set of gamma which we will find the top p-value at
	gamma_list<-seq(-ORmax,ORmax,by=ORgrid)
	n_gamma<-length(gamma_list)
	top_SNP<-matrix(-1,n_gamma,n_gamma)
	#from the frequencies, compute the object giving the probability of seeing each {X,W} genotype vector
	GenoProbList<-make_GenoProbList(snps,c(CV1,CV2),freq)
	#how many causal variants are there?
	m<-2
	for (ii in 1:n_gamma){
			for (jj in 1:n_gamma){
				try_gamma<-c(gamma_list[ii],gamma_list[jj])
				est_z_score<-abs(est_statistic(N0,N1,snps,c(CV1,CV2),try_gamma,freq,GenoProbList))
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

	#find the subset of points at which we want labels in our plot 
	subset<-which(abs(round(gamma_list,1)-gamma_list)<0.001)
	col.l <- colorRampPalette(c('white', 'firebrick1', 'firebrick4', 'gold'))

	#make the plot
	levelplot(top_SNP,scales=list(x=list(at=subset, labels=gamma_list[subset]),y=list(at=subset, labels=gamma_list[subset])),xlab=list("log OR for A",cex=1.5),ylab=list("log OR for D",cex=1.5),col.regions=col.l,colorkey=FALSE, panel = function(...){
		panel.levelplot(...)
		panel.abline(h = which(gamma_list==0))
            	panel.abline(v = which(gamma_list==0))
		panel.xyplot((gamma.true[1]+ORmax+ORgrid)/ORgrid,(gamma.true[2]+ORmax+ORgrid)/ORgrid, cex = 1, col = 1,pch=19)
        })




