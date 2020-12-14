## Author: Xue Zeng (Modified by Sheng Chih Jin)
## Latest modified date: 10-15-2018
## Start from aggregated exome_calls.vcf
## Regenotyping fixing PL
args=commandArgs(TRUE)
if(length(args)!=4){
       cat("Please enter the trio number you want to analyze.","\n")
        cat("Only one trio is allowed per analysis!","\n")
        q(save = "no")
} else{
        ##the ID of the proband of a trio
        Input = args[1]                # exome_calls.vcf
        ped = args[2]                   # /home/xz374/PCGC6.ped
        CohortName = args[3]
        Familyno = args[4]
}

system(paste("mkdir ",CohortName,sep = ''),intern = F)

# Read in the ped file
Fam = read.table(file=ped,header=TRUE,stringsAsFactors=FALSE)

#set your working directory
setwd(CohortName)

##create a new directory for each sample
command0=paste("mkdir -p ",Familyno,sep="")
system(command0,intern=F)

##parse data of each sample in the new directory with the name of the sample
##change the working directory
setwd(Familyno)

##generate the ped file containing pedigree of only one trio
##get sample name of proband, father and mother 
index <- which(Fam$FamID == Familyno)
#command8=paste("grep -w '",Familyno,"' ",ped," > Trio_",Familyno,".ped",sep = "")
#system(command8,intern=F)
write.table(Fam[index,],file = paste("Trio_",Familyno,".ped",sep = ""),col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)
#Fam = unlist(strsplit(try(system(paste("head -1 Trio_",Familyno,".ped",sep = ''),intern = T)),'\t'))
Fam = Fam[index,]
index = which(Fam$Father != 0)
Proband_ID = Fam[index,2]
Father_ID = Fam[index,3]
Mother_ID = Fam[index,4]

##generate VCF file of each sample (trio)
##get rid of all the record with missing genotype
##only extract lines where AC != 0
command1=paste("java -Xmx64g -jar GenomeAnalysisTK_3.5.jar -nt 16 -R /ref_data/h_sapiens/1000genomes/2.5/b37/human_g1k_v37_decoy.fasta -T SelectVariants --variant ",Input," -o Trio_",Familyno,".vcf -env -sn ",Proband_ID,' -sn ',Father_ID,' -sn ',Mother_ID,sep="")
system(command1,intern=F)

##Regenotype
system(paste('java -Xmx64g -jar GenomeAnalysisTK_3.5.jar -nt 16 -R /ref_data/h_sapiens/1000genomes/2.5/b37/human_g1k_v37_decoy.fasta -T RegenotypeVariants --variant Trio_',Familyno,'.vcf -o Trio_',Familyno,'_reGT.vcf',sep = ''),intern = F)

##Split multi-allelic sites
system(paste('bcftools norm -m-both -o Trio_',Familyno,'_reGT_step1.vcf Trio_',Familyno,'_reGT.vcf',sep = ''),intern = F)

##Left normalization
system(paste('bcftools norm -f /ref_data/h_sapiens/1000genomes/2.5/b37/human_g1k_v37.fasta -o Trio_',Familyno,'_reGT_step2.vcf Trio_',Familyno,'_reGT_step1.vcf',sep = ''),intern = F)

##Remove extra information PID,PGT
system(paste('vcfkeepgeno Trio_',Familyno,'_reGT_step2.vcf GT AD DP GQ PL > Trio_',Familyno,'_reGT_step2_modified.vcf',sep = ''),intern = F)

##Remove ./., unfavored PL (should be none), and AC != 0 (should be none)
##remember to check the file name in the script before using!!
system(paste('python ParseTrioVCF.py ',Familyno,sep = ''),intern = F)

##annovar annotate updated VCF file
command7=paste("perl table_annovar.pl --vcfinput Trio_",Familyno,"_updated.vcf /programs/annovar/humandb/ -buildver hg19 -out Trio_",Familyno," -remove -protocol refGene,genomicSuperDups,snp138,dbnsfp33a,esp6500siv2_all,1000g2015aug_all,exac03,gnomad_exome,gnomad_genome -operation g,r,f,f,f,f,f,f,f -nastring '.'",sep = "")
system(command7,intern=F)

##triodenovo v0.04
command6=paste("triodenovo --ped Trio_",Familyno,".ped --in_vcf Trio_",Familyno,"_updated.vcf --out Trio_",Familyno,".denovo.Bayfilter.vcf --mixed_vcf_records",sep = "")

system(command6,intern=F)

##delete all the intermediate files
system(paste("rm Trio_",Familyno,".avinput",sep = ""),intern = F)
system(paste("rm Trio_",Familyno,".hg19_multianno.txt",sep = ""),intern = F)

##Rearrangement
system(paste('python PrepareMerge.py ',Familyno,sep = ''),intern = F)

#get information of the order of members in vcf file
Order = unlist(strsplit(try(system(paste("grep -w '#CHROM' Trio_",Familyno,"_updated.vcf",sep = ""),intern = T)),'\t'))
col15 = Order[10]
col16 = Order[11]
col17 = Order[12]

##process triodenovo processed file
Bayfilter=readLines(paste("Trio_",Familyno,".denovo.Bayfilter.content.txt",sep = ""))
Bayfilter=sapply(1:length(Bayfilter),function(i) unlist(strsplit(Bayfilter[i],"\t")))
Bayfilter=data.frame(t(Bayfilter),stringsAsFactors = F)
#colnames(Bayfilter)=c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT",paste(Familyno,"-02",sep = ""),paste(Familyno,"-01",sep = ""),Familyno)
colnames(Bayfilter)=c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT",Father_ID,Mother_ID,Proband_ID)
Bayfilter$POSITION=paste(Bayfilter$CHROM,Bayfilter$POS,Bayfilter$REF,Bayfilter$ALT,sep=":")

##process ANNOVAR annotated file
Anno=readLines(paste("Trio_",Familyno,".hg19_multianno.content.txt",sep = ""))
Anno=sapply(1:length(Anno),function(i) unlist(strsplit(Anno[i],"\t")))
Anno=data.frame(t(Anno),stringsAsFactors = F)
colnames(Anno)=c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","Info","Format",paste("Anno.",col15,sep = ""),paste("Anno.",col16,sep = ""),paste("Anno.",col17,sep = ""))
Anno$POSITION=paste(Anno$CHROM,Anno$POS,Anno$REF,Anno$ALT,sep=":")

##delete all the intermediate files
system(paste("rm Trio_",Familyno,".hg19_multianno.content.txt",sep = ""),intern = F)
system(paste("rm Trio_",Familyno,".denovo.Bayfilter.content.txt",sep = ""),intern = F)

##merge two files according to the record in triodenovo filtered file
AnnoBayfilter=merge(Bayfilter,Anno,by="POSITION",all.x=TRUE)
AnnoBayfilter=AnnoBayfilter[,c(1:8,10:13,21:25)]
colnames(AnnoBayfilter)=c("POSITION","CHROM","POS","ID","REF","ALT","QUAL","FILTER","FORMAT",Father_ID,Mother_ID,Proband_ID,"INFO","AnnoFormat",paste("Anno.",col15,sep = ""),paste("Anno.",col16,sep = ""),paste("Anno.",col17,sep = ""))

Output=NULL
for(i in 1:nrow(AnnoBayfilter)){
        CHROM=AnnoBayfilter$CHROM[i]
        POS=AnnoBayfilter$POS[i]
        ID=AnnoBayfilter$ID[i]
        REF=AnnoBayfilter$REF[i]
        ALT=AnnoBayfilter$ALT[i]
        QUAL=AnnoBayfilter$QUAL[i]
        FILTER=AnnoBayfilter$FILTER[i]
        FORMAT=AnnoBayfilter$FORMAT[i]
        FATHER=AnnoBayfilter[i,10]
        MOTHER=AnnoBayfilter[i,11]
        PROBAND=AnnoBayfilter[i,12]
        INFO=AnnoBayfilter$INFO[i]
        ANNO.PROBAND=AnnoBayfilter[i,which(colnames(AnnoBayfilter) == paste("Anno.",Proband_ID,sep = ""))] 
        ANNO.MOTHER=AnnoBayfilter[i,which(colnames(AnnoBayfilter) == paste("Anno.",Mother_ID,sep = ""))] 
        ANNO.FATHER=AnnoBayfilter[i,which(colnames(AnnoBayfilter) == paste("Anno.",Father_ID,sep = ""))] 
        
        FATHER=unlist(strsplit(FATHER,":"))
        ANNO.FATHER=unlist(strsplit(ANNO.FATHER,":"))
        GT_FATHER=FATHER[1]
        DQ_FATHER=FATHER[2]
        DGQ_FATHER=FATHER[3]
        DP_FATHER=FATHER[4]
        PL_FATHER=unlist(strsplit(FATHER[5],","))
        PL_HOM1_FATHER=PL_FATHER[1]
        PL_HET_FATHER=PL_FATHER[2]
        PL_HOM2_FATHER=PL_FATHER[3]
        AD_FATHER=unlist(strsplit(ANNO.FATHER[2],","))
        REF_COV_FATHER=AD_FATHER[1]
        NONREF_COV_FATHER=max(as.numeric(AD_FATHER[2:length(AD_FATHER)]))
        GQ_FATHER=ANNO.FATHER[4]
        
        MOTHER=unlist(strsplit(MOTHER,":"))
        ANNO.MOTHER=unlist(strsplit(ANNO.MOTHER,":"))
        GT_MOTHER=MOTHER[1]
        DQ_MOTHER=MOTHER[2]
        DGQ_MOTHER=MOTHER[3]
        DP_MOTHER=MOTHER[4]
        PL_MOTHER=unlist(strsplit(MOTHER[5],","))
        PL_HOM1_MOTHER=PL_MOTHER[1]
        PL_HET_MOTHER=PL_MOTHER[2]
        PL_HOM2_MOTHER=PL_MOTHER[3]
        AD_MOTHER=unlist(strsplit(ANNO.MOTHER[2],","))
        REF_COV_MOTHER=AD_MOTHER[1]
        NONREF_COV_MOTHER=max(as.numeric(AD_MOTHER[2:length(AD_MOTHER)]))
        GQ_MOTHER=ANNO.MOTHER[4]
        
        PROBAND=unlist(strsplit(PROBAND,":"))
        ANNO.PROBAND=unlist(strsplit(ANNO.PROBAND,":"))
        GT_PROBAND=PROBAND[1]
        DQ_PROBAND=PROBAND[2]
        DGQ_PROBAND=PROBAND[3]
        DP_PROBAND=PROBAND[4]
        PL_PROBAND=unlist(strsplit(PROBAND[5],","))
        PL_HOM1_PROBAND=PL_PROBAND[1]
        PL_HET_PROBAND=PL_PROBAND[2]
        PL_HOM2_PROBAND=PL_PROBAND[3]
        AD_PROBAND=unlist(strsplit(ANNO.PROBAND[2],","))
        REF_COV_PROBAND=AD_PROBAND[1]
        NONREF_COV_PROBAND=max(as.numeric(AD_PROBAND[2:length(AD_PROBAND)]))
        GQ_PROBAND=ANNO.PROBAND[4]
        
        INFO=strsplit(INFO,";")[[1]]
        
        AC1 = grep("^AC=",INFO,value=TRUE)
        AC2 = gsub("^AC=","",AC1)
        AC3 = paste(AC2,collapse=";")
        
        AF1 = grep("^AF=",INFO,value=TRUE)
        AF2 = gsub("^AF=","",AF1)
        AF3 = paste(AF2,collapse=";")
        
        AN1 = grep("^AN=",INFO,value=TRUE)
        AN2 = gsub("^AN=","",AN1)
        AN3 = paste(AN2,collapse=";")
        
        DP1 = grep("^DP=",INFO,value=TRUE)
        DP2 = gsub("^DP=","",DP1)
        DP3 = paste(DP2,collapse=";")
        
        Func.refGene1 = grep("^Func.refGene=",INFO,value=TRUE)
        Func.refGene2 = gsub("^Func.refGene=","",Func.refGene1)
        Func.refGene3 = paste(Func.refGene2,collapse=";")
        
        Gene.refGene1 = grep("^Gene.refGene=",INFO,value=TRUE)
        Gene.refGene2 = gsub("^Gene.refGene=","",Gene.refGene1)
        Gene.refGene3 = paste(Gene.refGene2,collapse=";")
        
        GeneDetail.refGene1 = grep("^GeneDetail.refGene=",INFO,value=TRUE)
        GeneDetail.refGene2 = gsub("^GeneDetail.refGene=","",GeneDetail.refGene1)
        GeneDetail.refGene3 = paste(GeneDetail.refGene2,collapse=";")
        
        ExonicFunc.refGene1 = grep("^ExonicFunc.refGene=",INFO,value=TRUE)
        ExonicFunc.refGene2 = gsub("^ExonicFunc.refGene=","",ExonicFunc.refGene1)
        ExonicFunc.refGene3 = paste(ExonicFunc.refGene2,collapse=";")
        
        AAChange.refGene1 = grep("^AAChange.refGene=",INFO,value=TRUE)
        AAChange.refGene2 = gsub("^AAChange.refGene=","",AAChange.refGene1)
        AAChange.refGene3 = paste(AAChange.refGene2,collapse=";")
        
        genomicSuperDups1 = grep("^genomicSuperDups=",INFO,value=TRUE)
        genomicSuperDups2 = gsub("^genomicSuperDups=","",genomicSuperDups1)
        genomicSuperDups3 = paste(genomicSuperDups2,collapse=";")
        
        snp1381 = grep("^snp138=",INFO,value=TRUE)
        snp1382 = gsub("^snp138=","",snp1381)
        snp1383 = paste(snp1382,collapse=";")
        
        esp6500siv2_all1 = grep("^esp6500siv2_all=",INFO,value=TRUE)
        esp6500siv2_all2 = gsub("^esp6500siv2_all=","",esp6500siv2_all1)
        esp6500siv2_all3 = paste(esp6500siv2_all2,collapse=";")
        
        Thousg2015aug_all1 = grep("^1000g2015aug_all=",INFO,value=TRUE)
        Thousg2015aug_all2 = gsub("^1000g2015aug_all=","",Thousg2015aug_all1)
        Thousg2015aug_all3 = paste(Thousg2015aug_all2,collapse=";")
        
        SIFT_score1 = grep("^SIFT_score=",INFO,value=TRUE)
        SIFT_score2 = gsub("^SIFT_score=","",SIFT_score1)
        SIFT_score3 = paste(SIFT_score2,collapse=";")
        
        SIFT_rankscore1 = grep("^SIFT_converted_rankscore=",INFO,value=TRUE)
        SIFT_rankscore2 = gsub("^SIFT_converted_rankscore=","",SIFT_rankscore1)
        SIFT_rankscore3 = paste(SIFT_rankscore2,collapse=";")

	SIFT_pred1 = grep("^SIFT_pred=",INFO,value=TRUE)
        SIFT_pred2 = gsub("^SIFT_pred=","",SIFT_pred1)
        SIFT_pred3 = paste(SIFT_pred2,collapse=";")
        
        MetaSVM_score1 = grep("^MetaSVM_score=",INFO,value=TRUE)
        MetaSVM_score2 = gsub("^MetaSVM_score=","",MetaSVM_score1)
        MetaSVM_score3 = paste(MetaSVM_score2,collapse=";")
        
	MetaSVM_rankscore1 = grep("^MetaSVM_rankscore=",INFO,value=TRUE)
        MetaSVM_rankscore2 = gsub("^MetaSVM_rankscore=","",MetaSVM_rankscore1)
        MetaSVM_rankscore3 = paste(MetaSVM_rankscore2,collapse=";")	

        MetaSVM_pred1 = grep("^MetaSVM_pred=",INFO,value=TRUE)
        MetaSVM_pred2 = gsub("^MetaSVM_pred=","",MetaSVM_pred1)
        MetaSVM_pred3 = paste(MetaSVM_pred2,collapse=";")
        
        CADD_raw1 = grep("^CADD_raw=",INFO,value=TRUE)
        CADD_raw2 = gsub("^CADD_raw=","",CADD_raw1)
        CADD_raw3 = paste(CADD_raw2,collapse=";")
        
	CADD_rawrank1 = grep("^CADD_raw_rankscore=",INFO,value=TRUE)
        CADD_rawrank2 = gsub("^CADD_raw_rankscore=","",CADD_rawrank1)
        CADD_rawrank3 = paste(CADD_rawrank2,collapse=";")

        CADD_phred1 = grep("^CADD_phred=",INFO,value=TRUE)
        CADD_phred2 = gsub("^CADD_phred=","",CADD_phred1)
        CADD_phred3 = paste(CADD_phred2,collapse=";")
        
        phyloP100way_vertebrate1 = grep("^phyloP100way_vertebrate=",INFO,value=TRUE)
        phyloP100way_vertebrate2 = gsub("^phyloP100way_vertebrate=","",phyloP100way_vertebrate1)
        phyloP100way_vertebrate3 = paste(phyloP100way_vertebrate2,collapse=";")
       
	phyloP100way_vertebraterank1 = grep("^phyloP100way_vertebrate_rankscore=",INFO,value=TRUE)
        phyloP100way_vertebraterank2 = gsub("^phyloP100way_vertebrate_rankscore=","",phyloP100way_vertebraterank1)
        phyloP100way_vertebraterank3 = paste(phyloP100way_vertebraterank2,collapse=";")	
 
        SiPhy_29way_logOdds1 = grep("^SiPhy_29way_logOdds=",INFO,value=TRUE)
        SiPhy_29way_logOdds2 = gsub("^SiPhy_29way_logOdds=","",SiPhy_29way_logOdds1)
        SiPhy_29way_logOdds3 = paste(SiPhy_29way_logOdds2,collapse=";")
        
	SiPhy_29way_logOddsrank1 = grep("^SiPhy_29way_logOdds_rankscore=",INFO,value=TRUE)
        SiPhy_29way_logOddsrank2 = gsub("^SiPhy_29way_logOdds_rankscore=","",SiPhy_29way_logOddsrank1)
        SiPhy_29way_logOddsrank3 = paste(SiPhy_29way_logOddsrank2,collapse=";")

        ExAC_ALL1 = grep("^ExAC_ALL=",INFO,value=TRUE)
        ExAC_ALL2 = gsub("^ExAC_ALL=","",ExAC_ALL1)
        ExAC_ALL3 = paste(ExAC_ALL2,collapse=";")
        
        ExAC_AFR1 = grep("^ExAC_AFR=",INFO,value=TRUE)
        ExAC_AFR2 = gsub("^ExAC_AFR=","",ExAC_AFR1)
        ExAC_AFR3 = paste(ExAC_AFR2,collapse=";")
        
        ExAC_AMR1 = grep("^ExAC_AMR=",INFO,value=TRUE)
        ExAC_AMR2 = gsub("^ExAC_AMR=","",ExAC_AMR1)
        ExAC_AMR3 = paste(ExAC_AMR2,collapse=";")
        
        ExAC_EAS1 = grep("^ExAC_EAS=",INFO,value=TRUE)
        ExAC_EAS2 = gsub("^ExAC_EAS=","",ExAC_EAS1)
        ExAC_EAS3 = paste(ExAC_EAS2,collapse=";")
        
        ExAC_FIN1 = grep("^ExAC_FIN=",INFO,value=TRUE)
        ExAC_FIN2 = gsub("^ExAC_FIN=","",ExAC_FIN1)
        ExAC_FIN3 = paste(ExAC_FIN2,collapse=";")
        
        ExAC_NFE1 = grep("^ExAC_NFE=",INFO,value=TRUE)
        ExAC_NFE2 = gsub("^ExAC_NFE=","",ExAC_NFE1)
        ExAC_NFE3 = paste(ExAC_NFE2,collapse=";")
        
        ExAC_OTH1 = grep("^ExAC_OTH=",INFO,value=TRUE)
        ExAC_OTH2 = gsub("^ExAC_OTH=","",ExAC_OTH1)
        ExAC_OTH3 = paste(ExAC_OTH2,collapse=";")
        
        ExAC_SAS1 = grep("^ExAC_SAS=",INFO,value=TRUE)
        ExAC_SAS2 = gsub("^ExAC_SAS=","",ExAC_SAS1)
        ExAC_SAS3 = paste(ExAC_SAS2,collapse=";")
        
	gnomAD_exome_ALL1 = grep("^gnomAD_exome_ALL=",INFO,value=TRUE)
	gnomAD_exome_ALL2 = gsub("^gnomAD_exome_ALL=","",gnomAD_exome_ALL1)
	gnomAD_exome_ALL3 = paste(gnomAD_exome_ALL2,collapse=";") 
	
	gnomAD_genome_ALL1 = grep("^gnomAD_genome_ALL=",INFO,value=TRUE)
        gnomAD_genome_ALL2 = gsub("^gnomAD_genome_ALL=","",gnomAD_genome_ALL1)
        gnomAD_genome_ALL3 = paste(gnomAD_genome_ALL2,collapse=";")


        Output=rbind.data.frame(Output,cbind.data.frame(paste(Familyno,"-",CHROM,":",POS,sep = ""),Familyno,CHROM,POS,REF,REF,ALT,QUAL,FILTER,AC3,AF3,AN3,Func.refGene3,Gene.refGene3,ExonicFunc.refGene3,AAChange.refGene3,genomicSuperDups3,snp1383,esp6500siv2_all3,Thousg2015aug_all3,ExAC_ALL3,ExAC_AMR3,ExAC_FIN3,ExAC_NFE3,ExAC_OTH3,ExAC_AFR3,ExAC_EAS3,ExAC_SAS3,gnomAD_exome_ALL3,gnomAD_genome_ALL3,SIFT_score3,SIFT_rankscore3,SIFT_pred3,MetaSVM_score3,MetaSVM_rankscore3,MetaSVM_pred3,CADD_raw3,CADD_rawrank3,CADD_phred3,phyloP100way_vertebrate3,phyloP100way_vertebraterank3,SiPhy_29way_logOdds3,SiPhy_29way_logOddsrank3,GT_PROBAND,REF_COV_PROBAND,NONREF_COV_PROBAND,DP_PROBAND,GQ_PROBAND,DGQ_PROBAND,PL_HOM1_PROBAND,PL_HET_PROBAND,PL_HOM2_PROBAND,GT_FATHER,REF_COV_FATHER,NONREF_COV_FATHER,DP_FATHER,GQ_FATHER,DGQ_FATHER,PL_HOM1_FATHER,PL_HET_FATHER,PL_HOM2_FATHER,GT_MOTHER,REF_COV_MOTHER,NONREF_COV_MOTHER,DP_MOTHER,GQ_MOTHER,DGQ_MOTHER,PL_HOM1_MOTHER,PL_HET_MOTHER,PL_HOM2_MOTHER,DQ_PROBAND,stringsAsFactors=F))
}

filename=paste("Trio_",Familyno,".BayesianFilter.DenovoM",sep = "")
out.colnames=c("Denovo_Candidate","Proband_ID","chr","pos","ref","parents","alt","QS","PASS","AC","AF","AN","Func.refGene","Gene.refGene","ExonicFunc.refGene","AAChange.refGene","genomicSuperDups","snp138","esp6500siv2_all","1000g2015aug_all","ExAC_ALL","ExAC_AMR","ExAC_FIN","ExAC_NFE","ExAC_OTH","ExAC_AFR","ExAC_EAS","ExAC_SAS","gnomAD_exome_ALL","gnomAD_genome_ALL","SIFT_score","SIFT_converted_rankscore","SIFT_pred","MetaSVM_score","MetaSVM_rankscore","MetaSVM_pred","CADD_raw","CADD_raw_rankscore","CADD_phred","phyloP100way_vertebrate","phyloP100way_vertebrate_rankscore","SiPhy_29way_logOdds","SiPhy_29way_logOdds_rankscore","PROBAND_GT","PROBAND_ref_cov","PROBAND_nonref_cov","PROBAND_DP","PROBAND_GQ","PROBAND_DGQ","PROBAND_PL_hom1","PROBAND_PL_het","PROBAND_PL_hom2","FATHER_GT","FATHER_ref_cov","FATHER_nonref_cov","FATHER_DP","FATHER_GQ","FATHER_DGQ","FATHER_PL_hom1","FATHER_PL_het","FATHER_PL_hom2","MOTHER_GT","MOTHER_ref_cov","MOTHER_nonref_cov","MOTHER_DP","MOTHER_GQ","MOTHER_DGQ","MOTHER_PL_hom1","MOTHER_PL_het","MOTHER_PL_hom2","Bayes_Factor")
write.table(Output,file=filename,col.names=out.colnames,row.names=F,quote=F,sep="\t")


