library(devtools)
library(denovolyzeR)
library(reshape)
library(dplyr)

## Read in data to be analyzed
setwd("/home/Trigmeminal/DeNovo/70\ Trios_New/")
raw_case = read.table("TN_n70_Input.txt",sep="\t",header=TRUE,stringsAsFactors = FALSE)
cases = raw_case
case_num <- 70

## Read in and reformat probability tables for both cases and controls
## code for reformatting the table ##
cases10 = read.delim("0819_hs37d5_coding_idt_med_v2_spikein_padded_Mar2018_adj_modified.txt",sep="\t",stringsAsFactors = F)
colnames(cases10) <- replace(colnames(cases10), colnames(cases10) == "Syn","syn")
colnames(cases10) <- replace(colnames(cases10), colnames(cases10) == "LoF","lof")
colnames(cases10) <- replace(colnames(cases10), colnames(cases10) == "All","all")
colnames(cases10) <- replace(colnames(cases10), colnames(cases10) == "Mis","mis")

reformat_pDNM = function(x, mis_filter="Mis_MetaSVM_or_CADD30"){
  names(x)[names(x)==mis_filter] <- "misD"
  x$misT <- x$mis - x$misD
  x$NegSyn = -(x$syn)
  x$prot <- x[,names(x) %in% c("all","NegSyn")] %>% apply(MARGIN=1, sum)
  x$protD = x[,names(x) %in% c("lof","misD")] %>% apply(MARGIN=1,sum)
  x <- melt(x)                                          #use anything that is 'chr' as id variables
  names(x)[names(x)=="variable"] <- "class"             #change whichever column with the name 'variable' to 'class'
  return(x)
}

reformat_pDNM(cases10) %>% select(class) %>% unique
pDNM_cases10<-reformat_pDNM(cases10)

## ------------- MetaSVM"D" + CADD>=30
burden_cases<-denovolyzeByClass(genes=cases$Cases10GeneName,classes= cases$denovo_call_metasvm_cadd30,geneId="gene",nsamples=case_num, probTable=pDNM_cases10, roundExpected = 3, includeClasses=c("all", "syn", "mis","misD","lof", "protD","prot"))
rate_obs_cases<-burden_cases$observed/case_num
rate_exp_cases<-burden_cases$expected/case_num
burden_cases_final<-cbind(as.character(burden_cases$class),burden_cases$observed,rate_obs_cases,burden_cases$expected, rate_exp_cases,burden_cases$enrichment,burden_cases$pValue)
colnames(burden_cases_final) <- c("class","observed","rate_obs_cases","expected","rate_exp_cases","enrichment","pValue")
burden_cases_final <- as.data.frame(burden_cases_final,stringsAsFactors = FALSE)

## T-Mis
observed <- as.numeric(burden_cases_final[which(burden_cases_final$class == "mis"),]$observed) - as.numeric(burden_cases_final[which(burden_cases_final$class == "misD"),]$observed)
expected <- as.numeric(burden_cases_final[which(burden_cases_final$class == "mis"),]$expected) - as.numeric(burden_cases_final[which(burden_cases_final$class == "misD"),]$expected)
pvalue = ppois(observed-1 , lambda = expected, lower.tail = F)
enrichment <- observed/expected
misT.tmp <- c("misT",observed,observed/case_num,expected,expected/case_num,enrichment,pvalue)
misT.tmp <- as.data.frame(t(misT.tmp))
colnames(misT.tmp) <- c("class","observed","rate_obs_cases","expected","rate_exp_cases","enrichment","pValue")
burden_cases_final <- rbind(burden_cases_final[match(c("all","syn"),burden_cases_final$class),],misT.tmp,burden_cases_final[match(c("misD","lof","protD","prot"),burden_cases_final$class),])

write.table(burden_cases_final,file="70trios_ObserveExpect_metasvm_cadd30.txt",col.names=T,row.names=F,sep="\t",append=F,quote=F)
