#Environment setup
#Obtain administrative right of this computer 
#Install ncbi blast plus suite from: https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
#Copy the fasta file of reference sequence into dir where blast+ suite is installed
#Run command prompt as administrator
#To create HXB2 library in blast+ suite, use commend:
#cd "C:\Users\yliu457\" 
#makeblastdb -in HXB2.fasta -parse_seqids -dbtype nucl
#be sure to change the path to the dir where the blast+ suite is installed
#Change below pathways based on your need
RawfastaDir<- "C:/Users/yliu457/Downloads/sequence.fasta"
MyBlastnDir<-"C:/Users/yliu457/"
Output<-"C:/Users/yliu457/Desktop/blastresult.txt"
mywd<-"C:/Users/yliu457"
Primer2ndF_HXB2 <- "GCGCCCGAACAGGGACCTGAAAGCGAAAG"
Primer2ndR_HXB2 <- "TAAGCCTCAATAAAGCTTGCCTTGAGTGC"

#Read fasta file
setwd(mywd)
library(Biostrings)
rawseq<-readDNAStringSet(RawfastaDir)
master<-data.frame(contig=names(rawseq))
n<-nrow(master)
master$pid<-unlist(strsplit(master$contig,"_"))[2+seq(0,n-1)*9]
master$plate<-unlist(strsplit(master$contig,"_"))[3+seq(0,n-1)*9]
master$assay<-unlist(strsplit(master$contig,"_"))[4+seq(0,n-1)*9]
master$well<-unlist(strsplit(master$contig,"_"))[5+seq(0,n-1)*9]
master$len<-unlist(strsplit(master$contig,"_"))[8+seq(0,n-1)*9]
master$rawseq<-as.character(rawseq)

#ncbi blast
system(
  paste(
    "blastn -query ",RawfastaDir," -db "
    ,MyBlastnDir
    ,"HXB2.fasta -num_alignments 1 -reward 1 -penalty -1 -gapopen 2 -gapextend 1 -out ",Output," -outfmt \"6 qseqid qlen sseqid sgi slen qstart qend sstart send evalue bitscore length pident nident btop stitle sstrand\""
    ,sep="")
)
rawdf<-read.delim(Output,header = F)
colnames(rawdf)<-c("SEQID", "qlen", "sseqid", "sgi", "slen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "length", "pident", "nident", "btop", "stitle", "orientation")
master$seqwell<-paste(master$pid,master$plate,master$well,master$assay,sep="_")
master$sample<-paste(master$pid,master$plate,master$well,sep="_")

#QC.01
#fail if no results comes back from blast
master$qc01<-ifelse((master$contig %in% unique(rawdf$SEQID)),"Passed","NonHIV")

#QC. 02
#fail if less than 500bp in length
master$qc02<-ifelse(as.numeric(master$len)>700,"Passed", "Too short")


#filter master table based on the previous result
filtered<-dplyr::filter(master,qc01=="Passed" & qc02=="Passed")

#Orient
rawdf<-dplyr::filter(rawdf,!(sstart<send & send<638) & !(sstart<send & sstart>9632) & !(sstart>send & sstart<638) & !(sstart>send & send>9632))
for (i in unique(filtered$contig)) {
  temp<-subset(rawdf,SEQID==i)
  
  if (length(unique(temp$orientation))!=1) {
    filtered$orient[filtered$contig==i]<-"mix"
  } else {
    filtered$orient[filtered$contig==i]<-unique(temp$orientation)
  }
}

#QC. 03
#fail if orient is mix
filtered$qc03<-ifelse(filtered$orient!="plus" & filtered$orient!="minus","Mix of orientation", "Passed")
for (i in unique(master$contig)) {
  if (i %in% filtered$contig) {
    master$orient[master$contig==i]<-filtered[filtered$contig==i,]$orient
  } else {
    master$orient[master$contig==i]<-"NA"
  }
}
master$qc03<-ifelse(master$orient!="plus" & master$orient!="minus","Mix of orientation or failed in QC", "Passed")

length(unique(filtered[filtered$qc01=="Passed"& filtered$qc02=="Passed" & filtered$qc03=="Passed",]$seqwell))

#Corrected sequences
#filtered<-dplyr::filter(filtered,qc01=="Passed" & qc02=="Passed" & qc03=="Passed")
for (i in unique(filtered$contig)) {
  if (filtered[filtered$contig==i,]$orient!="minus") {
    filtered$seq[filtered$contig==i]<-filtered[filtered$contig==i,]$rawseq
  } else {
    filtered$seq[filtered$contig==i]<-as.character(reverseComplement(DNAString(filtered[filtered$contig==i,]$rawseq)))
  }
}

# QC 04. fail if primers presented in multiple contigs within one well
for (i in unique(filtered$seqwell)) {
  temp<-subset(filtered,seqwell==i)
  if (nrow(temp)==1){
    filtered$qc04[filtered$seqwell==i]<-"Passed"
  } else{
    k<-0
    for (j in c(1:nrow(temp))) {
      if (vcountPattern(Primer2ndF_HXB2,temp$seq[j],max.mismatch = 3)>=1 & vcountPattern(Primer2ndR_HXB2,temp$seq[j],max.mismatch = 3)>=1) {
        k<-k+1
        filtered$Primer[filtered$contig==temp$contig[j]]<-"Yes"
      } else {
        filtered$Primer[filtered$contig==temp$contig[j]]<-"No"
      }
      if (k>=2) {
        break
      }
    }
    if (k==0) {
      filtered$qc04[filtered$seqwell==i]<-"Uncertain"
    }
    if (k==1) {
      filtered$qc04[filtered$seqwell==i]<-"Failed due to no primer"
      filtered$qc04[filtered$seqwell==i & filtered$Primer=="Yes"]<-"Passed"
    }
    if (k>=2) {
      filtered$qc04[filtered$seqwell==i]<-"Multi template"
    }
  }
  
}
#Further check the uncertain seqeunces, include those algined with HXB2
library(tidyverse)
test<-dplyr::filter(filtered,qc04=="Uncertain")
for (i in unique(test$seqwell)) {
  test$alignedlen<-0
  for (j in unique(test[test$seqwell==i,]$contig)) {
    temp<-subset(rawdf,SEQID==j)
    if ((sum(as.numeric(temp$len))/as.numeric(test[test$contig==j,]$len))>0.1) {
      test$alignedlen[test$contig==j]<-sum(as.numeric(temp$len))
    } else  {test$alignedlen[test$contig==j]<-0}
  }
  if(max(test$alignedlen)!=0) {
    filtered[which(filtered$contig==test[test$alignedlen==max(test$alignedlen),]$contig),]$qc04<-"Keep"
  } else {filtered[filtered$seqwell==i,]$qc04<-"Failed"}
  
}
filtered$Primer<-ifelse(is.na(filtered$Primer),"Notchecked",filtered$Primer)
master<-left_join(master,select(filtered,c("contig","Primer","qc04")),by="contig")
filtered<-dplyr::filter(filtered,(qc04=="Passed" & Primer!="No") | qc04=="Keep")

#qc 05
#failed if duplicate in FA/FD and Inner
library(tidyverse)
if (length(filtered[which(duplicated(filtered$sample)),]$sample)!=0) {
  for (i in filtered[which(duplicated(filtered$sample)),]$sample) {
    target<-subset(filtered,sample==i)
    temp<-append(temp,target[which(target$len!=max(target$len))[1],]$contig)
  }
  filtered$qc05<-ifelse(filtered$contig %in% temp,"Shorter Duplication","Passed")
  
} else {
  filtered$qc05<-"No duplication"
}
master<-left_join(master,select(filtered,c("contig","qc05")),by="contig")
filtered<-dplyr::filter(filtered,(qc05=="Passed"| qc05=="No duplication"))

#Export Sequences
Correctseq<-DNAStringSet()
for (i in unique(filtered$contig)) {
  Correctseq[i]<-DNAStringSet(filtered[filtered$contig==i,]$seq)
}
writeXStringSet(Correctseq,"/Users/yliu457/p2008filtered.fasta",append = T)
