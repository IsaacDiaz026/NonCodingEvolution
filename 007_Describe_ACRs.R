#!/usr/bin/env Rscript

setwd("/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/ACRs_classified")
library(dplyr)
library(ggplot2)

library(vcd)

#install.packages("ggstatsplot")
library(ggstatsplot,attach.required=T)

set.seed(123)


acr_freq_gene_dist <- read.table("/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/ACRs_classified/ACRs_class.distance.txt")

acr_freq_gene_dist$peakid <- paste(acr_freq_gene_dist$V1,acr_freq_gene_dist$V2,acr_freq_gene_dist$V3,sep="_")

genic_acrsIntWexons <- read.delim("/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/ACRs_classified/genic_acrsIntWexons.txt", header=FALSE)

genic_acrsIntWexons$acr_length <- genic_acrsIntWexons$V3 - genic_acrsIntWexons$V2
genic_acrsIntWexons$peakid <- paste(genic_acrsIntWexons$V1,genic_acrsIntWexons$V2,genic_acrsIntWexons$V3,sep="_")
genic_acrsIntWexons <- unique(genic_acrsIntWexons)

genic <- unique(acr_freq_gene_dist %>% dplyr::filter(V8 == 0) %>% select(1:4,9))
proximal <- unique(acr_freq_gene_dist %>% dplyr::filter(V8 > 0 & V8 <= 2000) %>% select(1:4,9))
distal <- unique(acr_freq_gene_dist %>% dplyr::filter(V8 > 2000) %>% select(1:4,9))


print("separate intronic from exonic genic acrs")
acr_genic_summary <- genic_acrsIntWexons[,c(1,2,3,8,9,10)]


exon_or_intron <- acr_genic_summary %>% group_by(peakid) %>% summarise(V1=V1,V2=V2,V3=V3,Exon_overlap = sum(V8), ACR_length = mean(acr_length))

exon_or_intron$percent_ACR_inexon <- exon_or_intron$Exon_overlap / exon_or_intron$ACR_length
p <- hist(exon_or_intron$percent_ACR_inexon,breaks=50)

print("noncoding ACRs have less than 10% overlap with exons")
noncoding_genic <- unique(exon_or_intron %>% dplyr::filter(percent_ACR_inexon < 0.1) %>% select(V1,V2,V3,peakid))

noncoding_genic$Type = "noncoding_genic"
noncoding_genic <- noncoding_genic %>% select(V1,V2,V3,Type,peakid)

#add info on frequency class 

noncoding_genic <- merge(noncoding_genic,genic , by="peakid") %>% select(2:4,9,1,5)
colnames(noncoding_genic) = c("V1","V2","V3","V4","peakid","Type")

#filter out noncoding from coding genic ACRs
in_noncoding <- genic$peakid %in% noncoding_genic$peakid

coding_genic <- genic[-which(in_noncoding),]

##add type column
coding_genic$Type = "coding_genic"
genic$Type = "genic"
distal$Type= "distal"
proximal$Type = "proximal"


write.table(noncoding_genic, "/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/ACRs_classified/Noncoding_genic_ACRs.bed", sep = "\t",col.names=F, row.names=F, quote=F)

write.table(coding_genic, "/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/ACRs_classified//bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/ACRs_classified/Species_level/Coding_genic_ACRs.bed", sep = "\t",col.names=F, row.names=F, quote=F)


write.table(genic, "/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/ACRs_classified/Genic_ACRs.bed", sep ="\t", col.names=FALSE, row.names=FALSE,quote=FALSE)
write.table(proximal, "/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/ACRs_classified/Proximal_ACRs.bed",sep="\t",col.names = FALSE, row.names=FALSE,quote =FALSE)
write.table(distal, "/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/ACRs_classified/Distal_ACRs.bed",sep="\t",col.names = FALSE, row.names=FALSE,quote =FALSE)




ACRs_classified <- rbind(coding_genic,noncoding_genic,proximal,distal)

ACRs_classified$Size <- log(ACRs_classified$V3 - ACRs_classified$V2)




pdf("/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/ACRs_classified/ACR_size_byfreq.pdf")
ggplot(ACRs_classified, aes(x=Size, color=V4)) +
 geom_density(fill="white", alpha=0.5, position="identity") + theme_classic()

grouped_ggbetweenstats(ACRs_classified, y=Size, x =V4, type="nonparametric",p.adjust.method = "hochberg",xlab="ACR occurence class", ylab="log(ACR size)",point.args=list(position = ggplot2::position_jitterdodge(dodge.width = 0.6), alpha =0.3, size = 0.5, stroke = 0, na.rm = TRUE), violin.args = list(width=0) , ggtheme = theme_classic(), grouping.var = Type)

ggbetweenstats(ACRs_classified, y=Size, x =Type, type="nonparametric",p.adjust.method = "hochberg",xlab="ACR occurence class", ylab="ACR size", point.args=list(position = ggplot2::position_jitterdodge(dodge.width = 0.6), alpha =0.3, size = 0.5, stroke = 0, na.rm = TRUE), violin.args = list(width=0), ggtheme = theme_classic())
  
dev.off()



to_plot <- ACRs_classified[,c(4,6)]

table_to_plot <- as.data.frame(prop.table(table(to_plot$V4, to_plot$Type), margin = 1))


chi_to_plot <- table(to_plot$V4,to_plot$Type)

b <- chisq.test(to_plot$V4,to_plot$Type)

level_order <- c('High_frequency', 'Common', 'Unique') 

freq_table_to_plot <- as.data.frame(prop.table(table(to_plot$V4)))



pdf("/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/ACRs_classified/2023-10-19_proportion_acr_class.pdf")


ggbarstats(to_plot, x= Type, y= V4, type="nonparametric", ggtheme=theme_classic(), xlab= "ACR occurence class")
ggbarstats(to_plot, x= V4, y= Type, type="nonparametric", ggtheme=theme_classic(),xlab= c("Region"))


mosaic(V4 ~ Type, data=to_plot, gp = shading_max, 
       split_vertical = TRUE)

ggplot(table_to_plot, aes(x=factor(Var1, level = level_order),fill=Var2, y = Freq)) + geom_bar(position = "fill", stat = "identity") + ggtitle("Proportion of ACRs")+theme_classic()+
 theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 0,vjust=.5,face="bold",size = 12)) + scale_fill_brewer(palette = "Pastel2") + xlab("") +ylab("Proportion") + labs(fill="Type") + theme(text = element_text(size=15))


ggplot(freq_table_to_plot, aes(x=factor(Var1, levels = level_order),fill=Var1, y = Freq)) + geom_bar(position = "stack", stat = "identity") + ggtitle("Proportion of ACRs")+theme_classic()+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 0,vjust=.5,face="bold",size = 12)) + scale_fill_brewer(palette = "Pastel2") + xlab("") +ylab("Proportion") + labs(fill="Type") + theme(text = element_text(size=15))

dev.off()

colnames(ACRs_classified)  <- c("ACR_chr","ACR_start","ACR_end","Freq","peakid","Type","Size")
write.table(ACRs_classified, "/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/ACRs_classified/ACRs_classified_freq_type.txt",sep="\t",quote=F, row.names = F, col.names = T)





# Repeat for species level


setwd("/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/ACRs_classified/Species_level/")
library(dplyr)
library(ggplot2)
library(vcd)

acr_freq_gene_dist <- read.table("/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/ACRs_classified/Species_level/ACRs_class.distance.txt")

acr_freq_gene_dist$peakid <- paste(acr_freq_gene_dist$V1,acr_freq_gene_dist$V2,acr_freq_gene_dist$V3,sep="_")

genic_acrsIntWexons <- read.delim("/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/ACRs_classified/Species_level/genic_acrsIntWexons.txt", header=FALSE)

genic_acrsIntWexons$acr_length <- genic_acrsIntWexons$V3 - genic_acrsIntWexons$V2
genic_acrsIntWexons$peakid <- paste(genic_acrsIntWexons$V1,genic_acrsIntWexons$V2,genic_acrsIntWexons$V3,sep="_")
genic_acrsIntWexons <- unique(genic_acrsIntWexons)

genic <- unique(acr_freq_gene_dist %>% dplyr::filter(V8 == 0) %>% select(1:4,9))
proximal <- unique(acr_freq_gene_dist %>% dplyr::filter(V8 > 0 & V8 <= 2000) %>% select(1:4,9))
distal <- unique(acr_freq_gene_dist %>% dplyr::filter(V8 > 2000) %>% select(1:4,9))


print("separate intronic from exonic genic acrs")
acr_genic_summary <- genic_acrsIntWexons[,c(1,2,3,8,9,10)]


exon_or_intron <- acr_genic_summary %>% group_by(peakid) %>% summarise(V1=V1,V2=V2,V3=V3,Exon_overlap = sum(V8), ACR_length = mean(acr_length))

exon_or_intron$percent_ACR_inexon <- exon_or_intron$Exon_overlap / exon_or_intron$ACR_length
p <- hist(exon_or_intron$percent_ACR_inexon,breaks=50)

print("noncoding ACRs have less than 10% overlap with exons")
noncoding_genic <- unique(exon_or_intron %>% dplyr::filter(percent_ACR_inexon < 0.1) %>% select(V1,V2,V3,peakid))

noncoding_genic$Type = "noncoding_genic"
noncoding_genic <- noncoding_genic %>% select(V1,V2,V3,Type,peakid)

#add info on frequency class 

noncoding_genic <- merge(noncoding_genic,genic , by="peakid") %>% select(2:4,9,1,5)
colnames(noncoding_genic) = c("V1","V2","V3","V4","peakid","Type")

#filter out noncoding from coding genic ACRs
in_noncoding <- genic$peakid %in% noncoding_genic$peakid

coding_genic <- genic[-which(in_noncoding),]

##add type column
coding_genic$Type = "coding_genic"
genic$Type = "genic"
distal$Type= "distal"
proximal$Type = "proximal"



write.table(noncoding_genic, "/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/ACRs_classified/Species_level/Noncoding_genic_ACRs.bed", sep = "\t",col.names=F, row.names=F, quote=F)

write.table(coding_genic, "/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/ACRs_classified/Species_level/Coding_genic_ACRs.bed", sep = "\t",col.names=F, row.names=F, quote=F)



write.table(genic, "/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/ACRs_classified/Species_level/Genic_ACRs.bed", sep ="\t", col.names=FALSE, row.names=FALSE,quote=FALSE)
write.table(proximal, "/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/ACRs_classified/Species_level/Proximal_ACRs.bed",sep="\t",col.names = FALSE, row.names=FALSE,quote =FALSE)
write.table(distal, "/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/ACRs_classified/Species_level/Distal_ACRs.bed",sep="\t",col.names = FALSE, row.names=FALSE,quote =FALSE)



ACRs_classified <- rbind(coding_genic,noncoding_genic,proximal,distal)

ACRs_classified$Size <- log(ACRs_classified$V3 - ACRs_classified$V2)



pdf("/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/ACRs_classified/Species_level/ACR_size_byfreq.pdf")
ggplot(ACRs_classified, aes(x=Size, color=V4)) +
  geom_density(fill="white", alpha=0.5, position="identity") + theme_classic()

grouped_ggbetweenstats(ACRs_classified, y=Size, x =V4, type="nonparametric",p.adjust.method = "hochberg",xlab="ACR occurence class", ylab="log(ACR size)",point.args=list(position = ggplot2::position_jitterdodge(dodge.width = 0.6), alpha =0.3, size = 0.5, stroke = 0, na.rm = TRUE), violin.args = list(width=0) , ggtheme = theme_classic(), grouping.var = Type)

ggbetweenstats(ACRs_classified, y=Size, x =Type, type="nonparametric",p.adjust.method = "hochberg",xlab="ACR occurence class", ylab="ACR size", point.args=list(position = ggplot2::position_jitterdodge(dodge.width = 0.6), alpha =0.3, size = 0.5, stroke = 0, na.rm = TRUE), violin.args = list(width=0), ggtheme = theme_classic())

dev.off()



to_plot <- ACRs_classified[,c(4,6)]

table_to_plot <- as.data.frame(prop.table(table(to_plot$V4, to_plot$Type), margin = 1))


chi_to_plot <- table(to_plot$V4,to_plot$Type)

b <- chisq.test(to_plot$V4,to_plot$Type)

level_order <- c('High_frequency', 'Common', 'Unique') 

freq_table_to_plot <- as.data.frame(prop.table(table(to_plot$V4)))



pdf("/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/ACRs_classified/Species_level/2023-10-19_proportion_acr_class.pdf")


ggbarstats(to_plot, x= Type, y= V4, type="nonparametric", ggtheme=theme_classic(), xlab= "ACR occurence class")
ggbarstats(to_plot, x= V4, y= Type, type="nonparametric", ggtheme=theme_classic(),xlab= c("Region"))


mosaic(V4 ~ Type, data=to_plot, gp = shading_max, 
       split_vertical = TRUE)

ggplot(table_to_plot, aes(x=factor(Var1, level = level_order),fill=Var2, y = Freq)) + geom_bar(position = "fill", stat = "identity") + ggtitle("Proportion of ACRs")+theme_classic()+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 0,vjust=.5,face="bold",size = 12)) + scale_fill_brewer(palette = "Pastel2") + xlab("") +ylab("Proportion") + labs(fill="Type") + theme(text = element_text(size=15))


ggplot(freq_table_to_plot, aes(x=factor(Var1, levels = level_order),fill=Var1, y = Freq)) + geom_bar(position = "stack", stat = "identity") + ggtitle("Proportion of ACRs")+theme_classic()+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 0,vjust=.5,face="bold",size = 12)) + scale_fill_brewer(palette = "Pastel2") + xlab("") +ylab("Proportion") + labs(fill="Type") + theme(text = element_text(size=15))

dev.off()

colnames(ACRs_classified)  <- c("ACR_chr","ACR_start","ACR_end","Freq","peakid","Type","Size")

write.table(ACRs_classified, "/bigdata/seymourlab/idiaz026/Citrus_non-coding_evolution/2023-08-13_PWN_socrates/ACRs_classified/Species_level/ACRs_classified_freq_type.txt",sep="\t",quote=F, row.names = F, col.names = T)




