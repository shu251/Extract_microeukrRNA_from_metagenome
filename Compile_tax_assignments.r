
library(reshape2);library(ggplot2);library(UpSetR);library(dplyr)
options(repr.plot.width = 4, repr.plot.height = 3) #set plot size output 

#rm(count_table)
#rm(tax_import)
tax_files<-list.files(pattern=".txt") # This line finds all file ending with .txt
#
for (i in tax_files){
    tax_import<-read.delim(i, header=FALSE, sep="\t", fill=TRUE, na.strings=".", stringsAsFactors=FALSE) # reads in each file # test this out with your specific input files, as it may require different parameters.
    names<-unlist(strsplit(i,"_")) # parses file name
    tax_import$Sample<-names[1] # replace column 1 header with first part of file name
    tax_import$SILVA<-names[2]
    tax_import$Count<-1 # add a count of n=1 for each rRNA sequence
    if (!exists("count_table")){
        count_table<-tax_import # create the final count table
    } else {
        count_table<-rbind(count_table, tax_import)
    }
    rm(tax_import) # remove excess df
}
head(count_table)

unique(count_table$Sample) # check sample IDs
unique(count_table$SILVA) # check if count came from SILVA database LSU or SSU
save(count_table, file="rawcount_table.RData")

#Add count column and start getting information on taxonomy counts
load("rawcount_table.RData")
df2<-aggregate(count_table$Count, by=list(Taxonomy=df1$V2, Sample=df1$Sample, SILVA=df1$SILVA),sum)
#head(df2)

df2$Tax2 <- gsub("D_.__", "", df2$Taxonomy)
head(df2)

x<-colsplit(df2$Tax2, ";", c("Level1","Level2","Level3","Level4","Level5","Level6", "Level7","Level8","Level9", "Level10", "Level11", "Level12"))

# Sum counts to different taxonomy levels
x2<-data.frame(df2,x)
Domain<-aggregate(x2$x, by=list(domain=x2$Level1, Sample=x2$Sample),sum)
Lev2<-aggregate(x2$x, by=list(domain=x2$Level1, Lev2=x2$Level2, Sample=x2$Sample),sum)
Lev3<-aggregate(x2$x, by=list(domain=x2$Level1, Lev2=x2$Level2, Lev3=x2$Level3, Sample=x2$Sample),sum)

head(Domain);head(Lev2); head(Lev3)

## Option to alter sample names (i.e. re-name or compile replicate samples)
# Domain$Loihi<-Domain$Sample
# pele<-c("PP1", "PP2", "PP5")
# Domain$Loihi[Domain$Sample %in% pele]= "PP"
# Domain_loihi<-aggregate(Domain$x, by=list(domain=Domain$domain, Loihi=Domain$Loihi),sum)
# head(Domain_loihi)

# use Domain df for each sample separate
# or Domain_loihi df for CTD versus Pele Pit
domain<-ggplot(Domain, aes(x=Sample, y=x, fill=domain))+
    geom_bar(stat="identity", position="stack",color="white")+
    theme_minimal()+
    labs(x="", title="", y="Total rRNA reads")+ # add title
    scale_fill_manual(values=c("#66c2a5","#fc8d62","#8da0cb","#e78ac3"))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(color="black"))+
    scale_y_continuous(expand = c(0, 0))
#
domain
#svg("domain.svg", w=4.5, h=3.5); domain ; dev.off() # save figure

summaryDomain <- Domain %>%
    group_by(Sample) %>%
    mutate(RelAbun=100*(x/sum(x))) %>%
    as.data.frame
#summaryDomain
summary<-dcast(summaryDomain[c(1:2,4)], Sample~domain)
#head(summary)
write.csv(summary, file="tmp.csv") # get count table at domain level

lev2_barplot<-ggplot(Lev2, aes(x=Sample, y=x, fill=Lev2))+
    geom_bar(stat="identity", position="stack")+
    theme_minimal()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(color="black"))+
    scale_y_continuous(expand = c(0, 0))+
    labs(x="", title="", y="Total rRNA reads")+
#    scale_fill_manual(values=c("#66c2a5","#fc8d62","#8da0cb","#e78ac3"))+ # option to add colors
    NULL
#lev2_barplot
lev3_barplot<-ggplot(Lev3, aes(x=Sample, y=x, fill=Lev3))+
    geom_bar(stat="identity", position="stack")+
    theme_minimal()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(color="black"))+
    scale_y_continuous(expand = c(0, 0))+
    labs(x="", title="", y="Total rRNA reads")+
#    scale_fill_manual(values=c("#66c2a5","#fc8d62","#8da0cb","#e78ac3"))+ # option to add colors
    NULL

lev2_barplot %+% subset(Lev2, domain %in% "Eukaryota")

lev3_barplot %+% subset(Lev3, Lev2 %in% "SAR")

head(x2)
write.csv(unique(x2[6:12]),file="tmp.csv")

# Manually curate Level 2 and 3:
rename<-function(x){
    x$Euk_curation<-x$Level2 # Most names at Level 2 are suitable
    x$Euk_curation[x$Level1 != "Eukaryota"]="Non-euk"
    x$Euk_curation[x$Level2 == "Incertae Sedis"]="Unannotated"
    x$Euk_curation[x$Level2 == ""]="Unannotated"
    x$Euk_curation[x$Level3 == ""]="Unannotated"
    x$Euk_curation[x$Level2 == "Centrohelida"]="Other"
#Stramenopile
    x$Euk_curation[x$Level3 == "Stramenopiles"]="Stramenopiles-Other"
    x$Euk_curation[x$Level4 == "Diatomea"]="Stramenopiles-Diatomea"
    x$Euk_curation[x$Level4 ==" Dictyochophyceae"]="Stramenopiles-Dictyochophyceae"
#Alveolata
    x$Euk_curation[x$Level3 == "Alveolata"]="Alveolata-Other"
    x$Euk_curation[x$Level4 == "Apicomplexa"]="Alveolata-Apicomplexa"
    x$Euk_curation[x$Level4 == "Ciliophora"]="Alveolata-Ciliophora"
    x$Euk_curation[x$Level4 == "Dinoflagellata"]="Alveolata-Dinoflagellata"
    x$Euk_curation[x$Level5 ==" Syndiniales"]="Alveolata-Syndiniales"
#Rhizaria
    x$Euk_curation[x$Level3 == "Rhizaria"]="Rhizaria-Other"
    x$Euk_curation[x$Level4 == "Cercozoa"]="Rhizaria-Cercozoa"
    x$Euk_curation[x$Level5 ==" Polycystinea"]="Rhizaria-Polycystinea"
    return(x)
}
#Other (defined as those groups that are only noted as "uncultured euk" and were in VERY low abundance here: centroheilda

tax<-rename(x2)
euk_curated_tax<-data.frame(x2,tax)
euk_tax<-aggregate(euk_curated_tax$x, by=list(Euk_tax=euk_curated_tax$Euk_curation, Sample=euk_curated_tax$Sample),sum)
unique(euk_tax$Euk_tax)

# factor:
tax<-c("Stramenopiles-Diatomea","Stramenopiles-Other","Alveolata-Apicomplexa","Alveolata-Ciliophora","Alveolata-Dinoflagellata","Alveolata-Other","Rhizaria-Cercozoa","Rhizaria-Other","Amoebozoa","Archaeplastida","Cryptophyceae","Excavata","Haptophyta","Picozoa","Opisthokonta","Other","Non-euk","Unannotated")
color<-c("#800026","#e31a1c","#fc4e2a","#fd8d3c","#fed976","#ffffcc","#d9f0a3","#78c679","#238443","#004529","#016c59","#3690c0","#08519c","#8c6bb1","#810f7c","#bdbdbd","#737373","#252525")
euk_tax$tax_order<-factor(euk_tax$Euk_tax, levels=tax)
names(color)<-tax

euk<-ggplot(euk_tax, aes(x=Sample, y=x, fill=tax_order))+
    geom_bar(stat="identity", position="stack",color="white")+
    theme_minimal()+
    labs(x="", title="", y="Total rRNA reads")+
    scale_fill_manual(values=color)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(color="black"))+
    scale_y_continuous(expand = c(0, 0))
#
rm<-c("Non-euk", "Unannotated")
euk %+% subset(euk_tax, !(Euk_tax %in% rm)) # try graphing with out the unannotated tax IDs

svg("euk.svg", w=4.5, h=4.5);euk %+% subset(euk_tax, !(Euk_tax %in% rm)); dev.off() # save plot


