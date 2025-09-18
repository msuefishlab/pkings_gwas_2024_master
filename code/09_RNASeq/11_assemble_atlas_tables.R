library(tximport)

samples<-read.table('sample_table_for_kallisto_mt.csv',sep=",",header=T)

mouse_samples<-subset(samples,species=="Mus_musculus")
human_samples<-subset(samples,species=="Homo_sapiens")
zebrafish_samples<-subset(samples,species=="Danio_rerio")
paramormyrops_samples<-subset(samples,species=="Paramormyrops_kingsleyae")

mouse_tx2gene<-read.csv('mouse_mt_tx2gene.csv')
human_tx2gene<-read.csv('human_mt_tx2gene.csv')
zebrafish_tx2gene<-read.csv('zebrafish_mt_tx2gene.csv')
paramormyrops_tx2gene<-read.csv('paramormyrops_mt_and_somatic_tx2gene.csv')

mouse_tx2gene$X<-NULL
human_tx2gene$X<-NULL
zebrafish_tx2gene$X<-NULL
paramormyrops_tx2gene$X<-NULL

mouse_files<-file.path("kallisto_output_mt",mouse_samples$tissue,mouse_samples$species,mouse_samples$individual,"abundance.h5")
human_files<-file.path("kallisto_output_mt",human_samples$tissue,human_samples$species,human_samples$individual,"abundance.h5")
zebrafish_files<-file.path("kallisto_output_mt",zebrafish_samples$tissue,zebrafish_samples$species,zebrafish_samples$individual,"abundance.h5")
paramormyrops_files<-file.path("kallisto_output_mt",paramormyrops_samples$tissue,paramormyrops_samples$species,paramormyrops_samples$individual,"abundance.h5")

names(zebrafish_files)<-paste0(zebrafish_samples$tissue,"_",zebrafish_samples$individual)
names(human_files)<-paste0(human_samples$tissue,"_",human_samples$individual)
names(mouse_files)<-paste0(mouse_samples$tissue,"_",mouse_samples$individual)
names(paramormyrops_files)<-paste0(paramormyrops_samples$tissue,"_",paramormyrops_samples$individual)

zebrafish_data<-tximport(zebrafish_files,type="kallisto",tx2gene=zebrafish_tx2gene)
human_data<-tximport(human_files,type="kallisto",tx2gene=human_tx2gene)
mouse_data<-tximport(mouse_files,type="kallisto",tx2gene=mouse_tx2gene)
paramormyrops_data<-tximport(paramormyrops_files,type="kallisto",tx2gene=paramormyrops_tx2gene)

zebrafish_cfap_data<-subset(zebrafish_data$abundance, rownames(zebrafish_data$abundance)=="si:ch73-222h13.1")
human_cfap_data<-subset(human_data$abundance, rownames(human_data$abundance)=="CFAP221")
mouse_cfap_data<-subset(mouse_data$abundance, rownames(mouse_data$abundance)=="Cfap221")
paramormyrops_cfap_data<-subset(paramormyrops_data$abundance,rownames(paramormyrops_data$abundance)=="cfap221")

paramormyrops_amhr2a_data<-subset(paramormyrops_data$abundance,rownames(paramormyrops_data$abundance)=="LOC111835276")
paramormyrops_amhr2b_data<-subset(paramormyrops_data$abundance,rownames(paramormyrops_data$abundance)=="LOC111854052")
paramormyrops_amhr2_data<-rbind(paramormyrops_amhr2a_data,paramormyrops_amhr2b_data)

write.csv(zebrafish_cfap_data,"zebrafish_cfap_data.csv")
write.csv(mouse_cfap_data,"mouse_cfap_data.csv")
write.csv(human_cfap_data,"human_cfap_data.csv")
write.csv(paramormyrops_cfap_data,"paramormyrops_cfap_data.csv")

write.csv(paramormyrops_amhr2_data,"paramormyrops_amhr2_data.csv")
