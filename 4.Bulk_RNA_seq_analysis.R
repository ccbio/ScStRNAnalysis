##vocano plot


marker = read.csv("marker_CHolesterol.cell.csv",header=T,stringsAsFactors = F)


library(ggplot2)




chol = c("MSMO1","DHCR7","MVD","INSIG1","IDI1","CYP51A1","HMGCR","HMGCS1","HSD17B7","TM7SF2","FDFT1","FDPS","CES1","MVK",
         "SC5D","LSS","DHCR24","G6PD","SREBF2","ACLY","SQLE","ACAT2")

marker = subset(marker, rownames(marker) %in% c("FADS2","STXBP6","TESMIN","AP002387.2","AL033397.1") == FALSE)
data <- data.frame(
  log2FC = marker$avg_log2FC,
  specificity = marker$pct.1 / marker$pct.2,
  gene = rownames(marker)
)
data$color = ifelse(data$gene %in% chol, "#ff9a01","grey")
data$label = ifelse(data$gene %in% chol, data$gene, "")
# 
p2 = ggplot(data, aes(x = specificity, y = log2FC)) +
  geom_point(aes(size =  log2FC, colour  = I(data$color))) +
  geom_text(aes(label=label), vjust=-1) +
  #scale_fill_manual(values = c("gray", "blue")) +
  #scale_size_area(max_size = 10) +
  theme_minimal() +
  #theme(axis.text.x = element_blank()) +
  coord_flip() +
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "Specificity")
ggplot2::ggsave("vocano_CHolesterol.cells_markers.pdf",p2,height = 6,width=10)



##aurvival
clin$sigscore = calcSigscore(marker, TCGAesca, pos=i, neg=j)
clin$sigscore_group = ifelse(clin$sigscore > median(clin$sigscore),2,1)
select = which(clin$OS.month != "Unknown" & clin$OS.month != "Unk" & clin$OS.month <50)#
clinOnlyCancer0 = clin[select,]
OS_months = as.character(clin$OS.month[select])
OS_status = clin$OS.status[select]
OS_status = ifelse(grepl("Deceased",as.character(OS_status),perl=TRUE),1,0)
SurvClin <- data.frame(cbind(samples = as.character(clinOnlyCancer0$PatientId),OS_status,
                             OS_months = as.numeric(OS_months), group = factor(clinOnlyCancer0$sigscore_group)),
                       stringsAsFactors =FALSE)
SurvClin$OS_status = as.numeric(SurvClin$OS_status)
SurvClin$OS_months = as.numeric(SurvClin$OS_months)
surv_OS = Surv(SurvClin[,3], SurvClin[,2]) #Surv(time,status)
summary(coxph(surv_OS ~ group  , data = SurvClin))

g <- ggsurvplot(
  survfit(surv_OS ~ group, data = SurvClin), 
  conf.int = TRUE,                
  conf.int.style = "ribbon",      
  conf.int.alpha = 0.2,           
  palette = c("#1F77B4FF", "#FF3900"), 
  risk.table = TRUE,              
  risk.table.height = 0.32,       
  size = 0.7,                     
  censor.size = 2.5,              
  pval = TRUE,                   
  break.x.by = 20,                
  surv.median.line = "hv",        
  xlab = "Time (Months)",         
  ylab = "Overall survival probability", 
  legend = "right",               
  legend.title = "Group",         
  legend.labs = c( "Low","High")  
)
g
ggplot2::ggsave(file="sigscore_OS.pdf",arrange_ggsurvplots(list(g), ncol = 1, nrow = 1),height = 4, width = 5)



clin$sigscore = calcSigscore(marker, TCGAesca)
clin$sigscore_group = ifelse(clin$sigscore > median(clin$sigscore),2,1)
select = which(clin$PFS.month != "Unknown" & clin$PFS.month != "Unk" )
clinOnlyCancer0 = clin[select,]
OS_months = as.character(clin$PFS.month[select])
OS_status = clin$PFS.status[select]
OS_status = ifelse(grepl("Progressed",as.character(OS_status),perl=TRUE),1,0)
SurvClin <- data.frame(cbind(samples = as.character(clinOnlyCancer0$PatientId),OS_status,
                             OS_months = as.numeric(OS_months), group = factor(clinOnlyCancer0$sigscore_group)),
                       stringsAsFactors =FALSE)
SurvClin$OS_status = as.numeric(SurvClin$OS_status)
SurvClin$OS_months = as.numeric(SurvClin$OS_months)
surv_OS = Surv(SurvClin[,3], SurvClin[,2]) #Surv(time,status)
summary(coxph(surv_OS ~ group  , data = SurvClin))

g <- ggsurvplot(
  survfit(surv_OS ~ group, data = SurvClin), 
  conf.int = TRUE,                
  conf.int.style = "ribbon",     
  conf.int.alpha = 0.2,           
  palette = c("#1F77B4FF", "#FF3900"),  
  risk.table = TRUE,              
  risk.table.height = 0.32,       
  size = 0.7,                    
  censor.size = 2.5,             
  pval = TRUE,                   
  break.x.by = 20,                
  surv.median.line = "hv",       
  xlab = "Time (Months)",         
  ylab = "Overall survival probability",  
  legend = "right",               
  legend.title = "Group",         
  legend.labs = c( "Low","High")  
)
g
ggplot2::ggsave(file="sigscore_PFS.pdf",arrange_ggsurvplots(list(g), ncol = 1, nrow = 1),height = 4, width = 5)




