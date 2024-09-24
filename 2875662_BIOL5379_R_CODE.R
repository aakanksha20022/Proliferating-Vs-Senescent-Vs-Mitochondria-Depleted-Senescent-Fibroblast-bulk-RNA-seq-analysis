#loading the relevant libraries
library(reshape2)
library(ggplot2)
library(ggrepel)
library(clusterProfiler)
library(org.Hs.eg.db)
library(amap)
library(devEMF)

#Loading the files
ss = read.table("C:/Users/Aakanksha Choudhary/OneDrive/Desktop/data_assess/Dataset and Description-20240112/sample_sheet.csv", row.names= 1, header= TRUE,  sep="\t")
de_mp = read.table("C:/Users/Aakanksha Choudhary/OneDrive/Desktop/data_assess/Dataset and Description-20240112/DE_Senes_MtD_vs_Prolif.csv", header= TRUE, row.names =1, sep="\t")
de_ms = read.table("C:/Users/Aakanksha Choudhary/OneDrive/Desktop/data_assess/Dataset and Description-20240112/DE_Senes_MtD_vs_Senes.csv", header= TRUE, row.names =1, sep="\t")
de_sp = read.table("C:/Users/Aakanksha Choudhary/OneDrive/Desktop/data_assess/Dataset and Description-20240112/DE_Senes_vs_Prolif.csv", header= TRUE, row.names =1, sep="\t")
annot = read.table("C:/Users/Aakanksha Choudhary/OneDrive/Desktop/data_assess/Dataset and Description-20240112/Human_Background_GRCh38.p13.csv", header= TRUE, row.names =1, sep="\t")
em= read.table("C:/Users/Aakanksha Choudhary/OneDrive/Desktop/data_assess/Dataset and Description-20240112/EM.csv", header= TRUE, row.names =1, sep="\t")

#merge em with annotation
em_annotated= merge(em, annot, by.x=0, by.y=0)

#makes gene symbols the rownames
row.names(em_annotated)= em_annotated[,"SYMBOL"]
names(em_annotated)[1]= "GENE_ID"
  
#expression table with symbols
em_symbols= em_annotated[, as.vector(row.names(ss))]

#scaling the expression values
em_scaled = data.frame(scale(t(em_symbols)))
em_scaled= data.frame(t(em_scaled))
em_scaled = na.omit(em_scaled)

#creating a matrix of scaled expressions
em_matrix = as.matrix(sapply(em_scaled, as.numeric))
pca = prcomp(t(em_matrix))
pca_coordinates = data.frame(pca$x)
ss$SAMPLE = row.names(ss)

#calculating percentage variability of each principle component
vars = apply(pca$x, 2, var)
prop_x = round(vars["PC1"] / sum(vars),4) * 100
prop_y = round(vars["PC2"] / sum(vars),4) * 100

x_axis_label = paste("PC1 ", " (",prop_x, "%)",sep="")
y_axis_label = paste("PC2 ", " (",prop_y, "%)",sep="")

#PCA plot
pca_plot= ggplot(pca_coordinates, aes(x=PC1, y= PC2, color = ss$SAMPLE_GROUP)) + 
  geom_point()+
  scale_color_manual(values = c("lightsteelblue4", "palevioletred", "dodgerblue3"))+
  geom_text_repel(aes(label = ss$SAMPLE), segment.color = "transparent")+
  labs(title= "PCA plot", x = x_axis_label, y= y_axis_label, colour = "Groups" )+
  theme_bw()+
  theme(legend.text = element_text(face = "bold", size = 10),
        axis.text = element_text(face = "bold", size = 8),
        axis.title = element_text(face = "bold", size = 12),
        plot.title = element_text(face = "bold", size = 14))
print(pca_plot)
save_file("C:/Users/Aakanksha Choudhary/OneDrive/Desktop/data_assess/pca.emf",pca_plot, 5, 6)

#melting the em table and assigning groups to each sample 
em_symbols.m = melt(em_symbols)

em_symbols.m$group = ifelse(grepl("Prolif", em_symbols.m$variable), "Prolif",
                            ifelse(grepl("Senes_MtD", em_symbols.m$variable), "MTD Senes",
                                   ifelse(grepl("Senes", em_symbols.m$variable), "senes", NA)))


#creating a faceted density plot
density = ggplot(em_symbols.m, aes(x = log10(value), fill = factor(group))) +
  geom_density(size = 0.5, alpha = 0.5) +
  scale_fill_manual(values = c("lightsteelblue4", "palevioletred", "dodgerblue3"))+
  facet_wrap(~variable, ncol = 3, ) +
  theme_bw()+
  labs(title = "Density Plots", x = "Log10(Value)", y = "Density", color = "Groups")+
  theme(legend.text = element_text(face = "bold", size = 10),
        axis.text = element_text(face = "bold", size = 8),
        axis.title = element_text(face = "bold", size = 12),
        plot.title = element_text(face = "bold", size = 14))
  print(density)
  save_file("C:/Users/Aakanksha Choudhary/OneDrive/Desktop/data_assess/density.emf",density, 5, 6)
  
#ANALYSIS OF SENESCENT VS PROLIFERATING CELLS

#using the parse_de function to parse differential senescent v proliferation files
p_threshold = 0.001
fold_threshold = 2

master_sp = parse_de("C:/Users/Aakanksha Choudhary/OneDrive/Desktop/data_assess/Dataset and Description-20240112/DE_Senes_vs_Prolif.csv", p_threshold, fold_threshold)

#extracting the sig genes from master_sp
em_scaled_sp_sig = get_em_scaled_sig(master_sp, em_symbols, em_scaled)

#vocano plot for senescent vs proliferating cells
sp_volcano = plot_volcano(master_sp)
save_file("C:/Users/Aakanksha Choudhary/OneDrive/Desktop/data_assess/sp_volc.emf",sp_volcano, 5, 6)

#multi-boxplot of top 15 significant genes of senescent vs proliferating cells
sp_box= multi_boxplot(master_sp, em_scaled, ss)
save_file("C:/Users/Aakanksha Choudhary/OneDrive/Desktop/data_assess/sp_box.emf",sp_box, 3, 6)

#ma_plot
sp_ma = plot_ma(master_sp)
save_file("C:/Users/Aakanksha Choudhary/OneDrive/Desktop/data_assess/sp_ma.emf",sp_ma, 4, 5)

#heatmap for only senescent vs proliferating cells
em_scaled_sp_sig1 = em_scaled_sp_sig[,c(1:6)]
sp_heatmap = create_heatmap(em_scaled_sp_sig1[1:100, ])
save_file("C:/Users/Aakanksha Choudhary/OneDrive/Desktop/data_assess/sp_heat.emf", sp_heatmap,5,3)

#creating a rug
sp_rug = rug(ss)
save_file("C:/Users/Aakanksha Choudhary/OneDrive/Desktop/data_assess/sp_rug.emf",sp_rug, 5, 6)


#gse on sig genes
master_sp_sig= subset(master_sp, sig == T)
sp_gse= gsea_ridge(master_sp_sig)
save_file("C:/Users/Aakanksha Choudhary/OneDrive/Desktop/data_assess/sp_ridge",sp_gse, 7, 6)


#getting the significant up and down regulated genes of senescent Vs proliferating cells
master_sp_up = subset(master_sp, p.adj < 0.001 & log2fold > 2)
sp_up = row.names(master_sp_up)
master_sp_down = subset(master_sp, p.adj < 0.001 & log2fold < -2)
sp_down = row.names(master_sp_down)

#ORA on upregulated genes
sp_ora_up = ora_results(sp_up)
sp_ora_bar_up = barplot(sp_ora_up, showCategory=10)+
  theme_bw()+
  labs(title= "Over-representation analysis Bar-Plot", y= "Biological Pathway", x= "count")
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(legend.text = element_text(face = "bold", size = 10),
        axis.text = element_text(face = "bold", size = 8),
        axis.title = element_text(face = "bold", size = 12),
        plot.title = element_text(face = "bold", size = 14))
print(sp_ora_bar_up)
save_file("C:/Users/Aakanksha Choudhary/OneDrive/Desktop/data_assess/sp_ora_up.emf",sp_ora_bar_up, 4, 5)

#oRA on downregulated genes
sp_ora_down = ora_results(sp_down)
sp_ora_bar_down = barplot(sp_ora_down, showCategory=10)+
  theme_bw()+
  labs(title= "Over-representation analysis Bar-Plot", y= "Biological Pathway", x= "count")
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(legend.text = element_text(face = "bold", size = 10),
        axis.text = element_text(face = "bold", size = 8),
        axis.title = element_text(face = "bold", size = 12),
        plot.title = element_text(face = "bold", size = 14))
print(sp_ora_bar_down)
save_file("C:/Users/Aakanksha Choudhary/OneDrive/Desktop/data_assess/sp_ora_down.emf",sp_ora_bar_down, 4, 5)



#ANALYSIS OF MIT DEPLETED SENESCENT VS SENESCENT CELLS
#using the parse_de function to parse differential senescent mitochondrial depleted vs senescent files

master_ms = parse_de("C:/Users/Aakanksha Choudhary/OneDrive/Desktop/data_assess/Dataset and Description-20240112/DE_Senes_MtD_vs_Senes.csv", p_threshold, fold_threshold)

#extracting the sig genes from master_ms
em_scaled_ms_sig = get_em_scaled_sig(master_ms, em_symbols, em_scaled)

#vocano plot 
ms_volcano = plot_volcano(master_ms)
save_file("C:/Users/Aakanksha Choudhary/OneDrive/Desktop/data_assess/ms_volc.emf",ms_volcano, 5, 7)

#multi-boxplot for the top 15 significant genes
ms_box= multi_boxplot(master_ms, em_scaled, ss)
save_file("C:/Users/Aakanksha Choudhary/OneDrive/Desktop/data_assess/ms_box.emf",ms_box, 3, 6)

#ma_plot
ms_ma = plot_ma(master_ms)
save_file("C:/Users/Aakanksha Choudhary/OneDrive/Desktop/data_assess/ms_ma.emf",ms_volcano, 4, 5)

#heatmap for the top 100 significant genes
em_scaled_ms_sig1 = em_scaled_ms_sig[,c(4:9)]
ms_heatmap = create_heatmap(em_scaled_ms_sig1[1:100, ])
save_file("C:/Users/Aakanksha Choudhary/OneDrive/Desktop/data_assess/ms_heat.emf",ms_heatmap, 5, 3)

#gse on sig genes
master_ms_sig= subset(master_ms, sig == T)
ms_gse= gsea_ridge(master_ms_sig)
save_file("C:/Users/Aakanksha Choudhary/OneDrive/Desktop/data_assess/ms_gse.emf",ms_gse, 7, 6)


#getting the up and downregulated genes
master_ms_up = subset(master_ms, p.adj < 0.001 & log2fold > 2)
ms_up = row.names(master_ms_up)
master_ms_down = subset(master_ms, p.adj < 0.001 & log2fold < -2)
ms_down = row.names(master_ms_down)

#ora on the upregulated genes
ms_ora_up = ora_results(ms_up)
ms_ora_bar_up = barplot(ms_ora_up, showCategory=10) +
  labs(title= "Over-representation analysis Bar-Plot", y= "Biological Pathway", x= "count")
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(legend.text = element_text(face = "bold", size = 10),
        axis.text = element_text(face = "bold", size = 8),
        axis.title = element_text(face = "bold", size = 12),
        plot.title = element_text(face = "bold", size = 14))
print(ms_ora_bar_up)
save_file("C:/Users/Aakanksha Choudhary/OneDrive/Desktop/data_assess/ms_ora_up.emf",ms_ora_bar_up, 4, 5)


ms_ora_down = ora_results(ms_down)
ms_ora_bar_down = barplot(ms_ora_down, showCategory=10)+
  labs(title= "Over-representation analysis Bar-Plot", y= "Biological Pathway", x= "count")
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(legend.text = element_text(face = "bold", size = 10),
        axis.text = element_text(face = "bold", size = 8),
        axis.title = element_text(face = "bold", size = 12),
        plot.title = element_text(face = "bold", size = 14))
  save_file("C:/Users/Aakanksha Choudhary/OneDrive/Desktop/data_assess/ms_ora_down.emf",ms_ora_bar_down, 4, 5)
  
#ANALYSIS OF MTD DEPLETED SENESCENT VS PROLIFERATING CELLS
master_mp = parse_de("C:/Users/Aakanksha Choudhary/OneDrive/Desktop/data_assess/Dataset and Description-20240112/DE_Senes_MtD_vs_Prolif.csv", p_threshold, fold_threshold)

#extracting the sig genes from master_sp
em_scaled_mp_sig = get_em_scaled_sig(master_mp, em_symbols, em_scaled)


#volcano plot for 
mp_volcano = plot_volcano(master_mp)
save_file("C:/Users/Aakanksha Choudhary/OneDrive/Desktop/data_assess/mp_volc.emf", mp_volcano,5, 7)

#multi-boxplot of sp
mp_box= multi_boxplot(master_mp, em_scaled, ss)
save_file("C:/Users/Aakanksha Choudhary/OneDrive/Desktop/data_assess/mp_box.emf",mp_box,3, 6)

#ma_plot
mp_ma = plot_ma(master_mp)
save_file("C:/Users/Aakanksha Choudhary/OneDrive/Desktop/data_assess/mp_ma.emf",mp_ma, 4,5)


#heatmap 
em_scaled_mp_sig1 = em_scaled_mp_sig[,c(1,2,3,7,8,9)]
mp_heatmap = create_heatmap(em_scaled_mp_sig1[1:100, ])
save_file("C:/Users/Aakanksha Choudhary/OneDrive/Desktop/data_assess/mp_heat.emf", mp_heatmap, 5, 3)

#gse on sig genes
master_mp_sig= subset(master_mp, sig == T)
mp_gse= gsea_ridge(master_mp_sig)
save_file("C:/Users/Aakanksha Choudhary/OneDrive/Desktop/data_assess/mp_ridge.emf",mp_gse, 7, 5)


#ora on upreg and downregulated genes
master_mp_up = subset(master_mp, p.adj < 0.01 & log2fold > 2)
mp_up = row.names(master_mp_up)
master_mp_down = subset(master_mp, p.adj < 0.01 & log2fold < -2)
mp_down = row.names(master_mp_down)

mp_ora_up = ora_results(mp_up)
mp_ora_bar_up = barplot(mp_ora_up, showCategory=10)
print(mp_ora_bar_up)

mp_ora_down = ora_results(mp_down)
mp_ora_bar_down = barplot(mp_ora_up, showCategory=10)
print(mp_ora_bar_down)

#SIGNATURE ANALYSIS FOR SEN V PROL AND MTD SEN V SEN
#creating new de tables with sig columns
de_sp_sig = de_sp
de_ms_sig = de_ms
de_mp_sig = de_mp


de_sp_sig$sig =  factor(de_sp_sig$p.adj < 0.001 & abs(de_sp_sig$log2fold) > 2)
de_ms_sig$sig =  factor(de_ms_sig$p.adj < 0.001 & abs(de_ms_sig$log2fold) > 2)
de_mp_sig$sig =  factor(de_mp_sig$p.adj < 0.001 & abs(de_mp_sig$log2fold) > 2)

#merging all three differential tables containing only sig genes
master_temp= merge(de_sp_sig, de_ms_sig, by.x=0, by.y=0, suffixes = c(".sp", ".ms"))
master_all = merge(master_temp, de_mp_sig, by.x=1, by.y=0)
master_all= merge(master_all, em_annotated, by.x=1, by.y=1)
row.names(master_all)= master_all$SYMBOL
names(master_all)[1]= "GENE_ID"

master_all_sig = subset(master_all, sig.sp == TRUE | sig.ms == TRUE | sig == TRUE)
sig_genes = row.names(master_all_sig)

#getting the overlapping significant genes
ms_sig_count = nrow(subset(master_all, sig.ms == T))
sp_sig_count = nrow(subset(master_all, sig.sp == T))
mp_sig_count = nrow(subset(master_all, sig == T))

sp_ms_sig_count = nrow(subset(master_all, sig.ms == T & sig.sp == T))
sp_mp_sig_count = nrow(subset(master_all, sig.ms == T & sig == T))
ms_mp_sig_count = nrow(subset(master_all, sig.ms == T & sig == T))
all_sig_count = nrow(subset(master_all, sig.ms == T))

#heatmap of significant genes from all three groups
complex_heatmap = create_heatmap(em_scaled[sig_genes,])
save_file("C:/Users/Aakanksha Choudhary/OneDrive/Desktop/data_assess/complex_heat.emf",complex_heatmap, 5, 3)

#Subsetting the master_all table into signatures based on the expression patterns observed in the clustered heatmap
signature_1 = row.names(subset(master_all, (sig.sp== TRUE &log2fold.sp > 2) & (sig.ms == TRUE & log2fold.ms > 2) ))
signature_2 = row.names(subset(master_all, (sig.sp == TRUE & log2fold.sp < -2) & (sig.ms == TRUE & log2fold.ms < -2)))
signature_3 = row.names(subset(master_all, (sig.sp == TRUE & log2fold.sp < -2) & (sig.ms == TRUE & log2fold.ms > 2)))
signature_4 = row.names(subset(master_all, (sig.sp == TRUE & log2fold.sp > 2 ) & (sig.ms == TRUE & log2fold.ms < -2)))

#performing  the correlation on mtd senes vs proliferating and senescent vs proliferating log2fold change
cor_test = cor.test(master_all$log2fold, master_all$log2fold.sp)
print(cor_test)
#performing  the correlation on mtd senes vs proliferating and mtd senes vs proliferating log2fold change
cor_test = cor.test(master_all$log2fold, master_all$log2fold.ms)
print(cor_test)

#fold v fold plot for mtd senes vs proliferating and senescent vs proliferating log2fold change
fold1 = ggplot(master_all, aes(x = log2fold, y = log2fold.sp)) +
  geom_point(data = master_all, aes(colour = ifelse(as.logical(sig.sp) & as.logical(sig), "Both", ifelse(as.logical(sig.sp), "Sig", ifelse(as.logical(sig), "MP", "Other")))), size = 0.7) +
  labs(x = "log2 fold change (MTD Vs Prolif)", y = "log2 fold change (Senes Vs Prolif)", title = "Fold Change vs. Fold Change Plot ") +
  scale_color_manual(values = c("Sig" = "dodgerblue3", "MP" = "palevioletred", "Both" = "orange", "Other" = "grey"),
                    labels = c("SIG_SP","SIG_MP","SIG_BOTH","SIG_NONE"), guides(colour= "")) +
  theme_bw()
  theme(legend.text = element_text(face = "bold", size = 10),
        axis.text = element_text(face = "bold", size = 8),
        axis.title = element_text(face = "bold", size = 12),
        plot.title = element_text(face = "bold", size = 14))

      # Specify colors for significant and non-significant points
save_file("C:/Users/Aakanksha Choudhary/OneDrive/Desktop/data_assess/fold1.emf",fold1, 4, 5)


#fold v fold plot for mtd senes vs proliferating and mtd senes vs proliferating log2fold change
fold2= ggplot(master_all, aes(x = log2fold, y = log2fold.ms)) +
  geom_point(data = master_all, aes(colour = ifelse(as.logical(sig.ms) & as.logical(sig), "Both", ifelse(as.logical(sig.ms), "Sig", ifelse(as.logical(sig), "MP", "Other")))), size = 0.7) +
  labs(x = "log2 fold change (MTD Vs Prolif)", y = "log2 fold change (MTD Senes Vs Senes)", title = "Fold Change vs. Fold Change Plot ") +
  scale_color_manual(values = c("Sig" = "dodgerblue3", "MP" = "palevioletred", "Both" = "orange", "Other" = "grey"),
                     labels = c("SIG_MS","SIG_MP","SIG_BOTH","SIG_NONE"), guides(colour= "")) +
  theme_bw()+
  theme(legend.text = element_text(face = "bold", size = 10),
        axis.text = element_text(face = "bold", size = 8),
        axis.title = element_text(face = "bold", size = 12),
        plot.title = element_text(face = "bold", size = 14))# Specify colors for significant and non-significant points
        theme_bw()
save_file("C:/Users/Aakanksha Choudhary/OneDrive/Desktop/data_assess/fold2.emf",fold2, 4, 5)
        
  print(fold2)     

#performing the ORA analysis on each signature and visualising the results with a cnet plot
sig1_ora = ora_results(signature_1)
sig1_ora_cnet= cnetplot(sig1_ora, showCategory = 10)
print(sig1_ora_cnet)
save_file("C:/Users/Aakanksha Choudhary/OneDrive/Desktop/data_assess/sig1_ora.emf", sig1_ora_cnet, 5, 7)


sig2_ora = ora_results(signature_2)
sig2_ora_cnet= cnetplot(sig2_ora, showCategory = 10)
print(sig2_ora_cnet)
save_file("C:/Users/Aakanksha Choudhary/OneDrive/Desktop/data_assess/sig2_ora.emf", sig2_ora_cnet, 5, 7)


sig3_ora = ora_results(signature_3)
sig3_ora_cnet= cnetplot(sig3_ora, showCategory = 10)
print(sig3_ora_cnet)
save_file("C:/Users/Aakanksha Choudhary/OneDrive/Desktop/data_assess/sig3_ora.emf", sig3_ora_cnet, 5, 7)


sig4_ora = ora_results(signature_4)
sig4_ora_cnet= cnetplot(sig4_ora, showCategory = 10)
print(sig4_ora_cnet)
save_file("C:/Users/Aakanksha Choudhary/OneDrive/Desktop/data_assess/sig4_ora.emf", sig4_ora_cnet, 5, 7)



plot_data = data.frame(
  Group = rep(c("sp", "ms", "mp"), each = 2),
  Regulation = rep(c("Up", "Down"), times = 3),
  Count = c(862,429,466,1106,1158,1021)
)

#barplot to display the count of up and downregulated genes in all groups

bar_reg= ggplot(plot_data, aes(x = Group, y = Count, fill = Regulation)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Up and Down-regulated Genes",
       x = "Groups",
       y = "Count") +
  scale_fill_manual(values = c("Up" = "palevioletred", "Down" = "dodgerblue3")) +
  theme_bw()
  emf("C:/Users/Aakanksha Choudhary/OneDrive/Desktop/data_assess/bar_reg.emf", height = 3, width = 4)
  print(bar_reg)
  dev.off()








