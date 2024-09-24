#parse de function
#merges the de file with the annotation and expression files
parse_de = function(de_path, p_threshold, fold_threshold)
{
  de = read.table(de_path, header = TRUE, row.names = 1, sep = "\t")
  de_merged = merge(em_annotated, de, by.x=1, by.y=0)
  row.names(de_merged)= de_merged[,"SYMBOL"]
  de_merged= na.omit(de_merged)
  sorted_de = order(de_merged[,"p.adj"], decreasing = FALSE)
  de_merged = de_merged[sorted_de,]
  de_merged$mean_exp= rowMeans(de_merged[as.vector(row.names(ss))])
  de_merged$mlog10p = -log10(de_merged$p)
  de_merged$mlog10p[is.infinite(de_merged$mlog10p)] = 350
  de_merged$sig= as.factor(de_merged$p.adj < p_threshold & abs(de_merged$log2fold) > fold_threshold)
  return (de_merged)
}

#em_scaled_sig function
get_em_scaled_sig = function(de_merged, em_symbols, em_scaled)
{
  de_merged_sig=subset(de_merged, sig == T)
  sig_genes = row.names(de_merged_sig)
  em_symbols_merged_sig = em_symbols[sig_genes,]
  em_scaled_merged_sig = em_scaled[sig_genes,]
  em_scaled_merged_sig = na.omit(em_scaled_merged_sig)
  return(em_scaled_merged_sig)
}

#volcano plot function
plot_volcano = function(de_merged)
{
  de_sig_up = subset(de_merged, p.adj < 0.001 & log2fold > 2)
  de_sig_down = subset(de_merged, p.adj < 0.001 & log2fold < -2)
  de_sig_up_top5 = de_sig_up[1:5,]
  de_sig_down_top5 = de_sig_down[1:5,]
  ggp = ggplot(de_merged, aes(x = log2fold, y = mlog10p)) +
    geom_point(aes(colour = "a"), size= 0.5) +
    geom_point(data = de_sig_up, aes(colour = "b"), size= 0.5) +
    geom_point(data = de_sig_down, aes(colour = "c"), size= 0.5) +
    labs(title = "Volcano Plot", x = "Log2fold", y = "-Log(p,10)") +
    geom_vline(xintercept = -1, linetype = "dashed", colour= "grey") +
    geom_vline(xintercept = 1, linetype = "dashed", colour="grey") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour= "grey") +
    geom_text_repel(data = de_sig_up_top5, aes(label = SYMBOL), size= 3,fontface = "bold", show.legend = FALSE, segment.color= "transparent") +
    geom_text_repel(data = de_sig_down_top5, aes(label = SYMBOL),size= 3,fontface= "bold", show.legend = FALSE, segment.color= "transparent") +
    scale_colour_manual(values = c("lightsteelblue4", "palevioletred", "dodgerblue3"), labels = c("Non-significant", "Upregulated", "Downregulated"), name= "")+
    theme_bw()+
    theme(legend.text = element_text(face = "bold", size = 10),
          axis.text = element_text(face = "bold", size = 12),
          axis.title = element_text(face = "bold", size = 14),
          plot.title = element_text(face = "bold", size = 16))
  print(ggp)
}

#multi-boxplot function
multi_boxplot = function(de_merged, em_scaled, ss) 
  {
  de_sig15 = de_merged[1:15,]
  gene_data = em_scaled[row.names(de_sig15),]
  gene_data = data.frame(t(gene_data))
  gene_data$SAMPLE_GROUP = ss$SAMPLE_GROUP
  gene_data.m = melt(gene_data, id.vars = "SAMPLE_GROUP")
  
  ggp = ggplot(gene_data.m, aes(x = variable, y = value, fill = SAMPLE_GROUP)) +
    geom_boxplot() +
    labs(title = "Multi-boxplot for the most significant genes", y= "Expression values",x= "Genes")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    scale_fill_manual(values = c("lightsteelblue4", "palevioletred", "dodgerblue3"), 
                      labels = c("Prolif", "Senes", "MTD Senes"))+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    theme(legend.text = element_text(face = "bold", size = 10),
          axis.text = element_text(face = "bold", size = 8),
          axis.title = element_text(face = "bold", size = 12),
          plot.title = element_text(face = "bold", size = 14))
    print(ggp) 
}

#MA plot
plot_ma= function(de_merged)
{
  de_sig_up = subset(de_merged, p.adj < 0.01 & log2fold > 2)
  de_sig_down = subset(de_merged, p.adj < 0.01 & log2fold < -2)
  de_sig_up_top5 = de_sig_up[1:5,]
  de_sig_down_top5 = de_sig_down[1:5,]
  ggp = ggplot(de_merged, aes(x= log10(mean_exp), y = log2fold)) +
    geom_point(aes(colour = "a"), size= 0.5)+
    geom_point(data = de_sig_up, aes(colour= "b"), size= 0.5)+
    geom_point(data = de_sig_down, aes(colour = "c"), size= 0.5)+
    labs(title = "MA Plot", x= "Mean Expression (log10)", y = "Log2fold")+
    geom_hline(yintercept=-1,linetype="dashed") + 
    geom_hline(yintercept=1,linetype="dashed") +
    scale_colour_manual(values = c("lightsteelblue4", "palevioletred", "dodgerblue3"), labels = c("Non-significant", "Upregulated", "Downregulated"), name= "")+
    theme_bw()+
    theme(legend.text = element_text(face = "bold", size = 10),
          axis.text = element_text(face = "bold", size = 12),
          axis.title = element_text(face = "bold", size = 14),
          plot.title = element_text(face = "bold", size = 16))
  print(ggp)
  
}

#Save file function
save_file = function(path, plot, height, width)
{
  emf(path, height = height, width = width)
  print(plot)
  dev.off()
}

#rug function
rug= function(ss)
{
  group_data = as.matrix(as.numeric(as.factor(ss$SAMPLE_GROUP)))
  group_data = melt(group_data)
  
  
  group_colours= c("lightsteelblue4", "palevioletred", "dodgerblue3")
  rug_plot = ggplot(group_data, aes(x = Var1, y = Var2, fill = value)) +
    geom_tile(linetype = "blank") +
    scale_fill_gradientn(colors = group_colours) +
    labs(x = "", y = "") +
    theme(legend.position = "none", legend.title = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank())
  print(rug_plot)
}


#GSEA
gsea_ridge= function(de_merged)
{
  gsea_input = de_merged$log2fold
  names(gsea_input) = row.names(de_merged)
  gsea_input = na.omit(gsea_input)
  gsea_input = sort(gsea_input, decreasing = TRUE)
  
  gse_results = gseGO(geneList= gsea_input, 
                      ont ="BP", 
                      keyType = "SYMBOL", 
                      nPerm = 10000, 
                      minGSSize = 3, 
                      maxGSSize = 800, 
                      pvalueCutoff = 0.05, 
                      verbose = TRUE, 
                      OrgDb = org.Hs.eg.db, 
                      pAdjustMethod = "none")
  ggp = ridgeplot(gse_results) + theme(axis.text.y = element_text(size = 5))+
    labs(title = "Ridgeplot for upregulated genes",
         x= "NES",
         y= "Biological Pathway")
  print(ggp)
}

#ORA
ora_results = function(sig_genes)
{
  sig_genes_entrez = bitr(sig_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  ora_results = enrichGO(gene = sig_genes_entrez$ENTREZID, OrgDb = org.Hs.eg.db, readable = T, ont = 
                           "BP", pvalueCutoff = 0.05, qvalueCutoff = 0.10)
  return (ora_results)
  
}

#Heatmap function
create_heatmap = function(input_matrix) 
{
  hm.matrix = as.matrix(input_matrix[1:100,])
  hm.matrix = na.omit(hm.matrix) 
  y.dist = Dist(hm.matrix, method="spearman")
  y.cluster = hclust(y.dist, method="average")
  y.dd = as.dendrogram(y.cluster)
  y.dd.reorder = reorder(y.dd,0,FUN="average")
  y.order = order.dendrogram(y.dd.reorder)
  hm.matrix_clustered = hm.matrix[y.order,]
  hm.matrix_clustered = melt(hm.matrix_clustered)
  
  colours = c("dodgerblue3", "palevioletred")
  color_palette = colorRampPalette(colours)(100)
  
  #heatmap
  heatmap_plot = ggplot(hm.matrix_clustered, aes(x = Var2, y = Var1, fill = value)) +
    geom_tile() +
    ylab("")+
    xlab("")+
    theme(axis.text.y = element_blank(), axis.ticks=element_blank(), legend.title = element_blank(), legend.spacing.x = unit(0.25, 'cm'))+
    scale_fill_gradientn(colors = color_palette) +
    
    labs(title = "Heatmap of Expression Values",
         x = "Samples",
         y = "Genes")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title= element_text(size= 14, face= "bold"))
  print(heatmap_plot)
}


