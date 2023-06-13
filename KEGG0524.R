library(clusterProfiler)

# 导入基因列表
gene_all <- read.table("ExpressionMatrixFile",
                   header = T, sep="\t",quote="")

#gene <- subset(gene_all, FDR<0.05 & abs(Log2FC)>1)#筛选标准

gene <- subset(gene, Log2FC > 0)#上调
#gene <- subset(gene, Log2FC < 0)#下调

gene <- as.factor(gene$edgeR_res.gene)
# 导入注释文件
term2gene <- read.table("ko_fin.txt",
                        header=F,sep="\t")
term2name <- read.table("ko_anno.txt",
                        header=F,sep="\t")
term2name <- as.data.frame(cbind(term2name$V1, term2name$V4))

# 富集分析
x <- enricher(gene,TERM2GENE=term2gene,TERM2NAME=term2name,pvalueCutoff = 0.05, 
              pAdjustMethod = "BH",qvalueCutoff = 0.1);head(x)
# 绘制气泡图
dotplot(x,font.size = 12)
# 输出结果
write.table(x,"KEGGResultFile",
            quote = F,
            sep = "\t")




