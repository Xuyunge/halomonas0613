#----------PART I edgeR方法的差异表达分析----------

## 1.导入并处理RPKM matrix
m <- read.table("CountMatrixFile",
                header=T,
                sep="\t",
                quote="")

m <- m[m[,20]>2,] #过滤低表达量的基因
m1 <- m[,c(1:7)] #每列分别为：基因名1，对照组表达量3，实验组表达量3
Geneid <- m1$Gene.ID
df <- m1[,-1]
row.names(df) <- Geneid
a <- as.matrix(df)

## 2.画个热图先
data <- a
colnames(data) <- c("control-1",
                    "control-2", 
                    "control-3",
                    "case-1",
                    "case-2",
                    "case-3")
# 构建列注释信息（行名与表达矩阵的列名col保持一致）
annotation_col = data.frame(
  group = rep(c("control", "case"), each = 3),
  row.names = colnames(data))
#修改注释标签的颜色
ann_colors = list(
  group = c(control = "#EAEDED", case = "#D5DBDB"))
# 绘制热图
mpf <- pheatmap(data, 
                scale = "row", 
                #cluster_cols = FALSE,
                show_rownames = FALSE,  # 不显示行名称
                treeheight_row = 0,  # 不显示行树的高度
                treeheight_col = 25,
                color = colorRampPalette(c("#5DADE2", "white", "#EC7063"))(100),
                #fontsize_row = 8, 
                fontsize = 9,
                border_color = NA,
                angle_col = 45,
                cellwidth = 26.7,
                cellheight = 0.07,
                cutree_cols = 3,
                #cutree_rows = 4,
                annotation_col = annotation_col,
                annotation_colors = ann_colors,
                #main = "Control vs Case",
)


## 3.edgeR method
group <- factor(c(rep("control", 3), rep("case", 3)),
                levels=c("control","case"),order=F)
metaData <- data.frame(id = colnames(a), group)
#BiocManager::install("edgeR")
library(edgeR)

dgList <- DGEList(counts = a, genes = rownames(a), group=group)
dgList <- calcNormFactors(dgList, method="TMM")

design.mat <- model.matrix(~dgList$sample$group)
colnames(design.mat) <- levels(dgList$sample$group)

d2 <- estimateGLMCommonDisp(dgList, design=design.mat)
d2 <- estimateGLMTrendedDisp(d2, design=design.mat)
d2 <- estimateGLMTagwiseDisp(d2, design=design.mat)

fit <- glmFit(d2, design.mat)
lrt <- glmLRT(fit,coef=2)

edgeR_res <- topTags(lrt, n = nrow(dgList))$table;head(edgeR_res)
edgeR_res <- na.omit(edgeR_res) #去掉数据中有NA的行或列
edgeR_diff <- subset(edgeR_res, FDR<0.05 & abs(logFC)>1) #筛选标准
#省略了之后的添加列的步骤
#本人R语言小白，不太会用R的函数，所以写了for循环，代码很臃肿，不想放了
#添加的列：对照组平均表达量、实验组平均表达量、基因长度、是否essential、简要描述
#列名：mean1, mean2, gene_length, essential, description

write.table(edgeR_diff,
            file="ExpressionMatrixFile",
            row.names=F,
            sep="\t",
            quote=F)#导出表格

## 4.后续过滤
#后续的这步过滤是为了挑基因，后面画火山图用的不是这个文件，是上面的文件
edgeR_diff_filtered <- subset(edgeR_diff, FDR<0.01 & abs(Log2FC)>2
                              & (mean1+mean2)>1000) 
write.table(edgeR_diff_filtered,
            file="FilterExpressionMatrixFile",
            row.names=F,
            sep="\t",
            quote=F)

#----------PART II 差异基因的火山图绘制----------

#BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)
m_exp <- read.table("ExpressionMatrixFile", header=T, sep="\t", quote="")

gene <- m_exp[,1]
Log2FC <- as.numeric(m_exp$Log2FC)
FDR <- as.numeric(m_exp$FDR)

data <- cbind(Log2FC, FDR)

row.names(data) <- gene
data <- as.data.frame(data)

p1 <- EnhancedVolcano(data,
                      lab = rownames(data),
                      x = "Log2FC",
                      y = "FDR",
                      xlim = c(-10, 10),
                      ylim = c(0, 75 + 5),
                      xlab = bquote(~Log[2] ~ "fold change"),
                      ylab = bquote(~-Log[10] ~ italic(FDR)),
                      title = "Case vs Control",
                      subtitle = NULL,
                      pCutoff = 0.05,
                      FCcutoff = 1,
                      pointSize = 1.5,
                      labSize = 0,
                      labFace = "bold",
                      shape = 20,
                      border = "full",
                      borderWidth=1,
                      titleLabSize=15,
                      captionLabSize=11,
                      cutoffLineType="longdash",
                      cutoffLineCol = "black",
                      #colAlpha=1,
                      gridlines.major=F,
                      gridlines.minor=F,
                      legendPosition = "none",
                      #drawConnectors = TRUE,
                      #widthConnectors = 0.5,
                      #typeConnectors = "open",
                      #lengthConnectors = unit(0.01, 'npc'),
                      #min.segment.length = 5,
                      #max.overlaps = "Inf",
                      #arrowheads = F,
                      #colConnectors = "black",
                      #endsConnectors = "last",
                      #directionConnectors = "x",
                      #boxedLabels = TRUE
)

p2<-p1+theme(axis.text.x = element_text(color="black", size=12),
             axis.text.y = element_text(color="black", size=12),
             plot.title = element_text(hjust = 0.5));p2
