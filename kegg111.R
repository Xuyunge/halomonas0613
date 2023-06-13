
b <- read.table("/Users/xyg/Desktop/KEGG-ttpo/all-ko-raw-uniq.txt",
                sep = "\t", header = F)
# write.table(b, file = "/Users/xyg/Desktop/KEGG-ttpo/teat.txt", quote=F,
#             row.names = F, col.names = F)



m <- read.table("/Users/xyg/Desktop/KEGG-ttpo/gene_with_ko_large.txt",
                sep = "\t", header = T)


col_name <- rep(0, 20)

for (k in c(1:20)) {
  col_name[k] <- paste0("A", k)
}

colnames(m) <- col_name



vec_list <- list()  # 创建一个空列表

for (l in 1:225) {
  vec_name <- paste0("vector_", l)  # 定义向量的名称
  vec_list[[vec_name]] <- character(0)  # 创建一个空的character型向量，并将其存储到列表中
}

vec_list$vector_1  # 访问名称为vector_1的向量


for (k in c(2:20)) {
  col_name <- paste0("A", k)
  for (j in c(1:225)) {
    for (i in c(1:1473)) {
      if (m[[col_name]][i] == b[j,1]){
        vec_list[[j]] <- c(vec_list[[j]], m[[1]][i])
      }
    }
  }
}

# 将list中的向量转换为data.frame
vec_df <- data.frame(do.call(rbind, vec_list))

# 添加行名
c <- read.table("/Users/xyg/Desktop/KEGG-ttpo/ko-pure.txt",
                sep = "\t", header = F)


rownames(vec_df) <- c$V1

t_vec_df <- t(vec_df)

write.table(t_vec_df, file = "/Users/xyg/Desktop/KEGG-ttpo/ko-with-gene-rep.txt"
            , row.names = F,quote = F, sep = "\t")

t_vec_df_sort <- apply(t_vec_df, 2, function(x) sort(x))

write.table(t_vec_df_sort, file = "/Users/xyg/Desktop/KEGG-ttpo/ko-with-gene-rep.txt"
            , row.names = F,quote = F, sep = "\t")

