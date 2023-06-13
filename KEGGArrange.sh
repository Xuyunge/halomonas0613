#以分号为分隔符，打印333.txt的第1-19列
#每列为一个文件，文件存储在/seek1目录下，文件名为ko_列名.txt
for val in $(seq 1 19)
do 
    cat 333.txt | awk -v col="$val" -F';' '{print $col}' > seek1/ko_${val}.txt
    echo "$val"
done

#按行拼贴/seek1中的19个文件至all-ko-raw.txt
touch all-ko-raw.txt
for val in $(seq 1 19)
do 
    cat seek1/ko_${val}.txt >> all-ko-raw.txt
    echo "$val"
done
grep -v '^$' all-ko-raw.txt > all-ko-raw2.txt #过滤空白行

cat all-ko-raw2.txt | sort | uniq > all-ko-raw-uniq.txt
cat all-ko-raw-uniq.txt | awk  -F'|' '{print $3}' > teat.txt #以|为分隔的第三列
cat teat.txt | sed 's/.//' > level3.txt #删除每行的第一个字符
rm teat.txt

cat all-ko-raw-uniq.txt | awk  -F'|' '{print $2}' > teat.txt
cat teat.txt | sed 's/.//' > level2.txt
rm teat.txt

cat all-ko-raw-uniq.txt | awk  -F'|' '{print $1}' > teat.txt
cat teat.txt | awk  -F']' '{print $1}' > teat2.txt
cat teat.txt | awk  -F']' '{print $2}' > teat3.txt
cat teat3.txt | sed 's/.//' > level1.txt
cat teat2.txt | sed 's/......//' > ko-pure.txt #删除每行的前六个字符
rm teat*.txt

cut -c 9- < level1.txt | > level1_pure.txt #删除每行的前八个字符
cut -c 9- < level2.txt | > level2_pure.txt
cut -c 9- < level3.txt | > level3_pure.txt

paste ko-pure.txt level1.txt level2.txt level3.txt > ko-level.txt

paste 111.txt seek1/ko_?.txt seek1/ko_1?.txt > gene_with_ko_large.txt

#后面的ko-with-gene-rep.txt文件是通过R得到的
#用了for循环，跑得贼慢，代码文件名是kegg111.R，具体代码不放这里了
#作用就是提取每个ko下包含的基因放到一个表格里

#把每个ko及包含的基因分别提取为一个文件，存储在/seek2目录下
for val in $(seq 1 225)
do 
    cat ko-with-gene-rep.txt | cut -f $val | uniq > seek2/ko_gene_${val}.txt
    echo "$val"
done

paste seek2/ko_gene_?.txt seek2/ko_gene_??.txt seek2/ko_gene_???.txt > ko-with-gene-no-rep.txt

#提取每个ko文件的行数
#行数-1即为ko包含的基因数（第一行是ko名）
touch lines.txt
for val in $(seq 1 225)
do 
    cat seek2/ko_gene_${val}.txt | wc -l >> lines.txt
    echo "$val"
done

#后面的ko-55sum.txt文件是含基因数量最多的前20%的ko的名称，共计55个ko
#用excel直接拉的（能用excel的时候就是要用啊）

#把这55个ko及包含的基因分别提取为一个文件，存储在/seek3目录下
for val in $(seq 1 55)
do 
    a=`sed -n "${val}{p;q;}" ko-55sum.txt`#ko-55sum.txt的第val行
    for val2 in $(seq 1 225)
    do
        b=`head -n 1 seek2/ko_gene_${val2}.txt`
        if [ $a = $b ];then
            touch seek3/ko_gene_${val}.txt
            cat seek2/ko_gene_${val2}.txt > seek3/ko_gene_${val}.txt
        fi
    done
    echo "$val"
done

#按行拼贴/seek3中的55个文件中的基因名至ko_genes_long.txt
touch ko_genes_long.txt
for val in $(seq 1 55)
do 
    cat seek3/ko_gene_${val}.txt | sed '1d' >> ko_genes_long.txt #去除第一行
    echo "$val"
done

#按行拼贴/seek3中的55个文件中的基因名至ko_genes_long.txt
touch ko_process_long.txt
for val in $(seq 1 55)
do 
    cat seek3/ko_gene_${val}.txt | sed '$d' >> ko_process_long.txt #去除末行
    echo "$val"
done

paste ko_process_long.txt ko_genes_long.txt > ko_process2_long.txt
#ko_process2_long.txt之后把第一列的除了ko名以外的部分删除了
#忘了怎么实现的了，总之应该挺好实现的，可能excel就行吧
cut -f 1 ko_process2_long.txt > ko_process3_long.txt

#以下为一个单独的脚本，脚本名为ko_process.sh
#目的是把ko_process3_long.txt空行用ko名填满，以满足富集分析的格式要求
#!/bin/bash
current_ko_term=""
while read line
do
    # 判断是否为空行
    if [ -z "$line" ]
    then
        # 如果是空行，则将当前KO term添加到这一行前面
        echo "$current_ko_term"
    else
        # 如果不是空行，则将当前行的KO term保存，并输出整行
        current_ko_term=$(echo $line | awk '{print $1}')
        echo "$line"
    fi
done < ko_process3_long.txt

zsh ko_process.sh > ko_names_rep.txt
paste ko_names_rep.txt ko_genes_long.txt > ko_fin.txt




