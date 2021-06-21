
## 计算染色体长度
seqtk comp  genome.fasta |grep  "^LG" |awk '{print $1"\t"$2}' > genome.len

## 生成染色体文件 7列
awk '{print "chr\t-\t"$1"\t"$1"\t0\t"$2"\tchr"NR}' genome.len > Sind_karyotype.txt

## 生成窗口文件， 窗口大小50Kb
bedtools  makewindows -w 50000 -g  genome.len > genome.window.bed

## 计算每个窗口平均GC含量
seqtk subseq genome.fasta  genome.window.bed  > genome.window.fasta
seqtk comp  genome.window.fasta |awk '{print $1 "\t" ($4+$5)/($3+$4+$5+$6) } ' |awk -F ":|-" '{print $1"\t"$2"\t"$3"\t"$4}'> Sind_gc.txt

## 计算每个窗口基因条数
bedtools intersect  -a genome.window.bed -b Sind.bed -c -F 0.1  > Sind_genecount.txt

## 计算每个窗口重复序列含量
bedtools coverage -a  genome.window.bed -b  repeat.gff |awk '{print $1 "\t" $2 "\t" $3 "\t" $7}' > Sind_repeat.txt

## 生成共线性link文件
perl get_jcvi_block.pl  Sind.bed  Sind.Sind.anchors 8  >  Sind_links.txt
