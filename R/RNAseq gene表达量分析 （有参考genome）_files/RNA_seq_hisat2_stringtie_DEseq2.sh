#trimmomatic

for i in `ls *_R1.fq.gz`
do 
x=${i/1.fq.gz/}
echo java -jar  /home/wenhao/tools/Trimmomatic-0.33/trimmomatic-0.33.jar PE -phred33 \
${x}"1.fq.gz" ${x}"2.fq.gz" ${x}"1_clean.fq.gz" ${x}"1_unpaired_prepeocessing.fq.gz" \
${x}"2_clean.fq.gz" ${x}"2_unpaired_preprocessing.fq.gz" \
ILLUMINACLIP:/home/wenhao/tools/Trimmomatic-0.33/adapters/TruSeq3-PE.fa:2:30:10 HEADCROP:9 CROP:135 AVGQUAL:20 TRAILING:20 MINLEN:50
done > trimmomatic.sh


## 2---hisat2

#hisat2-build -p 15 --ss Danio_rerio.GRCz10.85.ss --exon Danio_rerio.GRCz10.85.exon  Danio_rerio.GRCz10.dna.toplevel.fa Danio_rerio.GRCz10.dna.toplevel

for i in `ls *1_clean.fq.gz`
do
x=${i/1_clean.fq.gz/}
echo hisat2 -p 16 --dta-cufflinks --no-softclip -x /home/wenhao/work/9_overy_temperature_BS-seq/10_lwh_temperature_RNA-seq/Danio_rerio.GRCz10.dna.toplevel -1 $i -2 ${x}"2_clean.fq.gz" -S ${x}".sam" 2\> ${x}".log"
done > total.hisat2.sh


#############################################################################################################################################################

## 3--- sam2bam and sort

for i in `ls *.sam`
do
x=${i/.sam/}
echo samtools sort -@ 8 -o ${x}".bam" $i
done > samtools.sh

nohup ParaFly -c samtools.sh -CPU 9 &

## 4--- stringtie

for i in `ls *.bam`
do
x=${i/.bam/}
echo stringtie -e -B -p 15 -G /home/ruiqin/data/zebrafish_genome/Danio_rerio.GRCz10.85.LD4.gtf        -o /home/ruiqin/data/transcriptom/CHG030585/$x/$x.gtf  $i
done > stringtie.genome_annotation.sh

##--------DEGseq 差异基因

mkdir ballgown

#cp -r 3_denovo_count ../ballgown
alias python2='/home/wenhao/tools/python27/bin/python'
python2 prepDE.py #与ballgown在同一个目录中，则自动计算count数




