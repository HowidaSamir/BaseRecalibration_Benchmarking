
Downloading the ref.
wget ftp://ftp.ensembl.org/pub/release-96/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.21.fa.gz

Building ref. index
bowtie2-build ~/BaseRecalibration_Benchmarking/Homo_sapiens.GRCh38.dna.chromosome.21.fa index_two_bowtie2/Homo_sapiens.fa

R1="$HOME/BaseRecalibration_Benchmarking/SRR8115017.fastq.gz”
RGID=$(cat $R1 | head -n1 | sed 's/:/_/g' |cut -d "." -f1)
PU=$RGID.$LB
LB="SRR8115017_same"
PL="Illumina"
Alignment step:

bowtie2 -p 20 -q --no-unal -x index_two_bowtie2/Homo_sapiens.fa -U SRR8115017.fastq.gz --rg-id $RGID --rg SM:$SM --rg PL:$PL --rg LB:$LB --rg PU:$PU 2> align_stats.txt| samtools view -Sb -o bowtie2.bam

Sorting:
samtools sort bowtie2.bam -o SRR8115017.sorted.bam

Mark-duplicates:
picard_path=$CONDA_PREFIX/share/picard-2.19.0-0

java -Xmx2g -jar $picard_path/picard.jar MarkDuplicates INPUT=SRR8115017.sorted.bam OUTPUT=SRR8115017.dedup.bam METRICS_FILE=SRR8115017.metrics.txt
Indexing

java -Xmx2g -jar $picard_path/picard.jar BuildBamIndex VALIDATION_STRINGENCY=LENIENT INPUT=SRR8115017.dedup.bam

java -Xmx2g -jar $picard_path/picard.jar CreateSequenceDictionary R=Homo_sapiens.GRCh38.dna.chromosome.21.fa O=Homo_sapiens.GRCh38.dna.chromosome.21.dict

samtools faidx Homo_sapiens.GRCh38.dna.chromosome.21.fa

Downloading known variants:
wget ftp://ftp.ensembl.org/pub/release-96/variation/vcf/homo_sapiens/homo_sapiens-chr21.vcf.gz

gunzip homo_sapiens-chr21.vcf.gz 

Indexing:
gatk IndexFeatureFile -F Homo_sapiens_chr21.vcf

head -1291601 Homo_sapiens_chr21.vcf | tail -1   #Then remove this record 

Base recalibration:

gatk --java-options "-Xmx2G" BaseRecalibrator -R Homo_sapiens.GRCh38.dna.chromosome.21.fa -I SRR8115017.dedup.bam --known-sites Homo_sapiens_chr21.vcf -O SRR8115017.report

gatk --java-options "-Xmx2G" ApplyBQSR -R Homo_sapiens.GRCh38.dna.chromosome.21.fa -I SRR8115017.dedup.bam -bqsr SRR8115017.report -O SRR8115017.bqsr.bam --add-output-sam-program-record --emit-original-quals

Without base rec.

gatk --java-options "-Xmx2G" HaplotypeCaller -R Homo_sapiens.GRCh38.dna.chromosome.21.fa -I SRR8115017.dedup.bam --emit-ref-confidence GVCF --pcr-indel-model NONE -O SRR8115017.gvcf

With base rec.

gatk --java-options "-Xmx2G" HaplotypeCaller -R Homo_sapiens.GRCh38.dna.chromosome.21.fa -I SRR8115017.bqsr.bam --emit-ref-confidence GVCF --pcr-indel-model NONE -O SRR8115017.bqsr.gvcf

