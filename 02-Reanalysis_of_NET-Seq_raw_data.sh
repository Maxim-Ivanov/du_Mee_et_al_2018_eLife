# This Bash pipeline was used to reanalyze the raw yeast NET-Seq
# data from Marquardt et al., 2014 (accession GSE55982);

# Extract Fastq files from SRA archives:
for file in *sra; do echo $file && fastq-dump.2.6.3 $file; done

# Quality and adapter trimming (observe the custom adapter
# sequence):
for file in *fastq.gz; do trim_galore --adapter \
"ATCTCGTATGCCG" --dont_gzip $file; done

# Align trimmed reads to sacCer3 using Bowtie2:
for file in *trimmed.fq; do echo $file && bowtie2 \
--very-sensitive-local -p 4 -x /index/bt2/saccer3 -U $file \
-S ${file/.fq/.sam}; done

# Remove unmapped reads (FLAG == 4):
for file in *sam; do samtools view -h -F 4 $file > \
${file/.sam/_F4.sam}; done

# Convert SAM to sorted BAM:
for file in *F4.sam; do echo $file && samtools view -hu $file | \
samtools sort - -o ${file/.sam/.bam}; done

# Filter out reads which were aligned to tRNA or rRNA genes
# (UCSC):
for file in *bam; do echo $file && bedtools intersect -v \
-abam $file -b /ann/SacCer3/tRNA_rRNA_sacCer3.bed > \
${file/.bam/_clean.bam}; done

# Remove reads with MAPQ values below 10:
for file in *clean.bam; do echo $file && samtools view -hbq 10 \
$file > ${file/.bam/_mapq10.bam}; done

# Make stranded Bedgraph files (observe the strand switch):
for str in "+" "-"; do
  echo $str
  [ "$str" = "+" ] && n="rev" || n="fw"
  for file in *mapq10.bam; do
    echo $file && bedtools genomecov -ibam $file -bg -5 \
-strand $str | sort -k 1,1 -k 2,2n > \
${file/_trimmed_F4_clean_mapq10.bam/}_${n}.bg
  done
done

# Merge forward and reverse Bedgraph files for the same sample:
for file1 in *fw.bg; do
  file2=${file1/fw/rev}
  outfile=${file1/fw.bg/fw_rev.bg.gz}
  echo $file1 "+" $file2 "=" $outfile
  awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"-"$4}' $file2 | \
cat $file1 - | sort -k1,1 -k2,2n | sed '1i track type=bedGraph \
color=0,100,200 altColor=200,100,0' | gzip > $outfile
done

# Normalize Bedgraph files to 1M tags:
for file in *fw_rev.bg.gz; do
  norm=$( zcat $file | sed '/^[#t]/d' | awk 'BEGIN{SUM=0}\
{SUM+=sqrt($4^2)*($3-$2)}END{print SUM / 1000000}' )
  echo $file $norm
  zcat $file | awk -v norm=$norm 'BEGIN{OFS="\t"}\
{if ($0~/^[#t]/) print $0; else print $1, $2, $3, $4 / norm}' | \
gzip > ${file/.bg.gz/_norm1M.bg.gz}
done
