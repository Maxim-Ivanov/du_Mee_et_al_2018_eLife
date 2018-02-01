# This Bash pipeline was used to reanalyze the raw yeast 4sU-Seq
# and Nrd1 PAR-CLIP data from Schulz et al., 2014 (accession
# E-MTAB-1766);

### PART I ###
# Remapping the raw 4sU-Seq and Nrd1 PAR-CLIP reads

# Quality and adapter trimming. Observe that 4sU-Seq samples from
# this study have Illumina Universal adapters, whereas the
# PAR-CLIP samples have Illumina Small RNA adapters:
for file in *fastq.gz; do echo $file && trim_galore --length 15 --no_report_file $file; done

# Align trimmed reads to sacCer3 using STAR:
for file in *trimmed_fq.gz; do echo $file && STAR --genomeDir /index/saccer3/star --readFilesIn $file --runThreadN 4 --outFileNamePrefix ${file/fq.gz/} --outSAMmultNmax 1 --alignEndsType Local --readFilesCommand zcat; done

rm *out *tab; rmdir *STARtmp

# Convert SAM to sorted BAM:
for file in *sam; do echo $file && samtools view -hu $file | samtools sort - -o ${file/.sam/_sorted.bam}; done

# PAR-CLIP libraries are overloaded with rRNA reads. Filter out
# reads which were aligned to tRNA and rRNA genes (UCSC):
for file in *sorted.bam; do echo $file && bedtools intersect -v -abam $file -b /ann/SacCer3/tRNA_rRNA_sacCer3.bed > ${file/.bam/_clean.bam}; done

# Remove reads with MAPQ values below 10:
for file in *clean.bam; do echo $file && samtools view -hbq 10 $file > ${file/.bam/_mapq10.bam}; done

# Generate 4sU-Seq BedGraph files (scaled to 1M reads);
# Observe that strand info was not switched (unlike NET-Seq data
# processed in the "02-Remapping_NET-Seq_raw_data.sh" pipeline);
# Also observe that the 4sU-Seq Bedgraph files represent the
# coverage of the whole 4sU-Seq reads (unlike NET-Seq where only
# 5' terminal bases of reads were used):
for file in *4sU*mapq10.bam; do
  sample=${file/_trimmed.Aligned.out_sorted_mapq10.bam/}
  outfile=${sample}_fw_rev_scaled1M.bg.gz
  reads=$(samtools flagstat $file | sed -n '1p' | awk '{print $1}')
  scaling=$( awk -v rds=$reads 'BEGIN{print 1000000 / rds}' )
  echo $sample $reads $scaling
  bedtools genomecov -ibam $file -bg -split -scale $scaling -strand "-" | awk 'BEGIN{OFS="\t"}{print $1, $2, $3, "-"$4}' > tempfile_rev && echo "  Reverse reads processed;"
  bedtools genomecov -ibam $file -bg -split -scale $scaling -strand "+" | cat - tempfile_rev | sort -k 1,1 -k 2,2n | sed "1i track type=bedGraph color=0,100,200 altColor=200,100,0" | gzip > $outfile
  rm tempfile_rev && echo "  Done!"
done

# Also generate unscaled 4sU-Seq Bedgraph file from the wild-type
# sample to be used as input for the
# "06-Quantification_of_Nrd1_occupancy_in_CUTs.R" pipeline.
# Observe that replicates of the wild-type sample were merged
# together:
file=4sU_wt_merged.bam
samtools merge $file 4sU_wt_rep1*mapq10.bam 4sU_wt_rep2*mapq10.bam
bedtools genomecov -ibam $file -bg -split -strand "-" | awk 'BEGIN{OFS="\t"}{print $1, $2, $3, "-"$4}' > tempfile_rev
bedtools genomecov -ibam $file -bg -split -strand "+" | cat - tempfile_rev | sort -k 1,1 -k 2,2n | sed "1i track type=bedGraph color=0,100,200 altColor=200,100,0" | gzip > ${file/.bam/.bg.gz}
rm tempfile_rev

### PART II ###
# Calling the crosslinking sites in the Nrd1 PAR-CLIP data

# Split Nrd1 PAR-CLIP BAM files into forward and reverse reads:
for file in *Nrd1_PAR-CLIP*mapq10.bam; do echo $file && samtools view -hb -F 16 $file > ${file/.bam/_fwReads.bam} && samtools view -hb -f 16 $file > ${file/.bam/_revReads.bam}; done

# Generate mpileup files:
samtools mpileup -C 50 -f /genomes/sacCer3.fa *fwReads.bam -o All_fwReads.mpileup
samtools mpileup -C 50 -f /genomes/sacCer3.fa *revReads.bam -o All_revReads.mpileup

# Call all SNPs from mpileup files:
java -Xmx8G -jar VarScan.v2.4.3.jar mpileup2snp All_fwReads.mpileup --min-coverage 6 --min-reads2 2 --min-var-freq 0.01 --min-avg-qual 5 --p-value 0.05 --strand-filter 0 --output-vcf 1 > All_fwReads.vcf
java -Xmx8G -jar VarScan.v2.4.3.jar mpileup2snp All_revReads.mpileup --min-coverage 6 --min-reads2 2 --min-var-freq 0.01 --min-avg-qual 5 --p-value 0.05 --strand-filter 0 --output-vcf 1 > All_revReads.vcf

# Extract the conversion events (T-to-C for the forward strand
# reads, A-to-G for the reverse strand reads):
sed '/^#/d' All_fwReads.vcf | awk '{if ($4=="T" && $5=="C") print}' > All_fwReads_TC_conv.vcf
sed '/^#/d' All_revReads.vcf | awk '{if ($4=="A" && $5=="G") print}' > All_revReads_AG_conv.vcf

# Extract the counts of reads supporting the conversion events.
# The resultant Bedgraph file represents the peaks of Nrd1
# occupancy with single base resolution:
python3 Extract_PAR-CLIP_signal_from_VarScan_VCF.py *conv.vcf

# Normalize Bedgraph files to 1M tags and Gzip compress:
for file in PAR-CLIP_signal*bg; do
  norm=$( sed '/^[#t]/d' $file | awk 'BEGIN{SUM=0}{SUM+=sqrt($4^2)*($3-$2)}END{print SUM / 1000000}' )
  echo $file $norm
  awk -v norm=$norm 'BEGIN{OFS="\t"}{if ($0~/^[#t]/) print $0; else print $1, $2, $3, $4 / norm}' $file | gzip > ${file/.bg/_norm1M.bg.gz}
done
