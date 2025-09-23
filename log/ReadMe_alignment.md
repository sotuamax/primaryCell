**Table of Content**
- [Data processing steps for the human primary cell project](#data-processing-steps-for-the-human-primary-cell-project)
  - [reference genome](#reference-genome)
  - [Primary celltypes](#primary-celltypes)
  - [Data organization](#data-organization)
  - [Seq alignment](#seq-alignment)
  - [Markdup](#markdup)
  - [Quality control (QC)](#quality-control-qc)
  - [Merge replicates (no longer required - June 05, 2025)](#merge-replicates-no-longer-required---june-05-2025)
  - [Evaluation](#evaluation)
  - [Downstream processed files](#downstream-processed-files)
  - [Peaks](#peaks)
  - [Note](#note)
  - [Important update](#important-update)

### Data processing steps for the human primary cell project 

#### reference genome
- hg38 (downloaded from GENCODE)
- nucleotide seq of GRCh38 primary genome assembly (/data/jim4/Reference/human/GRCh38.p14/fasta/GRCh38.primary_assembly.genome.fa); 


#### Primary celltypes
- five blood celltypes (B, CD4, CD8, Monocyte, NK, all selected from blood sample);
- twenty-four organ celltypes.

#### Data organization 
- Each dataset (RNA, DNase, Hi-TrAC) has separate folder; 
- Under each dataset, alignment file stored in "raw" (at "individual" sample level); 
- Remove duplicates file stored in "markdup" (at "individual" sample level, when necessay);
- QC alignment file stored in "QC" (at "individual" sample level); 
- For sample replicates, merged and stored in "merged" (after QC, at celltype level);

#### Seq alignment 
- "raw" fastq aligned against the reference genome using BWA and sorted using samtools;
- RNA-seq data processed using different programs (STAR, RSEM, featureCount);
- specifically for Hi-TrAC data, "-5SP" used in BWA alignment;
- sorted by position with index generated.
`bwa mem -t 48 $hg38 $r1 $r2 | samtools view -@ 48 -Su - | samtools sort -@ 48 - -o output.bam && samtools index -@ 24 output.bam`

#### Markdup 
- "markdup" picard MarkDuplicates used to remove duplicates by assigning "--REMOVE_DUPLICATES true";
`java -jar $picard MarkDuplicates -I input.bam -O output.bam -M output.metric --REMOVE_DUPLICATES true && samtools index -@ 24 output.bam`;

#### Quality control (QC)
- "QC" see `QC.py` for details about this step: 
    * MAPQ > 10; primary alignment for DNase-seq; 
    * MAPQ > 30; primary alignment for Hi-TrAC. Note: blood cells are not included as the quality is bad. 
    * MAPQ > 30; **("NH", 1)** in read tags, and on chromosome region for RNA-seq (downstream of STAR alignment); 
    * For transcriptome alignment, its BAM file was filtered based on genome alignment QC selection.

#### Merge replicates (no longer required - June 05, 2025)
- samples with multiple replicates, after QC the alignment files are merged (specifically for DNase-seq and Hi-TrAC):
    * see "name_sorted" (temporary): for read sorted by _readname_;
    * see "clean" for paired-read pass QC and sorted by _readname_;
    * see "pos_sorted" for position sorted BAM based on 'clean' BAM;
    * see "split" (temporary) for split sample bam files, used to re-generate merged-bam fils of good quality samples;
    * see "clean_w_good_sample" for new merged bam file with samples pass the criteria. If it is the same as in "clean" folder, a softlink created. 
    * see "pos_sorted_w_good_sample" for new merged bam file with good samples. If it is the same as in the "pos_sorted" folder, a softlink created.
    * see "split_good" for good samples split from the merged BAM file from "pos_sorted_w_good_sample". 
    * Note:  To solve the issue related to one-end of the paired reads passed QC, paired-end controlled BAM files (both reads have to pass QC) are added in "clean" folder.

- For record, RG (readGroup) is added when merging BAM files;
  `samtools merge -r -o merged.bam sampl1.bam sample2.bam` 

#### Evaluation 
- "eval" folder contains metrics for library quality evaluation. 
  * DNase-seq:
    * SPOT score; 
    all sample replicates are pooled to perform peak call. Sample 5 M reads to count the fraction of reads reside on peak enrichment region. 
    `spot_enrich.py -bed peak.bed -bam sample.bam -n 48`
    also output in file *SPOT.txt*
    * TSS enrichment score 
    selected TSS used for evaluate the enrichment strength of TSS in comparing to its surrounding background. 
    `b2bw_profile.py -bam sample.bam -bed TSS.bed -o sample`
    also see *TSS_enrich.txt*
    * About window size: 
    flanking region side of TSS can be 1000, 2000. There is no differences in flanking region size. Use 1000 bp to save process time. 

  * Hi-TrAC: 
    * SPOT score: 
    update 02/11/2025: "spot_enrich.py" updated considering read-pair reside within peak enrichment loci (by counting read query_name occurence twice). It is no need to rerun for DNase-seq, as reads are properly paired, which means when it occur, its mate is likely also within peak region. 
    `spot_enrich.py -bed peak.bed -bam sample.bam -n 48`
    After updating "spot_enrich.py", the enrichment score is noticeble decreased. 
    Based on experience, cutoff at 0.1 should be enough. 
    * TSS enrichment score:
    same as for DNase-seq. 
    `b2bw_profile.py -bam sample.bam -bed TSS.bed -o sample`
    see output in *TSS_enrich.txt"
    * cis distribution
    cis (intra-chromosomal PETs) ratio (against total)
    distant ratio (against cis)
    `bam_cis.py sample.bam -t 48`
    see output in *Cis_ratio.txt*
    It may be confusing use either cis ratio only, or distant cis ratio only, and to combine the two metrics, I used their multiplied value. At least 20% as distant cis versus total is set.

  * RNA-seq (non-stranded specific): 
    * To quantify gene expression, I used two strategies: 
      * read-count overlaying gene body: 
       `featureCounts -p --countReadPairs -a $hg38_gtf --extraAttributes gene_name -T 36 -s 0 -o featureCount_QC.txt QC/*bam`
        > Note that I filtered RNA-seq alignment (on 03/31/2025) and rerun the read-count, and featured-reads ratio is significantly improved (~85%). It indicates the recent QC on alignment file can help get cleaner gene expression;
        
        > See `count_norm.R` for normalization using DESeq. 
      * normalized expression by RSEM: 
        for rsem_build, see: "/data/jim4/Reference/human/GRCh38.p14/rsem_build" 
      `rsem-calculate-expression --star-gzipped-read-file --star --strandedness none -p 32 --paired-end r1.fastq.gz r2.fastq.gz /data/jim4/Reference/human/GRCh38.p14/rsem_build/hg38_rsem sample`
         > Note that RSEM is based on alignment file on transcriptome, so keep its BAM file for future debug. 

         > Two runs: run rsem first for transcriptome-alignment;
          then QC BAM file to get a cleaner version of alignment;
          rerun rsem with new alignment to estimate expression level. 

#### Downstream processed files 

- process pipeline for HiTrAC: 
  - `hic_align.py pre`; 
  - `hic_align.py align`;
  - `hic_align.py markdup`;
  - `hic_align.py qc`;
  - `hic_format.py <bam> -o <out> -n 12`;

- bed: 
  * stores "bedpe" file generated using bedtools from **clean BAM** 
  `bedtools bamtobed -i input.bam -bedpe > output.bedpe`
  * stores "bed" file generated using bedtools from clean BAM 
  `bedtools bamtobed -i input.bam > output.bed`
  * see "bed_w_good_sample" folder for bed file with bam file using good samples. When it is the same as in bed file, a soft-link generated.

- bigwig: 
  * generated using deeptools on *clean bam* for merged alignment (passing QC alignment filtering)
    `bamCoverage -b input.bam --blackListFileName hg38_blacklist.bed -o output.bw -of bigwig -p 24`
    * bam file must be position-sorted and indexed.
    * --normalizeUsing RPKM (default: None) and default is very similar. To simplify the method, using *default setting*. 
    * For downstream visualization purposes, RPKM normalization performed. 
    * BAM input: '/data/jim4/Seq/primary_cell_project/alignment/DNase/merged/pos_sorted_w_good_sample/*bam' most updated on 02/27/2025 
    `bamCoverage -b input.bam --blackListFileName hg38_blacklist.bed -o output.bw -of bigwig -p 24 --normalizeUsing RPKM`
    * RPKM (per bin) = number of reads per bin / (number of mapped reads (in millions) * bin length (kb))

- pairs 
  * this format is a intermediate format prepared for cooler and juicer to generate cool and hic file 
  Use [pairtools parse](https://pairtools.readthedocs.io/en/latest/examples/pairtools_walkthrough.html#pairtools-parse): 
  `pairtools parse input.bam -c $hg38_chrsize -o output.pairs --assembly hg38 --nproc-in 4 --nproc-out 4 --drop-sam`
    * columns: 
    readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type
    Pairs have to be sorted for HiC generation (output in sorted and compressed pairs), see: 
    `pairtools sort input.pairs -o output.pairs.gz --nproc 24 --memory 40g`

- mcool
  * specific for Hi-TrAC data 
  Use [cooler cload paris](https://cooler.readthedocs.io/en/latest/cli.html#cooler-cload-pairs): at 1 kb resolution cool file
  `cooler cload pairs $hg38_chrsize:1000 sample.pairs sample.cool --assembly hg38 -c1 2 -p1 3 -c2 4 -p2 5`
  explanation: 
  * columns info for the input pairs file need to be assigned. 
  -c1: chrom1 column num (1-based)
  -p1: pos1 column num
  -c2: chrom2 column num 
  -p2: pos2 column num 
  * at multiple resolution mcool file
  `cooler zoomify input.cool -r 1000,2000,5000,10000,25000,50000,100000 -o output.mcool`

- HiC
  * specific for Hi-TrAC data
  Use "`juicer_tools pre`" 
  `juicer_tools pre -r 1000,2000,5000,10000,25000,50000,100000 -k KR,ICE,VC input.pairs.gz sample.hic --threads 24 hg38`

#### Peaks
- specific for DNase-seq and Hi-TrAC data; 
- refer to [here](https://github.com/macs3-project/MACS/discussions/435) about its setting; 
- Generated using MACS3 with "--nomodel" setting 
  * For DNase-seq: 
  `macs3 callpeak -q 0.01 --keep-dup all -g hs -t input.bed -f BED --nomodel --shift -100 --extsize 200 --outdir . -n output_name`
    * see "raw" for peaks called using all samples combined BED; 
    * see "peak_w_good_sample" for peaks called using good samples; 
    * see "noblacklist" for peaks that filtered out peaks overlapping blacklist region (based on peaks generated using combined good sample BED files); 
    * see "QC" for peaks that pass QC filtering:
      "fc2.5_p0.01" - fc > 2.5 and -log10p greater than *quantile0.01*, [reference](https://github.com/crazyhottommy/ChIP-seq-analysis?tab=readme-ov-file#peak-calling)
      "fc2.5" fc > 2.5 only 

  * For Hi-TrAC: 
  `macs3 callpeak -q 0.01 --keep-dup all -g hs -t input.bed -f BED --nomodel --shift -100 --extsize 200 --outdir . -n output_name`
  Note: *.xls stores the command history, so **do not** delete it. 
  Note: BAM (for some reason) cannot successfully make peak call; use bed 
    * Peaks call for Hi-TrAC data: 
      Using all reads; 
      Using all cis reads; 
      Using cis reads with distance < 5000 bp;
      Using cis reads with distance > 5000 bp.
    * "raw" peaks: output from macs3 
    * "noblacklist" peaks: remove peaks overlapping blacklist region 
    * "QC" peak: clean peaks that pass quality control (fc & -log10p, see QC.py).

#### Note 
- about MAPQ: 
(DNase-seq) cutoff at 10 or 30. Between 10-30 MAPQ, thre read is mostly scattered randomly across the genome; however, I did find in one case, there is multiple reads enriched at one real peak site. To count for more data, I decided to use cutoff at MAPQ=10. 
(Hi-TrAC) test run indicates that reads with MAPQ 10-30 are mostly interchromosomal distributed, and random scattered and are mostly noises. To save re-run time-cost, I still use MAPQ=30.

- about reads-pairs with different MAPQ: 
the issue is that in some cases, one end of the read pass the cutoff, and the other end does not pass. It may cause issue for downstream analysis. To fix this problem, I need generate BEDPE format (MACS does not support pairs format). I also wrote a code to filter reads with both end pass the cutoff (in  the bam file).

- DNase-seq peakcall with MACS3: test result see /peak_test
from discussion here https://github.com/macs3-project/MACS/discussions/435, people change paired BAM to single-end BED for peak-calling using "--nomodel".
* Using BEDPE as input with default settings: 
`macs3 callpeak -q 0.01 --keep-dup all -g hs -t input.bedpe -f BEDPE --outdir . -n output_bedpe`
* Using BED as input with "--nomodel" combined with "--extsize" and "--shift": 
`macs3 callpeak -q 0.01 --keep-dup all -g hs -f BED --nomodel --extsize 200 --shift -100 -t input.bed --outdir . -n output_name`
    --shift = -1*half of extsize
    extsize = smaller value may give better resolution. As it is combined with Hi-TrAC data, high-resolution may not be accessible. 
* --nomodel vs. --model: 
model show: more distinguished peaks when tethered together (better resolution in distinguishing peaks closeby); 
nomodel show: unite peaks when closeby, and more sensitive for minor peaks; 
* well-resolutioned peaks may not be that useful here, as we need to connect DNase-seq with Hi-TrAC data; also, Hi-TrAC have better power in detect weak peaks but DNase-seq does not have sufficient coverage for this purpose. 
* In conclusion, --nomodel is a better choice. 

- about peak QC
To make method writing easy and less confusion, filtering peaks based on fc only (DHS fc > 2.5). It is also based on the observation that more stringent p-value may loss more information. 

- Normalize in bigwig:
by default, bamCoverage does not perform normalization.

- Peaks call for Hi-TrAC data: 
    Using all reads; all peaks but with many minor peaks 
    Using all cis reads; based on observation, cis data give peaks more scattered when closeby peaks
    Using cis reads with distance < 5000 bp; large background and less reads 
    Using cis reads with distance > 5000 bp. less background but local peaks without interaction cannot be located. 
Given above information, I decided to use all reads for Hi-TrAC peak detection. 

- About pairs generated from "pairtools parse" 
https://pairtools.readthedocs.io/en/latest/parsing.html 
It is noticable that "pairtools parse" behaves different from my default thought. They classified pair types given the read alignment condition. In the cases I found, when a read showing cigarstring "30M21S", which has 21 bp soft-clipped, it considers this read as containing ligations, or from at least two different positions of the genome. Although it can be uniquely aligned, or "rescued", its actual source of position is not sure. It is labeled as "null rescued" (NR). 

(I considered all of them as "unique-unique" UU previously when it passes MAPQ control)

Although this default setting of pairtools may cause some lossing of information, it helps retain high-quality PETs for downstream analysis. 
I thus proceeded with "pairtools parse".

- Trace the progress of Hi-TrAC libraries
- 02/13/2025: I used three metrics (see above), including SPOT score, TSS enrichment, and distant cis PETs ratio, to set cutoffs for good quality libraries. 
    Up to now, samples comparable for downstream analysis: 
    - nUVECp
    - nCAEC
    - nCASMC


#### Important update 
- correlation analysis indicated poor correlation of "nUVEC_1_GC4703735" versus other samples. 
I then dropped "nUVEC_1_GC4703735" (DNase-seq). 
BAM files updated under "/data/jim4/Seq/primary_cell_project/alignment/DNase/merged/pos_sorted_w_good_sample". 
`samtools merge -r -@ 48 -o - /data/jim4/Seq/primary_cell_project/alignment/DNase/merged/split_good/nUVEC_*.bam | samtools sort -@ 48 - -o uUVEC.bam && samtools index -@ 48 uUVEC.bam`

