**Table of content**

- [TO DO](#to-do)
- [Consensus](#consensus)
- [Loops](#loops)
- [Expression](#expression)
- [Compartment](#compartment)
- [Insulation](#insulation)
- [DHS](#dhs)
- [Updates](#updates)

#### TO DO
  - [x] Connect loops with target gene (build a model for loop classification based on gene); 
  - [x] update RSEM gene quantity by applying custom cutoffs on the BAM file (unique and contain MAPQ > 30); 
  - [x] gene filtering on expression data (featureCount and TPM), note that TPM is preferred in between sample comparisons; 
  - [x] loop frequency is correlated with gene expression (this is seen in high interaction loci, and its overlapping genes are constitutively highly expressed.)
  - [?] adapt ABC model for the analysis of loop contact frequencies and expression variation; 
  - [ ] start from differentially expressed genes, and see how the loop interactions vary; 
  - [ ] concerted analysis on variable DHS/loops and enriched motif and TF; 
  - [ ] examine the quality of loop annotation; 
  - [ ] add H3K27ac to loop understanding; 
  - [ ] TAD & its relationship with AG loops; 
  - [ ] classify genes into GG-prone, or IG-prone; 
  - [ ] For genes overweighted by GG-loops, check the Gene-gene expression relations; 
  - [ ] what is chromatin stripes in landscape ? stripe peak? 
  - [ ] **linear mixed effect model apply on DNase-seq for testing**
  - [x] genome-wide signal profile correlation between DNase and Hi-TrAC; 
  - [x] use DNase signal to find open TSS + gene expressed to filter gene for further consideration about its regulation by loop interactions (**no longer consider**);
  - [x] super interaction loci and SIL-associated genes; 
  - [x] Organize GWAS variants information from GWAS Catalog
  - [ ] Organize ClinVar variants information (also see dbSNP)
  - [ ] how to test enrichment of disease-associated SNPs in SILs (variable SILs, permutation test, chi-square test), including alleviating LD problem.
  - [x] if loop strength correlated with gene expression;  
  - [x] compare to published DNase-seq database, including DHS from publication of 2022 and SCREEN;
  - [x] conservative score profile around DNase-seq peak loci;
  - [x] read count for DNase-seq peak and its entropy score and [static score](https://www.nature.com/articles/s41586-020-2151-x)
  - [x] specific peak overlap with DHS index dataset
  - [x] annotate DHS peaks (using homer)
  - [x] from DNase-seq specific peaks, to link peak to its gene through loop data and see what pathways enriched.
  - [x] Try FitHChIP (~11.7k loop can be identified, ideal number of peaks should be around 100k)
  - [x] Try MANGO (not succeed in installing the program)
  * [x] run loop detection for each sample separately, and then compare to the joint loop detection result (combining all samples together). For celltype specific loops, check its detection differences in the two methods (separate and joint)
  - [x] CTCF orientation within loop connections (the role of thumb is consistent: +- convergence dominant, -+ divergent orientation less)
  - [x] compare loops identified in Hi-TrAC with 24 celltype unified loops (n=124830)
  - [x] combine all sample (3 celltypes for now) for pooled cool 
  - [x] collect blood DHS data & download them;
  - [ ] when all celltypes DNase-seq data available, integrate blood with all the other primary celltype DNase-seq data.

#### Consensus 
To compare among samples, consistent peaks based on all input peaks are generated. 

- Consensus DNase-seq: 
    - blood DNase-seq here: `/data/jim4/Seq/primary_cell_project/data/blood_DHS/peakcall`; 
    - evaluation for DNase-seq: `/data/jim4/Seq/primary_cell_project/alignment/DNase/eval/eval_score.xlsx`; 
    - DNase-seq resolution may depend on the sequencing coverage, for shallow sequencing, peak summit may not be accurate. if set very strict summit distance for clustering, it may inherit noises from the sequencing data into downstream analysis (use 100/150 bp for summit clustering). 
    - In conclusion, it may still useful to have the DNase-seq data report in a best-of-resolution for the data can possibly provide;
    - base peak file used for consensus peakcall, see `/data/jim4/Seq/primary_cell_project/analysis/DHS/peak`; 
    - fc2.5 (see `/data/jim4/Seq/primary_cell_project/alignment/DNase/peak/QC/fc2.5`)
        - I used "fc2.5" as the final cutoff to call consensus peaks, run: 
        - `peak_consensus.py -peaks /data/jim4/Seq/primary_cell_project/analysis/DHS/peak/*Peak -min_dist 100 -o DNase_consensus_dist100 -saf`
          - "-saf": generate saf for featureCount input (see below);
    - featureCount of DNase-seq consensus: 
        - `ml subread`
        - `featureCounts -F SAF -p --countReadPairs -B -T 21 --donotsort -a $saf -o DNase_consensus_dist100.featureCount -s 0 input1.bam input2.bam input3.bam`
          - "-p": paired-end reads; 
          - "-B": both ends assigned to the feature; 
          - "-s": strands ignored; 
    - normalization of read-count: 
        - I applied RPM normalization by dividing on the total-counts;
        - `count_norm.R DNase_consensus_dist100.featureCount -norm RPM -o DNase_consensus_dist100.featureCount.norm`
    - Correlation analysis: 
        - see "func.r/corr_func"; 
        - For correlation between replicates and among samples;
    - PCA based clustering:   
        - Based on the PCA results, I can then clustering celltypes based on their PCA distance [ref](https://doi.org/10.1038/s41586-020-2023-4)
          - Use cumulative PC distances weighted by its explained proportion of variance to cluster similar celltypes together; 
        - All blood celltypes form a separate branch; 
    - cluster-specific peaks; 
        - I considered t.test (default welch two samples), wilcox.test (non-parametric t.test); 
        - For small sample size, not follow normality assumption, t.test is not good, and I then used wilcox.test; 
        - for wilcox.test (ranked sum), real value ignored, to add the signal strength, I calculated its score using p-value and mean of signal strength of the cluster; 
        - currently, top 5 k peaks are applied; 

- Consensus Hi-TrAC peakcall
    - updated "peak_consensus.py" for peaks that overlapping each other within min_dist, then clustered. 
    - Use: `peak_consensus.py -peaks sample1.narrowPeak sample2.narrowPeak sample3.narrowPeak -min_dist 500 -o HiTrAC_consensens_peak`
    - update on 03/12/2025, bin-based clustering (when peaks overlapping the same bin at resolution 1 kb, peaks merged, no summits considered). To maximize the resolution power in peak interaction call. 
    - See `HiTrAC_cluster.bed`


- Annotation peaks 
    - Use [homer](http://homer.ucsd.edu/homer/ngs/annotation.html) to annotate peaks directly.  
      - although I have wrote a script to do so, it is quite slow.  
    - Use: `annotatePeaks.pl $bed $genome -cpu 4 -annStats annotate.stat > bed_annotated`  

- disease vs. normal DHS 
  - see analysis/disease; 
  - comparison between fLF vs. nLF; 
  - generated consensus DHS peaks: `peak_consensus.py -peaks *narrowPeak -min_dist 100 -o nfLF -saf` ; 
  - perform featureCount: `featureCounts -F SAF -p --countReadPairs -B -T 21 --donotsort -a $saf -o disease.featureCount -s 0 *bam`
  - normalization and differential analysis: see `diff_analysis.Rmd`
  - test for TF enrichment: `TF_enrich.py -bg nfLF.bed -query diff_nfLF.bed -TF /data/jim4/Seq/primary_cell_project/data/motif/motif_from_hg19/*bed -o diff_nfLF_TF_overlap`

#### Loops 
- minimum counts at 5 is reliable, and minimum dist > 5000 should be applied. 
- 
=======

1. Hi-TrAC 1-D profile comparing to DNase-seq profile
   - run `bw2bedGraph.py` to retrive the signal intensity using a specific step size; 
     - `bw2bedGraph.py input.bw -step 5000 -o output`
   - To combine all score into one dataframe: 
     - run `script/HiTrAC_DHS_cor.py`;
   - Use "Pearson's correlation coefficient" for its signal strength similarity across the genome (example see `expl.R`); 
   -  signal profile comparing between HiTrAC and DNase-seq:
     - DNase-seq profiles are the most similar to itself; Hi-TrAC profiles are the most similar to itself;
     - When comparing between Hi-TrAC & DNase-seq: Hi-TrAC is the most similar to the same celltype's DNase-seq profile;
     - results see "Loops/output/profile_cp_DNase"

2. loop detection 
    - Use "HiTrAC_cluster.bed" as anchor (see "anchor_processing.py" about its clustering strategy);
      - Specifically, when anchors have overlapped bin(s) at specific resolution (1 kb), they are merged together; 
    - Use: (in `discovery` mode)
        `mpiexec -n 12 peak_interaction.py discovery -cool sample_1k.cool -peak HiTrAC_cluster.bed -w 5 10 100 -fdr 0.01 -o HiTrAC_interaction`

3. Super Interactive loci 
   - methodology (see `SIL.py`): 
     - highly interactive anchors are list as super interactive anchor candidates;
     - use 10 kb distance as cutoff, tether these candidate loop anchors when distance less than the cutoff; 
     - for each celltype's loop data, run `SIL.py loop.bedpe -min_interaction 3 -tether_dist 10000 -report_loop -o output`
   - Get loop data for "GM12878" and identify SIL, to compare SIL and SE of its own celltype; 
     - the previously reported loops for GM12878 is not high-quality and complete, very few SIL reported, and the anchor/tether-anchor distribution is also different; 
     - I need to get the cool data for GM12878, and call loops using my method; 
     - TBC ... 
   - GM12878: Lymphoblastoid (transformed B cells)
   - celltype specific SIL should be comparing to SE of the same celltype; 
     - see `script/SIL_SE_compare.py`; 
     - 

4. Count read for each of the significant interactions 
    - Use `peak_interaction.py count` mode
     - `mpiexec -n 4 peak_interaction.py count -cool nCAEC.mcool -interaction interaction.txt -o count`
    - normalization: the total number of reads (on loop sites) is similarily agree with the genome-wide PETs. I can calculate "RPM" before enable between sample comparisons.
     - use `script/loop_count.py`
    - To elleviate potential biases imported by celltype specific loops, I filtered loops with RPM > 1 in all samples (which means in most cases, the loops are shared among the celltypes in study. And this is also true in considering that minimum loop > 5 k); 
    - 

5. CTCF orientation within loops (use joint loopcall)
    - identify loops overlap CTCF motif; 
    - there are multiple motifs known/discovered for CTCF binding. To make the searching of its motif be simple, I selected one motif specifically here (CTCF_known1, 'YDRCCASYAGRKGGCRSYV'); 
    - See `CTCF_orientation.py`.
    - based on data, I found that: 
        > set a score cutoff for CTCF motif can improve the "convergence orientation" ratio;
        > most of '+-' (convergence orientation) are long-distant, and instead, short-distant interactions are more likely be false positives.
    - when set motif-score cutoff (score > 15), and set distance cutoff (> 10 k), I found: 
        > convergent orientation dominant long-distance; 
        > divergent orientation dominant short-distance. 

6. joint master cool file: 
    - by combining all samples' pairs file, master cool file was generated. 
    `pairtools merge --memory 3G -o master.pairs.gz --nproc 4 nCASMC.pairs.gz nCAEC.pairs.gz nUVECp.pairs.gz`
    `cooler cload pairs $hg38_chrsize:1000 master.pairs.gz master.cool --assembly hg38 -c1 2 -p1 3 -c2 4 -p2 5`

7.  Loop annotation 
    - see `loop_annotation.py` for detail (report annotation on selected transcript level and its frequency for each type of loops classified);
      - `loop_annotation.py loop.bedpe -transcript ncbiRefSeqSelect.bed -o loop_ann`
    - there are 5 types of loops: SG, GG, IG, AG, NG;
    - Update RefSeqSelect on 04/16/2025 (original data was download from [UCSC](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/));
    - 

8. public loop/domain data 
   - https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63525 (Rao, 2014)


9. published Hi-TrAC data contains "XX" in pairs file 
  - I need to figure out how could "XX" present in pairs files (for example, K562);
  - for K562.n.bam, 442622878 total reads (222825581 + 0 read1, 219797297 + 0 read2); 
  - for K562.bam, 442622878 reads in total (read number no change between formats);
  - filter K562.n.bam, make it clean file K562.nn.bam; 
  - for K562.nn.bam, 348815682 total reads (174407841 PETs);
  - for K562.pairs, 174407841 UU PETs (match to clean BAM file); 
  - although I cannot understand why "XX" in pairs file, but all valid reads are retained in pairs file; 


#### Expression

1. Note that I updated the BAM alignment by applying stringent QC on RNA-seq (03/31/2025)
   - pleaset see "alignment/RNA/star_align/clean/" for detail (briefly, unique reads, MAPQ > 30, and reside on chromosome region);
   - for transcriptome alignment, reads have met the criteria of genome alignment file. 

2. Bad sample (final, 04/02/2025)
    - Use different strategies, correlation & PCA, indicated at least one sample of "nDFn" low quality; 
    - For record, **nDFn_1** and **nDFn_3** removed (for example, MYC gene expression is very similar between replicates, and __nDFn_3 is very different from other replicates__); 

3. Entropy 
    - See "script/entropy.py" [ref](https://www.nature.com/articles/s41586-020-2151-x), to identify consitutive expressed genes.
    - constitutive: active, and stably expressed in all celltype in study. 

4. T-statistics 
    - this is recommended when the sample size is relatively large, to compare the distribution of samples' mean to its population mean; 
    - T = (X – μ) / [ s/√(n) ] (see `t.r` for t-statistics function)
    - Used to identify gene specific to a cell-type.  
    - The higher value indicates a more biased expression in the specific celltype compared to other celltypes. 

5. linear mixed effects model (expression ~ celltype)
    - it functions similar as t-statistics, used to detect genes with celltype specific expression pattern; 
    - 

6. CP_enrich 
    - Given specific gene set, perform GO enrichment analysis based on canonical pathway terms. 
    > Example:  
    `Rscript /data/jim4/Seq/primary_cell_project/code/GSEA.R -gmt /data/jim4/Seq/primary_cell_project/data/MSigDB/c2.cp.v2023.2.Hs.symbols.gmt -gene gene_set.txt -O enrich_cp -bg bg_gene.txt`  
    **Note:**  
    given gene should have columns "symbol";  
    bg or background gene set is optional.

#### Compartment 
- dchic used for compartment call; 
  - `mamba activate dchic`
  - use `dchic_pipe.py`; 
    - `dchic_pipe.py -cool input.mcool -r 5000 -ref_path /data/jim4/Seq/primary_cell_project/analysis/compartment/data/hg38_goldenpathData/ -o output`
    - inputs: cool file, resolution, ref (hg38), golden path (if previously downloaded), output prefix 
- based on experiment run, 10k can be a good resolution (also try 5k); 
- Analysis between loops and compartment: 
  - loops are connecting A-A compartment (> 85%), and few percentage are from B-B interaction; 
- Correlation between Histone modification marks & compartment: 
  - see `compartment_mark.py`; 
  - 
- Gene expression comparison: 
  - "A" genes: gene completely reside within A; 
  - "AB" genes: gene across A and B compartment; 
  - "B" genes: genes completely reside within B; 
  - global analysis indicated: A genes' expression is significantly higher than AB genes, and all more than B genes (that is genes reside completely within B compartment)
- 

#### Insulation
- Consider insulation as the way to segment genome into domains: between domains, the insulation score is greater; 
- Method: slide a diamond-shape window along the genome (one corner on the diagonal of the matrix), and score for the sum of contacts within the window for each position. At certain location, the score is significantly lower than its surrounding region (reflecting lowered contact frequencies between upstream and downstream loci), this position is referred to as boundary position. 
- see domain paper [here](https://www.nature.com/articles/nature11082);
- 

#### DHS 
- This folder is to store all final processed data for DHS; 
  - `bigwig`: folder for all bigwig files generated (RPKM normalization by deeptools);
  - `bam`: folder for clean bam files per cell-type (required for featureCounts);
  - `peak`: DHS peaks that pass filtering; 
  - 
- visualize specific DNase-peak using deeptools: 
  - all bigwig is normalized using RPKM; 
  - example: `computematrix scale-regions -S sample1.bw sample2.bw -R input.bed -b 1000 -a 1000 -o region.mat.gz -bs 20 --skipZeros`
  - `dir="/data/jim4/Seq/primary_cell_project/analysis/DHS/bigwig"`
  `computeMatrix scale-regions -S $dir/{MC_,NK_,B_,CD4,CD8}* $dir/{nEK,nGK,nSAEC}* $dir/{nDF,nLF,ub}* $dir/*SMC*.bw $dir/{nAEC,nCAEC,nPAEC,nUVEC}* $dir/nEMn.bw -R /data/jim4/Seq/primary_cell_project/analysis/DHS/output/cell_specific/*bed -b 1000 -a 1000 -o DHS_cluster.mat.gz -bs 20 --skipZeros`
  - `plotHeatmap -m DHS_sp.mat.gz -out DHS_sp.png`
- TF motif enrichment for cluster-specific peaks: 
  - background use can affect the enrichment results; 
  - when I tried to combine all cluster-specific peaks (see `DNase_bg1.bed`), and then use it as background, I can see that smc is no much enriched for any TF, and PBMC enriched for many; this is because smc selected peaks are not contrasting to other selected peaks (not strongly), and instead, pbmc selected peaks are constrasting other peaks very strongly and its enriched TFs are the most; 
  - I tried to use peaks that not showing any significant differences as background (plus its own query peak), for test, I run smc specific peak with updated background (see `DNase_bg2.bed`). It is clear that more TF show up as significant enrichment. 
  - To avoid biases imported by different background applied, try use all as background (`DNase_bg3.bed`); 
  - based on test on different background datasets, I still decided using **whole dataset**. The enrichment is the most significant and most factor can came out, it also do not inherite potential biases imported by varying background choices. 
  - For results, see `TF_enrich/bg3` folder.
- Conservation scores on unified DHS set: 
  - for fast test: `computeMatrix scale-regions -S /data/jim4/Seq/primary_cell_project/data/Phylo/hg38.phyloP100way.bw -R DNase_consensus_dist100.bed -bs 10 --skipZeros -b 1000 -a 1000 -o DNase_consensus_PhyloP100way.mat.gz`
  - previously, I wrote script `peak_conservation.py` for this purpose; 
  - I integrated conservation examination into `bed_profile.py` (used to examine profile scores across all regions in bed); 
  - Update (05/01/2025), use `bed_profile.py`; (under development)
    - example: `bed_profile.py <region.bed> <score.bw> ` 


#### Updates
- `peak_consensus.py` 
  - Simplified the process by calling on peak cluster and summit cluster only. For peak cluster, multiple summits can be retained if their distance > min_dist cutoff. 
  - when multiple summits retained for the same peak interval, to improve the resolution in binding count, peak interval break at the midpoint between two neighboring summits. 

- DNase-seq min_dist used for consensus peak call 
I tried min_dist -> 150, 100, 50, 20. 
For downstream analysis in high-resolution, 20 bp used for now. 

- 02/12/2025: Many Hi-TrAC libraries are not ready. I generated consensus using CAEC, CASMC, UVECp. 
control fc2.5 to keep consistent with peak QC on DNase-seq data. 
To do on Hi-TrAC consensus:
    - loop static score and its association with gene expression; 
    - static score may be replaced by coefficient of variation (CV); 

- 04/04/2025: About gene filtering:
  - Initially, I wanted to filter genes based on in-house expression and DNase-seq data (on transcripts). this is not working well, expression is not easy to use as, for example, some genes may be tandemly duplicates, and the expression level for each is unclear; DHS site may be overlapping two gene TSS region when, for example, they are facing against each other; 
  - Also, if filtering on in-house data, the loop classification may be restricted and cannot have a broader application when different geneset utilized; 
  - third, the refSeq select gene set is per the efforts on NCBI team who build it with multiple factor in consideration, including expression level, evolutionary conservation, and functional significance. 
  - Importantly, it can simplify my work :)

- 