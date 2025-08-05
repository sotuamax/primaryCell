# primaryCell
scripts developed for the primary cell project which includes Hi-TrAC, DNase, RNA-seq data.

- To generate HiC/cool/pairs format based on BAM input, use: 
  `hic_format.py`
- To pipeline HiTrAC data, use (need tests): 
  `hic_pipe.py ` 
- To call compartment using dchic, use: 
  `dchic_pipe.py`
- Use loop data to call SIL: 
  `SIL.py `
- To find SIL target gene: 
  `SIL_gene.py`
- Annotate loops in bedpe format: 
  `loop_annotation.py`
