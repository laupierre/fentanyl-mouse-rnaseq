# Fentanyl mouse RNA-Seq

This is the RNA-Seq data from the hippocampus of 3 groups of mice treated with fentanyl (W, C, F groups)  

1- Create the mouse STAR index for 50 nt single-end reads (mouse_star_index.pbs)  
2- Run the STAR pipeline  (star_pipeline.pbs)  
3- Get the annotations from mouse Gencode (gene_annotation.R)  
4- Analyze the different groups with DESeq2  (deseq2_wald_two_groups.R)  

Four files are created including:  
star_gene_raw_counts.xlsx (This is the input file for the DESeq2 analysis)  

and the results of the 3 comparisons:  
hippocampus_deseq2_FentanylvsControl_differential_expression.xlsx  (F vs C)  
hippocampus_deseq2_WithdrawalvsControl_differential_expression.xlsx (W vs C)    
hippocampus_deseq2_WithdrawalvsFentanyl_differential_expression.xlsx  (W vs F)  
