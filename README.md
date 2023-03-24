# Fentanyl mouse RNA-Seq

This is the RNA-Seq data from the hippocampus of mice treated with fentanyl (W, C, F groups)  

1- Create the mouse STAR index for 50 nt single-end reads  
2- Run the STAR pipeline  
3- Analyze the different groups with DESeq2  

Four files are created:  
star_gene_raw_counts.xlsx (This is the input file for DESeq2)  

hippocampus_deseq2_FentanylvsControl_differential_expression.xlsx  (F vs C)  
hippocampus_deseq2_WithdrawalvsControl_differential_expression.xlsx (W vs C)    
hippocampus_deseq2_WithdrawalvsFentanyl_differential_expression.xlsx  (W vs F)  
