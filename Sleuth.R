#load package
library(sleuth)
R_cwd = getwd()
R_cwd
#read in the table you made describing samples and kallisto output, 
#assign to variable name stab 
stab = read.table("sample_table_final.txt",header=TRUE)
#initialize sleuth object using sleuth_prep function from sleuth library
so = sleuth_prep(stab)

#fit a model comparing the two conditions 
so = sleuth_fit(so, ~condition, 'full')

#fit the reduced model to compare in the likelihood ratio test so = sleuth_fit(so, ~1, 'reduced')
so = sleuth_fit(so, ~1, 'reduced')
#perform the likelihood ratio test for differential expression between conditions 
so = sleuth_lrt(so, 'reduced', 'full')

#load the dplyr package for data.frame filtering
library(dplyr)


#extract the test results from the sleuth object 
sleuth_table = sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE) 

#filter most significant results (FDR/qval < 0.05) and sort by pval
sleuth_significant = dplyr::filter(sleuth_table, qval <= 0.05) |> dplyr::arrange(pval) 

output_dir = paste0(R_cwd, "/PipelineProject.log")

#write FDR < 0.05 transcripts to file
write.table(sleuth_significant[, c('target_id', 'test_stat', 'pval', 'qval')], file = output_dir, quote = FALSE,row.names = FALSE, append = TRUE)

