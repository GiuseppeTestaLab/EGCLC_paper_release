library(DSS)

#to modify
Type_1 <- "hEGCLCs_p10"
Type_2 <- "hiPSCs"

#to not modify

bsseq_obj <- readRDS("/group/testa/Project/PGCLC/emseq_rebuttal/2.PreliminaryAnalysis/bsseq_obj_sharedby75ofall.rds")  #getting bsseq object where sample selection and CpG filtering were alredy performed 

group_1 <- pData(bsseq_obj)[pData(bsseq_obj)$Type %in% Type_1, "SampleID"]
group_2 <- pData(bsseq_obj)[pData(bsseq_obj)$Type %in% Type_2, "SampleID"]

dmlTest <- DSS::DMLtest(bsseq_obj, group1=group_1,
                        group2=group_2, smoothing = TRUE)  #performing differential methylation test

output_file_name <- paste0("/group/testa/Project/PGCLC/Rebuttal_paperEG2025/Output_DM/", Type_1, "vs", Type_2, ".rds")

saveRDS(dmlTest, file = output_file_name) #saving output of the test

