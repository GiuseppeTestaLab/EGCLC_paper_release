TopGO_markdown <- "TopGO_all.Rmd" # TopGO script
Folder_TopGO_Output <- "../../../../DataDir/5.DMRInterpretation/TopGO/" # Folder for TopGO outputs

rmarkdown::render(
  input = TopGO_markdown,
  params = list(
    OutputFolder = paste0(Folder_TopGO_Output, "hEGCLCs_vs_hiPSCs/All/"),
    Condition_1 = "hEGCLCs",
    Condition_2 = "hiPSCs"
  ),
  output_file = "TopGO_hEGCLCs_vs_hiPSCs_all"
)
rm(list = ls())

TopGO_markdown <- "TopGO_all.Rmd" # TopGO script
Folder_TopGO_Output <- "../../../../DataDir/5.DMRInterpretation/TopGO/" # Folder for TopGO outputs

rmarkdown::render(
  input = TopGO_markdown,
  params = list(
    OutputFolder = paste0(Folder_TopGO_Output, "hEGCLCs_vs_hPGCLCs/All/"),
    Condition_1 = "hEGCLCs",
    Condition_2 = "hPGCLCs"
  ),
  output_file = "TopGO_hEGCLCs_vs_hPGCLCs_all"
)
rm(list = ls())

TopGO_markdown <- "TopGO_all.Rmd" # TopGO script
Folder_TopGO_Output <- "../../../../DataDir/5.DMRInterpretation/TopGO/" # Folder for TopGO outputs

rmarkdown::render(
  input = TopGO_markdown,
  params = list(
    OutputFolder = paste0(Folder_TopGO_Output, "hiPSCs_vs_hPGCLCs/All/"),
    Condition_1 = "hiPSCs",
    Condition_2 = "hPGCLCs"
  ),
  output_file = "TopGO_hiPSCs_vs_hPGCLCs_all"
)
rm(list = ls())

TopGO_markdown <- "TopGO_all.Rmd" # TopGO script
Folder_TopGO_Output <- "../../../../DataDir/5.DMRInterpretation/TopGO/" # Folder for TopGO outputs

rmarkdown::render(
  input = TopGO_markdown,
  params = list(
    OutputFolder = paste0(Folder_TopGO_Output, "iMeLCs_vs_hPGCLCs/All/"),
    Condition_1 = "iMeLCs",
    Condition_2 = "hPGCLCs"
  ),
  output_file = "TopGO_iMeLCs_vs_hPGCLCs_all"
)
rm(list = ls())
