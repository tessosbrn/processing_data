
library(configr)

#'NOTE ABOUT FILE NAMES:
#'All files should be descriptively named as platenum_yy_mm_dd_timepoint(hr) 
#'The mm_yy descriptors are in variable name defined during processing
#'Further descriptors can be added to lab notebook/config file 
#'Further descriptors will then be associated to processed data by date
#'
cfg <- list(
  input_file_prefix = "2931_23_05_04_",
  time_point = 24,
  mutant_name = "WT",
  pos_ctrl_name = "MT",
  plot_title = "plot_title",
  pos_ctrl_legend_name = "Inactive 2A",
  test_cond_legend_name = "Active 2A",
  hit_output_file_name = "2931_hits_10_2023",
  output_file_name = "2931_filtered_10_2023",
  output_file_path = "10_12_2023_output"
)

write.config(cfg, "config_compoundscreening.yaml", write.type = "yaml")
