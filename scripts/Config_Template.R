library(configr)

cfg <- config(
  input_file_prefix = "input_prefix",
  time_point = 24,
  mutant_name = "WT",
  pos_ctrl_name = "MT",
  plot_title = "plot_title",
  pos_ctrl_legend_name = "Inactive 2A",
  test_cond_legend_name = "Active 2A",
  hit_output_file_name = "1243_hits_10_2023",
  output_file_name = "1243_filtered_10_2023",
  output_file_path = "10_10_2023_output"
)

write_config(cfg, "config.yaml")
