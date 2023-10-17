#load all the required libraries
library(tidyverse)
library(dplyr)
library(tibble)
library(readr)
library(here)
library(stringr)
library(ggplot2)
library(plotly)
library(ggrepel)
library(yaml)

#' next we will define our config param's

#' to start, we'll read the config file in

config <- read_yaml("config_compoundscreening.yaml")

#' now we'll actually specify the parameters

input_file_prefix <- config$input_file_prefix
time_point <- config$time_point
mutant_name <- config$mutant_name
pos_ctrl_name <- config$pos_ctrl_name
plot_title <- config$plot_title
pos_ctrl_legend_name <- config$pos_ctrl_legend_name
test_cond_legend_name <- config$test_cond_legend_name
hit_output_file_name <- config$hit_output_file_name
output_file_name <- config$output_file_name
output_file_path <- config$output_file_path

#'Below are all of the functions we will be using to process HTS data
#'
#' The first function will clean the data, we'll use it when reading the csv in
#' if this is erroring, double check the column name containing absorbance 
#' values, change as needed

csv_edits <- function(data, plate, time) {
  data %>%
    mutate(time = {{time}},
           plate_ID = {{plate}}) %>%
    select(-c(Content, X)) %>%
    rename(OD_600 = Raw.Data..600.) %>% 
    mutate(OD_600 = as.numeric(OD_600)) %>%
    na_if("") %>%
    na.omit
}

#'Function to read file in and process it
#'Avoids manually entering each file name and hard-coded file paths
#'Selects file based on config file prefix
#'All files should be descriptively named as platenum_yy_mm_dd_timepoint(hr) 
#'The mm_yy descriptors are in variable name defined during processing
#'Further descriptors can be added to lab notebook/config file 
#'Further descriptors will then be associated to processed data by date
#'Arguments: 
#'   input_file_prefix: defined by config file, prefix is file name excluding 
#'   hour that all files should contain
import_data_pattern_recog <- function(input_file_prefix) {
  
  file_list <- list.files(path = here("input"), 
                          pattern = str_flatten("^", {{input_file_prefix}}))
  processed_dfs <- list()
  
  for (file in file_list) {
    df <- read.csv(here("input", file))
    file_name <- str_remove(file, "\\.[cC][sS][vV]")
    
    components <- str_split(file_name, "_")[[1]]
    variable_name <- str_flatten(c("plate",components[1], components[2], 
                                   components[3], components[5]), "_")
    
    hour <- as.numeric(components[5])
    print(hour)
    plate_ID <- str_sub(variable_name, 7, 10)
    processed_df <- csv_edits(df, plate_ID, hour)
    processed_dfs[[variable_name]] <- processed_df
  }
  return(processed_dfs)
}

#'Function cleans up the names of wells
#'Arguments:
#'   data: name of processed dataframe, should contain all time points
#'   mutant_name: defined in config file, name of condition which all compounds 
#'      tested on
well_name_replace <- function(data, pos_ctrl_name, mutant_name) {
  data %>%
    mutate(
      Well = case_when(
        str_detect(Well, "^[a-zA-Z].*(01|02|23|24)$") ~ 
          str_c(Well, {{pos_ctrl_name}}, sep = "_"),
        str_detect(Well, "^[a-zA-Z].*[0123456789]$") ~ 
          str_c(Well, {{mutant_name}}, sep = "_"),
        TRUE ~ Well
      )
    )
}

#'Function widens existing df to extract string we added with well_name_replace
#'Arguments
#'   plate_df: the processed df with a clean 'Well' column
widen_df <- function(plate_df) {
  pivot_wider({{plate_df}},
              names_from = Well,
              values_from = OD_600)
}

#'Function lengthens df to add a 'Condition' or 'Mutant' column filled 
#'from string
#'Arguments
#'  wide_plate_df: the widened df created from the widen_df function
lengthen_df <- function(wide_plate_df) {
  pivot_longer({{wide_plate_df}},
               cols = -c(time, plate_ID),
               names_to = c("Well", "Condition"),
               names_pattern = "(.*)_(.*)")
}

#'Functions below take cleaned and processed data
#'
#'Function filters cleaned df for a specified time point of interest
#'then calculates the mean and sd OD's at that time point for the condition of 
#'interest (i.e. condition which has compounds added)
#'Arguments:
#'   plate_df: the df we are working with, cleaned and processed
#'   timepoint: defined in config file, the time point of interest
filter_df <- function(plate_df, timepoint) {
  plate_df %>%
    filter(time == {{timepoint}}) %>%
    group_by(Condition, time) %>%
    mutate(mean_OD = mean(value),
           sd_OD = sd(value)) %>%
    rowid_to_column() 
}

#'Function to filter to condition of interest only
#'Arguments
#'   filtered_df: filtered/cleaned df, has filter_df function applied
#'   mutant_name: specified in config file, mutant or condition with compounds
condition_filter_df <- function(filtered_df, mutant_name) {
  filtered_df %>%
    filter(Condition == {{mutant_name}}) 
} 

#'Function to return summary stats for condition of interest
#'Arguments
#'   mutant_df: data filtered for time point and condition
#'   of interest, used to calculate mean/sd
calculate_summary_stats <- function(mutant_df) {
  mean_value <- mean(mutant_df$mean_OD)
  threesigma_value <- 3 * (mutant_df$sd_OD[1]) + mean_value
  return(list(mean_value = mean_value, three_sigma_value = threesigma_value))
}

# Function to create a scatter plot for plate data at specified time point
# Arguments:
#   plate_df: Df cleaned and filtered for time point of interest
#   three_sigma_id: Three sigma value, defined in summary stats
#   mean_value_id: Mean value, definded in summary stats
#   title: Title for the plot
#   group1_label: Label for the first condition (pos control, no comps)
#   group2_label: Label for the second group (test condition, with comps)
plate_plotting_clean <- function(plate_df, three_sigma_id, mean_value_id, title,
                                 group1_label, group2_label) {
  ggplot({{ plate_df }}, aes(rowid, value,
                             colour = Condition,
                             alpha = Condition == "MT")) +
    geom_point() +
    scale_alpha_manual(values = c(1, 0.5), guide = "none") +
    geom_hline(yintercept = {{ three_sigma_id }},
               col = "red", linetype = "dashed") +
    geom_hline(yintercept = {{ mean_value_id }}, linetype = "dashed") +
    theme(legend.key.size = unit(1.5, "cm")) +
    theme_minimal(base_size = 30) +
    theme(text = element_text(size = 30, family = "Helvetica")) +
    annotate("text",
             x = 390, y = {{ three_sigma_id }} + 0.1,
             label = expression(~ 3 * sigma), colour = "red") +
    annotate("text",
             x = 390, y = {{ mean_value_id }} + 0.05,
             label = expression(~mu)) +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.key.height = unit(3, "lines"),
      legend.spacing.y = unit(3, "lines")
    ) +
    ylim(c(0, NA)) +
    xlim(c(0, NA)) +
    scale_color_manual(
      name = "Legend",
      labels = c({{ group1_label }}, {{ group2_label }}),
      values = c("#8575C4", "#5BBE6C", "#5C6BC0")
    ) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black")
    ) +
    ggtitle({{ title }}) +
    labs(x = "Compound", y = "Growth Reading")
}

#' Function to filter to desired info and write to csv file
#'Arguments:
#'   filtered_plate: The filtered df to use, filtered by time point of interest
#'   threshold_value: The threshold value to determine hits (3-sigma)
#'   output_file_path: The path that you would like to add saved files to
#'   output_file_name: The name to give your df with hits or filtered data
filter_and_save_data <- function(filtered_plate, threshold_value, 
                                 output_file_path, hit_output_file_name) {
  hits_filtered <- filtered_plate %>%
    filter(value > threshold_value,
           Condition != mutant_name)  
  
  # Check if the parent directory exists
  if (!dir.exists(dirname(output_file_path))) {
    stop("Parent directory does not exist.")
  }
  
  # Create the target directory
  if (!dir.exists(output_file_path)) {
    dir.create(output_file_path)
  }
  
  output_file <- file.path(output_file_path, {{hit_output_file_name}})
  
  write_csv(hits_filtered, output_file)
}

#'Function just to save data, if saving filtered data is preferred
#'Arguments:
#'   filtered_plate: Filtered for time point of interest
#'   output_file_path: Same output file path defined in config
#'   output_file_name: File output name defined in config
save_data <- function(filtered_plate, output_file_path, output_file_name) {
  # Check if the parent directory exists
  if (!dir.exists(dirname(output_file_path))) {
    stop("Parent directory does not exist.")
  }
  
  # Create the target directory
  if (!dir.exists(output_file_path)) {
    dir.create(output_file_path)
  }
  
  output_file <- file.path(output_file_path, {{output_file_name}})
  
  write_csv(filtered_plate, output_file)
}
#'Below is script execution
#'
#'Read in and process csv files
list_processed_dfs <- import_data_pattern_recog(input_file_prefix)

#'Bind the list into 1 df
grouped_df <- bind_rows(list_processed_dfs) %>%
  arrange(time)

#'Clean up 'Well' column
grouped_df <- well_name_replace(grouped_df, pos_ctrl_name, mutant_name)

#'Widen and lengthen df
wide_df <- widen_df(grouped_df)
long_df <- lengthen_df(wide_df)

#'Filter df for time point of interest, defined in config file
filtered_df <- filter_df(long_df, time_point)

#'Filter df to condition of interest only, defined in config file
condition_filtered_dfs <- condition_filter_df(filtered_df, mutant_name)

#'Calculate summary stats for conditionally filtered df
summary_stats <- calculate_summary_stats(condition_filtered_dfs)

#'Plot plate
plate_plotting_clean(filtered_df, summary_stats$three_sigma_value, 
                     summary_stats$mean_value, plot_title,
                     pos_ctrl_legend_name, test_cond_legend_name)

#'Filter and save hit data
filter_and_save_data(filtered_df, summary_stats$three_sigma_value,
                     output_file_path, hit_output_file_name)

#'Save filtered data
save_data(filtered_df, output_file_path, output_file_name)
