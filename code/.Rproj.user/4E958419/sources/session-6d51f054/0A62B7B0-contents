write_latex <- function(output_file, df, decimal_precision = 2, append = FALSE, section_title = NULL) {
  # Validate input
  if (!is.data.frame(df)) {
    stop("The input 'df' must be a data frame.")
  }
  
  # Open a connection to the output file
  file_mode <- ifelse(append, "a", "w")
  output_stream <- file(output_file, file_mode)
  
  # Write the section title if provided
  if (!is.null(section_title)) {
    cat(sprintf("\n%% %s\n", section_title), file = output_stream) # As a comment
    # Or as a custom LaTeX command
    # cat(sprintf("\\section*{%s}\n", section_title), file = output_stream)
  }
  
  # Helper function to write LaTeX macros to the file
  write_macro_line <- function(label, value) {
    formatted_value <- formatC(value, format = "f", digits = decimal_precision)
    cat(sprintf("\\newcommand{\\%s}{%s\\xspace}\n", label, formatted_value), file = output_stream)
  }
  # write_macro_line <- function(label, value) {
  #   formatted_value <- formatC(value, format = "f", digits = decimal_precision)
  #   cat(sprintf("\\newcommand{\\%s}{%s}\n", label, formatted_value), file = output_stream)
  # }
  
  # Iterate over each column in the data frame
  for (col_name in names(df)) {
    values <- df[[col_name]]
    if (length(values) > 1) {
      for (i in seq_along(values)) {
        roman_numeral <- as.roman(i)  # Corrected function to convert integer to Roman numeral
        write_macro_line(paste0(col_name, roman_numeral), values[i])
      }
    } else {
      write_macro_line(col_name, values)
    }
  }
  
  # Close the connection to the output file
  close(output_stream)
}


# # Example usage
# df_to_write <- data.frame(
#   stat = rnorm(1), # Single value column
#   pval = runif(1), # Single value column
#   a2tstat = rnorm(4) # Multiple values column
# )
# 
# write_latex_incremental_macros_df("incremental_macros_df.tex", df_to_write)
