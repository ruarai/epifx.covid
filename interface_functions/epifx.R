
venv_prefix <- ". venv_dev_local/bin/activate; "

create_forecast_file <- function(template_file, forecast_args, output_file) {
  
  template_file_text <- readr::read_file(template_file)
  
  
  for(row_i in 1:nrow(forecast_args)){
    template_file_text <- str_replace_all(template_file_text, fixed(forecast_args$key[row_i]), forecast_args$value[row_i])
  }
  
  
  readr::write_file(template_file_text, output_file)
}
