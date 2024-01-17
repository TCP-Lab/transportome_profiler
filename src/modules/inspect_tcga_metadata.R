suppressMessages({
  library(tidyverse)
})

args <- commandArgs(trailingOnly = TRUE)

data <- read_csv(args[1])

# Compute missing values for every tumor type for every variable
data |> group_by(study) |>
  summarise(across(everything(), ~ round(sum(is.na(.x)) / n() * 100, 2) )) -> missing_no

write_csv(missing_no, stdout())
