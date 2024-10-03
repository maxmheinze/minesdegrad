
message("Setting paths.")

DATA <- "/data/"
OUT <- "/home/max/minesdegrad/output/"
CODE <- "/home/max/minesdegrad/code/"


# Data Path ---------------------------------------------------------------

pd <- \(path, ..., pre = DATA) { 
  paste0(pre, path, ...)
}

if(!dir.exists(DATA)) {
  message("\tCouldn't find the directory `DATA`. ",
          "Attempting to create the directory...")
  dir.create(DATA, showWarnings = FALSE, recursive = TRUE) # Make sure this exists
}


# Output Path -------------------------------------------------------------

po <- \(path, ..., pre = OUT) { 
  paste0(pre, path, ...)
}

if(!dir.exists(OUT)) {
  message("\tCouldn't find the directory `OUT`. ",
          "Attempting to create the directory...")
  dir.create(OUT, showWarnings = FALSE, recursive = TRUE) # Make sure this exists
}

# Code Path ---------------------------------------------------------------

pc <- \(path, ..., pre = CODE) { 
  paste0(pre, path, ...)
}

if(!dir.exists(CODE)) {
  message("\tCouldn't find the directory `CODE`. ",
          "Attempting to create the directory...")
  dir.create(CODE, showWarnings = FALSE, recursive = TRUE) # Make sure this exists
}

