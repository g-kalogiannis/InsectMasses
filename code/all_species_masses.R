# Load Packages
library(readxl)
library(tidyverse)
library(taxizedb)
library(parallel)

# Read configuration file - handle quotes and NA properly
config <- read.csv("../data/dataset.csv", 
                   stringsAsFactors = FALSE, 
                   na.strings = c("", "NA"),
                   quote = "\"",           # Handle quoted fields
                   strip.white = TRUE)     # Remove extra whitespace

# Function to load data based on configuration
load_dataset <- function(row) {
  dataset_name <- row$dataset_name
  file_path <- row$file_path
  file_type <- row$file_type
  
  # Skip manual entries
  if (file_type == "manual") {
    return(NULL)
  }
  
  cat(paste("Loading:", dataset_name, "\n"))
  
  # Load file based on type
  if (file_type == "csv") {
    sep_char <- switch(as.character(row$separator),
                      "comma" = ",",
                      "tab" = "\t",
                      ",")
    if (!is.na(row$encoding)) {
      data <- read.csv(file_path, sep = sep_char, encoding = row$encoding)
    } else {
      data <- read.csv(file_path, sep = sep_char)
    }
  } else if (file_type == "xlsx") {
    sheet_num <- ifelse(is.na(row$sheet), 1, as.numeric(row$sheet))
    data <- read_xlsx(file_path, sheet = sheet_num)
  } else if (file_type == "rdata") {
    env <- new.env()
    load(file_path, envir = env)
    # Get the first object from the loaded environment
    obj_name <- ls(env)[1]
    data <- get(obj_name, envir = env)
  }
  
  # Apply post-processing if specified
  if (!is.na(row$post_processing)) {
    if (row$post_processing == "calculate_mean_row") {
      columns <- names(data)[4:224]
      data[data == "NA"] <- NA
      data[columns] <- lapply(data[columns], as.numeric)
      data$mean_row <- rowMeans(data[, c("FoL_var_male_average", "FoL_var_female_average", 
                                          "FoL_HR_average", "FoL_Ten_average")], na.rm = TRUE)
    } else if (row$post_processing == "fix_header") {
      colnames(data) <- data[2,]
      data <- data[-c(1:2),]
    }
  }
  
  # Apply filter if specified
  if (!is.na(row$filter_condition)) {
    data <- data %>% filter(eval(parse(text = row$filter_condition)))
  }
  
  return(list(
    data = data,
    config = row
  ))
}

# Data Standardisation & Organisation Functions
check_data <- function(data){
  data <- data %>%
    filter(complete.cases(data$Species),
           !grepl("[0-9]|(sp|spp)$|(sp.|spp.)", data$Species),
           Value > 0)
  data$Species <- str_replace_all(data$Species, "\\s*\\([^\\)]+\\)", "")
  return(data)
}

bodymass <- function(file, species, value, trait, metric, state, estimate, doi){
  # Convert NA to string "NA" for validation
  state <- ifelse(is.na(state), "NA", state)
  estimate <- ifelse(is.na(estimate), "NA", estimate)
  doi <- ifelse(is.na(doi), "NA", doi)
  
  if (!state %in% c("Dry", "Live", "NA")){
    stop("Invalid collection type.")
  }
  if (!estimate %in% c("Yes", "No", "NA")){
    stop("Invalid estimate type.")
  }
  if (!metric %in% c("mg", "g", "kg", "cm", "mm")){
    stop("Invalid metric.")
  }
  if (!trait %in% c("mass", "size")){
    stop("Invalid trait")
  }
  
  # Fix encoding issues before processing
  species <- iconv(species, from = "latin1", to = "UTF-8", sub = "")
  
  tmp  <- data.frame(Species = species,
                      Trait = rep(trait, length(species)),
                      Value = as.numeric(value),
                      Metric = rep(metric, length(species)),
                      Collection = rep(state, length(species)), 
                      Estimate = rep(estimate, length(species)),
                      doi = rep(doi, length(species)),
                      stringsAsFactors = FALSE
                      ) 
  tmp$Species <- sub("_", " ", tmp$Species)
  tmp <- na.omit(check_data(data = tmp))
  output <- rbind(file, tmp)
  return(output)
}

# Process all datasets
file <- c()

for (i in 1:nrow(config)) {
  row <- config[i, ]
  
  # Skip manual entries
  if (row$file_type == "manual") {
    # Handle manual entries separately
    if (row$dataset_name == "diorhabda") {
      diorhabda <- data.frame(
        Species = "Diorhabda carinulata", 
        Value = 0.01089, 
        Trait = "mass", 
        Metric = "g", 
        State = "Live", 
        Estimate = "No", 
        doi = "https://hdl.handle.net/10539/25023"
      )
      file <- bodymass(file, diorhabda$Species, diorhabda$Value, diorhabda$Trait, 
                      diorhabda$Metric, diorhabda$State, diorhabda$Estimate, diorhabda$doi)
    }
    next
  }
  
  # Load dataset
  tryCatch({
    dataset_info <- load_dataset(row)
    if (is.null(dataset_info)) next
    
    data <- dataset_info$data
    cfg <- dataset_info$config
    
    # Extract columns
    species_col <- cfg$species_col
    value_col <- cfg$value_col
    
    # Add to file
    file <- bodymass(
      file, 
      data[[species_col]], 
      data[[value_col]], 
      cfg$trait, 
      cfg$metric, 
      cfg$state, 
      cfg$estimate, 
      cfg$doi
    )
    
  }, error = function(e) {
    cat(paste("Error loading", row$dataset_name, ":", e$message, "\n"))
  })
}

print(paste("Data entry complete. Final dataset dimensions:", nrow(file), "x", ncol(file)))
print(paste("Number of unique species in dataset:", n_distinct(file$Species)))

# Check species with size and measured mass
species_count <- file %>%
  filter(Estimate %in% c("Yes", "No")) %>%
  group_by(Species) %>%
  filter(
    any(Trait == "size" & Estimate == "No") &
    any(Trait == "mass" & Estimate == "No")
  ) %>%
  distinct(Species) %>%
  nrow()

print(paste("Number of species with size and measured (non-estimated) mass:", species_count))


# Link taxonomy
insect_mass = file; colnames(insect_mass)[c(1,7)] = c("species", "source_doi") # Create a duplicate file 
ids <- name2taxid(unique(insect_mass$species), out_type="summary", db = "gbif") # Get ncbi id for taxonomy from name
classes = classification(ids$id) # output taxonomy from id

process_df <- function(df) { 
  df %>%
    t() %>%
    as_tibble() %>%
    slice(-c(2, 3)) %>%
    set_names(make.unique(as.character(df$rank))) %>%
    as.data.frame(stringsAsFactors = FALSE)
} # Function to reorder classification data
processed_data <- mclapply(classes, process_df, mc.cores = detectCores()) # Apply function using multi-core lapply (it's faster)
taxonomy <- bind_rows(processed_data) # Bind output rows

taxonomy = taxonomy[, c("class", "order", "suborder", "family", "genus", "species")] # subset taxonomy
insect_mass = merge(taxonomy, insect_mass, by = "species") # merge taxonomy with data
insect_mass = insect_mass %>% filter(class == "Insecta") # reorder columns
