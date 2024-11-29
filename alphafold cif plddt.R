dir <- here::here("alphafold cifs")
files <- list.files(dir)

# Read the already processed data
plddt_scores_cif <- readRDS("plddt_scores_cif.Rds")

# Read the data new from the cif files, unnecessary if readRDS worked
plddt_scores_cif <- c()
i <- 1
for (file in files) {
  print(i)
  lines <- readLines(paste0(dir, "/", file))
  plddt_line <- grep("_ma_qa_metric_global.metric_value", lines, value = TRUE)
  accession <- grep("_ma_target_ref_db_details.db_accession", lines, value = TRUE)
  accession <- strsplit(accession, "                 ")[[1]][2]

  val <- strsplit(plddt_line, " ")[[1]][2]
  plddt_scores_cif[accession] <- as.numeric(val)
  i = i + 1
}

idps <- which(plddt_scores_cif < 50)
names(idps)



