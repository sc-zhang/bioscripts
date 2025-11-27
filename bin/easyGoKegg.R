#!/usr/bin/env Rscript
library(optparse)
library(clusterProfiler)
library(KEGGREST)
library(dplyr)
library(stringr)
library(AnnotationForge)
library(jsonlite)
library(purrr)
library(RCurl)
library(ggplot2)

download_plant_path <- function(db_file) {
  ### Download plant pathway
  org <- data.frame(keggList("organism"))
  plants <- org[grep("Plants", org$phylogeny), ]
  plant_pathways_total <- vector()

  for (i in seq_along(length(plants$organism))) {
    try({
      pathways <- keggLink("pathway", plants[i, 2])
      pathways <- sub(paste(".*", plants[i, 2], sep = ""), "", pathways)
      pathways <- unique(pathways)
      plant_pathways_total <- append(plant_pathways_total, pathways)
      plant_pathways_total <- unique(plant_pathways_total)
    })
  }

  plant_pathways_total <- paste0("ko", plant_pathways_total)
  write.table(plant_pathways_total, file = db_file, quote = FALSE, row.names = FALSE, col.names = FALSE)
}

load_plant_path <- function(db_file) {
  ### Load plant pathway
  plant_pathways_total <- read.table(db_file, sep = "", header = FALSE)
  colnames(plant_pathways_total) <- "Pathway"
  plant_pathways_total
}

generate_go_db <- function(annotation_file, db_path, genus, species, tax_id) {
  cat("Reading annotation file\n")
  egg <- read.table(annotation_file, header = TRUE, sep = "\t", quote = "")
  egg[egg == ""] <- NA

  gterms <- egg %>%
    dplyr::select("query", "GOs") %>%
    na.omit()

  gene2go <- data.frame(
    GID = character(),
    GO = character(),
    EVIDENCE = character()
  )

  gene_ids <- egg$query
  eggnog_lines_with_go <- egg$GOs != "-" & egg$GOs != ""
  eggnog_annoations_go <- strsplit(egg$GOs[eggnog_lines_with_go], ",")
  gene2go <- data.frame(
    GID = rep(gene_ids[eggnog_lines_with_go],
      times = sapply(eggnog_annoations_go, length)
    ),
    GO = unlist(eggnog_annoations_go),
    EVIDENCE = "IEA"
  )

  gene_info <- egg %>%
    dplyr::select(GID = "query", GENENAME = "Preferred_name") %>%
    na.omit()

  gene2ko <- egg %>%
    dplyr::select(GID = "query", Ko = "KEGG_ko") %>%
    na.omit()

  gene2ko$Ko <- gsub("ko:", "", gene2ko$Ko)

  cat("Saving GO database\n")
  genus <- "Custom genus"
  species <- "CUSTOM"
  tax_id <- "0000"
  makeOrgPackage(
    gene_info = gene_info,
    go = gene2go,
    ko = gene2ko,
    version = "0.1",
    maintainer = "maintainer <tmp@tmp>",
    author = "author <tmp@tmp>",
    outputDir = db_path,
    tax_id = tax_id,
    genus = genus,
    species = species,
    goTable = "go"
  )

  anno_db <- list(go = gene2go, ko = gene2ko)
  return(anno_db)
}

update_kegg_db <- function(db_path, kegg_json, kegg_db_file) {
  url <- "https://www.kegg.jp/kegg-bin/download_htext?htext=ko00001&format=json&filedir="
  if (!file.exists(kegg_json)) {
    json <- paste(db_path, "ko00001.json", sep = "/")
    download.file(url, json)
  } else {
    json <- kegg_json
  }
  pathway2name <- tibble(Pathway = character(), Name = character())
  ko2pathway <- tibble(Ko = character(), Pathway = character())
  kegg <- fromJSON(json)
  for (a in seq_along(kegg[["children"]][["children"]])) {
    A <- kegg[["children"]][["name"]][[a]]
    for (b in seq_along(kegg[["children"]][["children"]][[a]][["children"]])) {
      B <- kegg[["children"]][["children"]][[a]][["name"]][[b]]
      for (c in seq_along(kegg[["children"]][["children"]][[a]][["children"]][[b]][["children"]])) {
        pathway_info <- kegg[["children"]][["children"]][[a]][["children"]][[b]][["name"]][[c]]
        pathway_id <- str_match(pathway_info, "ko[0-9]{5}")[1]
        pathway_name <- str_replace(pathway_info, " \\[PATH:ko[0-9]{5}\\]", "") %>% str_replace("[0-9]{5} ", "")
        pathway2name <- rbind(pathway2name, tibble(Pathway = pathway_id, Name = pathway_name))
        kos_info <- kegg[["children"]][["children"]][[a]][["children"]][[b]][["children"]][[c]][["name"]]
        kos <- str_match(kos_info, "K[0-9]*")[, 1]
        ko2pathway <- rbind(ko2pathway, tibble(Ko = kos, Pathway = rep(pathway_id, length(kos))))
      }
    }
  }
  save(pathway2name, ko2pathway, file = kegg_db_file)
}

load_kegg_db <- function(kegg_db_file, gene2ko) {
  kegg_db <- new.env()
  load(file = kegg_db_file, envir = kegg_db)
  gene2pathway <- gene2ko %>%
    left_join(kegg_db$ko2pathway, by = "Ko", relationship = "many-to-many") %>%
    dplyr::select("GID", "Pathway") %>%
    na.omit()
  kegg_db <- list(pathway2name = kegg_db$pathway2name, gene2pathway = gene2pathway)
  return(kegg_db)
}

keep_plant_path_only <- function(pathway2name, gene2pathway, plant_pathways_total) {
  pathway2name <- pathway2name %>% filter(pathway2name$Pathway %in% plant_pathways_total$Pathway)
  gene2pathway <- gene2pathway %>% filter(gene2pathway$Pathway %in% pathway2name$Pathway)
  filtered_kegg_db <- list(pathway2name = pathway2name, gene2pathway = gene2pathway)
  return(filtered_kegg_db)
}

run_go_kegg <- function(db_path, gene2pathway, pathway2name, gene_file,
                        pvalue_cutoff, qvalue_cutoff, padjust_method, ontology) {
  gene <- read.table(gene_file, header = FALSE)
  gene_list <- gene[, 1]
  GO <- enrichGO(
    gene = gene_list,
    OrgDb = "org.CCUSTOM.eg.db",
    keyType = "GID",
    ont = "ALL",
    pAdjustMethod = padjust_method,
    pvalueCutoff = pvalue_cutoff,
    qvalueCutoff = qvalue_cutoff
  )

  cat("Saving GO result\n")
  GO_df <- as.data.frame(GO)
  write.table(GO_df, file = "GO.results.tsv", sep = "\t", quote = FALSE)

  cat("Saving GO barplot\n")
  pdf(file = "GO_barplot.pdf", width = 15, height = 20)
  print(barplot(GO, drop = TRUE, showCategory = 10, split = "ONTOLOGY") + facet_grid(ONTOLOGY ~ ., scale = "free"))
  dev.off()

  cat("Saving GO bubble\n")
  pdf(file = "GO_bubble.pdf", width = 15, height = 20)
  print(dotplot(GO, showCategory = 10, split = "ONTOLOGY") + facet_grid(ONTOLOGY ~ ., scale = "free"))
  dev.off()

  KEGG <- enricher(
    gene = gene_list,
    TERM2GENE = gene2pathway[c("Pathway", "GID")],
    TERM2NAME = pathway2name[c("Pathway", "Name")],
    pAdjustMethod = padjust_method,
    pvalueCutoff = pvalue_cutoff,
    qvalueCutoff = qvalue_cutoff
  )


  cat("Saving KEGG result\n")
  KEGG_df <- as.data.frame(KEGG)
  write.table(KEGG_df, file = "KEGG.results.tsv", sep = "\t", quote = FALSE)

  cat("Saving KEGG barplot\n")
  pdf(file = "KEGG_barplot.pdf", width = 15, height = 20)
  print(barplot(KEGG, drop = TRUE, showCategory = 10))
  dev.off()

  cat("Saving KEGG bubble\n")
  pdf(file = "KEGG_bubble.pdf", width = 15, height = 20)
  print(dotplot(KEGG))
  dev.off()
}

main_pipe <- function() {
  opt <- parse_args(opt_parser)
  opt_names <- names(opt)
  if ("input" %in% opt_names && "anno" %in% opt_names && "db" %in% opt_names) {
    anno_file <- opt$anno
    db_path <- opt$db
    if (!dir.exists(db_path)) {
      dir.create(db_path)
    }
    go_db_name <- paste("org.", substr(opt$genus, 1, 1), opt$species, ".eg.db", sep = "")
    db_pack <- paste(db_path, go_db_name, sep = "/")

    pvalue_cutoff <- opt$pvalue
    qvalue_cutoff <- opt$qvalue
    padjust_method <- opt$padjust
    ontology <- opt$ontology
    cat("Generating GO database\n")
    if (!file.exists(anno_file)) {
      cat("Fatal: annotation file not exists\n")
      return()
    }

    if (file.exists(db_pack)) {
      unlink(db_pack, recursive = TRUE)
    }
    anno_db <- generate_go_db(anno_file, db_path, opt$genus, opt$species, opt$tax_id)

    cat("Loading GO database\n")
    install.packages(db_pack, repos = NULL, type = "sources")
    do.call(library, list(go_db_name))
    kegg_json <- ""
    if ("kegg_json" %in% opt_names) {
      kegg_json <- opt$kegg_json
    }
    kegg_db_file <- paste(db_path, "KEGG_db.RData", sep = "/")
    if (!file.exists(kegg_db_file) || "update" %in% opt_names) {
      cat("Generating KEGG database\n")
      update_kegg_db(db_path, kegg_json, kegg_db_file)
    }

    cat("Loading KEGG database\n")
    kegg_db <- load_kegg_db(kegg_db_file, anno_db$ko)

    if ("plant" %in% opt_names) {
      cat("Keeping plants pathway\n")
      plants_kegg_file <- opt$plant_kegg
      if (!file.exists(plants_kegg_file)) {
        plants_kegg <- "plants.kegg.txt"
        plants_kegg_file <- paste(db_path, plants_kegg, sep = "/")
        download_plant_path(plants_kegg_file)
      }
      plant_pathways_total <- load_plant_path(plants_kegg_file)
      kegg_db <- keep_plant_path_only(kegg_db$pathway2name, kegg_db$gene2pathway, plant_pathways_total)
    }

    gene_file <- opt$input
    cat("Running GO and KEGG\n")
    run_go_kegg(
      db_path, kegg_db$gene2pathway, kegg_db$pathway2name, gene_file,
      pvalue_cutoff, qvalue_cutoff, padjust_method, ontology
    )

    cat("Finished\n")
  } else {
    print_help(opt_parser)
  }
}

opt_list <- list(
  make_option(c("-i", "--input"), type = "character", help = "Input gene list file"),
  make_option(c("-a", "--anno"), type = "character", help = "Functional annotation file"),
  make_option(c("-d", "--db"), type = "character", help = "Database path"),
  make_option(c("--kegg_json"), type = "character", help = "Pre-downloaded kegg json file"),
  make_option(c("--genus"),
    type = "character", help = "Genus name for creating GO database, default=\"Custom genus\"",
    default = "Custom genus"
  ),
  make_option(c("--pvalue"), type = "numeric", help = "P value cutoff for GO and KEGG, default=0.05", default = 0.05),
  make_option(c("--qvalue"), type = "numeric", help = "Q value cutoff for GO and KEGG, default=0.05", default = 0.05),
  make_option(c("--padjust"),
    type = "character", help = "P adjust method for GO and KEGG, default=\"BH\"",
    default = "BH"
  ),
  make_option(c("--ontology"),
    type = "character", help = "Ontology for GO, default=\"ALL\"",
    default = "ALL"
  ),
  make_option(c("--species"),
    type = "character", help = "Species name for creating GO database, default=\"CUSTOM\"",
    default = "CUSTOM"
  ),
  make_option(c("--tax_id"),
    type = "character", help = "Tax id for creating GO database, default=\"0000\"",
    default = "0000"
  ),
  make_option(c("--update"), action = "store_true", help = "Update databases"),
  make_option(c("--plant"), action = "store_true", help = "enrich with plant pathway only"),
  make_option(c("--plant_kegg"), type = "character", help = "Pre-generated plant kegg db file")
)

opt_parser <- OptionParser(option_list = opt_list, usage = "This Script is used for running GO and KEGG")
main_pipe()
