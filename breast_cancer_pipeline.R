# Modulo 1: Setup e librerie.
# Funzione che accetta una lista di dipendenze da installare.
 # Se il pacchetto NON è installato, lo installa usando BiocManager.
 # 'update = FALSE': evita di aggiornare le dipendenze già installate.
setup_libraries <- function(required_packages) {
  install_if_missing <- function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      BiocManager::install(pkg, update = FALSE, ask = FALSE)
    }
  }
  invisible(lapply(required_packages, install_if_missing))
  lapply(required_packages, library, character.only = TRUE)
}

# Modulo 2: Scaricamento e preparazione dei dati.
# Funzione che scarica i dati da un progetto e li prepara per l'elaborazione successiva.
# Crea una query per il database GDC usando i parametri forniti nel metodo:
  # - 'project': nome del progetto (es. "TCGA-BRCA").
  # - 'category': categoria di dati (es. "Transcriptome Profiling").
  # - 'data_type': tipo di dati (es. "Gene Expression Quantification").
  # - 'workflow': workflow specifico (es. "STAR - Counts").
# Scarica i dati corrispondenti alla query, dal database GDC.
# Output: Una lista contenente:
  # expr_data: matrice di dati espressione genica.
  # metadata: informazioni aggiuntive sui campioni.
download_and_prepare_data <- function(project, category, data_type, workflow) {
  print("Scaricamento dati dal progetto...")
  query <- GDCquery(
    project = project,
    data.category = category,
    data.type = data_type,
    workflow.type = workflow
  )
  GDCdownload(query)
  brca_data <- GDCprepare(query)
  list(expr_data = assay(brca_data), metadata = colData(brca_data))
}

# Modulo 3: Filtraggio dei campioni e dei geni.
# Funzione per filtrare i dati di espressione genica e i metadati in base a dei criteri.
#Filtraggio:
  # Geni: Mantiene solo quelli con espressione media maggiore di una soglia.
  # Campioni: Seleziona un numero equilibrato di campioni per ogni condizione (opzionale).
  # Variabilità genica: Mantiene i geni con maggiore varianza (opzionale).
# Output: Una lista con i dati di espressione filtrati e i metadati aggiornati.
filter_data <- function(expr_data, metadata, filter_expression = 10, num_samples_per_condition = NULL, num_genes = NULL) {
  print("Filtraggio dei geni e dei campioni...")
  
  expr_data <- expr_data[rowMeans(expr_data) > filter_expression, ]
  metadata$condition <- ifelse(metadata$definition == "Primary solid Tumor", "Tumor", "Normal")
  
  if (!is.null(num_samples_per_condition)) {
    print(paste("Campionamento di", num_samples_per_condition, "campioni per condizione"))
    sampled_indices <- unlist(lapply(unique(metadata$condition), function(cond) {
      cond_indices <- which(metadata$condition == cond)
      sample(cond_indices, size = min(num_samples_per_condition, length(cond_indices)))
    }))
    expr_data <- expr_data[, sampled_indices]
    metadata <- metadata[sampled_indices, ]
  }
  
  if (!is.null(num_genes)) {
    print(paste("Campionamento dei top", num_genes, "geni per variabilità"))
    gene_variances <- apply(expr_data, 1, var)
    top_genes <- names(sort(gene_variances, decreasing = TRUE))[1:num_genes]
    expr_data <- expr_data[top_genes, ]
  }
  
  list(expr_data = expr_data, metadata = metadata)
}

# Modulo 4: Suddivisione in blocchi
# Funzione che permette di suddividere in blocchi il dataset e utilizzare l'elaborazione parallela
# Per ciascun blocco:
  # Estrae i dati relativi ai geni 
  # Esegue l'analisi differenziale utilizzando DESeq2
  # Restituisce i risultati come dataframe
# Una lista di dataframe, ciascuno contenente i risultati dell'analisi differenziale per un blocco di dati.
# Sono combinati successivamente.
process_in_chunks <- function(expr_data, metadata, chunk_size, design_formula, contrast) {
  num_chunks <- ceiling(nrow(expr_data) / chunk_size)

  process_chunk <- function(i) {
    start <- (i - 1) * chunk_size + 1
    end <- min(i * chunk_size, nrow(expr_data))
    subset_genes <- rownames(expr_data)[start:end]
    
    expr_data_chunk <- expr_data[subset_genes, ]
    print(paste("Processamento blocco", i, "di", num_chunks, ":", start, "-", end, "geni"))
  
    dds_chunk <- DESeqDataSetFromMatrix(
      countData = expr_data_chunk,
      colData = metadata,
      design = design_formula
    )
  
    dds_chunk <- DESeq(dds_chunk)
    results_chunk <- results(dds_chunk, contrast = contrast)
    as.data.frame(results_chunk)
  }
  
  print("Esecuzione dell'analisi in parallelo...")
  mclapply(seq_len(num_chunks), process_chunk, mc.cores = detectCores() - 1)
}

# Modulo 5: Analisi quantitativa
# Funzione che analizza l'espressione genica media per ciascuna condizione (Tumor e Normal)
# e indentifica i geni più differenziati tra le due condizioni.
# Passaggi:
 # Calcola la media dell'espressione genica per ciascun gene e condizione.
 # Determina la differenza tra le medie delle condizioni "Tumor" e "Normal".
 # Ordina i geni in base alla differenza di espressione in ordine decrescente.
# Output:
#Un dataframe contenente:
 # Gene: i nomi dei geni.
 # Tumor_Mean: media dell'espressione nei campioni "Tumor".
 # Normal_Mean: media dell'espressione nei campioni "Normal".
 # Difference: differenza tra le due medie.
# Stampa i 10 geni con la maggiore differenza positiva di espressione.
analyze_gene_expression_by_condition <- function(expr_data, metadata) {

  print("Analisi dell'espressione genica per condizione...")
  
  condition_means <- data.frame(
    Gene = rownames(expr_data),
    Tumor_Mean = rowMeans(expr_data[, metadata$condition == "Tumor"]),
    Normal_Mean = rowMeans(expr_data[, metadata$condition == "Normal"])
  )
  
  condition_means$Difference <- condition_means$Tumor_Mean - condition_means$Normal_Mean
  condition_means <- condition_means[order(condition_means$Difference, decreasing = TRUE), ]
  print("Top 10 geni con espressione più alta nei campioni Tumor:")
  print(head(condition_means, 10))
  
  return(condition_means)
}


# Modulo 6: Visualizzazioni
# Funzione che genera:
 # Volcano Plot: per identificare i geni significativi in base a p-value e log2FoldChange.
 # Heatmap: per visualizzare i pattern di espressione genica dei geni più significativi.
generate_visualizations <- function(results_combined, expr_data, metadata) {
  print("Creazione del volcano plot...")

  required_columns <- c("log2FoldChange", "pvalue", "padj")
  if (!all(required_columns %in% colnames(results_combined))) {
    stop("Mancano colonne necessarie in results_combined: ",
         paste(setdiff(required_columns, colnames(results_combined)), collapse = ", "))
  }
  
  results_combined <- results_combined[!is.na(results_combined$log2FoldChange) & 
                                         !is.na(results_combined$pvalue) & 
                                         !is.na(results_combined$padj), ]
  
  
  print(paste("Numero di righe valide:", nrow(results_combined)))
  
  results_combined$Significance <- "Not Significant"
  results_combined$Significance[results_combined$padj < 0.05] <- "Significant"
  results_combined$Significance[results_combined$padj < 0.01 & abs(results_combined$log2FoldChange) > 2] <- "Highly Significant"
  
  if (nrow(results_combined) > 0) {
    dev.new()
   
    volcano <- ggplot(results_combined, aes(x = log2FoldChange, y = -log10(pvalue))) +
      geom_point(aes(color = Significance), alpha = 0.5) +
      scale_color_manual(values = c("Not Significant" = "gray", 
                                    "Significant" = "blue", 
                                    "Highly Significant" = "red")) +
      theme_minimal() +
      xlab("Log2 Fold Change") +
      ylab("-Log10 P-value") +
      geom_hline(yintercept = -log10(0.05), col = "blue", linetype = "dashed") +
      ggtitle("Volcano Plot")
    
    
    print(volcano)
  } else {
    warning("Non ci sono dati validi per generare il volcano plot.")
  }
  
  print("Creazione della heatmap...")
  
  annotation_col <- data.frame(
    Condition = metadata$condition
  )
  rownames(annotation_col) <- colnames(expr_data)
  top_genes <- head(results_combined[order(results_combined$padj, na.last = NA), ], 30)
  
  
  if (nrow(top_genes) > 0) {
    dev.new()
    
    pheatmap(
      expr_data[rownames(top_genes), ],
      scale = "row",
      annotation_col = annotation_col,
      show_rownames = TRUE,
      show_colnames = TRUE
    )
    
  } else {
    warning("Non ci sono abbastanza geni significativi per generare una heatmap.")
  }
}



# Modulo main
main <- function() {
 
  # Step 1: Setup
  # Specifica i pacchetti richiesti per la pipeline.
  print("Step 1 avviato...")
  required_packages <- c(
    "TCGAbiolinks", "DESeq2", "pheatmap", "ggplot2", 
    "clusterProfiler", "org.Hs.eg.db", "parallel"
  )
  # Esegui la la funzione `setup_libraries` per installare e caricare i pacchetti richiesti.
  setup_libraries(required_packages)
  print("Step 1 concluso.")

  
  # Step 2: Scarica e prepara i dati
  # Scarica i dati dal progetto TCGA-BRCA, categoria "Transcriptome Profiling",
  # con tipo di dati "Gene Expression Quantification" e workflow "STAR - Counts".
  print("Step 2 avviato...")
  data <- download_and_prepare_data(
    project = "TCGA-BRCA",
    category = "Transcriptome Profiling",
    data_type = "Gene Expression Quantification",
    workflow = "STAR - Counts"
  )
  print("Step 2 concluso.")

  # Step 3: Filtra i dati
  # Filtra i dati:
  # - Mantieni i geni con espressione media > 10.
  # - Campiona 100 campioni per condizione (Tumor e Normal).
  # - Seleziona i 5000 geni con maggiore variabilità.
  print("Step 3 avviato...")
  filtered_data <- filter_data(
    expr_data = data$expr_data,
    metadata = data$metadata,
    filter_expression = 10,
    num_samples_per_condition = 100, 
    num_genes = 5000
  )
  print("Step 3 concluso")

  
  # Step 4: Analizza in blocchi
  # Esegue l'analisi differenziale in blocchi di 500 geni:
  # - Usa il design sperimentale definito da `~ condition`.
  # - Confronta le condizioni "Tumor" e "Normal".
  print("Step 4 avviato...")
  results_all <- process_in_chunks(
    expr_data = filtered_data$expr_data,
    metadata = filtered_data$metadata,
    chunk_size = 500,
    design_formula = ~ condition,
    contrast = c("condition", "Tumor", "Normal")
  )
  print("Step 4 concluso.")

  
  # Step 5: Combina i risultati
  # Combina i risultati di tutti i blocchi in un unico dataframe.
  print("Step 5 avviato...")
  results_combined <- do.call(rbind, results_all)
  # Crea una directory chiamata "results" per salvare i risultati.
  dir.create("results", showWarnings = FALSE, recursive = TRUE)
  # Esporta i risultati combinati in un file CSV.
  write.csv(results_combined, "results/differential_expression_results_combined.csv")
  print("Step 5 concluso")
 
  
  # Step 6: Visualizzazioni
  # Genera le visualizzazioni:
  # - Volcano plot per i risultati combinati.
  # - Heatmap per i top geni più significativi.
  print("Step 6 avviato...")
  generate_visualizations(results_combined, filtered_data$expr_data, filtered_data$metadata)
  print("Step 6 concluso")
  
  # Step 7: Analisi dell'espressione genica per condizione
  # Calcola le medie di espressione per ciascun gene nelle condizioni Tumor e Normal.
  print("Step 7 avviato...")
  condition_means <- analyze_gene_expression_by_condition(
    expr_data = filtered_data$expr_data,
    metadata = filtered_data$metadata
  )
  # Esporta i risultati dell'analisi per condizione in un file CSV.
  write.csv(condition_means, "results/gene_expression_by_condition.csv")
  print("Step 7 concluso")

  print("Analisi completata.")
}

main()
