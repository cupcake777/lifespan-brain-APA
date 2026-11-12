# =============================================================================
# src/utils_report.R
# =============================================================================
# Interactive HTML Report Generation for EasyColoc Results
# =============================================================================
# Generates publication-ready, interactive HTML reports with:
# - Summary statistics dashboard
# - Interactive Manhattan/volcano plots
# - Gene ontology enrichment tables
# - Downloadable result tables
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(jsonlite)
  library(glue)
})

report_value <- function(row, primary, fallback = NULL, default = NA_character_) {
  candidates <- c(primary, fallback)
  candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
  candidates <- unique(candidates)

  for (col in candidates) {
    if (col %in% names(row)) return(row[[col]])
  }

  default
}

read_report_results <- function(results_dir) {
  merged_file <- file.path(results_dir, "all_colocalization_results.csv")

  if (file.exists(merged_file)) {
    return(fread(merged_file))
  }

  abf_files <- list.files(file.path(results_dir, "abf"),
                          pattern = "_results\\.csv$",
                          full.names = TRUE)
  merge_results(abf_files)
}

# =============================================================================
# generate_html_report: Main report generation function
# =============================================================================
generate_html_report <- function(results_dir,
                                  output_file = "coloc_report.html",
                                  project_name = "EasyColoc Analysis") {

  message("[Report] Generating interactive HTML report...")

  susie_files <- list.files(file.path(results_dir, "susie"),
                            pattern = "_susie\\.csv$",
                            full.names = TRUE)

  all_results <- read_report_results(results_dir)

  if (is.null(all_results) || nrow(all_results) == 0) {
    warning("[Report] No results found to report")
    return(FALSE)
  }

  # Generate report HTML
  html_content <- build_report_html(
    all_results = all_results,
    susie_files = susie_files,
    project_name = project_name,
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  )

  # Write HTML file
  writeLines(html_content, output_file)
  message(sprintf("[Report] Saved: %s", output_file))

  # Generate embedded JSON data for interactivity
  json_data <- list(
    summary = get_summary_stats(all_results),
    results = all_results,
    timestamp = Sys.time()
  )

  json_file <- gsub("\\.html$", "_data.json", output_file)
  write_json(json_data, json_file, pretty = TRUE, auto_unbox = TRUE)
  message(sprintf("[Report] Saved data: %s", json_file))

  return(TRUE)
}

# =============================================================================
# merge_results: Combine all locus results
# =============================================================================
merge_results <- function(abf_files) {

  if (length(abf_files) == 0) return(NULL)

  results_list <- lapply(abf_files, function(f) {
    tryCatch({
      dt <- fread(f)
      dt$source_file <- basename(f)
      dt
    }, error = function(e) {
      warning(sprintf("Failed to read %s: %s", basename(f), e$message))
      NULL
    })
  })

  valid_results <- Filter(Negate(is.null), results_list)

  if (length(valid_results) == 0) return(NULL)

  rbindlist(valid_results, fill = TRUE)
}

# =============================================================================
# get_summary_stats: Calculate summary statistics
# =============================================================================
get_summary_stats <- function(results) {

  pp4_vals <- results$PP4[!is.na(results$PP4)]
  gene_labels <- if ("Gene" %in% names(results)) {
    results$Gene
  } else if ("Phenotype" %in% names(results)) {
    results$Phenotype
  } else {
    character(0)
  }

  phenotype_locus_count <- if ("phenotype_locus_count" %in% names(results)) {
    results$phenotype_locus_count
  } else if ("lead_locus_count" %in% names(results)) {
    results$lead_locus_count
  } else {
    rep(NA_integer_, nrow(results))
  }

  list(
    total_tests = nrow(results),
    significant_pp4_08 = sum(results$PP4 >= 0.8, na.rm = TRUE),
    significant_pp4_07 = sum(results$PP4 >= 0.7, na.rm = TRUE),
    significant_pp4_05 = sum(results$PP4 >= 0.5, na.rm = TRUE),
    mean_pp4 = round(mean(pp4_vals, na.rm = TRUE), 4),
    median_pp4 = round(median(pp4_vals, na.rm = TRUE), 4),
    max_pp4 = round(max(pp4_vals, na.rm = TRUE), 4),
    min_pp4 = round(min(pp4_vals, na.rm = TRUE), 4),
    mean_n_snps = round(mean(results$n_snps, na.rm = TRUE), 1),
    unique_genes = length(unique(gene_labels)),
    unique_loci = length(unique(results$Locus)),
    multi_locus_phenotype_rows = sum(!is.na(phenotype_locus_count) & phenotype_locus_count > 1),
    max_phenotype_locus_count = if (all(is.na(phenotype_locus_count))) NA_integer_ else max(phenotype_locus_count, na.rm = TRUE)
  )
}

# =============================================================================
# build_report_html: Generate HTML content
# =============================================================================
build_report_html <- function(all_results, susie_files, project_name, timestamp) {

  stats <- get_summary_stats(all_results)

  # Top results table (by PP4)
  top_results <- all_results[order(-all_results$PP4), ][seq_len(min(50, nrow(all_results))), ]

  # PP4 values for JSON/plotting
  pp4_vals <- all_results$PP4[!is.na(all_results$PP4)]

  sprintf('
<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>%s - Colocalization Report</title>
  <script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
  <script src="https://cdn.datatables.net/1.13.7/js/jquery.dataTables.min.js"></script>
  <link rel="stylesheet" href="https://cdn.datatables.net/1.13.7/css/jquery.dataTables.min.css">
  <style>
    :root {
      --primary-color: #2c3e50;
      --secondary-color: #3498db;
      --success-color: #27ae60;
      --warning-color: #f39c12;
      --danger-color: #e74c3c;
      --light-bg: #ecf0f1;
    }

    * { box-sizing: border-box; margin: 0; padding: 0; }

    body {
      font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
      line-height: 1.6;
      color: #333;
      background: var(--light-bg);
    }

    header {
      background: linear-gradient(135deg, var(--primary-color), var(--secondary-color));
      color: white;
      padding: 2rem;
      text-align: center;
    }

    header h1 { margin-bottom: 0.5rem; }
    header p { opacity: 0.9; }

    .container {
      max-width: 1400px;
      margin: 0 auto;
      padding: 2rem;
    }

    .summary-cards {
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
      gap: 1rem;
      margin: 2rem 0;
    }

    .card {
      background: white;
      border-radius: 8px;
      padding: 1.5rem;
      box-shadow: 0 2px 10px rgba(0,0,0,0.1);
    }

    .card h3 {
      color: var(--primary-color);
      font-size: 0.9rem;
      margin-bottom: 0.5rem;
    }

    .card .value {
      font-size: 2rem;
      font-weight: bold;
      color: var(--secondary-color);
    }

    .card.highlight {
      background: linear-gradient(135deg, var(--secondary-color), var(--primary-color));
      color: white;
    }

    .card.highlight h3, .card.highlight .value {
      color: white;
    }

    .section {
      background: white;
      border-radius: 8px;
      padding: 1.5rem;
      margin: 2rem 0;
      box-shadow: 0 2px 10px rgba(0,0,0,0.1);
    }

    .section h2 {
      color: var(--primary-color);
      margin-bottom: 1rem;
      padding-bottom: 0.5rem;
      border-bottom: 2px solid var(--secondary-color);
    }

    .pp4-badge {
      display: inline-block;
      padding: 0.25rem 0.75rem;
      border-radius: 20px;
      font-weight: bold;
      font-size: 0.85rem;
    }

    .pp4-strong { background: #d4edda; color: #155724; }
    .pp4-suggestive { background: #fff3cd; color: #856404; }
    .pp4-weak { background: #e2e3e5; color: #383d44; }

    table.dataTable {
      width: 100%% !important;
      border-collapse: collapse;
    }

    table.dataTable thead th {
      background: var(--primary-color);
      color: white;
      padding: 0.75rem;
    }

    table.dataTable tbody tr:hover {
      background: var(--light-bg);
    }

    table.dataTable tbody tr:nth-child(even) {
      background: #f8f9fa;
    }

    .plot-container {
      width: 100%%;
      height: 600px;
    }

    footer {
      text-align: center;
      padding: 2rem;
      color: #666;
      font-size: 0.9rem;
    }

    .download-btn {
      display: inline-block;
      background: var(--secondary-color);
      color: white;
      padding: 0.75rem 1.5rem;
      border-radius: 5px;
      text-decoration: none;
      margin: 1rem 0;
    }

    .download-btn:hover {
      background: var(--primary-color);
    }
  </style>
</head>
<body>
  <header>
    <h1>%s</h1>
    <p>Colocalization Analysis Report</p>
    <p>Generated: %s</p>
  </header>

  <div class="container">

    <!-- Summary Cards -->
    <div class="summary-cards">
      <div class="card highlight">
        <h3>Total Tests</h3>
        <div class="value">%d</div>
      </div>
      <div class="card">
        <h3>Strong (PP4 ≥ 0.8)</h3>
        <div class="value">%d</div>
      </div>
      <div class="card">
        <h3>Suggestive (0.5-0.8)</h3>
        <div class="value">%d</div>
      </div>
      <div class="card">
        <h3>Mean PP4</h3>
        <div class="value">%.3f</div>
      </div>
      <div class="card">
        <h3>Max PP4</h3>
        <div class="value">%.3f</div>
      </div>
      <div class="card">
        <h3>Unique Genes</h3>
        <div class="value">%d</div>
      </div>
    </div>

    <!-- PP4 Distribution Plot -->
    <div class="section">
      <h2>PP4 Distribution</h2>
      <div id="pp4-plot" class="plot-container"></div>
    </div>

    <!-- Results Table -->
    <div class="section">
      <h2>All Results</h2>
      <a href="#" class="download-btn" onclick="downloadTable()">Download CSV</a>
      <table id="results-table" class="display">
        <thead>
          <tr>
            <th>Locus</th>
            <th>Gene</th>
            <th>PP4</th>
            <th>PP3</th>
            <th>PP.H0</th>
            <th>n_snps</th>
            <th>Phenotype Loci</th>
            <th>Evidence</th>
          </tr>
        </thead>
        <tbody>
          %s
        </tbody>
      </table>
    </div>

    <!-- SuSiE Results -->
    <div class="section">
      <h2>SuSiE Fine-Mapping Results</h2>
      <p>SuSiE analysis was performed for loci with PP4 ≥ %.2f</p>
      %s
    </div>

  </div>

  <footer>
    <p>Generated by EasyColoc v1.2</p>
    <p>Powered by R, coloc, and SuSiE</p>
  </footer>

  <script>
    // Initialize DataTable
    $(document).ready(function() {
      $("#results-table").DataTable({
        pageLength: 25,
        order: [[2, "desc"]],  // Sort by PP4 descending
        responsive: true
      });
    });

    // PP4 Distribution Plot
    var pp4Data = %s;

    var trace1 = {
      x: pp4Data,
      type: "histogram",
      marker: {
        color: "rgba(52, 152, 219, 0.7)"
      }
    };

    var layout = {
      title: "Distribution of Colocalization Probabilities (PP4)",
      xaxis: {
        title: "PP4",
        range: [0, 1]
      },
      yaxis: {
        title: "Count"
      },
      shapes: [
        {
          type: "line",
          x0: 0.8, x1: 0.8, y0: 0, y1: 1,
          line: { color: "green", dash: "dash" },
          annotation: { text: "Strong", y: 0.95, x: 0.82 }
        },
        {
          type: "line",
          x0: 0.5, x1: 0.5, y0: 0, y1: 1,
          line: { color: "orange", dash: "dash" },
          annotation: { text: "Suggestive", y: 0.85, x: 0.52 }
        }
      ]
    };

    Plotly.newPlot("pp4-plot", [trace1], layout, {responsive: true});

    function downloadTable() {
      var table = document.getElementById("results-table");
      var csv = [];
      for (var i = 0; i < table.rows.length; i++) {
        var row = table.rows[i];
        var cells = [];
        for (var j = 0; j < row.cells.length; j++) {
          cells.push(row.cells[j].innerText);
        }
        csv.push(cells.join(","));
      }
      var blob = new Blob([csv.join("\\n")], { type: "text/csv" });
      var url = URL.createObjectURL(blob);
      var a = document.createElement("a");
      a.href = url;
      a.download = "coloc_results.csv";
      a.click();
    }
  </script>
</body>
</html>
',
    project_name,  # Title
    project_name,  # Header h1
    timestamp,
    stats$total_tests,
    stats$significant_pp4_08,
    stats$significant_pp4_05 - stats$significant_pp4_08,  # Suggestive only
    stats$mean_pp4,
    stats$max_pp4,
    stats$unique_genes,
    # Table rows
    paste(apply(top_results, 1, function(row) {
      row_list <- as.list(row)
      pp4 <- as.numeric(report_value(row_list, "PP4", default = NA_real_))
      pp3 <- suppressWarnings(as.numeric(report_value(row_list, "PP3", default = NA_real_)))
      pp0 <- suppressWarnings(as.numeric(report_value(row_list, "PP.H0.abf", default = NA_real_)))
      n_snps <- suppressWarnings(as.numeric(report_value(row_list, "n_snps", default = NA_real_)))
      phenotype_locus_count <- suppressWarnings(as.integer(report_value(row_list, "phenotype_locus_count", fallback = "lead_locus_count", default = NA_integer_)))
      locus <- report_value(row_list, "Locus", default = "")
      gene_label <- report_value(row_list, "Gene", fallback = "Phenotype", default = "")
      badge_class <- if(pp4 >= 0.8) "pp4-strong" else if(pp4 >= 0.5) "pp4-suggestive" else "pp4-weak"
      sprintf('
          <tr>
            <td>%s</td>
            <td>%s</td>
            <td><span class="pp4-badge %s">%.3f</span></td>
            <td>%.3f</td>
            <td>%.3f</td>
            <td>%d</td>
            <td>%s</td>
            <td>%s</td>
          </tr>',
        locus,
        gene_label,
        badge_class,
        pp4,
        ifelse(is.na(pp3), NaN, pp3),
        ifelse(is.na(pp0), NaN, pp0),
        ifelse(is.na(n_snps), 0L, as.integer(round(n_snps))),
        ifelse(is.na(phenotype_locus_count), "", as.character(phenotype_locus_count)),
        ifelse(pp4 >= 0.8, "Strong", ifelse(pp4 >= 0.5, "Suggestive", "Weak"))
      )
    }), collapse = "\n"),
    # SuSiE threshold
    0.3,  # Default threshold
    # SuSiE summary
    if(length(susie_files) > 0) {
      sprintf("<p>Found %d SuSiE result files</p>", length(susie_files))
    } else {
      "<p>No SuSiE results available</p>"
    },
    # JSON data for plot
    toJSON(pp4_vals, auto_unbox = TRUE)
  )
}

# =============================================================================
# download_json_data: Helper for interactive features
# =============================================================================
download_json_data <- function(results, output_file) {
  jsonlite::write_json(results, output_file, pretty = TRUE, auto_unbox = TRUE)
  message(sprintf("[Report] Saved JSON data: %s", output_file))
}
