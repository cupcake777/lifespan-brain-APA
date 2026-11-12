suppressPackageStartupMessages({
  library(glue)
})

easycoloc_command_exists <- function(command) {
  nzchar(Sys.which(command))
}

easycoloc_require_command <- function(command, install_hint = NULL, fatal = TRUE) {
  command_path <- Sys.which(command)
  if (nzchar(command_path)) {
    return(command_path)
  }

  message_text <- if (!is.null(install_hint) && nzchar(install_hint)) {
    glue("Required command '{command}' not found in PATH. Install hint: {install_hint}")
  } else {
    glue("Required command '{command}' not found in PATH.")
  }

  if (isTRUE(fatal)) {
    stop(message_text, call. = FALSE)
  }
  warning(message_text, call. = FALSE)
  ""
}

easycoloc_pick_downloader <- function() {
  if (easycoloc_command_exists("wget")) {
    return("wget")
  }
  if (easycoloc_command_exists("curl")) {
    return("curl")
  }
  stop(
    "Neither 'wget' nor 'curl' is available. Install one of them to enable resumable downloads.",
    call. = FALSE
  )
}

easycoloc_download_file <- function(url, dest, force = FALSE, quiet = FALSE) {
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)

  if (!isTRUE(force) && file.exists(dest) && file.info(dest)$size > 0) {
    return("exists")
  }

  if (startsWith(url, "file://")) {
    src <- sub("^file://", "", url)
    if (!file.exists(src)) {
      stop("Local file URL does not exist: ", url, call. = FALSE)
    }
    ok <- file.copy(src, dest, overwrite = TRUE)
    if (!isTRUE(ok)) {
      stop("Failed to copy local URL source to ", dest, call. = FALSE)
    }
    return("downloaded")
  }

  if (!grepl("^(https?|ftp)://", url)) {
    if (!file.exists(url)) {
      stop("Source file does not exist: ", url, call. = FALSE)
    }
    ok <- file.copy(url, dest, overwrite = TRUE)
    if (!isTRUE(ok)) {
      stop("Failed to copy local source to ", dest, call. = FALSE)
    }
    return("downloaded")
  }

  downloader <- easycoloc_pick_downloader()
  status <- if (identical(downloader, "wget")) {
    args <- c("-c", "-O", dest, url)
    if (isTRUE(quiet)) {
      args <- c("-q", args)
    }
    system2("wget", args)
  } else {
    args <- c("-L", "-C", "-", "-o", dest, url)
    if (isTRUE(quiet)) {
      args <- c("-sS", args)
    }
    system2("curl", args)
  }

  if (!identical(status, 0L) || !file.exists(dest) || file.info(dest)$size <= 0) {
    stop(glue("Download failed for {url} -> {dest}"), call. = FALSE)
  }
  "downloaded"
}

easycoloc_run_command <- function(command, args, fail_message = NULL, stdout = "", stderr = "") {
  status <- system2(command, args = args, stdout = stdout, stderr = stderr)
  if (!identical(status, 0L)) {
    if (is.null(fail_message)) {
      fail_message <- glue("Command failed: {command} {paste(args, collapse = ' ')}")
    }
    stop(fail_message, call. = FALSE)
  }
  invisible(status)
}

easycoloc_parse_chromosomes <- function(chromosomes) {
  if (length(chromosomes) == 0 || is.null(chromosomes) || is.na(chromosomes)) {
    return(as.character(1:22))
  }
  if (length(chromosomes) > 1) {
    return(as.character(chromosomes))
  }
  chr_str <- gsub("\\s+", "", as.character(chromosomes))
  parts <- unlist(strsplit(chr_str, ",", fixed = TRUE))
  out <- character()
  for (part in parts) {
    if (grepl("^[0-9]+-[0-9]+$", part)) {
      bounds <- as.integer(strsplit(part, "-", fixed = TRUE)[[1]])
      out <- c(out, as.character(seq(bounds[1], bounds[2])))
    } else {
      out <- c(out, part)
    }
  }
  unique(out)
}
