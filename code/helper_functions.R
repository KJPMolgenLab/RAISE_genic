#' load & install packages
#' ...


#' Get a glimpse of your data, including variable labels.
#'
#' This is like a transposed version of print: columns run down the page,
#' and data runs across. This makes it possible to see every column in
#' a data frame. It's a little like \code{str()} applied to a data frame
#' but it tries to show you as much data as possible. (And it always shows
#' the underlying data, even when applied to a remote data source.)
#'
#' This is a modified version of the tibble package's 'glimpse.tbl' S3 method
#' that includes the contents of the 'label' variable attribute, if present
#' (e.g., if you imported a Stata dataset with variable labels using
#' \code{haven::read_stata()}).
#'
#' @param x An object to glimpse at.
#' @param width Width of output: defaults to the setting of the option
#'   `tibble.width` (if finite) or the width of the console.
#' @param labels Include variable labels in the printout. If `TRUE` and no
#'   variables have the `label` attribute, labels will be printed as 10
#'   spaces to demonstrate their absence.
#' @return x original x is (invisibly) returned, allowing \code{glimpse_labels()}
#'   to be used within a data pipe line.
#' @importFrom tibble type_sum
#' @importFrom stringr str_wrap
#' @importFrom utils combn head
#' @export
#' @examples
#'  dat1 <- tibble::data_frame(x = 1:20, y = 21:40)
#'  attributes(dat1[["x"]])$label <- "this is my label"
#'  glimpse_labels(dat1)
#'
#'  dat2 <- dat1
#'  attributes(dat2[["x"]])$label <- paste0(rep("my label", 12), collapse = " ")
#'  glimpse_labels(dat2)
#'
#'  dat3 <- dat1
#'  attributes(dat3[["x"]])$label <- paste0(rep("my label", 15), collapse = " ")
#'  glimpse_labels(dat3)
glimpse_labels <- function (x, width = NULL, labels = TRUE) {
  # Same as glimpse
  width <- tibble:::tibble_glimpse_width(width)
  stopifnot(is.finite(width))
  cat("Observations: ", tibble:::big_mark(nrow(x)), "\n", sep = "")
  if (ncol(x) == 0)
    return(invisible())
  cat("Variables: ", tibble:::big_mark(ncol(x)), "\n", sep = "")
  rows <- as.integer(width/3)
  df <- as.data.frame(head(x, rows))
  var_types <- vapply(df, tibble::type_sum, character(1))

  # Edits to add labels
  if (!labels) {
    var_names <- paste0("$ ", format(names(df)), " <", var_types,
                        "> ")
  } else {
    if (!("label" %in% names(attributes(x[[1]])))) {
      tbl_labels <- rep("", length(names(x)))
    } else {
      tbl_labels <- vapply(x, function(x) {
        ifelse(is.null(attributes(x)[["label"]]), paste0(rep(" ", 10), collapse = ""),
               attributes(x)[["label"]])
      }, character(1))
    }

    tbl_labels_fmt <- format(tbl_labels)

    # max label width is longest the label can be while leaving room for data ellipses on one line
    names_width <- max(nchar(names(df)))
    types_width <- max(nchar(var_types))
    max_label_width <- width - names_width - types_width - 7 - 6 # 7 chars for $ and <>
    # 6 chars for shortest data allowed

    # wrap long lines to the leftmost label position
    max_line_width <- width - 3 # glimpse() leaves 2 chars off by default, +1 to force wrap
    left_label_pos <- names_width + 3 # 3 chars for $ and 2 spaces
    label_wrap_width <- max_line_width - left_label_pos

    tbl_labels <- ifelse(nchar(tbl_labels) > max_label_width,
                         stringr::str_wrap(tbl_labels, width = label_wrap_width,
                                           exdent = left_label_pos),
                         tbl_labels)

    # TO DO: deal with labels that don't wrap across lines because of single words that are too long
    #wrap_prep <- (nchar(tbl_labels) > max_label_width) & !grepl(" ", tbl_labels)

    # length of longest label that isn't longer than the max width
    tbl_labels_width_nowrap <- max(nchar(tbl_labels)[nchar(tbl_labels) < max_label_width])

    # if last wrapped line is short, pad end so var types line up
    pad_to <- names_width + tbl_labels_width_nowrap + 4

    last_line_width <- function(tbl_labels) {
      newline_chars <- vapply(gregexpr("\n", tbl_labels), function(x) {
        out <- x[length(x)]
        out
      }, numeric(1))
      newline_chars <- ifelse(newline_chars == -1, NA, newline_chars)

      last_line_width <- ifelse(!is.na(newline_chars),
                                nchar(substr(tbl_labels, newline_chars, nchar(tbl_labels))),
                                NA)
      last_line_width
    }

    last_line_widths <- last_line_width(tbl_labels)

    to_pad <- which(!is.na(last_line_widths) & last_line_widths < pad_to)

    if (length(to_pad) > 0) {
      for (i in seq_along(to_pad)) {
        tbl_labels[to_pad[i]] <- paste0(tbl_labels[to_pad[i]],
                                        paste0(rep(" ", pad_to - last_line_widths[to_pad[i]]),
                                               collapse = ""))
      }

      last_line_widths <- last_line_width(tbl_labels)
    }

    tbl_labels <- ifelse(nchar(tbl_labels) > max_label_width,
                         tbl_labels,
                         substr(tbl_labels_fmt, 1, tbl_labels_width_nowrap))

    var_names <- paste0("$ ", format(names(df)), " ",
                        tbl_labels, " <", var_types, "> ")
  }

  data_width <- width - nchar(var_names) - 2
  data_width <- ifelse(data_width < 0, width - last_line_widths - 8, data_width)

  formatted <- vapply(df, function(x) paste0(tibble:::format_n(x), collapse = ", "),
                      character(1), USE.NAMES = FALSE)

  truncated <- stringr::str_trunc(formatted, data_width)
  cat(paste0(var_names, truncated, collapse = "\n"), "\n",
      sep = "")
  invisible(x)
}


## package installation & loading ##############################################

# libraries
load_inst_pkgs <- function(..., silent = FALSE){
  pkgs <- list(...) # arguments must be of type char

  # load/install unixtools if in pkgs, needed for set.tempdir()
  if ("unixtools" %in% pkgs){
    if (!require("unixtools")){
      load_inst_pkgs("remotes") # for install_github()
      install_github("s-u/unixtools")
      library("unixtools")
    }
    pkgs <- pkgs[pkgs != "unixtools"]
  }

  for (pkg in pkgs){
    if (!require(pkg, character.only = TRUE)){
      install.packages(pkg,
                       dependencies = TRUE,
                       repos = "https://ftp.fau.de/cran/"
      )
      library(pkg, character.only = TRUE)
    }
  }
  if (!silent) print("Done: Load libraries")
}

# Bioconductor packages
load_inst_bioc.pkgs <- function(...){ # arguments must be of type char
  bioc.pkgs <- list(...)

  if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")
  }

  for (bioc.pkg in bioc.pkgs){
    if (!require(bioc.pkg, character.only = TRUE)){
      BiocManager::install(bioc.pkg, lib = .libPaths()[1])
      library(bioc.pkg, character.only = TRUE)
    }
  }
  print("Done: Load Bioconductor packages")
}

