# ______________________________________________________________________________
#
# CIRCOS (circos.R)
#
# Collection of secondary functions to get Circos plots from genomic data.
# Created by José R. Hernández (joseramon.hernandez@crg.eu) on 2020/06/16.
# This script is licensed under a GPLv3 (GPL-3) license.
#
# Last modifications on 2020/06/22.
# ______________________________________________________________________________

#LOR ANALYSIS
chromatin_circos("LOR_CLL_vs_normal.csv",
                 run_name = NULL,
                 na_imputation = FALSE,
                 normalization = NULL,
                 smoothing = 0.12,
                 ladder_values = c(0.5,0, -0.5),
                 inner_closest = TRUE,
                 chr_colors = rep("blue",20)) 

#You can change the distances argument to either LOR_CLL_vs_normal.csv" or "LOR_MCL_vs_normal.csv"

# CHROMATIN CIRCOS -------------------------------------------------------------
chromatin_circos <- function(distances,
                             run_name = NULL,
                             na_imputation = FALSE,
                             normalization = NULL,
                             smoothing = NULL,
                             ladder_values = NULL,
                             inner_closest = FALSE,
                             chr_colors = NULL) {
  #
  # ... generates Circos plots from a supplied table of 3D chromatin distances
  #   (distances to a single genomic region obtained by 4C, 5C or Hi-C methods).
  #
  # Args.:
  #   distances (dataframe or chr): table or pathway to table with distances
  #     (one of the column names should be 'chrom' or 'chr'; coverage or GC
  #     content data can be optionally added in 'cov' and 'gc' columns; the rest
  #     of columns -one column for each experiment- should contain distances).
  #   run_name (chr): name for the new folder where plots will be saved (this
  #     folder is created inside the distances file parent folder, when a file
  #     is supplied, or inside current folder; plots are not saved if 'NULL').
  #   na_imputation (logi): permission to carry out imputation of missing data
  #     (i.e.: NA's replacement with values in the middle of supplied ones).
  #   normalization (chr vector): list of methods to correct skeweness of the
  #     distance values (i.e.: "none", "squaroot", "logarith", "logistic",
  #     "atangent", "gompertz" and "hyperbol").
  #   smoothing (int or double): number of points for smoothing windows (smaller
  #     or equal than 1 to define proportion of points for each section, higher
  #     than 1 when defining the exact number of points for every section).
  #   ladder_values (int or double vector): number of reference points or list
  #     of reference values to be present in plots (e.g.: 'c(0.5, 1, 2, 4, 8)').
  #   inner_closest (logi): permission to invert representation of distances,
  #     plotting closest peaks in the inner area of chromosome sections.
  #   chr_colors (chr vector): colors to be used for each section/chromosome.
  #
  # Returns (recordedplot list) plots in a list of 'recordPlot()' objects.


  # Loading needed packages:
  packages <- c("circlize", "sigmoid")
  invisible(
    suppressPackageStartupMessages(
      sapply(packages, function(x) {
        install.packages(x[!c(require(x, character.only = TRUE))])
        library(x, character.only = TRUE)
      })
    )
  )

  # Loading data:
  opt <- 1 + is.character(distances)
  setwd(c(getwd(), dirname(dirname(distances)))[opt])
  distances <- list(distances, read.csv(distances))[[opt]]
  colnames(distances)[colnames(distances) == "chrom"] <- "chr"

  # Folder to save results:
  opt <- 1 + !is.null(run_name)
  now <- format(Sys.time(), "%Y.%m.%d.%H.%M.%OS")
  folder <- c(getwd(), paste0(run_name, "_", now))[opt]
  dir.create(
    path = folder,
    showWarnings = FALSE,
    recursive = TRUE
  )

  # Transforming genome coordinates:
  genome <- sapply(unique(distances$chr), function(chr) {
    min_val <- min(as.integer(rownames(distances)[distances$chr == chr]))
    max_val <- max(as.integer(rownames(distances)[distances$chr == chr]))
    c(chr, 0, max_val - min_val)
  })
  genome <- as.data.frame(t(genome))
  colnames(genome) <- cols_coor <- c("chr", "start", "end")
  genome$chr <- as.factor(genome$chr)

  # Transforming data values to bed format:
  cols_nochr <- c(colnames(distances) != "chr")
  bed <- lapply(unique(distances$chr), function(chr) {
    cases <- distances$chr == chr
    ncases <- sum(cases)
    cbind(
      rep(chr, ncases),
      seq(0, ncases - 1, 1),
      seq(0, ncases - 1, 1),
      distances[cases, colnames(distances)[cols_nochr]]
    )
  })
  bed <- do.call(rbind, bed)
  bed <- as.data.frame(bed)
  colnames(bed)[1:3] <- colnames(genome)
  bed$chr <- as.factor(bed$chr)
  cols_data <- colnames(bed)[!colnames(bed) %in% c(cols_coor, "cov", "gc")]
  range_global <- range(bed[, cols_data], na.rm = T)
  width <- range_global[2] - range_global[1]

  # Pre-processing:
  inversion <- c(1, -1)[1 + inner_closest]
  normalization <- list("none", normalization)[[1 + !is.null(normalization)]]
  smoothing <- c(0, smoothing)[1 + !is.null(smoothing)]
  ladder_values <- list(
    c(range_global[1] - width, range_global[2] + width),
    ladder_values
  )[[1 + !is.null(ladder_values)]]

  # Function to normalize values distribution to a less skewed one:
  norm_funcs <- list(
    none = function(x) {
      x
    },
    logarith = function(x) {
      x <- x + 0.0625
      x[!is.na(x)] <- log2(x[!is.na(x)])
      x[is.infinite(x)] <- min(x[!is.infinite(x)], na.rm = TRUE)
      x
    },
    squaroot = function(x) {
      x[!is.na(x)] <- sqrt(x[!is.na(x)])
      x
    },
    logistic = function(x) {
      1 / (1 + exp(-x))
    },
    hyperbol = function(x) {
      tanh(x)
    },
    atangent = function(x) {
      atan(x)
    },
    gompertz = function(x) {
      sigmoid(x, "Gompertz")
    }
  )

  # Function to impute missing values taking intermediate values:
  impu_fun <- function(distance_vec, chr_from) {

    # Working each chromosome independently:
    imputed_vec <- sapply(unique(chr_from), function(chr) {
      chr_impu <- distance_vec[which(chr_from == chr)]
      chr_lim <- c(1, length(chr_impu))
      na_rows <- which(is.na(chr_impu))
      opt <- 1 + (length(na_rows) > 0)

      # Imputing for each group of contiguous NAs:
      starting <- list(0, sapply(na_rows, function(x) !((x - 1) %in% na_rows)))
      starting <- starting[[opt]]
      starting <- na_rows[starting]
      ending <- list(0, sapply(na_rows, function(x) !((x + 1) %in% na_rows)))
      ending <- ending[[opt]]
      ending <- na_rows[ending]
      for (na_group in seq_len(length(starting))) {

        # Min. value for the imputation:
        opt <- (starting[na_group] - 1) >= chr_lim[1]
        min_val <- c(
          chr_impu[min(ending[na_group] + 1, length(chr_impu))],
          chr_impu[starting[na_group] - 1]
        )
        min_val <- min_val[1 + opt]

        # Max. value for the imputation:
        opt <- (ending[na_group] + 1) <= chr_lim[2]
        max_val <- c(
          chr_impu[max(starting[na_group] - 1, 1)],
          chr_impu[ending[na_group] + 1]
        )
        max_val <- max_val[1 + opt]

        # Contiguous imputed values (first and last values not used):
        n <- ending[na_group] - starting[na_group] + 3
        contiguous <- head(seq(min_val, max_val, length.out = n)[-1], -1)
        chr_impu[starting[na_group]:ending[na_group]] <- contiguous
      }
      chr_impu
    })
    do.call(c, imputed_vec)
  }

  # Function to smooth values using weighted points windows:
  smooth_fun <- function(distance_vec, chr_from, smooth) {
    smooth_vec <- sapply(unique(chr_from), function(chr) {
      chr_rows <- chr_from == chr
      chr_smoothed <- distance_vec[chr_rows]
      span_percen <- smooth
      span_points <- smooth / sum(chr_rows)
      final_span <- c(span_percen, span_points)[1 + (smooth > 1)]
      na_rows <- is.na(chr_smoothed)
      smooth_matrix <- cbind(chr_smoothed[!na_rows], 1:sum(!na_rows))
      smooth_matrix <- data.frame(smooth_matrix)
      if (final_span != 0) {
        chr_smoothed[!na_rows] <- loess(
          formula = "X1 ~ X2",
          data = smooth_matrix,
          degree = 1,
          span = final_span
        )$fitted
      }
      chr_smoothed
    })
    do.call(c, smooth_vec)
  }

  # Function to generate random colors for chromosomes/sectors:
  color_fun <- function(n_colors) {
    chr_colors <- rainbow(n_colors)
    chr_colors <- chr_colors[c(
      seq(1, length(chr_colors), 4), seq(2, length(chr_colors), 4),
      seq(3, length(chr_colors), 4), seq(4, length(chr_colors), 4)
    )]
  }

  # Function to check if a command line is available:
  command_available <- function(command) {
    check_result <- suppressWarnings(system2(
      command = command,
      args = "--version",
      stdout = FALSE,
      stderr = FALSE
    ))
    check_result == 0
  }

  # Function to join saved plots in a single PDF file:
  join_pdf <- function(input, output) {
    parameters <- c(
      "-q",
      "-dBATCH",
      "-dNOPAUSE",
      "-sDEVICE=pdfwrite",
      "-dPDFSETTINGS=/prepress"
    )
    possible_actions <- list(
      nothing = function(...) {

      },
      joining = function(input, output) {
        system2(
          command = "gs",
          args = c(
            parameters,
            paste0("-sOutputFile=", output),
            input
          ),
          stdout = FALSE,
          stderr = FALSE,
          wait = TRUE
        )
        invisible(file.remove(input))
      }
    )
    opt <- 1 + (command_available("gs") && (length(input) > 0))
    possible_actions[[opt]](input, output)
  }

  # Useful parameters for the plots:
  free_espace <- 0.1 * sum(c("cov" %in% colnames(bed), "gc" %in% colnames(bed)))
  free_espace <- 0.88 - free_espace
  ladder_pos <- round(nrow(genome) * 2.5 / 4)
  opt <- 1 + (length(chr_colors) != nrow(genome))
  chr_colors <- list(chr_colors, color_fun(nrow(genome)))[[opt]]
  text_color <- "black"
  back_color <- "grey"

  # Processing data (working each normalization independently):
  all_plots <- sapply(normalization, function(norm) {
    n_norm <- which(normalization == norm)

    # Data normalization and smoothing (each experiment/column independently):
    norm_values <- sapply(cols_data, function(column) {

      # Transformation/normalization to correct skewness:
      col_values <- inversion * norm_funcs[[norm]](bed[, column])

      # Imputation, replacing NAs for intermediate values (unmixed chrs.):
      opt <- 1 + na_imputation
      col_values <- list(col_values, impu_fun(col_values, bed$chr))
      col_values <- col_values[[opt]]

      # Smoothing data (without mixing chrs. again):
      opt <- 1 + (smoothing > 0)
      col_values <- list(col_values, smooth_fun(col_values, bed$chr, smoothing))
      col_values <- col_values[[opt]]

      # Returning finally ready values:
      col_values
    })
    range_global_norm <- range(norm_values, na.rm = TRUE)

    # Ladder normalization (without smoothing!):
    opt <- 1 + (length(ladder_values) == 1 && ((ladder_values[1] %% 1) == 0))
    int_ladder <- seq(
      from = range_global_norm[1],
      to = range_global_norm[2],
      length.out = round(abs(ladder_values[1])) + 1
    )[-1]
    vec_ladder <- inversion * norm_funcs[[norm]](ladder_values)
    norm_ladder <- list(vec_ladder, int_ladder)[[opt]]
    norm_ladder <- norm_ladder[c(norm_ladder >= range_global_norm[1] &
      norm_ladder <= range_global_norm[2])]

    # Plotting each experiment/column independently (after transforming all of
    #   them to know global max. and min. -same scale for all columns/cases-):
    norm_plots <- lapply(cols_data, function(column) {
      range_column <- range(bed[, column], na.rm = TRUE)
      range_column_norm <- range(norm_values[, column], na.rm = TRUE)

      # Initializing the plot:
      circos.clear()
      circos.par(
        track.height = 0.8,
        gap.degree = 2,
        cell.padding = c(0, 0, 0, 0)
      )
      circos.initialize(
        factors = genome$chr,
        xlim = genome[, c("start", "end")]
      )
      text(
        x = 0,
        y = 0,
        col = "#f8f8ff20",
        labels = paste(rep(paste0(rep(paste("Beekman Lab @ CRG -", now), 2),
          collapse = "   "
        ), 24), collapse = "\n")
      )
      draw.sector(
        rou1 = 0.02,
        col = back_color,
        border = back_color
      )

      # Genome structure:
      circos.track(
        ylim = c(0, 1),
        bg.col = rep(c("azure4","azure3"),20),
        bg.border = FALSE,
        track.height = 0.06,
        panel.fun = function(x, y) {
          chr <- CELL_META$sector.index
          xlim <- CELL_META$xlim
          ylim <- CELL_META$ylim
          circos.text(
            x = mean(xlim),
            y = mean(ylim),
            labels = chr,
            cex = 0.55,
            col = text_color,
            facing = "bending.inside",
            niceFacing = TRUE
          )
        }
      )

      # Genomes x axis:
      brk <- seq(0, max(genome$end), length.out = 6)
      circos.track(
        track.index = get.current.track.index(),
        bg.border = FALSE,
        panel.fun = function(x, y) {
          circos.axis(
            h = "top",
            major.at = brk,
            labels = "",
            # labels      = seq(0, 250, by = 50),
            # labels      = round(brk),
            labels.cex = 0.4,
            col = text_color,
            labels.col = text_color,
            lwd = 0.7,
            labels.facing = "clockwise"
          )
        }
      )

      # Coverage:
      if ("cov" %in% colnames(bed)) {
        circos.genomicTrack(
          data = bed[, c(cols_coor, "cov")],
          numeric.column = 4,
          track.height = 0.08,
          bg.border = FALSE,
          panel.fun = function(region, value, ...) {
            circos.genomicLines(
              region = region,
              value = value,
              type = "l",
              col = "grey50",
              lwd = 0.6
            )
            circos.segments(
              x0 = 0,
              x1 = max(bed$end),
              y0 = min(bed$cov, na.rm = TRUE),
              y1 = min(bed$cov, na.rm = TRUE),
              lwd = 0.6,
              lty = "11",
              col = back_color
            )
            circos.segments(
              x0 = 0,
              x1 = max(bed$end),
              y0 = max(bed$cov, na.rm = TRUE),
              y1 = max(bed$cov, na.rm = TRUE),
              lwd = 0.6,
              lty = "11",
              col = back_color
            )
          }
        )
      }

      # GC content:
      if ("gc" %in% colnames(bed)) {
        circos.genomicTrack(
          data = bed[, c(cols_coor, "gc")],
          numeric.column = 4,
          track.height = 0.08,
          bg.border = FALSE,
          ylim = range(bed$gc, na.rm = TRUE),
          panel.fun = function(region, value, ...) {
            circos.genomicLines(
              region = region,
              value = value,
              type = "l",
              col = "grey50",
              lwd = 0.6
            )
            circos.segments(
              x0 = 0,
              x1 = max(bed$end),
              y0 = 0.3,
              y1 = 0.3,
              lwd = 0.6,
              lty = "11",
              col = back_color
            )
            circos.segments(
              x0 = 0,
              x1 = max(bed$end),
              y0 = 0.5,
              y1 = 0.5,
              lwd = 0.6,
              lty = "11",
              col = back_color
            )
            circos.segments(
              x0 = 0,
              x1 = max(bed$end),
              y0 = 0.7,
              y1 = 0.7,
              lwd = 0.6,
              lty = "11",
              col = back_color
            )
          }
        )
      }

      # 3D distances:
      circos.track(
        factors = bed$chr,
        x = bed$start,
        y = norm_values[, column],
        ylim = range_global_norm,
        track.height = free_espace,
        bg.border = FALSE,
        panel.fun = function(x, y) {
          chr <- CELL_META$sector.index

          # Adding the ladder:
          sapply(norm_ladder, function(z) {
            circos.lines(
              x = c(1, length(x) - 1),
              y = c(z, z),
              col = back_color,
              lty = 2,
              lwd = 1.2
            )
            txt <- c("", as.character(round(z * inversion, 1)))
            txt <- paste0("\n", txt[1 + (chr == ladder_pos)])
            suppressMessages(
              circos.text(
                x = mean(CELL_META$xlim),
                y = z,
                labels = txt,
                cex = 0.55,
                col = back_color, # text_color,
                facing = "downward",
                adj = c(0.3, 0.7),
                niceFacing = TRUE
              )
            )
          })

          # Plotting distance lines:
          circos.lines(
            x = x,
            y = y,
            type = "l",
            col = chr_colors[genome$chr == chr],
            cex = 0.4,
            lwd = 1.5
          )

          # # Adding a final reference for the inner value:
          width <- range_column_norm[1] - range_column_norm[2]
          last_ladder <- c(
            suppressWarnings(min(norm_ladder)),
            range_column_norm[2]
          )
          last_ladder <- last_ladder[1 + (length(norm_ladder) == 0)]
          center_dis <- range_column_norm[1] - last_ladder
          enough_dis <- (center_dis / width) >= 0.1
          if (chr == ladder_pos && inner_closest && enough_dis) {
            txt <- round(range_global_norm[1] * inversion, 1)
            pos <- min(norm_values, na.rm = TRUE)
            suppressMessages(
              circos.text(
                x = mean(CELL_META$xlim),
                y = pos,
                labels = as.character(txt),
                cex = 0.55,
                col = back_color, # text_color,
                facing = "downward",
                adj = c(-0.2, 0.1),
                niceFacing = TRUE
              )
            )
          }
        }
      )
      title(paste0("'", column, "' data"))
      mtext(paste0(
        "(range for the raw data: ",
        paste0(round(range_column, 2), collapse = " to "),
        "; range for the data after normalization and smoothing: ",
        min(round(range_column_norm * inversion, 2)),
        " to ",
        max(round(range_column_norm * inversion, 2)),
        "; plot range: ",
        min(round(range_global_norm * inversion, 2)),
        " to ",
        max(round(range_global_norm * inversion, 2)),
        ")\n(plot obtained by using '",
        norm,
        "' data normalization, smoothing with '",
        smoothing,
        "' points and closest regions in the ",
        c("out", "inn")[1 + inner_closest],
        "er part of the section)\n"
      ),
      side = 1,
      cex = 0.7
      )

      # Saving the plot:
      column_plot <- recordPlot()
      if (!is.null(run_name)) {
        pdf(file.path(
          folder,
          paste0(
            now, "_",
            n_norm, "_", norm, "_",
            which(cols_data == column), "_", column, ".pdf"
          )
        ))
        print(column_plot)
        dev.off()
      }

      # Returning the plot for the each 'sapply' (each experiment/column):
      column_plot
    })

    # Returning the plots for the each transformation/normalization:
    names(norm_plots) <- cols_data
    norm_plots
  })

  # Joining plots in a single PDF:
  pattern <- c("an_impossible_pattern_as_it_is", paste0("^", now, "_.*\\.pdf$"))
  pattern <- pattern[1 + !is.null(run_name)]
  input_files <- list.files(
    path = folder,
    pattern = pattern,
    full.names = TRUE
  )
  output_file <- file.path(folder, paste0(run_name, ".pdf"))
  invisible(join_pdf(input_files, output_file))

  # Finishing the function and returning all the final plots:
  all_plots_names <- as.vector(outer(rownames(all_plots),
    colnames(all_plots),
    paste,
    sep = "."
  ))
  all_plots <- c(all_plots)
  names(all_plots) <- all_plots_names
  return(all_plots)

  #
  # TO-DO: Change the GC content track (and probably the coverage one) by tracks
  #   (kind of heatmaps) with ChIP-seq data. Allow multiple 'smoothings' in the
  #   same run.
  #
}


chromatin_circos <- function(distances,
                             run_name = NULL,
                             na_imputation = FALSE,
                             normalization = NULL,
                             smoothing = NULL,
                             ladder_values = NULL,
                             inner_closest = FALSE,
                             chr_colors = NULL)# LICENSE NOTICE ---------------------------------------------------------------
#
# Copyright(c) 2020, José R. Hernández and the Beekman Lab (Single Cell
# Epigenomics and Cancer Development group at Centre for Genomic Regulation).
#
# This script and its functions are free software: you can redistribute them
# and/or modify them under the terms of the GNU General Public License as
# published by the Free Software Foundation, either version 3 of the License,
# or any later version.
#
# This script and its functions are distributed in the hope that they will be
# useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# this program. If not, see <https://www.gnu.org/licenses/>.
#
#