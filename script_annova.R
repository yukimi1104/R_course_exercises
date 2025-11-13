#Two-way ANOVA: Maternal and larval host effects
#This script reads butterflies.csv,fits two-way ANOVAs (Type II) for DevelopmentTime,
#AdultWeight and GrowthRate,prints ANOVA tables in a "three-line" style,
#computes means ± SE and produces black-and-white line plots with error bars, 

#packages
## If needed, install packages once:
## install.packages(c("car", "dplyr", "ggplot2", "readr"))
library(car)        # Type-II ANOVA (car::Anova)
library(dplyr)      # Data manipulation (group_by, summarise)
library(ggplot2)    # Figures
library(readr)      # Robust CSV import

#path and directory
# Relative path to the data file.
data_path  <- "butterflies.csv"
# Directory where figures will be exported.
output_dir <- "output"
# Create output directory if it does not exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
#Import data and read the csv file
dat <- read_csv(data_path, show_col_types = FALSE)
# Convert host variables to categorical factor
dat$MaternalHost <- factor(dat$MaternalHost, levels = c("Barbarea", "Berteroa"))
dat$LarvalHost   <- factor(dat$LarvalHost,   levels = c("Barbarea", "Berteroa"))
# Two-way ANOVA (Type II)
run_anova <- function(resp) {
  # resp: character, name of the response variable
  #       e.g. "DevelopmentTime", "AdultWeight", "GrowthRate"
  
  fm <- lm(
    reformulate(
      c("MaternalHost", "LarvalHost", "MaternalHost:LarvalHost"),
      response = resp
    ),
    data = dat
  )
  an <- car::Anova(fm, type = 2)  # Type II sums of squares
  list(model = fm, anova = an)
}

# Compute group means and standard errors
mean_se <- function(resp) {
  # Compute mean and standard error by MaternalHost × LarvalHost
  dat %>%
    group_by(MaternalHost, LarvalHost) %>%
    summarise(
      mean = mean(.data[[resp]], na.rm = TRUE),
      se   = sd(.data[[resp]],  na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    )
}

#Print ANOVA table in a "three-line" style in the console
format_anova_console <- function(an, title = NULL) {
  df <- as.data.frame(an)
  df$Source <- rownames(df); rownames(df) <- NULL
  names(df) <- c("Sum Sq", "Df", "F value", "Pr(>F)", "Source")
  df <- df[, c("Source", "Sum Sq", "Df", "F value", "Pr(>F)")]
  # Rounding digits
  df$`Sum Sq`  <- round(df$`Sum Sq`, 3)
  df$`F value` <- ifelse(is.na(df$`F value`), NA, round(df$`F value`, 2))
  # Format p-values as scientific notation strings
  df$`Pr(>F)`  <- ifelse(
    is.na(df$`Pr(>F)`),
    "",
    format(df$`Pr(>F)`, scientific = TRUE, digits = 2)
  )
  # Print in a "three-line" style
  if (!is.null(title)) {
    cat("\n", title, "\n", sep = "")
  } else {
    cat("\nANOVA table\n")
  }
  cat(strrep("-", 60), "\n")
  cat(sprintf("%-20s %10s %5s %10s %10s\n",
              "Source", "Sum Sq", "Df", "F value", "Pr(>F)"))
  cat(strrep("-", 60), "\n")
  for (i in seq_len(nrow(df))) {
    cat(sprintf("%-20s %10.3f %5d %10s %10s\n",
                df$Source[i],
                df$`Sum Sq`[i],
                df$Df[i],
                ifelse(is.na(df$`F value`[i]), "",
                       sprintf("%.2f", df$`F value`[i])),
                df$`Pr(>F)`[i]))
  }
  cat(strrep("-", 60), "\n\n")
  
  invisible(df)
}
# Black-and-white line plot with error bars
plot_line_bw <- function(df, ylab, title) {
  ggplot(
    df,
    aes(x = LarvalHost,
        y = mean,
        group = MaternalHost,
        shape = MaternalHost,
        linetype = MaternalHost)
  ) +
    geom_line(size = 0.7, colour = "black") +
    geom_errorbar(
      aes(ymin = mean - se, ymax = mean + se),
      width = 0.05,
      colour = "black"
    ) +
    geom_point(size = 3, colour = "black") +
    # Barbarea: dashed line + open circle; Berteroa: solid line + filled circle
    scale_shape_manual(values = c("Barbarea" = 1, "Berteroa" = 16)) +
    scale_linetype_manual(values = c("Barbarea" = "dashed",
                                     "Berteroa" = "solid")) +
    labs(
      x = "Larval host",
      y = ylab,
      shape    = "Maternal host",
      linetype = "Maternal host",
      title = title
    ) +
    theme_bw() +
    theme(
      text = element_text(size = 12),
      legend.position = "right"
    )
}

#Run with each outcome
res_time   <- run_anova("DevelopmentTime")
res_weight <- run_anova("AdultWeight")
res_growth <- run_anova("GrowthRate")

# Adjusted R^2 values 
adjR2_time   <- summary(res_time$model)$adj.r.squared
adjR2_weight <- summary(res_weight$model)$adj.r.squared
adjR2_growth <- summary(res_growth$model)$adj.r.squared

#Print table

format_anova_console(res_time$anova,
                     title = "Type-II ANOVA for DevelopmentTime")
format_anova_console(res_weight$anova,
                     title = "Type-II ANOVA for AdultWeight")
format_anova_console(res_growth$anova,
                     title = "Type-II ANOVA for GrowthRate")

#Means and SE for plotting 
ms_time   <- mean_se("DevelopmentTime")
ms_weight <- mean_se("AdultWeight")
ms_growth <- mean_se("GrowthRate")

#Black-and-white line plots with error bars
#Development time
p_time_line <- plot_line_bw(
  ms_time,
  ylab  = "Development time (days)",
  title = "Development time by maternal and larval host"
)
print(p_time_line)
ggsave(
  filename = file.path(output_dir, "fig_devtime_line.png"),
  plot     = p_time_line,
  width    = 6,
  height   = 4,
  dpi      = 300
)
# Adult weight
p_weight_line <- plot_line_bw(
  ms_weight,
  ylab  = "Adult weight (mg)",
  title = "Adult weight by maternal and larval host"
)
print(p_weight_line)
ggsave(
  filename = file.path(output_dir, "fig_weight_line.png"),
  plot     = p_weight_line,
  width    = 6,
  height   = 4,
  dpi      = 300
)
#  Growth rate
p_growth_line <- plot_line_bw(
  ms_growth,
  ylab  = "Growth rate (1/day)",
  title = "Growth rate by maternal and larval host"
)
print(p_growth_line)
ggsave(
  filename = file.path(output_dir, "fig_growth_rate_line.png"),
  plot     = p_growth_line,
  width    = 6,
  height   = 4,
  dpi      = 300
)



