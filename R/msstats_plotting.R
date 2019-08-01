# Plotting ----------------------------------------------------------------
#' Draw condition plots
#'
#' Draw an MSstats condition plot for a specific protein. Will perform boxplot for n >= 5, or else a dotplot.
#'
#' @param df A dataframe with processed/normalized MSstats results
#' @param protein A protein identifier
#'
#' @return a ggplot object
#'
#' @examples
#' draw.Condition(rd.proposed$ProcessedData, "S100A12-P80511")
#'
#' @export
draw.Condition <- function(df, protein) {
  # Given a dataframe with processed/normalized MSstats results and a protein
  # Draw the condition plot

  if(!require(tidyverse)) library(tidyverse)
  if(!require(viridisLite)) library(viridisLite)

  prot.data <- df %>% filter(PROTEIN == protein & censored != TRUE)

  # Calculate mimum n in order to set plot type
  n.features <- length(unique(prot.data$FEATURE))

  prot.data %>%
    filter(LABEL == 'L') %>%
    group_by(GROUP_ORIGINAL) %>%
    summarize(n_points = n()) %>%
    mutate(group_n = n_points/n.features) -> n.table

  min.n <- min(n.table$group_n)

  # Calculate L/H ratio
  # Then get the mean abundance by run across features (which are pre-normalized)
  # Save to an object for plotting
  prot.data %>%
    mutate(RUNGROUPFEATURE = paste(RUN, GROUP_ORIGINAL, FEATURE, sep='-')) %>%
    select(RUNGROUPFEATURE, LABEL, ABUNDANCE) %>%
    spread(key=LABEL, value=ABUNDANCE) %>%
    mutate(norm.abundance = L/H) %>%
    separate(RUNGROUPFEATURE, into=c('RUN','GROUP_ORIGINAL','FEATURE'), sep='-') %>%
    group_by(RUN, GROUP_ORIGINAL) %>%
    summarize(run.mean = mean(norm.abundance, na.rm=T)) %>%
    filter(run.mean != Inf) -> plot.df

  # Plot the data depending on how many points there are
  if (min.n >= 5) { # Boxplot
    plot.df %>%
      ggplot(aes(x=GROUP_ORIGINAL, y=run.mean)) +
      geom_boxplot(fill=viridis(1), alpha=0.66) +
      theme_minimal(18) + ggtitle(paste(protein, 'Condition Plot')) +
      xlab('Group') + ylab('Normalized Abundance') -> g
  } else { # Dotplot due to low n
    plot.df %>%
      ggplot(aes(x=GROUP_ORIGINAL, y=run.mean)) +
      geom_point(size=3, color=viridis(1), alpha=0.66) +
      theme_minimal(18) + ggtitle(paste(protein, 'Condition Plot')) +
      xlab('Group') + ylab('Normalized Abundance') -> g
  }

  return(g)
}

#' Print conditon plots for all proteins to the screen
#'
#' Given a dataframe of normalized MSstats data, prints condition plots for all proteins
#'
#' @param df A dataframe with processed/normalized MSstats results
#'
#' @return NULL
#'
#' @examples
#' print.Conditions(rd.proposed$ProcessedData)
#'
#' @export print.Conditions
print.Conditions <- function(df){
  # Given a processed/normalized MSstats dataframe,
  # Print the condition plot for every protein

  for(prot in unique(df$PROTEIN)) {
    print(draw.Condition(df, prot))
  }
}

#' Draw a Volcano Plot
#'
#' Draw a Volcano plot from an MSstats comparison
#'
#' @param df A dataframe containing MSstats comparison results
#' @param comparison The comparison to plot
#' @param sig.p The adjusted p-value for significance. Default is 0.05
#'
#' @return A ggplot object
#'
#' @examples
#' draw.Volcano(pd.comparisons$ComparisonResult, 'A3-A1')
#'
#' @export
draw.Volcano <- function(df, comparison, sig.p=0.05) {
  # Given a dataframe containing MSstats comparisons,
  # draw a nice-looking volcano plot
  if(!require(tidyverse)) library(tidyverse)
  if(!require(ggrepel)) library(ggrepel)
  if(!require(viridisLite)) library(viridisLite)

  df.filt <- df %>% filter(Label == comparison) %>% drop_na(pvalue)

  df.filt$threshold = as.factor(df.filt$adj.pvalue < sig.p)

  if(!length(df.filt$log2FC) > 0) {
    stop("Error: comparison not found in data, check spelling")
  }

  # Construct the plot object
  # Automatically sets limits of x-axis based on data
  # Draws visual guide lines at +/- 2 (real) FC and sig adjusted p-value
  # Labels any protein > 2 (real) FC regardless of adjusted p-value
  g <- ggplot(data=df.filt,
              aes(x=log2FC, y =-log10(adj.pvalue),
                  colour=threshold)) +
    geom_point(alpha=0.66, size=3) +
    xlab("log2 fold change") + ylab("-log10 p-value") +
    xlim(-1*max(abs(df.filt$log2FC), na.rm=T)-0.5,
         max(abs(df.filt$log2FC), na.rm=T)+0.5) +
    geom_hline(yintercept=-log10(sig.p), linetype='dashed', color='grey') +
    geom_vline(xintercept=-1, linetype='dashed', color='grey') +
    geom_vline(xintercept=1, linetype='dashed', color='grey') +
    theme_bw(18) +
    geom_text_repel(aes(label=ifelse(abs(log2FC)>1, as.character(Protein),''), size=3)) +
    theme(legend.position="none", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    scale_color_viridis_d(direction=-1) +
    ggtitle(comparison)

  return(g)
}

#' Print all Volcano Plots
#'
#' Prints Volcano plots for all performed comparisons to the screen
#'
#' @param df A dataframe with MSstats comparison results
#'
#' @return NULL
#'
#' @examples
#' print.Volcanos(pd.comparisons$ComparisonResult)
#'
#' @export print.Volcanos
print.Volcanos <- function(df) {
  # Given a dataframe with >1 MSstats comparison
  # Display the volcano plot for each of them
  for(comparison in unique(df$Label)) {
    print(draw.Volcano(df, comparison))
  }
}
