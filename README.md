This package has some functions that make MSstats a little nicer to use.

## Data processing

* `preprocess.MSstats` takes a CSV file generated by Skyline using the MSstats report functionality and a data frame with "Condition", "Run", and "BioReplicate" columns, does some cleanup and annotation, and produces an object suitable for processing.
* `process.MSstats` takes the output of `preprocess.MSstats` and performs normalization using "best practice" parameters for MSstats. It will write the results to a file (according to the `filepath` argument), and will reload them from disk if they are present. This saves time when compiling Rmarkdown documents. The argument `recompute` should be set to `TRUE` if you want to reprocess the data even if a file object has been previously saved.
* `compare.MSstats` works similarly to `process.MSstats`, but its input is a processed MSstats object and a comparison matrix defined according to the MSstats user manual.

## Plotting

The built-in MSstats plots aren't terrible, but they don't return the ggplot object for further customization. These functions make some nice default plots, but importantly return the object so it can be customized with additional options or geoms.

* `draw.Condition` takes a normalized MSstats object and a specific protein ID and generates a plot showing the quantification for each experimental condition. It performs an internal simple Light/Heavy normalization. The output plot is automatically a dot plot when *n* < 5, or a boxplot when *n* >= 5
* `draw.Volcano` takes an MSstats comparison result object and a specific comparison name and generates a nice-looking volcano plot (adj.pvalue vs. log2FC). Points are colored by significance threshold, and any proteins with log2FC > 1 are labeled.
* `print.Conditions` and `print.Volcanos` are simple for loops that will iterate over an entire normalized MSstats object or comparison object, respectively, and generate those plots for all proteins or comparisons. These are just printed to the session and should be saved separately if desired, or piped to a file object.