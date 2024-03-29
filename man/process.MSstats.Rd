% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/msstats_processing.R
\name{process.MSstats}
\alias{process.MSstats}
\title{Process feature intensities with MSstats default parameters}
\usage{
process.MSstats(raw.features, filepath = "MSstats.rds", recompute = F)
}
\arguments{
\item{raw.features}{MSstats feature data. Can be generated by preprocess.MSstats()}

\item{recompute}{Boolean on whether to recompute normalization. Default is False}

\item{file.path}{Path to save the resulting normalized data object}
}
\value{
MSstats object containing normalized feature data
}
\description{
Normalizes and processes MRM data using default/recommended parameters for MSstats.
Saves result to disk and reloads on subsequent knits unless recompute is True
}
\examples{
rd.proposed <- process.MSstats(rd.msstats, filepath=paste0(proj.dir,'Analysis/MSstats.rds'), recompute = F)

}
