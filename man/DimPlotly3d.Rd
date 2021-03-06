% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DimPlotly3D.R
\name{DimPlotly3d}
\alias{DimPlotly3d}
\title{DimPlotly3D}
\usage{
DimPlotly3d(seuratObj, group.by = "ident", do.label = FALSE,
  label.size = 12, label.color = "000000", show.arrow = FALSE,
  do.return = FALSE, pt.size = 2, pt.shape = "circle", opacity = 1,
  reduction.use = "tsne", dim.1 = 1, dim.2 = 2, dim.3 = 3,
  palette.use = "Set1", plot.height = 900, plot.width = 900,
  plot.title = NULL, pt.info = NULL, legend = TRUE,
  legend.font.size = 12, plot.grid = FALSE, plot.axes = FALSE)
}
\arguments{
\item{seuratObj}{Seurat object}

\item{group.by}{Variable by which to group cells. Currently only works with the current ident and column names from meta.data (default: ident)}

\item{do.label}{Add a label showing thr group name to the graph (default: FALSE)}

\item{label.size}{Label font size (default: 12)}

\item{label.color}{Color for label border and arrow.  Need hex value. (default = '000000')}

\item{show.arrow}{Offset the position of the labels and instead point to each group with an arrow (default: FALSE)}

\item{do.return}{Return the plot object instead of displaying it (default: FALSE)}

\item{pt.size}{Size of the points in pixels (default: 2)}

\item{pt.shape}{Shape to use for the points (default: circle)}

\item{opacity}{Transparency level to use for the points on a 0-1 scale (default: 1)}

\item{reduction.use}{Dimensional reduction to display (default: tsne)}

\item{dim.1}{Dimension to display on the x-axis (default: 1)}

\item{dim.2}{Dimension to display on the y-axis (default: 2)}

\item{dim.3}{Dimension to display on the z-axis (default: 3)}

\item{palette.use}{Color palette to use.  Must be a palette available in the Paletteer package.  (default: 'Set1')}

\item{plot.height}{Plot height in pixels (default: 900)}

\item{plot.width}{Plot width in pixels (default: 900)}

\item{plot.title}{Plot title (default: reduction.use)}

\item{pt.info}{Meta.data columns to add to the hoverinfo popup. (default: ident)}

\item{legend}{Display legend? (default: TRUE)}

\item{legend.font.size}{Legend font size (default: 12)}

\item{plot.grid}{Display the major unit tick marks? (default: FALSE)}

\item{plot.axes}{Display the major x, y, and z axes? (default: FALSE)}
}
\value{
plotly object
}
\description{
DimPlotly3D
}
\examples{
DimPlotly3D(seuratObj, group.by = "res.0.6", do.label = TRUE, show.arrow = FALSE)
}
