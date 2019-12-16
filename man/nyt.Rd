\name{nyt}
\alias{nyt}
\docType{data}
\title{
The New York Times advertisement data
}
\description{
The data pertain to the Turning Point Project (TPP), where \eqn{n = 108} organizations signed one or more of \eqn{m = 25} environmental activist-oriented, full-page advertisements in the New York Times (NYT) during 1999-2000.
}
\usage{data("nyt")}
\format{
The data set is a binary network matrix with 108 rows (organizations) and 25 columns (full-page advertisements). An element of the binary network matrix assumes a value of 1 if organization \eqn{i} signed advertisement \eqn{j} and 0 otherwise.
}
%%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%%}
%%\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
%%}
\references{
Brusco, M., Doreian, P., Steinley, D., & Satornino, C. B. (2013). Multiobjective blockmodeling for social network analysis. \emph{Psychometrika}, 78 (3), 498-525

Brusco, M., & Doreian, P. (2015). A real-coded genetic algorithm for two-mode KL-means partitioning with application to homogeneity blockmodeling. \emph{Social Networks}, 41, 26-35.
}
\examples{
data(nyt)

# Plot the network.
# (The function plotMat is from blockmodeling package.)
# plotMat(nyt)
}
\keyword{datasets}
