\name{sampsonT3}
\alias{sampsonT3}
\docType{data}
\title{
Sampson's monastery data (time period 3)
}
\description{
Sampson (1968) collected network data among \eqn{n = 18} trainee monks at multiple time periods. Network ties were collected with respect to several different relations including affect, esteem, influence, and sanction.

With respect to affect, Sampson (1968) asked each of the 18 trainees (egos) to identify and rank the three other members (alters) of the cohort they liked the most, as well the three members they liked the least.

The data \code{sampsonT3} refers to affect network data, collected at the third time point (Doreian et al., 2005, p. 33).
}
\usage{data("sampsonT3")}
\format{
The data set is a 18 x 18 (signed weighted) network matrix. The egos are in lines and the alters are in columns. The entries are egos' answers. The elements of the network matrix are discrete values between -3 and 3.

The edge from a given ego to the most-liked alter (by this ego) is assigned with an edge weight of \eqn{+3}, the edge to the second-most-liked alter is assigned with a weight of \eqn{+2}, and the edge to the third-most-liked alter is assigned with \eqn{+1}. Likewise, the edge from a given ego to the most-disliked alter is assigned with an edge weight of \eqn{-3}, the edge to the second-most-disliked alter is assigned with a weight of \eqn{-2}, and the edge to the third-most-liked alter with a value of \eqn{-1}.
}
%%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%%}
%%\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
%%}
\references{
Sampson, S. F. (1968). A novitiate in a period of change: An experimental case study of relationships. \emph{Unpublished Ph.D. dissertation}, Department of Sociology, Cornell University, Ithaca, NY.

Doreian, P., Batagelj, V., & Ferligoj, A. (2005). \emph{Generalized blockmodeling}, Cambridge, UK: Cambridge University Press.
}
\examples{
data(sampsonT3)

# Plot the network.
# (The function plotMat is from blockmodeling package.)
# plotMat(sampsonT3)
}
\keyword{datasets}
