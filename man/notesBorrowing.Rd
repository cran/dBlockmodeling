\name{notesBorrowing}
\alias{notesBorrowing}
\docType{data}
\title{
The notes borrowing network between social-informatics students
}
\description{
The data come from a survey conducted in May 1993 on 13 social-informatics
students (Hlebec, 1996). The network was constructed from answers to the question, “How often did you borrow notes from this person?” for each of the fellow students. The respondents indicated the frequency of borrowing by choosing (on a computer) a line of length 1–20, where 1 meant no borrowing. 1 was deducted from all answers, so that 0 now means no borrowing. The data was first used for blockmodeling in Žiberna (2007).
}
\usage{data("notesBorrowing")}
\format{
The data set is a valued matrix with 13 rows and columns.
}
%%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%%}
%%\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
%%}
\references{
Hlebec, V., (1996). \emph{Metodološke značilnosti anketnega zbiranja podatkov v analizi omrežji: Magistersko delo}. FDV, Ljubljana.

Žiberna, A. (2007). Generalized blockmodeling of valued networks. \emph{Social Networks}, 29, 105–126. https://doi.org/10.1016/j.socnet.2006.04.002
}
\examples{
data(notesBorrowing)

# Plot the network.
# (The function plotMat is from blockmodeling package.)
# plotMat(notesBorrowing)
}
\keyword{datasets}
