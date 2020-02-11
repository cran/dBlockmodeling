#' One-Mode K-Means Heuristic
#'
#' @description This function runs one-mode K-means for an \eqn{RO x RO} network matrix.
#' @param A An \eqn{RO x RO} one-mode network matrix.
#' @param RC The number of clusters for row objects (\eqn{1 < RC < RO}).
#' @param TLIMIT A desired time limit.
#' @param IDIAG 0 if main diagonal to be ignored, any other value it will be included. Default is 0.
#' @return The function returns the following:
#' \itemize{
#' \item \code{sse} - the sum of the within-block sum-of-squared deviations from the block means;
#' \item \code{vaf} - the variance-accounted-for;
#' \item \code{RP} - an \eqn{RO}-dimensional vector of row cluser assignements;
#' \item \code{restarts} - the number of restarts within the time limit.
#' }
#' @examples
#' # Load the notes borrowing data..
#' data("notesBorrowing")
#'
#' #Run one-mode K-means procedure.
#' res <- omkm(notesBorrowing,RC = 3, TLIMIT = 1, IDIAG = 0)
#'
#' # See the results.
#' res
#' @author Michael Brusco
#' @references
#' Brusco, M. J., Doreian, P., & Steinley, D. (2019). Deterministic blockmodeling of signed and two-mode networks: a tutorial with psychological examples. \emph{British Journal of Mathematical and Statistical Psychology}.
#'
#' Baier, D., Gaul, W., & Schader, M. (1997). Two-mode overlapping clustering with applications in simultaneous benefit segmentation and market structuring. In R. Klar & O. Opitz (Eds), \emph{Classification and knowledge organization} (pp. 557-566), Heidelberg: Springer.
#'
#' Brusco, M., & Doreian, P. (2015). A real-coded genetic algorithm for two-mode KL-means partitioning with application to homogeneity blockmodeling. \emph{Social Networks}, 41, 26-35. http://dx.doi.org/10.1016/j.socnet.2014.11.007
#' Žiberna, A. (2020). K-means-based algorithm for blockmodeling linked networks. \emph{Social Networks}, 61, 153–169. https://doi.org/10.1016/j.socnet.2019.10.006
#' @export

omkm = function(A,RC,TLIMIT,IDIAG=0) {

	RO = dim(A)[1]
	VAF = 0
  ZBEST = 0
  NREPS = 0
	RBEST <- matrix(0, nrow = RO, ncol = 1)
#	dyn.load("c:/RBlockmodeling/omklm.dll")
	res =.Fortran("omkmf",as.integer(RO),as.integer(RC),as.integer(IDIAG),as.double(TLIMIT),as.double(A),as.integer(RBEST),as.double(ZBEST),as.double(VAF),as.integer(NREPS))
#	res =.Fortran("tmklm",as.integer(RO),as.integer(CO),as.integer(RC),as.integer(CC),as.double(TLIMIT),as.double(A))
	RP <- res[[6]]
	sse <- res[[7]]
	vaf <- res[[8]]
	restarts <- res[[9]]
	return(list(RP=RP, sse=sse, vaf=vaf, restarts=restarts))

}
