#' Two-Mode KL-Means Heuristic
#'
#' @description This function runs two-mode K-means for an \eqn{RO x CO} network matrix.
#' @param A An \eqn{RO x CO} two-mode network matrix.
#' @param RC The number of clusters for row objects (\eqn{1 < RC < RO}).
#' @param CC The number of clusters for column objects (\eqn{1 < CC < CO}).
#' @param TLIMIT A desired time limit.
#' @return The function returns the following:
#' \itemize{
#' \item \code{vaf} - the variance-accounted-for;
#' \item \code{RP} - an \eqn{RO}-dimensional vector of row cluser assignements;
#' \item \code{RC} - an \eqn{RC}-dimensional vector of column cluser assignements;
#' \item \code{restarts} - the number of restarts within the time limit.
#' }
#' @examples
#' # Load the Turning Point Project network (Brusco & Doreian, 2015) data.
#' data("nyt")
#'
#' # Run two-mode K-means procedure.
#' res <- tmklm(nyt,RC = 9,CC = 5,TLIMIT = 1)
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
#' @export

tmklm = function(A,RC,CC,TLIMIT) {
	RO = dim(A)[1]
	CO = dim(A)[2]
	VAF = 0
  NREPS = 0
	RBEST <- matrix(0, nrow = RO, ncol = 1)
	CBEST <- matrix(0, nrow = CO, ncol = 1)
	res =.Fortran("tmklmf",as.integer(RO),as.integer(CO),as.integer(RC),as.integer(CC),as.double(TLIMIT),as.double(A),as.integer(RBEST),as.integer(CBEST),as.double(VAF),as.integer(NREPS))
#	res =.Fortran("tmklmf",as.integer(RO),as.integer(CO),as.integer(RC),as.integer(CC),as.double(TLIMIT),as.double(A))
	RP <- res[[7]]
	CP <- res[[8]]
	vaf <- res[[9]]
	restarts <- res[[10]]
	return(list(RP=RP, CP=CP, vaf=vaf, restarts=restarts))
}
