#' Two-Mode Blockmodeling (Structural Equivalence) Heuristic
#'
#' @description This function runs two-mode KL-medians for an \eqn{RO x CO} two-mode binary network matrix.
#' @param A An \eqn{RO x CO} two-mode binary network matrix.
#' @param RC The number of clusters for row objects (\eqn{1 < RC < RO}).
#' @param CC The number of clusters for column objects (\eqn{1 < CC < CO}).
#' @param TLIMIT A desired time limit.
#' @return The function returns the following:
#' \itemize{
#' \item \code{objval} - total number of inconsistencies;
#' \item \code{RP} - an \eqn{RO}-dimensional vector of row cluser assignements;
#' \item \code{RC} - an \eqn{RC}-dimensional vector of column cluser assignements;
#' \item \code{restarts} - the number of restarts within the time limit.
#' }
#' @examples
#' # Load the Turning Point Project network (Brusco & Doreian, 2015) data.
#' data("nyt")
#'
#' # Run the two-mode blockmodeling heuristic procedure.
#' res <- tmklmed(nyt, RC = 9, CC = 5, TLIMIT = 1)
#'
#' # See the results.
#' res
#' @author Michael Brusco
#' @references
#' Brusco, M. J., Doreian, P., & Steinley, D. (2019). Deterministic blockmodeling of signed and two-mode networks: a tutorial with psychological examples. \emph{British Journal of Mathematical and Statistical Psychology}.
#'
#' Doreian, P., Batagelj, V., & Ferligoj, A. (2004). Generalized blockmodeling of two-mode network data. \emph{Social Networks}, 26, 29-53. doi:10.1016/j.socnet.2004.01.002
#'
#' Brusco, M., Stolze, H. J., Hoffman, M., Steinley, D., & Doreian, P. (2018). Deterministic blockmodeling of two-mode binary network data using two-mode KL-median partitioning. \emph{Journal of Social Structure}, 19, 1-21. Retrieved from: https://www.exeley.com/exeley/journals/journal_of_social_structure/19/1/pdf/10.21307_joss-2018-007.pdf
#' @export

tmklmed = function(A,RC,CC,TLIMIT) {
	RO = dim(A)[1]
	CO = dim(A)[2]
	GBEST = 0
      NREPS = 0
	GR <- matrix(0, nrow = RO, ncol = 1)
	GC <- matrix(0, nrow = CO, ncol = 1)
	res =.Fortran("tmklmedf",as.integer(RO),as.integer(CO),as.integer(RC),as.integer(CC),as.double(TLIMIT),as.integer(A),as.integer(GR),as.integer(GC),as.integer(GBEST),as.integer(NREPS))
	RP <- unlist(res[[7]])
	CP <- res[[8]]
	objval <- res[[9]]
	restarts <- res[[10]]
	return(list(RP=RP, CP=CP, objval=objval, restarts=restarts))
}
