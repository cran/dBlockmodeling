#' Relocation Heuristic for Relaxed Structural Balance
#'
#' @description This function runs relocation heuristic for relaxed structural balance on an \eqn{M x M} asymmetric matrix. The main diagonal is ignored.
#' @param A An \eqn{N x N} signed network matrix.
#' @param C The number of clusters (\eqn{1 < C < N}, where \eqn{N} is the number of nodes).
#' @param TLIMIT A desired time limit.
#' @return The function returns the following:
#' \itemize{
#' \item \code{obj} - the Doreian & Mrvar objective value;
#' \item \code{P} - \eqn{N}-dimensional vector of cluser assignements; and
#' \item \code{restarts} - the number of restarts within the time limit.
#' }
#' @examples
#' # Load the Sampson (1968) monastery network (3rd time point).
#' data("sampsonT3")
#'
#' # Run relocation heuristic for relaxed structural balance.
#' res <- rhrsbt(A = sampsonT3, C = 3, TLIMIT = 1)
#'
#'# See the results.
#'res
#' @author Michael Brusco
#' @references
#' Brusco, M. J., Doreian, P., & Steinley, D. (2019). Deterministic blockmodeling of signed and two-mode networks: a tutorial with psychological examples. \emph{British Journal of Mathematical and Statistical Psychology}.
#'
#'Doreian, P., & Mrvar, A. (2009). Partitioning signed social networks. \emph{Social Networks}, 31, 1-11. http://dx.doi.org/10.1016/j.socnet.2008.08.001
#' @export
#' @useDynLib dBlockmodeling, .registration = TRUE

rhrsbt = function(A,C,TLIMIT) {
	N = dim(A)[1]
	OBJVAL = 0
  NREPS = 0
	EBEST <- matrix(0, nrow = N, ncol = 1)
	res =.Fortran("rhrsbtf",as.integer(N),as.integer(C),as.double(TLIMIT),as.double(OBJVAL),as.integer(A),as.integer(EBEST),as.integer(NREPS))
	P <- res[[6]]
	obj <- res[[4]]
	restarts <- res[[7]]
	return(list(P=P, obj=obj, restarts=restarts))
}
