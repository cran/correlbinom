#' Correlated Binomial Probabilities
#'
#' This function reports the likelihoods of 0, 1, ..., n
#' successes given n trials of a binomial with a specified
#' correlation or association between trials and success
#' probability
#' @param rho The level of correlation or association between trials. In the Witt (2014) model, this parameter is the level of correlation between trials. In the Kuk (2004) model, it is the equivalent of one minus gamma from that paper, where a value of zero indicates independence. In both cases, this parameter must fall within the unit interval.
#' @param successprob The likelihood of success in one trial.
#' @param trials The number of trials.
#' @param precision Number of bits of precision. Defaults to 1024.
#' @param model Specify whether the 'kuk' or 'witt' model is to be used for calculation. Defaults to 'witt'.
#' @import Rmpfr methods
#' @export
#' @examples correlbinom(0.5,0.1,5)
#' @examples correlbinom(0.9,0.3,12,256)
#' @examples correlbinom(0.9,0.6,12,model="kuk")
#' @references Kuk, Anthony Y. C., 2004. A litter-based approach to risk assessment in developmental toxicity via a power family of completely monotone functions. \emph{Journal of the Royal Statistical Society, Series C (Applied Statistics)}, \emph{53}(2): 369-86. 
#' @references Witt, Gary, 2014. A simple distribution for the sum of correlated, exchangeable binary data. \emph{Communications in Statistics - Theory and Methods}, \emph{43}(20): 4265-80.

correlbinom <- function(rho,successprob,trials,precision=1024,model="witt"){
	if(is.na(rho) | is.na(successprob) | is.na(trials) | is.na(precision) | is.na(model)){
		stop("NA input detected.")
	}
	if(!is.numeric(rho) | !is.numeric(successprob) | !is.numeric(trials) | !is.numeric(precision)){
		stop("Non-numeric input detected.")
	}
	if(length(rho)!=1 | length(successprob)!=1 | length(trials)!=1 | length(precision)!=1){
		stop("Non-scalar input detected.")
	}
	if(min(rho,successprob)<0 | max(rho,successprob)>1){
		stop("Rho and success probability must both fall within the unit interval.")
	}
	if(trials<=0){
		stop("Number of trials must be positive.")
	}
	if(precision<2){
		stop("Precision must be at least two.")
	}
	if(model!="witt" & model!="kuk"){
		stop("Model must equal 'witt' or 'kuk'.")
	}
	if(as.integer(trials)!=trials | as.integer(precision)!=precision){
		stop("Trials and precision must both be integer-valued.")
	}
	rho <- Rmpfr::mpfr(rho,precision)
	successprob <- Rmpfr::mpfr(successprob,precision)
	if(model=="witt"){
		pows <- methods::new("mpfr",mapply(`^`,1-rho,0:(trials-1)))
		pow.sum <- cumsum(pows[1:(trials-1)])*rho
		cond.prob <- c(Rmpfr::mpfr(1,precision),successprob,pow.sum+pows[2:trials]*successprob)
		prob.diag <- cumprod(cond.prob)
	} else if(model=="kuk"){
		prob.diag <- successprob^((1-rho)*(0:trials))
	}
	gen.elem <- function(i){
		# pascal's triangle row of coefficients w/ alternating signs
		coeffs <- Rmpfr::chooseMpfr(trials+1-i,0:(trials+1-i))*(-1)^(0:(trials+1-i))
		return(sum(prob.diag[i:(trials+1)]*coeffs)*Rmpfr::chooseMpfr(trials,i-1))
	}
	return(Rmpfr::asNumeric(methods::new("mpfr",sapply(1:(trials+1),gen.elem))))
}