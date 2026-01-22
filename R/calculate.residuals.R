#' Calculate residuals
#'
#' Calculate residuals between quantiles of the input and quantiles of a Gaussian mixture model fitted to the data.
#'
#' @param x A numeric vector.
#'
#' @return A numeric vector of the same length as `x`.  Names are not retained.
#' @export
#' @examples
#' # Generate fake data.
#' set.seed(1234);
#' x <- rgamma(
#'     n = 20,
#'     shape = 2,
#'     scale = 2
#'     );
#' names(x) <- paste(
#'     'Sample',
#'     seq_along(x),
#'     sep = '.'
#'     );
#' calculate.residuals(
#'     x = x
#'     );
#'
calculate.residuals <- function(x) {
    # Apply 1% trimming to the data for parameter estimation
    x.trim <- trim.sample(
        x = x,
        trim = 0.01
        );

    # Generate the percentiles at which theoretical quantiles will be
    # computed for the fitted distribution.
    p <- stats::ppoints(
        n = x
        );
    # Get the quantiles of `x`.
    x.quantiles <- stats::quantile(
        x = x,
        probs = p
        );
    # Fit Gaussian mixture model to trimmed data
    mix <- mclust::densityMclust(
        x.trim,
        G = 1:10,
        modelNames = "V"
        );
    props <- mix$parameters$pro;  # mixture proportions
    mu <- mix$parameters$mean;    # component means
    sig2 <- mix$parameters$variance$sigmasq;  # component variances
    # Define mixture CDF
    pmix <- function(x, props, mu, sig2) {
        result <- 0;
        for (i in 1:length(props)) {
            result <- result + props[i] * stats::pnorm(x, mean = mu[i], sd = sqrt(sig2[i]));
            }
        return(result);
        }
    # Define quantile function
    qmix <- function(p, props, mu, sig2, 
        lower = min(mu) - 5 * max(sqrt(sig2)),
        upper = max(mu) + 5 * max(sqrt(sig2))) {
        sapply(p, function(pi) {
            stats::uniroot(function(x) pmix(x, props, mu, sig2) - pi,
                lower = lower, upper = upper)$root;
                }
            );
        }
    # Calculate theoretical quantiles
    theoretical.quantiles <- qmix(p, props, mu, sig2);
    # Calculate residuals between the quantiles of `x` and the
    # quantiles of the fitted distribution.
    residuals <- x.quantiles - theoretical.quantiles;
    residuals;
    }
