#' Simulate from a null distribution
#'
#' Simulate transcripts from a Gaussian mixture model fitted to the observed data.
#'
#' @param x A numeric vector of transcripts.
#' @param r A numeric vector of residuals calculated for this transcript.
#'
#' @return A numeric vector of the same length as `x`.  Names are not retained.
#' @export simulate.null
#' @examples
#' # Prepare fake data.
#' set.seed(1234);
#' x <- rgamma(
#'     n = 20,
#'     shape = 2,
#'     scale = 2
#'     );
#' names(x) <- paste('Sample', seq_along(x), sep = '.');
#' r <- calculate.residuals(
#'     x = x
#'     );
#' null <- simulate.null(
#'     x = x,
#'     r = r
#'     );
simulate.null <- function(
    x,
    r
    ) {
    #
    # Simulate transcripts
    #
    # Apply 1% trimming to the data
    x.nozero.trim <- trim.sample(x, trim = 0.01);
    # Fit Gaussian mixture model to trimmed data
    mix <- mclust::densityMclust(
        x.nozero.trim,
        G = 1:10,
        modelNames = "V"
        );
    # Set sample size (same as original data)
    n.samples <- length(x);
    # Sample directly from the fitted model
    simulated.null <- mclust::simV(
        modelName = mix$modelName,
        parameters = mix$parameters,
        n = n.samples
        );
    # Extract the simulated values
    simulated.null <- simulated.null[,2];
    #
    # Simulate noise
    #
    r.nozero.trim <- trim.sample(r, trim = 0.01);
    # Fit Gaussian mixture model to residuals
    mix.noise <- mclust::densityMclust(
        r.nozero.trim,
        G = 1:10,
        modelNames = "V"
        );
    # Sample noise from the fitted model
    simulated.noise <- mclust::simV(
        modelName = mix.noise$modelName,
        parameters = mix.noise$parameters,
        n = n.samples
        );
    simulated.noise <- simulated.noise[,2];
    #
    # Add the simulated noise to the simulated transcript.
    abs(simulated.null + simulated.noise)
    }
