#' Compute quantities for outlier detection
#'
#' Compute quantities for use in the detection of outliers.  Specifically, compute z-scores based on the mean / standard deviation, the trimmed mean / trimmed standard deviation, or the median / median absolute deviation, or the cluster assignment from k-means with two clusters.
#'
#' @param x A numeric vector.
#' @param method A string indicating the quantities to be computed.  Possible values are
#' * `'mean'` : z-scores based on mean and standard deviation or trimmed mean and trimmed standard deviation if `trim > 0`,
#' * `'median'` : z-scores based on median and median absolute deviation, or
#' * `'kmeans'` : cluster assignment from k-means with two clusters.
#' The default is z-scores based on the mean and standard deviation.
#' @param trim A number, the fraction of observations to be trimmed from each end of `x`.  Default is no trimming.
#' @param nstart A number, for k-means clustering, the number of random initial centers for the clusters.  Default is `1`.  See [stats::kmeans()] for further information.
#' @param exclude.zero A logical, whether zeros should be excluded (`TRUE`) or not excluded (`FALSE`, the default) from computations.  For `method = 'mean'` and `method = 'median'`, this means zeros will not be included in computing the summary statistics; for `method = 'kmeans'`, this means zeros will be placed in their own cluster, coded `0`.
#'
#' @return A numeric vector the same size as `x` whose values are the requested quantities computed on the corresponding elements of `x`.
#' @export
#' @examples
#' # Generate fake data.
#' set.seed(1234);
#' x <- rgamma(
#'     n = 20,
#'     shape = 2,
#'     scale = 2
#'     );
#' # Add missing values and zeros for demonstration.  Missing values are
#' # ignored, and zeros can be ignored with `exclude.zeros = TRUE`.
#' x[1:5] <- NA;
#' x[6:10] <- 0;
#'
#' # Compute z-scores based on mean and standard deviation.
#' quantify.outliers(
#'     x = x,
#'     method = 'mean',
#'     trim = 0
#'     );
#' # Exclude zeros from the calculation of the mean and standard
#' # deviation.
#' quantify.outliers(
#'     x = x,
#'     method = 'mean',
#'     trim = 0,
#'     exclude.zero = TRUE
#'     );
#'
#' # Compute z-scores based on the 5% trimmed mean and 5% trimmed
#' # standard deviation.
#' quantify.outliers(
#'     x = x,
#'     method = 'mean',
#'     trim = 0.05
#'     );
#'
#' # Compute z-scores based on the median and median absolute deviation.
#' quantify.outliers(
#'     x = x,
#'     method = 'median'
#'     );
#'
#' # Compute cluster assignments using k-means with k = 2.
#' quantify.outliers(
#'     x = x,
#'     method = 'kmeans'
#'     );
#' # Try different initial cluster assignments.
#' quantify.outliers(
#'     x = x,
#'     method = 'kmeans',
#'     nstart = 10
#'     );
#' # Assign zeros to their own cluster.
#' quantify.outliers(
#'     x = x,
#'     method = 'kmeans',
#'     exclude.zero = TRUE
#'     );
quantify.outliers <- function(x, method = 'mean', trim = 0, nstart = 1, exclude.zero = FALSE) {
    x.na <- stats::na.omit(x);
    if ('median' == method) {
        if (exclude.zero) {
            x.nonzero <- x.na[0 != x.na];
            data.median <- stats::median(x.nonzero);
            data.mad <- stats::mad(x.nonzero);
            }
        else {
            data.median <- stats::median(x.na);
            data.mad <- stats::mad(x.na);
            }
        if (0 == data.mad || is.na(data.mad) || is.na(data.median)) {
            return(rep(NA, length(x)));
            }
        result.na <- (x.na - data.median) / data.mad;
        x[which(!is.na(x))] <- result.na;
        return(x);
        }
    else if ('kmeans' == method) {
        if (exclude.zero) {
            if (1 == length(unique(x.na))) {
                kmeans.matrix <- rep(NA, length(x.na));
                names(kmeans.matrix) <- names(x.na);
                }
            else {
                data.order <- sort(x.na, decreasing = TRUE);
                non.zero <- data.order[data.order > 0];
                if (length(unique(non.zero)) <= 2) {
                    na.matrix <- rep(NA, length(non.zero));
                    cluster.zero <- c(na.matrix, rep(0, length(x.na[0 == x.na])));
                    kmeans.matrix <- cluster.zero[match(x.na, data.order)];
                    names(kmeans.matrix) <- names(x.na);
                    }
                else {
                    kmeans <- stats::kmeans(
                        x = non.zero,
                        centers = 2,
                        nstart = nstart
                        );
                    cluster <- kmeans$cluster;
                    cluster.zero <- c(cluster, rep(0, length(x[0 == x])));
                    kmeans.matrix <- cluster.zero[match(x.na, data.order)];
                    names(kmeans.matrix) <- names(x.na);
                    }
                }
            }

        else {
            if (1 == length(unique(x.na))) {
                kmeans.matrix <- rep(NA, length(x.na));
                names(kmeans.matrix) <- names(x.na);
                }
            else {
                kmeans <- stats::kmeans(
                    x = x.na,
                    centers = 2,
                    nstart = nstart
                    );
                cluster <- kmeans$cluster;
                kmeans.matrix <- cluster;
                names(kmeans.matrix) <- names(x.na);
                }
            }
        result.na <- kmeans.matrix;
        x[which(!is.na(x))] <- result.na;
        return(x);
        }
    else if ('mean' == method) {
        gene.order <- x.na[order(x.na, decreasing = TRUE)];
        if (exclude.zero) {
            gene.order.nonzero <- gene.order[0 != gene.order];
            if (trim == 0) {
                data.mean <- mean(gene.order.nonzero);
                data.sd <- stats::sd(gene.order.nonzero);
                    }
            else {
                trim.sample.number <- length(gene.order.nonzero) * trim;
                trim.sample.number.integer <- round(trim.sample.number);
                trim.count <- max(1, trim.sample.number.integer);
                trimmed.data <- gene.order.nonzero[(trim.count + 1):(length(gene.order.nonzero) - trim.count)];
                data.mean <- mean(trimmed.data);
                data.sd <- stats::sd(trimmed.data);
            }
        }
        else {
            if (trim == 0) {
                data.mean <- mean(gene.order);
                data.sd <- stats::sd(gene.order);
                }
            else {
                trim.sample.number <- length(x.na) * trim;
                trim.sample.number.integer <- round(trim.sample.number);
                trim.count <- max(1, trim.sample.number.integer);
                trimmed.data <- gene.order[(trim.count + 1):(length(gene.order) - trim.count)];
                data.mean <- mean(trimmed.data);
                data.sd <- stats::sd(trimmed.data);
                }
            }
        if (0 == data.sd || is.na(data.sd) || is.na(data.mean)) {
            return(rep(NA, length(x)));
            }
        result.na <- (x.na - data.mean) / data.sd;
        x[which(!is.na(x))] <- result.na;
        return(x);
        }
    }

#' Cosine similarity
#'
#' Compute cosine similarity for detection of outliers.  Generate theoretical quantiles based on the Gaussian mixture model of the data, and compute cosine similarity between a point made up of the largest observed quantile and the largest theoretical quantile and a point on the line y = x.
#' .
#' @param x A numeric vector.
#'
#' @return A number.
#' @export
#' @examples
#' # Generate fake data.
#' set.seed(1234);
#' x <- rgamma(
#'     n = 20,
#'     shape = 2,
#'     scale = 2
#'     );
#' outlier.detection.cosine(
#'     x = x
#'     );
outlier.detection.cosine <- function(x) {
    # Apply 1% trimming to the data
    x.trim <- trim.sample(
        x = x,
        trim = 0.01
        );
    # Generate the percentiles at which theoretical quantiles will be
    # computed for the optimal distribution of the data.
    p <- stats::ppoints(
        n = x
        );
    observed.quantiles <- stats::quantile(
        x = x,
        probs = p
        );
    # Fit Gaussian mixture model to trimmed data
    mix <- mclust::densityMclust(
        x.trim,
        G = 1:10,
        modelNames = "V"
        );
    # Extract GMM parameters
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
    # Calculate the cosine similarity
    cosine.similarity <- lsa::cosine(
        x = c(theoretical.quantiles[length(p)], observed.quantiles[length(p)]),
        y = c(1, 1)
        );
    cosine.similarity[1, 1];
    }
