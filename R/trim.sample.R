#' Trim a vector of numbers
#'
#' Symmetrically trim a vector of numbers after sorting it.
#'
#' @details a minimum of 1 observation trimmed from each end (enforced by `max(1, ...)`).
#'
#' @param x A numeric vector.
#' @param trim A number, the fraction of observations to be trimmed from each end of `x`.
#' @return A sorted, trimmed copy of `x`.
#' @export trim.sample
#' @examples
#' trim.sample(
#'     x = 1:20,
#'     trim = 0.01
#'     );
trim.sample <- function(x, trim = 0.01) {
    x <- sort(x);
    trim.sample.number <- length(x) * trim;
    trim.sample.number.integer <- round(trim.sample.number);
    trim.count <- max(1, trim.sample.number.integer);
    patient.trim.value <- (trim.count + 1):(length(x) - trim.count);
    x[patient.trim.value];
    }
