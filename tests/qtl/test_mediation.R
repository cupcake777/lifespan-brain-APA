# Synthetic test: mediation core of 04_qtl/eqtl/04-1_mediate.R (ols_fit, ols_coef, mediate_pair).
# Only these functions are parsed and evaluated, since the script needs HPC inputs.
# Usage: Rscript test_mediation.R <path/to/04-1_mediate.R>
# --- helpers ---
fails <- character(0)
check <- function(ok, msg) if (!isTRUE(ok)) fails <<- c(fails, msg)

# --- load functions ---
args <- commandArgs(trailingOnly = TRUE)
exprs <- parse(args[1], keep.source = FALSE)
env <- new.env()
# constants from the top of the script; N_BOOT is 300 here (script value 2000) to keep the test fast
assign("MIN_N", 50L, envir = env)
assign("MIN_MINOR", 10L, envir = env)
assign("N_BOOT", 300L, envir = env)
want <- c("ols_fit", "ols_coef", "mediate_pair")
found <- character(0)
for (e in exprs) {
  if (is.call(e) && (identical(e[[1]], as.name("<-")) || identical(e[[1]], as.name("="))) &&
      is.name(e[[2]]) && as.character(e[[2]]) %in% want) {
    eval(e, env)
    found <- c(found, as.character(e[[2]]))
  }
}
stopifnot(all(want %in% found))

# --- planted mediation: x -> m (a = 0.5), m -> y (b = 0.6), no direct effect ---
set.seed(20261019)
n <- 400
x <- rbinom(n, 2, 0.3)
cv <- matrix(rnorm(2 * n), ncol = 2)
m <- 0.5 * x + 0.3 * cv[, 1] + rnorm(n)
y <- 0.6 * m + 0.2 * cv[, 2] + rnorm(n)
r1 <- env$mediate_pair(x, m, y, cv)
check(identical(r1$status, "ok"), paste("status", r1$status))
check(abs(r1$a - 0.5) < 0.1, sprintf("a=%.3f", r1$a))
check(abs(r1$b - 0.6) < 0.1, sprintf("b=%.3f", r1$b))
check(r1$p_c_prime > 0.01, sprintf("direct effect p=%.3g should be non-significant", r1$p_c_prime))
check(r1$ab_fwd_lo > 0 && r1$p_ab_fwd < 1e-4, sprintf("ab CI [%.3f, %.3f], p=%.3g",
      r1$ab_fwd_lo, r1$ab_fwd_hi, r1$p_ab_fwd))
check(!is.na(r1$prop_med) && r1$prop_med > 0.6 && r1$prop_med < 1.4, sprintf("prop_med=%.3f", r1$prop_med))

# --- null: m unrelated to x, y depends on x directly (no mediation through m) ---
set.seed(7)
m0 <- rnorm(n) + 0.3 * cv[, 1]
y0 <- 0.5 * x + 0.2 * cv[, 2] + rnorm(n)
r0 <- env$mediate_pair(x, m0, y0, cv)
check(identical(r0$status, "ok"), paste("null status", r0$status))
check(r0$ab_fwd_lo < 0 && r0$ab_fwd_hi > 0, sprintf("null ab CI [%.3f, %.3f] should include 0",
      r0$ab_fwd_lo, r0$ab_fwd_hi))
check(r0$p_ab_fwd > 0.01, sprintf("null ab p=%.3g", r0$p_ab_fwd))

# --- report ---
if (length(fails)) {
  cat("FAIL", paste(fails, collapse = " | "), "\n")
  quit(status = 1)
}
cat(sprintf("PASS mediation_core: planted a=%.2f b=%.2f recovered (ab CI excludes 0, direct effect p=%.2g); null CI includes 0\n",
            r1$a, r1$b, r1$p_c_prime))
