#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

source("src/utils_runtime.R")

assert_true <- function(cond, msg) {
  if (!isTRUE(cond)) stop(msg, call. = FALSE)
}

should_process_shard_item <- function(item_index, shard_index, shard_count) {
  ((as.integer(item_index) - 1L) %% as.integer(shard_count)) + 1L == as.integer(shard_index)
}

assignments <- lapply(seq_len(3), function(shard_idx) {
  which(vapply(seq_len(10), should_process_shard_item, logical(1), shard_index = shard_idx, shard_count = 3L))
})

assigned_items <- unlist(assignments)
assert_true(identical(sort(assigned_items), seq_len(10)), "shard planner should cover each item once")
assert_true(!any(duplicated(assigned_items)), "shard planner should not duplicate items")
assert_true(identical(assignments[[1]], c(1L, 4L, 7L, 10L)), "shard 1 assignment mismatch")

tmp_dir <- tempfile("easycoloc_flush_")
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
initialize_runtime_tracker(tmp_dir, enabled = TRUE, config_fingerprint = "smoke")
assert_true(!runtime_should_flush_task_state(3L), "first flush interval tick should not flush")
assert_true(!runtime_should_flush_task_state(3L), "second flush interval tick should not flush")
assert_true(runtime_should_flush_task_state(3L), "third flush interval tick should flush")

cat("[SMOKE] shard runtime smoke test passed\n")
