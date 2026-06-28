PYTHON_PATH <- Sys.getenv("CANOPY_PYTHON",
  unset = "/home/s38leste_hpc/.conda/envs/canopy_rgee/bin/python3.12")
reticulate::use_python(PYTHON_PATH, required = TRUE)

source("scripts/patches.R")

pkgs <- c("rgee", "readr", "mcera5", "microclimf", "microclimdata",
          "terra", "luna", "reticulate", "ecmwfr", "parallel")

for (p in pkgs) {
  result <- tryCatch({ library(p, character.only = TRUE); "OK" },
                     error = function(e) paste("FAILED:", e$message))
  cat(sprintf("%-20s %s\n", p, result))
}

cat("\n--- reticulate Python ---\n")

cat(sprintf("ee available: %s\n", reticulate::py_module_available("ee")))
tryCatch(reticulate::import("ee"), error = function(e) cat("ee import error:", e$message, "\n"))

cat("\n--- terra GDAL ---\n")
cat(sprintf("terra version: %s\n", as.character(packageVersion("terra"))))
cat(sprintf("GDAL version:  %s\n", terra::gdal()))

cat("\nAll done.\n")