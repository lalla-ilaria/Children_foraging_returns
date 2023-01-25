
if (dir.exists("plots")) unlink("plots", recursive = TRUE)
if (dir.exists("3_outcomes/model_fit")) unlink("3_outcomes/model_fit", recursive = TRUE)

model_binaries <- list.files("models", full.names = TRUE)
model_binaries <- model_binaries[-grep("\\.stan$", model_binaries)]
file.remove(model_binaries)

dir.create("plots")
dir.create("plots/validate_model")
dir.create("3_outcomes/model_fit")
