###########################################################
### Script to explore the height component of the joint NDVI + LiDAR height
### model (Model_height.R), fitted by fit_model_height2026.R
### Based on plotting_landsat.R
###########################################################

libs <- c("tidyverse", "coda", "rjags", "cowplot", "ggridges", "viridis")
lapply(libs, require, character.only = TRUE)

#file locations and names (dates are appended by fit_model_height2026.R)
mdatwd <- "data/"
mname <- "peninsulaHeight2026"
fitdate <- "2026-10-05" #date the model was fitted

foutput <- paste0(mdatwd, mname, fitdate, "_modeloutput.Rdata")
envdata <- paste0(mdatwd, mname, "_envdata.Rdata", fitdate)
inputdata <- paste0(mdatwd, mname, fitdate, "_inputdata_small.Rdata")

load(foutput)    # m: mcmc.list
load(envdata)    # env: data.frame (intercept, covariates, UIJ, jag_id)
load(inputdata)  # tdat, hdat

dir.create("figures", showWarnings = FALSE)

#covariate names, in the order used for the beta coefficients (see fit_model_height2026.R)
envars <- c("elevation", "slope", "tpi", "northness", "eastness")
covnames <- c("intercept", envars)

#posterior draws (all chains combined)
draws <- as.matrix(m)

#summarise draws of an indexed parameter (e.g. "Hmax[3]") as one row per index
summ_par <- function(par) {
  cols <- grep(paste0("^", par, "\\["), colnames(draws), value = TRUE)
  tibble(par = cols,
         jag_id = as.integer(str_extract(cols, "(?<=\\[)\\d+(?=\\])")),
         mean = colMeans(draws[, cols, drop = FALSE]),
         sd = apply(draws[, cols, drop = FALSE], 2, sd),
         lower = apply(draws[, cols, drop = FALSE], 2, quantile, 0.025),
         upper = apply(draws[, cols, drop = FALSE], 2, quantile, 0.975))
}

site_ids <- env |> select(UIJ, jag_id)

###########################################################
### Convergence diagnostics for key height parameters
###########################################################

hpars <- c("delta.mu", "delta.sigma", "H0.mu", "H0.sigma", "Hmax.sigma",
           "sigma.height", "sigma.ndvi", paste0("Hmax.beta[", seq_along(covnames), "]"))
hpars <- intersect(hpars, colnames(draws))

gd <- gelman.diag(m[, hpars], multivariate = FALSE)$psrf
print(round(gd, 3))
print(summary(m[, hpars]))

pdf("figures/height_traceplots.pdf", width = 10, height = 12)
plot(m[, hpars])
dev.off()

###########################################################
### Regression coefficients on asymptotic height (Hmax.beta)
###########################################################

hb <- draws[, grep("^Hmax\\.beta\\[", colnames(draws)), drop = FALSE]
colnames(hb) <- covnames[as.integer(str_extract(colnames(hb), "\\d+"))]

hbeta <- as_tibble(hb) |>
  pivot_longer(everything(), names_to = "covariate", values_to = "value") |>
  filter(covariate != "intercept")

p_beta <- ggplot(hbeta, aes(x = covariate, y = value)) +
  geom_boxplot(outlier.size = 0.5) +
  geom_hline(yintercept = 0, colour = "gray50") +
  coord_flip() +
  labs(y = "Effect on log(Hmax) (standardised covariates)", x = NULL) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggsave("figures/height_Hmax_coefficients.png", p_beta,
       width = 12, height = 10, units = "cm", dpi = 300)

###########################################################
### Relative recovery timescale (height vs NDVI)
###########################################################

# delta.mu = 0: same timescale; > 0: height recovers more slowly
p_delta <- ggplot(tibble(delta = draws[, "delta.mu"]), aes(delta)) +
  geom_density(fill = "grey70") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = expression(delta[mu] ~ "= mean log(" * lambda[height] / lambda[NDVI] * ")"),
       y = "Density") +
  theme_bw()

ggsave("figures/height_delta_posterior.png", p_delta,
       width = 12, height = 8, units = "cm", dpi = 300)

cat("P(height recovers more slowly than NDVI) =",
    mean(draws[, "delta.mu"] > 0), "\n")

###########################################################
### Site-level parameters: H0, Hmax, lambda.height, lambda.ndvi
###########################################################

hpar <- bind_rows(
  summ_par("H0") |> mutate(parameter = "H0 (m)"),
  summ_par("Hmax") |> mutate(parameter = "Hmax (m)"),
  summ_par("lambda.height") |> mutate(parameter = "lambda height (yr)"),
  summ_par("lambda.ndvi") |> mutate(parameter = "lambda NDVI (yr)")
) |>
  left_join(site_ids, by = "jag_id")

p_par <- ggplot(hpar, aes(x = mean)) +
  geom_histogram(bins = 20, fill = "grey60", colour = "white") +
  facet_wrap(~parameter, scales = "free") +
  labs(x = "Posterior mean", y = "Number of sites") +
  theme_bw()

ggsave("figures/height_site_parameters.png", p_par,
       width = 18, height = 14, units = "cm", dpi = 300)

# Height vs NDVI recovery timescale per site
lam <- hpar |>
  filter(str_detect(parameter, "lambda")) |>
  select(jag_id, UIJ, parameter, mean) |>
  pivot_wider(names_from = parameter, values_from = mean)

p_lam <- ggplot(lam, aes(`lambda NDVI (yr)`, `lambda height (yr)`)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_point() +
  scale_x_log10() + scale_y_log10() +
  theme_bw()

ggsave("figures/height_lambda_height_vs_ndvi.png", p_lam,
       width = 12, height = 12, units = "cm", dpi = 300)

# Site-level Hmax against covariates
hmax_env <- summ_par("Hmax") |>
  left_join(env |> select(jag_id, all_of(envars)), by = "jag_id") |>
  pivot_longer(all_of(envars), names_to = "covariate", values_to = "value")

p_hmax_env <- ggplot(hmax_env, aes(value, mean)) +
  geom_pointrange(aes(ymin = lower, ymax = upper), size = 0.2, colour = "grey40") +
  geom_smooth(method = "lm", se = FALSE, colour = "firebrick") +
  facet_wrap(~covariate, scales = "free_x") +
  labs(x = "Scaled covariate", y = "Hmax (m)") +
  theme_bw()

ggsave("figures/height_Hmax_vs_covariates.png", p_hmax_env,
       width = 20, height = 14, units = "cm", dpi = 300)

###########################################################
### Posterior height trajectories with observed LiDAR heights
###########################################################

# pred.age must match the vector used in fit_model_height2026.R
pred.age <- seq(0, 25, by = 0.5)

pcols <- grep("^height\\.pred\\[", colnames(draws), value = TRUE)
idx <- str_match(pcols, "\\[(\\d+),(\\d+)\\]")
hp <- tibble(col = pcols, jag_id = as.integer(idx[, 2]), t = as.integer(idx[, 3]),
             age = pred.age[as.integer(idx[, 3])],
             mean = colMeans(draws[, pcols]),
             lower = apply(draws[, pcols], 2, quantile, 0.025),
             upper = apply(draws[, pcols], 2, quantile, 0.975)) |>
  left_join(site_ids, by = "jag_id")

# All sites
p_all <- ggplot() +
  geom_line(data = hp, aes(age, mean, group = UIJ), colour = "grey50", alpha = 0.5) +
  geom_point(data = hdat, aes(DA, Height, colour = factor(Year)), size = 1.5) +
  scale_colour_viridis_d(name = "LiDAR year") +
  labs(x = "Years since fire", y = "Height (m)") +
  theme_bw()

ggsave("figures/height_trajectories_all.png", p_all,
       width = 16, height = 12, units = "cm", dpi = 300)

# Individual sites with 95% credible intervals
sites <- unique(hdat$UIJ)
nshow <- min(12, length(sites))
show <- sites[round(seq(1, length(sites), length.out = nshow))]

p_sites <- ggplot() +
  geom_ribbon(data = filter(hp, UIJ %in% show),
              aes(age, ymin = lower, ymax = upper), alpha = 0.3) +
  geom_line(data = filter(hp, UIJ %in% show), aes(age, mean)) +
  geom_point(data = filter(hdat, UIJ %in% show), aes(DA, Height),
             colour = "firebrick", size = 1.5) +
  facet_wrap(~UIJ) +
  labs(x = "Years since fire", y = "Height (m)") +
  theme_bw()

ggsave("figures/height_trajectories_sites.png", p_sites,
       width = 22, height = 16, units = "cm", dpi = 300)

###########################################################
### Model fit: observed vs predicted height at LiDAR observations
###########################################################

hp_par <- function(p) summ_par(p) |> select(jag_id, mean) |> rename(!!p := mean)

fit <- hdat |>
  left_join(hp_par("H0"), by = "jag_id") |>
  left_join(hp_par("Hmax"), by = "jag_id") |>
  left_join(hp_par("lambda.height"), by = "jag_id") |>
  mutate(Predicted = H0 + (Hmax - H0) * (1 - exp(-DA / lambda.height)),
         Residual = Height - Predicted)

p_fit <- ggplot(fit, aes(Predicted, Height)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_point(aes(colour = factor(Year))) +
  scale_colour_viridis_d(name = "LiDAR year") +
  labs(x = "Predicted height (m)", y = "Observed height (m)") +
  theme_bw()

p_res <- ggplot(fit, aes(DA, Residual)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(aes(colour = factor(Year))) +
  scale_colour_viridis_d(name = "LiDAR year") +
  labs(x = "Years since fire", y = "Residual (m)") +
  theme_bw()

ggsave("figures/height_model_fit.png",
       plot_grid(p_fit + theme(legend.position = "none"), p_res, rel_widths = c(1, 1.3)),
       width = 22, height = 10, units = "cm", dpi = 300)

cat("RMSE (m):", sqrt(mean(fit$Residual^2)), "\n")
cat("R2 (obs vs pred):", cor(fit$Height, fit$Predicted)^2, "\n")

###########################################################
### Predicted height at selected ages across sites
###########################################################

p_ridge <- hp |>
  filter(age %in% c(2, 5, 10, 15, 20, 25)) |>
  ggplot(aes(x = mean, y = factor(age), fill = factor(age))) +
  geom_density_ridges(alpha = 0.7, show.legend = FALSE) +
  scale_fill_viridis_d() +
  labs(x = "Predicted height (m)", y = "Years since fire") +
  theme_bw()

ggsave("figures/height_predicted_by_age_ridges.png", p_ridge,
       width = 14, height = 12, units = "cm", dpi = 300)
