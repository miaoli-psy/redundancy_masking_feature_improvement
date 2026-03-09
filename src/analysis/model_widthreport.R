library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(lme4)
library(lmerTest)
library(emmeans)


setwd("D:/OneDrive/projects/RM_feature_improvement/src/analysis/")
data <- readr::read_csv("D:/OneDrive/projects/RM_feature_improvement/data/processed_exp1.csv")
# data <- data %>%
#   filter(correct_num == 3)


data <- data %>%
  mutate(
    Participant = as.character(subID),
    N_act     = correct_num,
    W_act     = correct_width,
    G_act     = actual_edge_to_edge_spacing,       
    N_rep     = response_num,
    W_rep     = response_width,
    G_rep     = response_edge_to_edge_spacing,
    RM        = as.integer(N_rep < N_act),
    RM_f      = factor(RM, levels = c(0, 1), labels = c("noRM", "RM")),
    W_dev     = W_rep - W_act,
    Extent    = N_act * W_act + (N_act - 1) * G_act   # physical extent
  )

# check data
nrow(data)
length(unique(data$Participant))
round(mean(data$RM) * 100, 1) # percentage RM trials


# --- M0: Null ---
# Width_pred = W_actual
predict_independent <- function(W_act, ...) {
  W_act
}

# --- M1: fixed extent ---
# Extent is perfectly preserved.
# W_predicted = (Extent - (N_rep - 1) × G_rep) / N_rep
# No parameters

predict_fixed_extent <- function(Extent, N_rep, G_rep, ...) {
  (Extent - (N_rep - 1) * G_rep) / N_rep
}


# ---- Model 2: Compressed Extent ----
# Extent compresses in RM trials.
# Extent_perceived = Extent × (N_rep / N_act)^α
# W_predicted = (Extent_perceived - (N_rep - 1) × G_rep) / N_rep
# parameter: α
predict_compressed_extent <- function(Extent, N_rep, N_act, G_rep, alpha, ...) {
  ratio <- N_rep / N_act
  Extent_perceived <- Extent * (ratio ^ alpha)
  (Extent_perceived - (N_rep - 1) * G_rep) / N_rep
}

# ---- Model 3: Adaptive Integration ----
# Width is a blend of actual width and the compressed-extent prediction.
# The blend weight differs for RM and non-RM trials.

# W_predicted = (1 - c) × W_actual + c × Compressed Extent_prediction
#   where c = c_noRM on non-RM trials
#         c = c_RM on RM trials
# parameters: α, c_noRM, c_RM

predict_adaptive <- function(W_act, Extent, N_rep, N_act, G_rep, RM,
                             alpha, c_noRM, c_RM, ...) {
  extent_pred <- predict_compressed_extent(Extent, N_rep, N_act, G_rep, alpha)
  c_val <- ifelse(RM == 1, c_RM, c_noRM)
  (1 - c_val) * W_act + c_val * extent_pred
}


# fitting

compute_fit_stats <- function(sse, k, n) {
  sigma <- sqrt(sse / n)
  ll <- -n / 2 * log(2 * pi) - n * log(max(sigma, 1e-10)) - n / 2
  list(LL = ll,
       AIC = -2 * ll + 2 * (k + 1),
       BIC = -2 * ll + log(n) * (k + 1),
       RMSE = sqrt(sse / n))
}


fit_participant <- function(sub) {
  
  n  <- nrow(sub)
  W  <- sub$W_rep
  results <- list()
  
  # ---- Independent (0 params) ----
  pred <- predict_independent(sub$W_act)
  sse <- sum((W - pred)^2)
  s <- compute_fit_stats(sse, 0, n)
  results[[1]] <- data.frame(
    Model = "Independent", k = 0,
    alpha = NA, c_noRM = NA, c_RM = NA,
    SSE = sse, LL = s$LL, AIC = s$AIC, BIC = s$BIC, RMSE = s$RMSE)
  
  # ---- Fixed Extent (0 params) ----
  pred <- predict_fixed_extent(sub$Extent, sub$N_rep, sub$G_rep)
  sse <- sum((W - pred)^2)
  s <- compute_fit_stats(sse, 0, n)
  results[[2]] <- data.frame(
    Model = "Fixed Extent", k = 0,
    alpha = NA, c_noRM = NA, c_RM = NA,
    SSE = sse, LL = s$LL, AIC = s$AIC, BIC = s$BIC, RMSE = s$RMSE)
  
  # ---- Compressed Extent (1 param: α) ----
  obj <- function(alpha) {
    pred <- predict_compressed_extent(sub$Extent, sub$N_rep, sub$N_act,
                                      sub$G_rep, alpha)
    sum((W - pred)^2)
  }
  fit <- optim(0.5, obj, method = "Brent", lower = -2, upper = 5)
  s <- compute_fit_stats(fit$value, 1, n)
  results[[3]] <- data.frame(
    Model = "Compressed Extent", k = 1,
    alpha = fit$par, c_noRM = NA, c_RM = NA,
    SSE = fit$value, LL = s$LL, AIC = s$AIC, BIC = s$BIC, RMSE = s$RMSE)
  
  # ---- Adaptive Integration (3 params: α, c_noRM, c_RM) ----
  obj <- function(par) {
    pred <- predict_adaptive(sub$W_act, sub$Extent, sub$N_rep, sub$N_act,
                             sub$G_rep, sub$RM, par[1], par[2], par[3])
    sum((W - pred)^2)
  }
  fit <- optim(c(0.5, 0.2, 0.4), obj, method = "L-BFGS-B",
               lower = c(-2, 0, 0), upper = c(5, 1, 1))
  s <- compute_fit_stats(fit$value, 3, n)
  results[[4]] <- data.frame(
    Model = "Adaptive Integration", k = 3,
    alpha = fit$par[1], c_noRM = fit$par[2], c_RM = fit$par[3],
    SSE = fit$value, LL = s$LL, AIC = s$AIC, BIC = s$BIC, RMSE = s$RMSE)
  
  out <- do.call(rbind, results)
  out$Participant <- sub$Participant[1]
  out$n <- n
  out$n_rm <- sum(sub$RM == 1)
  rownames(out) <- NULL
  out
}



all_fits <- do.call(rbind, lapply(split(data, data$Participant), fit_participant))


model_order <- c("Independent", "Fixed Extent", "Compressed Extent", "Adaptive Integration")
all_fits$Model <- factor(all_fits$Model, levels = model_order)


# model results - model comparision

model_summary <- all_fits %>%
  group_by(Model) %>%
  summarise(
    Parameters = first(k),
    mean_AIC   = mean(AIC),
    mean_BIC   = mean(BIC),
    mean_RMSE  = mean(RMSE),
    .groups    = "drop"
  ) %>%
  arrange(mean_AIC) %>%
  mutate(across(where(is.numeric) & !matches("Parameters"), ~ round(., 3)))

model_summary

best_pp <- all_fits %>%
  group_by(Participant) %>%
  slice_min(AIC, n = 1, with_ties = FALSE) %>%
  ungroup()

print(table(best_pp$Model))


# results: likelihood ratio test
run_lrt <- function(model_null, model_alt, df_diff) {
  ll_0 <- all_fits %>% filter(Model == model_null) %>% arrange(Participant) %>% pull(LL)
  ll_1 <- all_fits %>% filter(Model == model_alt) %>% arrange(Participant) %>% pull(LL)
  chi2 <- -2 * (ll_0 - ll_1)
  p_vals <- pchisq(chi2, df = df_diff, lower.tail = FALSE)
  
  # Group-level combined test
  total_chi2 <- sum(chi2)
  total_df <- length(chi2) * df_diff
  group_p <- pchisq(total_chi2, df = total_df, lower.tail = FALSE)
  
  invisible(list(chi2 = chi2, p = p_vals, group_p = group_p))
}

lrt1 <- run_lrt("Fixed Extent", "Compressed Extent", 1)
lrt1

lrt2 <- run_lrt("Independent", "Adaptive Integration", 3)
lrt2


# parameters:

params <- all_fits %>%
  filter(Model == "Adaptive Integration") %>%
  select(Participant, alpha, c_noRM, c_RM)

# --- α: extent compression ---
# Extent compresses significantly under RM
mean_alpha = round(mean(params$alpha), 4)
mean_alpha

sd_alpha = round(sd(params$alpha), 4)
sd_alpha

t_alpha <- t.test(params$alpha, mu = 0)
t_alpha

#95% CI of alpha
round(t_alpha$conf.int[1], 3)
round(t_alpha$conf.int[2], 3)

# --- c_RM and c_noRM: extent weight  ---

round(mean(params$c_noRM),4)
round(sd(params$c_noRM), 4)

# when no rm, does the spatial extent influence width perception?
t_cnoRM <- t.test(params$c_noRM, mu = 0)
t_cnoRM

round(mean(params$c_RM),4)
round(sd(params$c_RM), 4)

# on RM trials, does the spatial extent influence width perception?
t_cRM <- t.test(params$c_RM, mu = 0)
t_cRM

# c_RM vs. c_noRM: does RM change how much the extent influence width?
t_cr <- t.test(params$c_RM, mu = 0)
t_cr

# =========plot parameters============
my_theme <- theme(
  axis.title.x = element_text(color = "black", size = 14, face = "bold", margin = margin(t = 10)),
  axis.title.y = element_text(color = "black", size = 14, face = "bold", margin = margin(r = 10)),
  axis.text.x  = element_text(size = 12, face = "bold", color = "black"),
  axis.text.y  = element_text(size = 12, face = "bold", color = "black"),
  axis.line    = element_line(colour = "black", linewidth = 0.8),
  panel.border     = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  strip.text       = element_text(size = 12, face = "bold"),
  legend.title     = element_text(size = 12, face = "bold"),
  legend.text      = element_text(size = 10),
  plot.title       = element_text(size = 16, face = "bold"),
  plot.subtitle    = element_text(size = 12, color = "grey30"),
  panel.spacing    = unit(1.5, "lines")
)

# plot alpha from M3

t_alpha <- t.test(params$alpha, mu = 0)
t_cnoRM <- t.test(params$c_noRM, mu = 0)
t_cRM <- t.test(params$c_RM, mu = 0)

make_star <- function(p) {
  if (p < 0.001) return("***")
  if (p < 0.01) return("**")
  if (p < 0.05) return("*")
  return("ns")
}

param_long <- params %>%
  pivot_longer(cols = c(alpha, c_noRM, c_RM),
               names_to = "Parameter", 
               values_to = "Value") %>%
  mutate(Parameter = factor(Parameter,
                            levels = c("alpha", "c_noRM", "c_RM"),
                            labels = c(
                              paste0("α (compression)\n", make_star(t_alpha$p.value)),
                              paste0("c (no RM)\n", make_star(t_cnoRM$p.value)),
                              paste0("c (RM)\n", make_star(t_cRM$p.value))
                            )))


c_long <- params %>%
  select(Participant, c_noRM, c_RM) %>%
  pivot_longer(
    cols = c(c_noRM, c_RM),
    names_to = "Condition",
    values_to = "c"
  ) %>%
  mutate(
    Condition = factor(
      Condition,
      levels = c("c_noRM", "c_RM"),
      labels = c("no RM", "RM")
    )
  )

c_sum <- c_long %>%
  group_by(Condition) %>%
  summarise(
    mean = mean(c, na.rm = TRUE),
    n = sum(!is.na(c)),
    sd = sd(c, na.rm = TRUE),
    se = sd / sqrt(n),
    ci = qt(0.975, df = n - 1) * se,
    .groups = "drop"
  )


plot_c <- ggplot(c_long, aes(x = Condition, y = c)) +
  geom_line(aes(group = Participant), alpha = 0.2, colour = "grey50") +
  geom_point(aes(colour = Condition), size = 3, alpha = 0.2) +
  geom_line(
    data = c_sum,
    aes(x = Condition,
        y = mean, 
        group = 1),
    inherit.aes = FALSE,
    linewidth = 1,
    colour = "black"
  ) +
  geom_point(
    data = c_sum,
    aes(x = Condition, 
        y = mean,
        colour = Condition),
    inherit.aes = FALSE,
    size = 4,
    shape = 19
  ) +
  geom_errorbar(
    data = c_sum,
    aes(x = Condition, 
        y = mean, 
        ymin = mean - ci, 
        ymax = mean + ci,
        colour = Condition),
    inherit.aes = FALSE,
    width = 0.00,
    alpha = 0.8,
    linewidth = 1
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey70") +
  scale_colour_manual(values = c("no RM" = "#298c8c", "RM" = "#800074")) +
  labs(x = "", y = "Extent integration weight (c)") +
  my_theme

plot_c


alpha_sum <- params %>%
  summarise(
    mean = mean(alpha, na.rm = TRUE),
    n = sum(!is.na(alpha)),
    sd = sd(alpha, na.rm = TRUE),
    se = sd / sqrt(n),
    ci = qt(0.975, df = n - 1) * se
  )

plot_alpha <- ggplot(params, aes(x = "alpha", y = alpha)) +
  geom_point(size = 3, 
             alpha = 0.1, 
             colour = "steelblue") +
  
  stat_summary(fun = mean, 
               geom = "point",
               colour = "steelblue",
               size = 4, 
               shape = 19) +
  geom_errorbar(
    data = alpha_sum,
    aes(x = "alpha", y = mean, ymin = mean - ci, ymax = mean + ci),
    inherit.aes = FALSE,
    width = 0.00,
    colour = "steelblue",
    linewidth = 1
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey70") +
  labs(
    x = expression(alpha),
    y = "Alpha Values",
  ) +
  my_theme

plot_alpha


plot <- plot_c + plot_alpha + plot_layout(widths = c(2, 1))
plot



# results: does RM improve width estiamtion?

# Compute predictions
data <- data %>%
  left_join(params, by = "Participant") %>%
  mutate(
    pred_adaptive = predict_adaptive(W_act, Extent, N_rep, N_act, G_rep, RM,
                                     alpha, c_noRM, c_RM),
    resid_adaptive = W_rep - pred_adaptive # add trial level residual
  )


# ---- figures model comparison ----

null_aic <- all_fits %>%
  filter(Model == "Independent") %>%
  select(Participant, AIC_null = AIC)

delta_df <- all_fits %>%
  filter(Model != "Independent") %>%
  left_join(null_aic, by = "Participant") %>%
  mutate(dAIC = AIC - AIC_null)


fig_model_comp <- ggplot(delta_df, 
               aes(x = Model,
                   y = dAIC, 
                   fill = Model)) +
  geom_hline(yintercept = 0,
             linetype = "dashed", 
             linewidth = 0.8) +
  geom_boxplot(alpha = 0.4,
               width = 0.6) +
  geom_jitter(width = 0.12, 
              alpha = 0.5, 
              size = 2) +
  scale_fill_manual(values = c("Fixed Extent" = "#5e4c5f",
                               "Compressed Extent" = "#999999",
                               "Adaptive Integration" = "#ffbb6f")) +
  labs(x = "", y = "ΔAIC (model − Independent)",
       title = "Model Comparison")+
  my_theme +
  theme(legend.position = "none")+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1))
    
fig_model_comp

# ----------

# Get per-participant parameters for models that need them
params_compressed <- all_fits %>%
  filter(Model == "Compressed Extent") %>%
  select(Participant, alpha_ce = alpha)

params_adaptive <- all_fits %>%
  filter(Model == "Adaptive Integration") %>%
  select(Participant, alpha_ai = alpha, c_noRM = c_noRM, c_RM = c_RM)

# Remove any old parameter columns before joining fresh ones
df_pred <- data %>%
  select(-any_of(c("alpha", "alpha_ce", "alpha_ai", 
                   "c_correct", "c_noRM", "c_RM", 
                   "c_corr", "c_rm", "c1", "c2",
                   "pred_adaptive", "resid_adaptive",
                   "pred_M0", "resid_M0"))) %>%
  left_join(params_compressed, by = "Participant") %>%
  left_join(params_adaptive, by = "Participant")

# Compute predictions step by step (no curly braces inside mutate)
df_pred$Pred_Independent <- df_pred$W_act

df_pred$Pred_Fixed_Extent <- (df_pred$Extent - (df_pred$N_rep - 1) * df_pred$G_rep) /
  df_pred$N_rep

ratio_ce <- df_pred$N_rep / df_pred$N_act
ext_perc_ce <- df_pred$Extent * (ratio_ce ^ df_pred$alpha_ce)
df_pred$Pred_Compressed <- (ext_perc_ce - (df_pred$N_rep - 1) * df_pred$G_rep) /
  df_pred$N_rep

ratio_ai <- df_pred$N_rep / df_pred$N_act
ext_perc_ai <- df_pred$Extent * (ratio_ai ^ df_pred$alpha_ai)
extent_pred_ai <- (ext_perc_ai - (df_pred$N_rep - 1) * df_pred$G_rep) / df_pred$N_rep
c_val <- ifelse(df_pred$RM == 1, df_pred$c_RM, df_pred$c_noRM)
df_pred$Pred_Adaptive <- (1 - c_val) * df_pred$W_act + c_val * extent_pred_ai




model_names <- c("Independent", "Fixed Extent", "Compressed Extent", "Adaptive Integration")
pred_cols <- c("Pred_Independent", "Pred_Fixed_Extent", "Pred_Compressed", "Pred_Adaptive")

stats_list <- lapply(seq_along(model_names), function(i) {
  pred <- df_pred[[pred_cols[i]]]
  obs <- df_pred$W_rep
  
  r_val <- cor(pred, obs, use = "complete.obs")
  r2_val <- r_val^2
  rmse <- sqrt(mean((obs - pred)^2, na.rm = TRUE))
  mae <- mean(abs(obs - pred), na.rm = TRUE)
  
  # Separate for RM and non-RM
  rm_idx <- df_pred$RM == 1
  r_rm <- cor(pred[rm_idx], obs[rm_idx], use = "complete.obs")
  r_nrm <- cor(pred[!rm_idx], obs[!rm_idx], use = "complete.obs")
  rmse_rm <- sqrt(mean((obs[rm_idx] - pred[rm_idx])^2, na.rm = TRUE))
  rmse_nrm <- sqrt(mean((obs[!rm_idx] - pred[!rm_idx])^2, na.rm = TRUE))
  
  data.frame(
    Model = model_names[i],
    r = r_val, R2 = r2_val, RMSE = rmse, MAE = mae,
    r_RM = r_rm, r_noRM = r_nrm,
    RMSE_RM = rmse_rm, RMSE_noRM = rmse_nrm
  )
})
stats_df <- do.call(rbind, stats_list)

# data to plot
df_long <- df_pred %>%
  select(Participant, W_rep, RM_f, all_of(pred_cols)) %>%
  pivot_longer(
    cols = all_of(pred_cols),
    names_to = "Model",
    names_prefix = "Pred_",
    values_to = "Predicted"
  ) %>%
  mutate(
    Model = case_when(
      Model == "Independent"  ~ "Independent",
      Model == "Fixed_Extent" ~ "Fixed Extent",
      Model == "Compressed"   ~ "Compressed Extent",
      Model == "Adaptive"     ~ "Adaptive Integration"
    ),
    Model = factor(Model, levels = model_names)
  )

#  annotation labels (r and RMSE per panel)
anno_rm_df <- stats_df %>%
  select(Model, r_noRM, r_RM, RMSE_noRM, RMSE_RM) %>%
  mutate(
    label_norm = paste0("r = ", round(r_noRM, 3), ", RMSE = ", round(RMSE_noRM, 3)),
    label_rm      = paste0("r = ", round(r_RM, 3), ", RMSE = ", round(RMSE_RM, 3)),
    Model = factor(Model, levels = model_names)
  )

# calculate axis range for consistent scales
all_vals <- c(df_long$Predicted, df_long$W_rep)
axis_min <- min(all_vals, na.rm = TRUE) - 0.05
axis_max <- max(all_vals, na.rm = TRUE) + 0.05


fig_scatter <- ggplot(df_long, aes(x = Predicted, 
                                   y = W_rep, 
                                   colour = RM_f)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(alpha = 0.1, size = 1.8) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1) +
  
  geom_text(data = anno_rm_df,
            aes(x = axis_min + 0.03, y = axis_max - 0.03, label = label_norm),
            hjust = 0, vjust = 1, 
            size = 3.2, colour = "#298c8c", 
            inherit.aes = FALSE) +
  
  geom_text(data = anno_rm_df,
            aes(x = axis_min + 0.03, y = axis_max - 0.20, label = label_rm),
            hjust = 0, vjust = 1, 
            size = 3.2, colour = "#800074", inherit.aes = FALSE) +

  facet_wrap(~ Model, nrow = 2) +
  
  scale_colour_manual(values = c("noRM" = "#298c8c",
                                 "RM" = "#800074")) +

  labs(x = "Model-predicted width (deg)",
       y = "Reported width (deg)",
       colour = "RM status")+
  
  my_theme+
  
  theme(legend.position = "bottom") +
  coord_cartesian(xlim = c(axis_min, axis_max), ylim = c(axis_min, axis_max))

fig_scatter

table_out <- stats_df %>%
  mutate(across(where(is.numeric), ~round(., 4))) %>%
  select(Model,
         `r (all)` = r,
         `r (correct)` = r_noRM,
         `r (RM)` = r_RM,
         `RMSE (all)` = RMSE,
         `RMSE (correct)` = RMSE_noRM,
         `RMSE (RM)` = RMSE_RM)
print(as.data.frame(table_out), row.names = FALSE)
