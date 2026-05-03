library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)


exp <- "exp1"  # change to "exp2" to run on experiment 2
data <- readr::read_csv(paste0("data/processed_", exp, ".csv"))


data <- data %>%
  mutate(
    Participant = as.character(subID),
    N_act = correct_num,
    N_rep = response_num,
    W_act = correct_width,
    W_rep = response_width,
    RM    = as.integer(N_rep < N_act),
    RM_f  = factor(RM, levels = c(0, 1), labels = c("noRM", "RM")),
    pool  = N_act / N_rep
  )

nrow(data)
length(unique(data$Participant))
round(mean(data$RM) * 100, 1)


# --- M0: Null ---
# W_pred = W_act
predict_null <- function(W_act, ...) W_act

# --- M1: Perfect Pooling ---
# lost bar's width perfectly redistributed
# W_pred = (N_act / N_rep) * W_act
predict_perfect <- function(W_act, pool, ...) pool * W_act

# --- M2: Partial Pooling (no compression) ---
# W_pred = (N_act / N_rep)^gamma * W_act
# gamma = 1 reduces to M1
predict_partial_nc <- function(W_act, pool, gamma, ...) pool^gamma * W_act

# --- M3: Compression + Pooling ---
# W_pred = alpha * (N_act / N_rep) * W_act
# alpha < 1 = compression
predict_compression <- function(W_act, pool, alpha, ...) alpha * pool * W_act

# --- M4: Compression + Partial Pooling ---
# W_pred = alpha * (N_act / N_rep)^gamma * W_act
# alpha = compression, gamma = pooling efficiency
predict_partial <- function(W_act, pool, alpha, gamma, ...) alpha * pool^gamma * W_act

# --- M5: Linear Compression + Partial Pooling ---
# W_pred = (alpha0 + alpha1 * W_act) * (N_act / N_rep)^gamma * W_act
# alpha0, alpha1 = width-dependent compression; gamma = pooling efficiency
predict_linear <- function(W_act, pool, alpha0, alpha1, gamma, ...)
  (alpha0 + alpha1 * W_act) * pool^gamma * W_act


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

  n <- nrow(sub)
  W <- sub$W_rep
  results <- list()

  # M0: Null (0 params)
  pred <- predict_null(sub$W_act)
  sse <- sum((W - pred)^2)
  s <- compute_fit_stats(sse, 0, n)
  results[[1]] <- data.frame(
    Model = "Null", k = 0, alpha = NA, alpha1 = NA, gamma = NA,
    SSE = sse, LL = s$LL, AIC = s$AIC, BIC = s$BIC, RMSE = s$RMSE)

  # M1: Perfect Pooling (0 params)
  pred <- predict_perfect(sub$W_act, sub$pool)
  sse <- sum((W - pred)^2)
  s <- compute_fit_stats(sse, 0, n)
  results[[2]] <- data.frame(
    Model = "Perfect Pooling", k = 0, alpha = NA, alpha1 = NA, gamma = NA,
    SSE = sse, LL = s$LL, AIC = s$AIC, BIC = s$BIC, RMSE = s$RMSE)

  # M2: Partial Pooling no compression (1 param: gamma)
  obj <- function(gamma) {
    pred <- predict_partial_nc(sub$W_act, sub$pool, gamma)
    sum((W - pred)^2)
  }
  fit <- optim(0.5, obj, method = "Brent", lower = -1, upper = 3)
  s <- compute_fit_stats(fit$value, 1, n)
  results[[3]] <- data.frame(
    Model = "Partial Pooling", k = 1,
    alpha = NA, alpha1 = NA, gamma = fit$par,
    SSE = fit$value, LL = s$LL, AIC = s$AIC, BIC = s$BIC, RMSE = s$RMSE)

  # M3: Compression + Pooling (1 param: alpha)
  obj <- function(alpha) {
    pred <- predict_compression(sub$W_act, sub$pool, alpha)
    sum((W - pred)^2)
  }
  fit <- optim(0.95, obj, method = "Brent", lower = 0.3, upper = 1.5)
  s <- compute_fit_stats(fit$value, 1, n)
  results[[4]] <- data.frame(
    Model = "Compression + Pooling", k = 1,
    alpha = fit$par, alpha1 = NA, gamma = NA,
    SSE = fit$value, LL = s$LL, AIC = s$AIC, BIC = s$BIC, RMSE = s$RMSE)

  # M4: Compression + Partial Pooling (2 params: alpha, gamma)
  obj <- function(par) {
    pred <- predict_partial(sub$W_act, sub$pool, par[1], par[2])
    sum((W - pred)^2)
  }
  fit <- optim(c(0.9, 0.5), obj, method = "L-BFGS-B",
               lower = c(0.3, -1), upper = c(1.5, 3))
  s <- compute_fit_stats(fit$value, 2, n)
  results[[5]] <- data.frame(
    Model = "Compression + Partial Pooling", k = 2,
    alpha = fit$par[1], alpha1 = NA, gamma = fit$par[2],
    SSE = fit$value, LL = s$LL, AIC = s$AIC, BIC = s$BIC, RMSE = s$RMSE)

  # M5: Linear Compression + Partial Pooling (3 params: alpha0, alpha1, gamma)
  obj <- function(par) {
    pred <- predict_linear(sub$W_act, sub$pool, par[1], par[2], par[3])
    sum((W - pred)^2)
  }
  fit <- optim(c(0.9, 0, 0.5), obj, method = "L-BFGS-B",
               lower = c(0.1, -3, -1), upper = c(2, 3, 3))
  s <- compute_fit_stats(fit$value, 3, n)
  results[[6]] <- data.frame(
    Model = "Linear Compression + Partial Pooling", k = 3,
    alpha = fit$par[1], alpha1 = fit$par[2], gamma = fit$par[3],
    SSE = fit$value, LL = s$LL, AIC = s$AIC, BIC = s$BIC, RMSE = s$RMSE)

  out <- do.call(rbind, results)
  out$Participant <- sub$Participant[1]
  out$n <- n
  rownames(out) <- NULL
  out
}


all_fits <- do.call(rbind, lapply(split(data, data$Participant), fit_participant))

model_order <- c("Null", "Perfect Pooling", "Partial Pooling",
                 "Compression + Pooling", "Compression + Partial Pooling",
                 "Linear Compression + Partial Pooling")
all_fits$Model <- factor(all_fits$Model, levels = model_order)


# model comparison
model_summary <- all_fits %>%
  group_by(Model) %>%
  summarise(Parameters = first(k),
            mean_AIC = mean(AIC), mean_BIC = mean(BIC),
            mean_RMSE = mean(RMSE), .groups = "drop") %>%
  arrange(mean_AIC) %>%
  mutate(across(where(is.numeric) & !matches("Parameters"), ~ round(., 3)))

model_summary

best_pp <- all_fits %>%
  group_by(Participant) %>%
  slice_min(AIC, n = 1, with_ties = FALSE) %>%
  ungroup()

print(table(best_pp$Model))


# likelihood ratio tests
run_lrt <- function(model_null, model_alt, df_diff) {
  ll_0 <- all_fits %>% filter(Model == model_null) %>% arrange(Participant) %>% pull(LL)
  ll_1 <- all_fits %>% filter(Model == model_alt) %>% arrange(Participant) %>% pull(LL)
  chi2 <- -2 * (ll_0 - ll_1)
  p_vals <- pchisq(chi2, df = df_diff, lower.tail = FALSE)
  total_chi2 <- sum(chi2)
  total_df <- length(chi2) * df_diff
  group_p <- pchisq(total_chi2, df = total_df, lower.tail = FALSE)
  invisible(list(chi2 = chi2, p = p_vals, group_p = group_p))
}

lrt1 <- run_lrt("Compression + Pooling", "Compression + Partial Pooling", 1)
lrt1

lrt2 <- run_lrt("Partial Pooling", "Compression + Partial Pooling", 1)
lrt2

# M5 vs M4: does width-dependent compression improve fit?
lrt3 <- run_lrt("Compression + Partial Pooling",
                "Linear Compression + Partial Pooling", 1)
lrt3


# parameters (M4: Compression + Partial Pooling)
params <- all_fits %>%
  filter(Model == "Compression + Partial Pooling") %>%
  select(Participant, alpha, gamma)

round(mean(params$alpha), 4)
round(sd(params$alpha), 4)
t.test(params$alpha, mu = 1)

round(mean(params$gamma), 4)
round(sd(params$gamma), 4)
t.test(params$gamma, mu = 0)
t.test(params$gamma, mu = 1)


# parameters (M5: Linear Compression + Partial Pooling)
params_m5 <- all_fits %>%
  filter(Model == "Linear Compression + Partial Pooling") %>%
  select(Participant, alpha, alpha1, gamma)

round(mean(params_m5$alpha),  4); round(sd(params_m5$alpha),  4)
round(mean(params_m5$alpha1), 4); round(sd(params_m5$alpha1), 4)
round(mean(params_m5$gamma),  4); round(sd(params_m5$gamma),  4)

t.test(params_m5$alpha,  mu = 1)
t.test(params_m5$alpha1, mu = 0)
t.test(params_m5$gamma,  mu = 0)
t.test(params_m5$gamma,  mu = 1)


# =========plots============

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
  panel.spacing    = unit(1.5, "lines")
)


# delta-AIC
null_aic <- all_fits %>%
  filter(Model == "Null") %>%
  select(Participant, AIC_null = AIC)

delta_df <- all_fits %>%
  filter(Model != "Null") %>%
  left_join(null_aic, by = "Participant") %>%
  mutate(dAIC = AIC - AIC_null)

fig_aic <- ggplot(delta_df, aes(x = Model, y = dAIC, fill = Model)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.8) +
  geom_boxplot(alpha = 0.4, width = 0.6) +
  geom_jitter(width = 0.12, alpha = 0.5, size = 2) +
  labs(x = "", y = expression(Delta * "AIC (model " %~% " Null)")) +
  my_theme +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

fig_aic


# parameters
param_long <- params %>%
  pivot_longer(cols = c(alpha, gamma),
               names_to = "Parameter", values_to = "Value") %>%
  mutate(Parameter = factor(Parameter,
                            levels = c("alpha", "gamma"),
                            labels = c("α (compression)", "γ (pooling)")))

param_sum <- param_long %>%
  group_by(Parameter) %>%
  summarise(mean = mean(Value), n = n(),
            sd = sd(Value), ci = qt(0.975, n - 1) * sd / sqrt(n),
            .groups = "drop")

fig_params <- ggplot(param_long, aes(x = Parameter, y = Value)) +
  geom_point(size = 3, alpha = 0.15, colour = "steelblue") +
  geom_point(data = param_sum, aes(y = mean),
             size = 4, colour = "steelblue") +
  geom_errorbar(data = param_sum,
                aes(y = mean, ymin = mean - ci, ymax = mean + ci),
                width = 0, colour = "steelblue", linewidth = 1) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey70") +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey70") +
  labs(x = "", y = "Parameter value") +
  my_theme

fig_params
