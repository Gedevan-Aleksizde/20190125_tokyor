require(conflicted)
require(ggdag)
require(dagitty)
require(latex2exp)
require(skimr)
require(np)
require(plotly)
require(coefplot)
require(plm)
require(gmm)
require(tidyverse)
require(ggthemes)
require(Metrics)
require(SemiPar)
require(mgcv)
require(ranger)
# require(bkmr)
# require(gplm)
require(boot)
require(estprod)
require(patchwork)

conflict_prefer("lag", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("TeX", "latex2exp")
conflict_prefer("spm", "SemiPar")

theme_set(new = theme_bw())
thm <- theme_classic(base_size = 20) + theme(legend.position = "bottom", axis.title.y = element_text(angle = 0, vjust = .5))

##### draw Cobb-Douglas 3D image ####
cobb_douglas <- function(k, l, a, b_k, b_l){
  a * k^b_k * l^b_l
}

# CD
d <- expand_grid(
  k=seq(0, 1, by = .01),
  l=seq(0, 1, by = .01)
) %>% mutate(y=cobb_douglas(k, l, a = 1, b_k = .5, b_l = .5)) %>%
  pivot_wider(id_cols = k, names_from = l, values_from = y) %>% column_to_rownames("k") %>% as.matrix 

draw_cdfunc <- function(d){
  plot_ly(
    colorscale = "Blues",
    contours = list(
      z = list(show = TRUE, start = 0, end = 1, size = 0.05, color = "white")
    )
  ) %>% add_surface(z = ~d) %>%
    hide_colorbar() %>% 
    layout(scene = list(
      xaxis = list(title = "L", showticklabels = T, tickvals = c(0), ticktext=c(0)),
      yaxis = list(title = "K", showticklabels = F),
      zaxis = list(title = "Y", showticklabels = F),
      camera=list(
        eye = list(x=-1.8, y=-1.8, z=.5)
      )
    )) %>%
    config(mathjax = 'cdn')
}

draw_cdfunc(d)
# maybe plotly cannot save image

graph <- dagify(omega_f ~ omega_b,
                inv ~ omega_b,
                k_f ~ k_b + inv + omega_f,
                x ~ k_f + omega_f,
                Exit ~ x,
                Func ~ x + k_f,
                Func ~ l,
                y ~ Func)
coordinates(graph) <- list(
  x = set_names(c(2, 2, 0, 0, 1, 1, 0, 1, 1.5, 2.5),
                c("Exit", "Func", "inv", "k_b", "k_f", "l", "omega_b", "omega_f", "x", "y")),
  y = set_names(
    c(.5, 1.5, 1, 2, 2, 3, 0, 0, 1, 1.5),
    c("Exit", "Func", "inv", "k_b", "k_f", "l", "omega_b", "omega_f", "x", "y")
  )
)
df_graph <- graph %>% tidy_dagitty() %>%
  arrange(name)
df_labels <- tibble(name = df_graph$data$name %>% unique) %>% mutate(
  param = c("Exit", "f(k, l)", "\\mathit{inv}_t", "k_t", "k_{t+1}", "l_{t+1}", "\\omega_t", "\\omega_{t+1}", "x_{t+1}", "y_{t+1}"),
  time = as.factor(c(1, 1, 0, 0, 1, 1, 0, 1, 1, 1))
)
df_graph$data <- df_graph$data %>% left_join(df_labels, by = "name") %>% as.data.frame()
df_graph$data <- df_graph$data %>% mutate(edge_value = if_else(name == "x" & to == "Func", "1", ""))
df_graph$data <- df_graph$data %>% mutate(edge_value = if_else(name == "x" & to == "Exit", "0", edge_value))
df_graph %>%
  ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_dag_point(aes(color = time), show.legend = F) +
  # geom_dag_edges(aes(label = edge_value)) +
  geom_dag_edges() +
  geom_dag_text(parse = T, aes(label = latex2exp::TeX(paste0("$", param, "$"), output = "character"))) +
  annotate(geom = "text", x =  c(1.7, 1.7), y = c(.7, 1.1), label = c("0", "1")) +
  scale_color_colorblind() +
  scale_x_continuous(breaks = c(0, 1), labels = c("t", "t+1")) +
  theme_dag_blank() + theme(
    plot.background = element_rect(fill = "transparent"),
    panel.background = element_rect(fill = "transparent",color = NA),
    axis.line.x = element_line(),
    axis.text.x = element_text(),
    axis.ticks.x = element_line())
ggsave(filename = "doc/img/op_graph.pdf", device = cairo_pdf, width = 16, height = 9, units = "cm")


##### computing #####

# estprod の乱数データはあやしい
df_estprod <- estprod_data %>% ungroup %>%
  rename(y = var1, l = var2, k = var3, inv = var4) %>% group_by(id)
fit_estprod <- lm(y ~ l + poly(k, inv, degree = 4), data = df_estprod)
ggplot(tibble(y = df_estprod$y, p = predict(fit_estprod))) + geom_point(aes(x = y, y = p))
rm(df_estprod, fit_estprod)

#### estimation ####

df_sample <- read_csv("data/df_sample.csv", col_types = cols(i = col_factor(), x = col_integer()))
df_ground_truth <- read_csv("data/df_ground_truth.csv", col_types = cols(i = col_factor(), x = col_integer()))


##### 1st stage #####

# kernel-based partial linear approximation
# time-consuming.

set.seed(42)
system.time({
  bw_1st <- npplregbw(formula = y ~ 1 + l | k + inv, data = df_sample %>% drop_na(y, k, l, inv)) # somehow `na.action` doesn't work so we need `drop_na()`
  fit_1st_np <- npplreg(bw_1st)
  write_rds(fit_1st_np, path = "tmp/fit_1st.rds")
})
summary(fit_1st_np)
# check model fitting
ggplot(df_sample %>% mutate(p = predict(fit_1st_np, newdata = df_sample))) + geom_point(aes(x = y, y = p))
summary(fit_1st_np)

# another trial: simply polynimial approximation (much faster)
set.seed(42)
df_fit_1st_poly <- tibble(q= 1:10) %>%
  mutate(model = map(q, function(q){
    if(q == 1){
      poly_term <- "+ k"
    } else{
      poly_term <-  paste("+ poly(k, inv, degree =", q, ")")
    }
    lm(as.formula(paste("y ~ 1 + l + k", poly_term )), data = df_sample)
  })) %>%
  mutate(
    rmse = map_dbl(model, function(x) sqrt(mean(x$residuals^2))),
    coef = map(model, function(x) enframe(coef(x)) %>% filter(name %in% c("k", "l")))
  ) %>% unnest(coef) %>% pivot_wider(names_from = name, values_from = value)
df_fit_1st_poly 
# we select 5th order, the estimate of beta_l = .161
fit_1st_poly <- df_fit_1st_poly$model[[5]]
ggplot(df_sample %>% mutate(pred = predict(fit_1st_poly, newdata = df_sample)) %>% drop_na(y, pred), aes(x = y, y = pred)) +
  geom_point() + stat_smooth()

# estimate survival anaylsis

std_Brier <- function(actual, prediction){
  sqrt(mse(actual, prediction) / mse(actual, mean(actual)))
}

# peeping
ggplot(df_ground_truth, aes(x = k_lag, y = boot::logit(p_x))) + geom_point() + stat_smooth() +
  thm + labs(x = TeX("$k_t$"), y = TeX("\\mathit{logit}(p_{t+1})"))
ggsave("doc/img/inv_logitp.pdf", device = cairo_pdf)
ggplot(df_ground_truth, aes(x = inv_lag, y = boot::logit(p_x))) + geom_point() + stat_smooth() +
  thm + labs(x = TeX("$\\mathit{inv}_t$"), y = TeX("\\mathit{logit}(p_{t+1})"))
ggsave("doc/img/k_logitp.pdf", device = cairo_pdf)

set.seed(42)
df_fit_survival <-tibble(q = 1:10) %>%
  mutate(model = map(q, function(q){
    if(q == 1){
      poly_term <- "k_lag + inv_lag"
    } else{
      poly_term <-  paste("poly(k_lag, inv_lag, degree =", q, ")")
    }
    glm(as.formula(paste("x ~", poly_term)),
        data = df_sample %>% drop_na(x, k_lag, inv_lag),
        family = binomial(link = "probit"), na.action = na.omit)
  }))
df_fit_survival <- df_fit_survival %>%
  mutate(
    actual = map(1:n(), function(x) drop_na(df_sample, x, k_lag, inv_lag)$x),
    pred_p = map(model, ~fitted(.x)),
    pred_label = map(pred_p, ~.x > .5),
    rate = map_dbl(actual, mean),
    acc = map2_dbl(actual, pred_label, accuracy),
    precision = map2_dbl(actual, pred_label, precision),
    recall = map2_dbl(actual, pred_label, recall),
    f1 = map2_dbl(actual, pred_label, f1),
    log_loss = map2_dbl(actual, pred_p, logLoss),
    neg_RIG = log_loss / (-rate * log(rate) - (1 - rate) * log(1 - rate)) - 1,
    std_Brier = map2_dbl(actual, pred_p, ~std_Brier(.x, .y)),
    cor = map2_dbl(actual, pred_p, ~cor(.x, .y)),
    unique_rate = map_dbl(pred_p, ~n_distinct(.x)/length(.x))
  ) %>% select(-actual, -pred_p, -pred_label)
df_fit_survival
fit_survival <- df_fit_survival$model[[2]] # we select 2nd order.


##### 2nd-stage procedure #####
# 2nd-stage
make_2nd_sample <- function(data, beta_l, phi, p_x){
  data %>%
    mutate(
      bl = beta_l * l,
      y_bl = y - bl,
      phi = phi,
      p_x = p_x
    ) %>%
    arrange(i, t) %>% group_by(i) %>% mutate(
      phi_lag = lag(phi),
    ) %>% ungroup
}

obj_op_2nd <- function(param, data, degree){
  # param = [beta_a, log(beta_k), gamma1, gamma2, ...]
  # note: param[2] constrains more than zero.
  h <- with(data, phi_lag - exp(param[2]) * k_lag)
  poly <- map_dfc(
    1:degree,
    function(d) {reduce(0:d, ~.x + param[2 + sum(1:d + 1) - d - 1 + 1 + .y] * with(data, h^-(d - .y) * p_x^-(.y)), .init = 0)}
  ) %>% rowSums
  # return MSE of the model
  return(mean((with(data, y_bl - param[1] - exp(param[2]) * k) - poly)^2, na.rm = T))
}

get_beta_a <- function(data, beta_k, beta_l){
  mean(with(data, y - beta_k * k -beta_l * l), na.rm = T)
}
get_fitted <- function(data, beta_k, beta_l, beta_a = NULL){
  if(is.null(beta_a)){
    beta_a <- get_beta_a(data, beta_k, beta_l)
  }
  with(data, beta_a + beta_k * k + beta_l * l)
}

plot_iterated_convergence <- function(data, max_degree = 6, size = 20){
  # plot iterated results
  data %>% filter(d <= max_degree, n <= size) %>% mutate(
    beta_a_start = map_dbl(start, ~.x[1]),
    beta_k_start = map_dbl(start, ~.x[2])
  ) %>% ggplot(aes(x = beta_a_start, y = beta_k_start, xend = beta_a, yend = beta_k, color = log(rmse))) +
    geom_segment(arrow = arrow(length = unit(.05, "npc")), size = .5) + geom_point() +
    facet_wrap(~d) + thm + theme(legend.position = "right") +
    labs(x = TeX("$\\beta_A$"), y = TeX("$\\beta_K$"))
}
plot_iterated_errorbar <- function(data){
  # plot summarised result
  skim(data %>% group_by(d)) %>% select(d, skim_variable, numeric.mean, numeric.sd) %>%
    filter(skim_variable %in% c("beta_a", "beta_k")) %>%
    mutate(variable = factor(skim_variable, labels = c(TeX("$\\beta_A$"), TeX("$\\beta_K$")))) %>%
    ggplot(aes(x = d, y =numeric.mean)) + geom_line() +
    geom_errorbar(aes(ymin = numeric.mean - numeric.sd, ymax = numeric.mean + numeric.sd)) +
    facet_wrap(~variable, scales = "free_y", labeller = label_parsed) +
    scale_x_continuous(breaks = 1:8) + thm + theme(axis.title.y = element_text(angle = 90)) + labs(x = TeX("$q$"), y = "estimates")
}


# make 2nd-stage dataset
phi_hat <- predict(fit_1st_np, newdata = df_sample %>% select(y, k, l, inv)) %>% as.numeric - df_sample$l * coef(fit_1st_np)["l"]
df_sample_2nd_np <- make_2nd_sample(
  df_sample, coef(fit_1st_np)["l"], phi_hat,
  predict(fit_survival, type = "response", newdata = df_sample)
)

df_sample_2nd_poly <- make_2nd_sample(
  df_sample, coef(fit_1st_poly)["l"],
  predict(fit_1st_poly, newdata = df_sample) - coef(fit_1st_poly)["l"] * df_sample$l,
  predict(fit_survival, type = "response", newdata = df_sample))

df_ground_truth_2nd <- make_2nd_sample(
  data = df_ground_truth, beta_l = beta_l.,
  phi =  beta_a. + beta_k. * df_ground_truth$k + df_ground_truth$omega,
  p_x = df_ground_truth$p_x
)


# check fitting by peeping the grount truth data 
lm(y_bl - omega ~ k, data = left_join(df_ground_truth_2nd, select(df_ground_truth, i, t, omega))) %>% summary
nlme::gls(y_bl - omega ~ k, data = left_join(df_ground_truth_2nd, select(df_ground_truth, i, t, omega)), na.action = na.omit) %>% summary

# peeping true correlation
df_sample_2nd_poly %>% left_join(select(df_ground_truth, i, t, omega)) %>%
  ggplot(aes(x = phi_lag, y = omega)) + geom_point() + stat_smooth() +
  labs(x = TeX("$\\hat{\\varphi}_{t}$"), y = TeX("$\\omega_{t+1}$")) + thm
ggsave(filename = "doc/img/phi_hat_omega.pdf", device = cairo_pdf)
df_sample_2nd_poly %>% left_join(select(df_ground_truth, i, t, omega)) %>%
  ggplot(aes(x = p_x, y = lag(omega))) + geom_point() + stat_smooth() +
  labs(x = TeX("$\\hat{p}_{t+1}$"), y = TeX("$\\omega_{t+1}$")) + thm
ggsave(filename = "doc/img/phat_omega.pdf", device = cairo_pdf)

# estimate with the result from np package
set.seed(42)
df_fit_2nd_np <- expand_grid(
  d = 1:8, n = 1:100
) %>% mutate(
  start = map(d, ~c(rnorm(n = 1, sd = .5),
                    log(1 - unname(coef(fit_1st_np)["l"])) + rnorm(n = 1, sd = .3),
                    rnorm(n = sum(1:.x + 1), sd = .5)))
) %>% mutate(
  model = map2(start, d, ~optim(par = .x, fn = obj_op_2nd, data = df_sample_2nd_np, degree = .y))
) %>% mutate(
  beta_a = map_dbl(model, ~.x$par[1]),
  beta_k = map_dbl(model, ~exp(.x$par[2])),
  rmse = map_dbl(model, ~sqrt(.x$value))
)
plot_iterated_convergence(df_fit_2nd_np)
ggsave("doc/img/beta_convergence_np.pdf", device = cairo_pdf)
plot_iterated_errorbar(df_fit_2nd_np)
ggsave("doc/img/beta_error_np.pdf", device = cairo_pdf)
fit_2nd_np <- filter(df_fit_2nd_np, d == 2)$model[[1]]


# polynomial version
set.seed(42)
df_fit_2nd_poly <- expand_grid(
  d = 1:8, n = 1:100
) %>% mutate(
  start = map(d, ~c(rnorm(n = 1, sd = .5),
                    log(1 - unname(coef(fit_1st_poly)["l"])) + rnorm(n = 1, sd = .3),
                    rnorm(n = sum(1:.x + 1), sd = .5)))
) %>% mutate(
  model = map2(start, d, ~optim(par = .x, fn = obj_op_2nd, data = df_sample_2nd_poly, degree = .y))
) %>% mutate(
  beta_a = map_dbl(model, ~.x$par[1]),
  beta_k = map_dbl(model, ~exp(.x$par[2])),
  rmse = map_dbl(model, ~sqrt(.x$value))
)

plot_iterated_convergence(df_fit_2nd_poly)
ggsave("doc/img/beta_convergence_poly.pdf", device = cairo_pdf)
plot_iterated_errorbar(df_fit_2nd_poly)
ggsave("doc/img/beta_error_poly.pdf", device = cairo_pdf)
fit_2nd_poly <- filter(df_fit_2nd_poly, d == 2)$model[[1]]


# addendum: estimate by GAM or random forest
fit_gam <- mgcv::gam(x ~ s(k_lag) + s(inv_lag), family = binomial(link = "probit"), data = df_sample %>% drop_na(x, k_lag, inv_lag))
sqrt(mse(fit_gam$y, fitted(fit_gam)) / mse(fit_gam$y, mean(fit_gam$y)))
fit_rf <- ranger::ranger(formula = x ~ k_lag + inv_lag, data = drop_na(df_sample, x, inv_lag, k_lag))
sqrt(mse(drop_na(df_sample, x, inv_lag, k_lag)$x, fit_rf$predictions) / mse(drop_na(df_sample, x, inv_lag, k_lag)$x, mean(drop_na(df_sample, x, inv_lag, k_lag)$x)))

# GAM
df_sample_2nd_gam <- make_2nd_sample(
  df_sample, coef(fit_1st_poly)["l"],
  predict(fit_1st_poly, newdata = df_sample) - coef(fit_1st_poly)["l"] * df_sample$l,
  predict(fit_gam, type = "response", newdata = df_sample))
df_fit_2nd_gam <- expand_grid(
  d = 1:8, n = 1:100
) %>% mutate(
  start = map(d, ~c(rnorm(n = 1, sd = .5), log(1 - coef(fit_1st_poly)["l"]) + rnorm(n = 1, sd = .3), rnorm(n = sum(1:.x + 1), sd = .5)))
) %>% mutate(
  model = map2(start, d, ~optim(par = .x, fn = obj_op_2nd, data = df_sample_2nd_gam, degree = .y))
) %>% mutate(
  beta_a = map_dbl(model, ~.x$par[1]),
  beta_k = map_dbl(model, ~exp(.x$par[2])),
  rmse = map_dbl(model, ~sqrt(.x$value))
)

plot_iterated_convergence(df_fit_2nd_gam)
plot_iterated_errorbar(df_fit_2nd_gam)
fit_2nd_gam <- filter(df_fit_2nd_gam, d == 2)$model[[1]]


# random forest
p_rf <- drop_na(df_sample, x, inv_lag, k_lag) %>% select(i, t) %>%
  mutate(p = predict(fit_rf, data = drop_na(df_sample, x, k_lag, inv_lag))$prediction) %>%
  left_join(df_sample, ., by = c("i", "t")) # somewhat inflexible...
df_sample_2nd_rf <- make_2nd_sample(
  df_sample, coef(fit_1st_poly)["l"],
  predict(fit_1st_poly, newdata = df_sample) - coef(fit_1st_poly)["l"] * df_sample$l,
  p_rf$p)
df_fit_2nd_rf <- expand_grid(
  d = 1:8, n = 1:100
) %>% mutate(
  start = map(d, ~c(rnorm(n = 1, sd = .5), log(1 - coef(fit_1st_poly)["l"]) + rnorm(n = 1, sd = .3), rnorm(n = sum(1:.x + 1), sd = .5)))
) %>% mutate(
  model = map2(start, d, ~optim(par = .x, fn = obj_op_2nd, data = df_sample_2nd_rf, degree = .y))
) %>% mutate(
  beta_a = map_dbl(model, ~.x$par[1]),
  beta_k = map_dbl(model, ~exp(.x$par[2])),
  rmse = map_dbl(model, ~sqrt(.x$value))
)

plot_iterated_convergence(df_fit_2nd_rf)
plot_iterated_errorbar(df_fit_2nd_rf)
fit_2nd_rf <- filter(df_fit_2nd_rf, d == 2)$model[[2]]


##### compute booststrap SE #####
estimate_op_poly <- function(data, degree_1, degree_p, degree_2nd){
  fit_1st <- lm(y ~ l + poly(k, inv, degree = degree_1), data = data, na.action = na.omit)
  beta_l <- unname(coef(fit_1st)["l"])
  rmse_1st <- sqrt(mean(resid(fit_1st)^2))
  fit_exit <- glm(x ~ l + poly(k_lag, inv_lag, degree = degree_p), data = data %>% drop_na(k_lag, inv_lag), family = binomial(link = "probit"))
  std_Brier <- std_Brier(fit_exit$data$x, fitted(fit_exit))
  df_2nd <- make_2nd_sample(
    data, beta_l,
    phi = predict(fit_1st, newdata = data) - beta_l * data$l,
    p_x = predict(fit_exit, newdata = data, type = "response")
  )
  fit_2nd <- optim(
    par = c(rnorm(n = 1, sd = .5), log(ifelse(1 - beta_l <=0, .5, 1 -beta_l)) + rnorm(n = 1, sd = .3), rnorm(n = sum(1:degree_2nd + 1), sd = .5)),
    fn = obj_op_2nd, data = df_2nd, degree = degree_2nd)
  beta_a <- mean(with(df_2nd, y_bl - exp(fit_2nd$par[2]) * k), na.rm = T)
  return(c(beta_a = beta_a, beta_l = beta_l, beta_k = unname(exp(fit_2nd$par[2])), rmse_1st = rmse_1st, std_Brier = std_Brier, rmse_2nd = sqrt(fit_2nd$value)))
}

estimate_op_gam <- function(data, degree_1, degree_2nd){
  fit_1st <- lm(y ~ l + poly(k, inv, degree = degree_1), data = data, na.action = na.omit)
  beta_l <- unname(coef(fit_1st)["l"])
  rmse_1st <- sqrt(mean(resid(fit_1st)^2))
  fit_exit_gam <- mgcv::gam(x ~ s(k_lag) + s(inv_lag), family = binomial(link = "probit"), data = data)
  std_Brier <- std_Brier(fit_exit_gam$y, fitted(fit_exit_gam))
  df_2nd <- make_2nd_sample(
    data, beta_l,
    phi = predict(fit_1st, newdata = data) - beta_l * data$l,
    p_x = predict(fit_exit_gam, newdata = data, type = "response")
  )
  fit_2nd <- optim(
    par = c(rnorm(n = 1, sd = .5), log(ifelse(1 - beta_l <=0, .5, 1 -beta_l)) + rnorm(n = 1, sd = .3), rnorm(n = sum(1:degree_2nd + 1), sd = .5)),
    fn = obj_op_2nd, data = df_2nd, degree = degree_2nd)
  beta_a <- mean(with(df_2nd, y_bl - exp(fit_2nd$par[2]) * k), na.rm = T)
  return(c(beta_a = beta_a, beta_l = beta_l, beta_k = unname(exp(fit_2nd$par[2])), rmse_1st = rmse_1st, std_Brier = std_Brier, rmse_2nd = sqrt(fit_2nd$value)))
}

estimate_op_rf <- function(data, degree_1, degree_2nd, ...){
  fit_1st <- lm(y ~ l + poly(k, inv, degree = degree_1), data = data, na.action = na.omit)
  beta_l <- unname(coef(fit_1st)["l"])
  rmse_1st <- sqrt(mean(resid(fit_1st)^2))
  fit_exit_rf <- ranger::ranger(formula = x ~ k_lag + inv_lag, data = drop_na(df_sample, x, inv_lag, k_lag))
  p <- drop_na(data, x, inv_lag, k_lag) %>% select(i, t) %>%
    mutate(p = predict(fit_rf, data = drop_na(data, x, k_lag, inv_lag))$prediction) %>%
    left_join(data, ., by = c("i", "t"))
  p <- p$p
  std_Brier <- std_Brier(drop_na(df_sample, x, k_lag, inv_lag)$x, fit_rf$predictions)
  df_2nd <- make_2nd_sample(
    data, beta_l,
    phi = predict(fit_1st, newdata = data) - beta_l * data$l,
    p_x = p
  )
  fit_2nd <- optim(
    par = c(rnorm(n = 1, sd = .5), log(ifelse(1 - beta_l <=0, .5, 1 -beta_l)) + rnorm(n = 1, sd = .3), rnorm(n = sum(1:degree_2nd + 1), sd = .5)),
    fn = obj_op_2nd, data = df_2nd, degree = degree_2nd)
  beta_a <- mean(with(df_2nd, y_bl - exp(fit_2nd$par[2]) * k), na.rm = T)
  return(c(beta_a = beta_a, beta_l = beta_l, beta_k = unname(exp(fit_2nd$par[2])), rmse_1st = rmse_1st, std_Brier = std_Brier, rmse_2nd = sqrt(fit_2nd$value)))
}

# test
# estimate_op_poly(df_sample, 4, 2, 3)
# estimate_op_gam(df_sample, 4, 3)
# estimate_op_rf(df_sample, 4, 3)


get_stat_boot <- function(data, indices, indices_, estimate, params){
  # estimate = estimate function
  # params = listed model hyper parameters (e.g., degrees of polynomials)
  data$y <- with(data, fitted + resid[indices])
  do.call(estimate, args = c(list(data = data), params))
}

op_params_np <- list(beta_k = exp(fit_2nd_np$par[2]), beta_l = unname(coef(fit_1st_np)))
op_params_poly <- list(beta_k = exp(fit_2nd_poly$par[2]), beta_l = unname(coef(fit_1st_poly)["l"]))
op_params_gam <- list(beta_k = unname(exp(fit_2nd_gam$par[2])), beta_l = unname(coef(fit_1st_poly)["l"]))
op_params_rf <- list(beta_k = unname(exp(fit_2nd_rf$par[2])), beta_l = unname(coef(fit_1st_np)["l"]))

set.seed(42)
boot_poly <- df_sample_2nd_poly %>%
  mutate(., fitted = if_else(is.na(y), NA_real_, do.call(get_fitted, args = c(list(data = .), op_params_poly))),
         resid = y - fitted) %>%
  boot(., get_stat_boot, R = 100, strata = .$x, estimate = estimate_op_poly, params = list(degree_1 = 3, degree_p = 2, degree_2nd = 2))
set.seed(42)
boot_gam <- df_sample_2nd_gam %>%
  mutate(., fitted = if_else(is.na(y), NA_real_, do.call(get_fitted, args = c(list(data = .), op_params_gam))),
         resid = y - fitted) %>%
  boot(., get_stat_boot, R = 100, strata = .$x, estimate = estimate_op_gam, params = list(degree_1 = 3, degree_2nd = 2))
set.seed(42)
boot_rf <- df_sample_2nd_rf %>%
  mutate(., fitted = if_else(is.na(y), NA_real_, do.call(get_fitted, args = c(list(data = .), op_params_rf))),
         resid = y - fitted) %>%
  boot(., get_stat_boot, R = 100, strata = .$x, estimate = estimate_op_rf, params = list(degree_1 = 3, degree_2nd = 2))


plot_boot_hist <- function(boot, name_params = c("beta_a", "beta_l", "beta_k", "RMSE_1st", "std_Brier", "RMSE_2nd")){
  as_tibble(boot$t) %>% set_names(name_params) %>%
    pivot_longer(cols = everything()) %>% filter(str_detect(name, "^beta")) %>%
    mutate(name = factor(name, labels = c(TeX("$\\beta_A$"), TeX("$\\beta_K$"), TeX("$\\beta_L$")))) %>%
    ggplot(aes(x = value)) + geom_histogram(bins = as.integer(max(boot$R/10), 10)) +
    facet_wrap(~name, scales="free", ncol = 1, labeller = label_parsed)
}

plot_boot_hist(boot_poly) + thm + theme(axis.title = element_blank()) + labs(title = "Polynomial Aprrox.")
plot_boot_hist(boot_gam) + thm + theme(axis.title = element_blank()) + labs(title = "GAM Approx.")
plot_boot_hist(boot_rf) + thm + theme(axis.title = element_blank()) + labs(title = "RF Approx.")


##### OP法の最終結果 #####
df_result_table_op <- list(OP_poly = boot_poly, OP_GAM = boot_gam, OP_RF = boot_rf) %>%
  map2_dfr(., names(.), ~skim(.x$t) %>% select(skim_variable, numeric.mean, numeric.sd) %>%
             mutate(name = names(.x$t0), model = .y) %>%
             rename(value = numeric.mean, std_err = numeric.sd) %>% select(name, value, std_err, model))


##### Assignment by Kawaguchi #####
# Kawaguchi's Assignment (GMM procedure without exit behavior)
moment_op <- function(theta, df){
  # theta = [beta_a, beta_k, alpha]
  z <- with(df, y_bl - theta[1] - theta[2] * k - theta[3] * (phi_lag - theta[1] - theta[2] * k))
  x <- df %>% ungroup %>% select(k, k_lag, inv_lag) %>% mutate(uni = 1) %>% as.matrix
  return(z * x)
}

set.seed(42)
fit_gmm <- gmm(g = moment_op, x = df_sample_2nd_np %>% drop_na(y, y_bl, k, k_lag, inv_lag, phi_lag),
               t0 = c(.1, .1, .1))
summary(fit_gmm)

# computing boostrap SE
fitted_gmm <- function(model, data, beta_l = unname(fit_1st_np$xcoef[1])){
  params <- model$coefficients
  with(data, params[1] + params[2] * k + beta_l * l)
}
resid_gmm <- function(model, data, beta_l = unname(fit_1st_np$xcoef[1])){
  params <- model$coefficients
  with(data, y - params[1] + params[2] * k + beta_l * l)
}
stat_gmm <- function(data, indices, indices_, bw){
  data$y <- with(data, fitted + resid[indices])
  fit_1st_np <- npplreg(bws = bw, txdat = as.data.frame(data[, "l"]), tydat = data$y, tzdat = as.data.frame(data[, c("k", "inv")]))
  beta_l <- unname(coef(fit_1st_np)["l"])
  phi_hat <- predict(fit_1st_np, newdata = data %>% drop_na(y, k, l, inv)) %>% as.numeric - drop_na(data, y, k, l, inv)$l * beta_l
  data <- make_2nd_sample(data, beta_l, phi_hat, 1)
  fit_gmm <- gmm(g = moment_op, x = data %>% drop_na(y, y_bl, k, k_lag, inv_lag, phi_lag), t0 = c(.1, .1, .1))
  params <- unname(fit_gmm$coefficients)
  beta_a <- params[1]; beta_k <- params[2]
  rmse = sqrt(mean(with(data, y - beta_a - beta_k * k -beta_l * l)^2, na.rm = T))
  return(c(beta_a = beta_a, beta_k = beta_k, beta_l = beta_l, rmse = rmse))
}

# compute boostrap SE
set.seed(42)
boot_gmm <- boot(df_sample %>% drop_na(y, k, l) %>% mutate(., fitted = fitted_gmm(fit_gmm, .), resid = resid_gmm(fit_gmm, .)),
                 stat_gmm, R = 100, bw = bw_1st)

df_sample %>% drop_na(y, k, l) %>% mutate(., fitted = fitted_gmm(fit_gmm, .), resid = resid_gmm(fit_gmm, .))

# GMMでの最終結果
as_tibble(boot_gmm$t) %>% set_names(names(boot_gmm$t0)) %>% pivot_longer(cols = everything()) %>% ggplot(aes(x = value)) +
  geom_histogram(bins = 10) + facet_wrap(~name, scales = "free", ncol = 1)
df_result_table_kawaguchi <- tibble(
  name = c("beta_a", "beta_k", "beta_l"),
  value = c(fit_gmm$coefficients[1:2], coef(fit_1st_np)) ,
  std_err = apply(boot_gmm$t, 2, sd)[1:3]
  ) %>% mutate(model = "GMM (Kawaguchi)")

#####  Comparison to standard panel data methods ######

# OLS
summary(fit_ols <- lm(y ~ k + l, data = df_sample))
fit_lsdv <- lm(y ~ k + l + i, data = df_sample)
# LSDV
summary(fit_lsdv)$coefficients %>% as_tibble(rownames = "Coefficient") %>% filter(!str_detect(Coefficient, "^i"))
# FE(within)
# lm() でも工夫すれば係数を求められるが標準誤差が異なる
fit_within_pseudo <- lm(y ~ -1 + k + l,
                        data = left_join(df_sample, df_sample %>% group_by(i) %>%
                                           summarise_all(list(mean = mean)), by = 'i') %>%
                          mutate(y = y - y_mean, k = k - k_mean, l = l - l_mean)
)
summary(fit_within_pseudo)
# by plm
fit_within <- plm(y ~ k + l, data = df_sample %>% arrange(i, t),
                  index = c("i", "t"),
                  model = "within")
summary(fit_within) # 係数は変わらないが標準誤差が異なる
summary(fit_between <- update(fit_within, model = "between"))
summary(fit_twoway <- update(fit_within, model = "within", effect = "twoways"))
summary(fit_FD <- update(fit_within, model = "fd"))
summary(fit_RE <- update(fit_within, model = "random"))
# Diff-GMM by Arellano and Bond, and System GMM by Blundel and Bond
fit_pgmm <- list(
  AB_twoway = pgmm(y ~ plm::lag(k) + k + l | plm::lag(k, 2:99) + plm::lag(l, 2:99), data = df_sample, effect = "twoways"),
  AB_indivi = pgmm(y ~ plm::lag(k) + k + l | plm::lag(k, 2:99) + plm::lag(l, 2:99), data = df_sample, effect = "individual"),
  BB_twoway = pgmm(y ~ plm::lag(k) + k + l | plm::lag(k, 2:99) + plm::lag(l, 2:99), data = df_sample, effect = "twoways", transformation = "ld"),
  BB_indivi = pgmm(y ~ plm::lag(k) + k + l | plm::lag(k, 2:99) + plm::lag(l, 2:99), data = df_sample, effect = "individual", transformation = "ld")
)
fit_pgmm %>% map(., ~summary(.x, robust = T))

fit_estprod <- estprod::olley_pakes(
  y ~ l | k | k + inv,
  data = df_ground_truth %>% mutate(y = y_true) %>% select(y, k, l, inv, i, t, p_x),
  id = "i", time = "t")

df_result_table_plm <- c(list(OLS = fit_ols, Within = fit_within, Between = fit_between, Twoway = fit_twoway, FD = fit_FD, RE = fit_RE),
                         fit_pgmm) %>% map2_dfr(., names(.), ~tibble(model = .y, fit = list(.x))) %>%
  mutate(
    beta_a_value = map_dbl(fit, ~coef(.x)["(Intercept)"]),
    beta_a_value = if_else(
      is.na(beta_a_value),
      map_dbl(fit, ~with(df_sample, mean(y - coef(.x)["k"] * k + coef(.x)["l"] * l, na.rm = T))),
      beta_a_value),
    beta_a_se = map_dbl(fit, function(x){
      coefs <- summary(x, robust = T)$coefficients
      if("(Intercept)" %in% rownames(coefs)){
        coefs["(Intercept)", "Std. Error"]
      } else{
        NA
      }
    }),
    beta_l_value = map_dbl(fit, ~coef(.x)["l"]),
    beta_l_se = map_dbl(fit, ~summary(.x)$coefficients["l", "Std. Error"]),
    beta_k_value = map_dbl(fit, ~coef(.x)["k"]),
    beta_k_se = map_dbl(fit, ~summary(.x)$coefficients["k", "Std. Error"])
  )  %>% pivot_longer(cols = starts_with("beta_"),
                      names_to = c("name", ".value"), names_pattern = "(beta_[alk])_(.+)") %>% rename(std_err = se) 

bind_rows(df_result_table_op, df_result_table_plm) %>%
  left_join(tibble(name = c("beta_a", "beta_k", "beta_l"), true = c(beta_a. ,beta_k., beta_l.)), by = c("name")) %>%
  mutate(bias = value - true) %>% arrange(desc(abs(bias)))

bind_rows(df_result_table_op, df_result_table_plm) %>%
  left_join(tibble(name = c("beta_a", "beta_k", "beta_l"), true = c(beta_a. ,beta_k., beta_l.)), by = c("name")) %>%
  mutate(bias = value - true) %>%
  ggplot(aes(x = name, y = bias, group = model, color = model)) +
  geom_point(shape = "x", size = 10) +
  thm + theme(legend.title = element_blank(), axis.title = element_blank()) + scale_color_colorblind() +
  scale_x_discrete(labels = c(expression(beta[A]), expression(beta[K]), expression(beta[L])))

df_result_table_op %>%
  ggplot(aes(x = name, y = value, ymin = value - std_err, ymax = value + std_err, group = model, color = model)) +
  geom_point(shape = "x", size = 10, position=position_dodge(width=0.5)) +
  geom_errorbar(width = .1, position=position_dodge(width=0.5)) +
  thm + theme(legend.title = element_blank(), axis.title = element_blank()) + scale_color_colorblind() +
  scale_x_discrete(labels = c(expression(beta[A]), expression(beta[K]), expression(beta[L])))

df_result_table_plm %>%
  ggplot(aes(x = name, y = value, ymin = value - std_err, ymax = value + std_err, group = model, color = model)) +
  geom_point(shape = "x", size = 10, position=position_dodge(width=0.5)) +
  geom_errorbar(width = .1, position=position_dodge(width=0.5)) +
  thm + theme(legend.title = element_blank(), axis.title = element_blank()) + scale_color_colorblind() +
  scale_x_discrete(labels = c(expression(beta[A]), expression(beta[K]), expression(beta[L])))

bind_rows(df_result_table_op, df_result_table_plm) %>%
  filter(name %in% c("beta_a", "beta_k", "beta_l")) %>%
  filter(!str_detect(model, "^BB_"), !str_detect(model, "^AB_"), model != "Twoway", !model %in% c("OP_RF", "OP_GAM")) %>%
  mutate(model = str_replace(model, "OP_poly", "OP (polynomial)"),
         name = factor(name, labels = c(TeX("$\\beta_A$"), TeX("$\\beta_K$"), TeX("$\\beta_L$")))
         ) %>%
  ggplot(aes(x = 0, y = value, ymin = value - std_err, ymax = value + std_err, group = model, color = model)) +
  geom_point(size = 4, position=position_dodge(width=0.5)) +
  geom_errorbar(width = .3, size = 1, position=position_dodge(width=0.5)) +
  geom_hline(aes(yintercept = val),
             data = tibble(name = factor(c("beta_a", "beta_k", "beta_l"),
                                         labels = c(TeX("$\\beta_A$"), TeX("$\\beta_K$"), TeX("$\\beta_L$"))),
                           val = c(beta_a., beta_k., beta_l.)), linetype = 2) +
  thm + theme(legend.title = element_blank(), axis.title = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  scale_color_colorblind() +
  facet_wrap(~name, labeller = label_parsed) +
  labs(caption = "dashed lines denote the true value.")
ggsave(filename = "doc/img/result.pdf", device = cairo_pdf)

bind_rows(df_result_table_op, df_result_table_plm, df_result_table_kawaguchi) %>% select(name, value, model) %>%
  filter(name %in% c("beta_a", "beta_k", "beta_l")) %>%
  pivot_wider(id_cols = model, names_from = name, values_from = value) %>%
  mutate(beta_a = beta_a - beta_a., beta_k = beta_k - beta_k., beta_l = beta_l - beta_l.) %>%
  set_names(c("model", "beta_A", "beta_K", "beta_L"))
