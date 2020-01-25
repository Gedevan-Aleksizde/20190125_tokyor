if(!"pacman" %in% installed.packages()) install.packages(pacman)
pacman::p_load(tidyverse, ggthemes, estprod, skimr, stargazer, estimatr, plotly)
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

#### Demonstration ####

df <- estprod_data %>% ungroup %>%
  rename(y = var1, l = var2, k = var3, inv = var4) %>% group_by(id)
skim(df %>% ungroup)
df_summary <- skim(df) %>% as_tibble(.name_repair = "unique") %>%
  filter(skim_type != "factor") %>% dplyr::select(starts_with("numeric"), skim_variable, -numeric.hist) %>%
  rename(variable = skim_variable)
df_summary %>% group_by(variable) %>% skim() %>% arrange(variable)

fit_ols <- lm(formula = y ~ k + l, data = df)
fit_lsdv <- lm(formula = y ~ k + l + id, data =df)
fit_op <- olley_pakes(
  y ~ l | k | inv, data = df,
  exit = ~exit, id = "id", time = "year", bootstrap = T)
fit_lp <- levinsohn_petrin(
  formula = y ~ l | k | inv, data = df,
  exit = ~exit, id = "id", time = "year", bootstrap = T)

coef_estprod <- function(model){
  b_a <- model$data %>% ungroup %>% summarise(b_a = mean(y - model$t0[1] * l + model$t0[2] * k)) %>% as.numeric
  return(c(c(`(Intercept)`=b_a), set_names(model$t0, model$varnames)))
}

predict_estprod <- function(model, newdata=NULL){
  b_a <- coef_estprod(model)[1]
  if(is.null(newdata)){
    mutate(model$data, yhat = b_a + model$t0[1] * l + model$t0[2] * k)$yhat
  } else{
    mutate(newdata, yhat = b_a + model$t0[1] * l + model$t0[2] * k)$yhat
  }
}

summary(fit_ols)
summary(fit_op)
summary(fit_lp)
coef_estprod(fit_op)

df_result <- left_join(
  data.frame(OLS = coef(fit_ols)) %>% rownames_to_column(),
  data.frame(OP = coef_estprod(fit_op)) %>% rownames_to_column(),
  by = "rowname"
  ) %>% left_join(
    y = data.frame(LP = coef_estprod(fit_lp)) %>% rownames_to_column(), by = "rowname"
    ) %>%
  left_join(
    data.frame(LSDV = coef(fit_lsdv)) %>% rownames_to_column(),
    by = "rowname"
  )  %>% rename(variable = rowname)

df_result %>% pivot_longer(cols = -variable, names_to = "model") %>%
  ggplot(aes(x = variable, y = value, group = model, color = model)) + geom_line(size = 2) + geom_point(aes(color = model), size = 4) +
  scale_color_colorblind() +
  theme_classic(base_size = 30) +
  theme(axis.title = element_blank(), legend.position = "bottom", legend.title = element_blank())

ggsave(filename = "img/result.pdf")

Hmisc::latex(df_result, file = "tab.tex")
