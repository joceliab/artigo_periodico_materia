##### ==========================================================
##### Pacotes 
##### ==========================================================
suppressPackageStartupMessages({
  library(rsm)
  library(ggplot2)
  library(patchwork)
  library(plotly)
  library(viridis)      # escala para ggplot/contornos
})

##### ==========================================================
##### Dados (FCCD com ponto central)
##### ==========================================================
dados <- data.frame(
  Concentracao = c(0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.12,0.12,0.12,0.12,0.12,0.12,0.12,0.12,0.08,0.12,0.10,0.10,0.10,0.10,0.10,0.10,0.10),
  Fluxo        = c(0.2,0.2,0.2,0.2,0.8,0.8,0.8,0.8,0.2,0.2,0.2,0.2,0.8,0.8,0.8,0.8,0.5,0.5,0.2,0.8,0.5,0.5,0.5,0.5,0.5),
  Distancia    = c(15,15,20,20,15,15,20,20,15,15,20,20,15,15,20,20,17.5,17.5,17.5,17.5,15,20,17.5,17.5,17.5),
  Voltagem     = c(15,25,15,25,15,25,15,25,15,25,15,25,15,25,15,25,20,20,20,20,20,20,15,25,20),
  Diametro     = c(3.067,3.755,3.217,1.712,2.53,3.002,2.887,2.561,3.333,2.24,3.855,2.925,2.884,3.098,3.44,3.405,2.037,2.66,2.698,3.196,2.661,2.699,2.661,2.565,2.108)
)

##### ==========================================================
#####  Codificação
##### ==========================================================
centros    <- list(Concentracao = 0.10, Fluxo = 0.50, Distancia = 17.5, Voltagem = 20)
intervalos <- list(Concentracao = 0.02, Fluxo = 0.30, Distancia = 2.50, Voltagem = 5.0)

codificados <- coded.data(
  dados,
  C_cod ~ (Concentracao - centros$Concentracao) / intervalos$Concentracao,
  F_cod ~ (Fluxo        - centros$Fluxo)        / intervalos$Fluxo,
  D_cod ~ (Distancia    - centros$Distancia)    / intervalos$Distancia,
  V_cod ~ (Voltagem     - centros$Voltagem)     / intervalos$Voltagem
)

##### ==========================================================
##### Ajuste do modelo quadrático completo (RSM)
##### ==========================================================
modelo_full <- rsm(Diametro ~ SO(C_cod, F_cod, D_cod, V_cod), data = codificados)

cat("\n=== SUMMARY (modelo_full) ===\n"); print(summary(modelo_full))
cat("\n=== ANOVA (modelo_full) ===\n");  print(anova(modelo_full))

cat("\n=== Análise canônica ===\n")
can <- summary(modelo_full)$canonical
print(can)

# Ponto estacionário (reais)
xs      <- can$xs
xs_real <- c(
  Concentracao = centros$Concentracao + intervalos$Concentracao * xs["C_cod"],
  Fluxo        = centros$Fluxo        + intervalos$Fluxo        * xs["F_cod"],
  Distancia    = centros$Distancia    + intervalos$Distancia    * xs["D_cod"],
  Voltagem     = centros$Voltagem     + intervalos$Voltagem     * xs["V_cod"]
)
cat("\n--- Ponto estacionário (variáveis REAIS) ---\n"); print(xs_real)

#### ==========================================================
##### Métricas de erro 
##### =========================================================
preditos <- predict(modelo_full)          
reais    <- dados$Diametro
erro     <- preditos - reais
mape     <- mean(abs(erro / reais)) * 100
mse      <- mean(erro^2)
rmse     <- sqrt(mse)
cat("\n=== Métricas (modelo_full) ===\n")
cat(sprintf("MAPE (%%): %.3f\nMSE: %.5f\nRMSE: %.5f\n", mape, mse, rmse))

##### ==========================================================
##### Superfície de resposta
##### ==========================================================
persp(modelo_full, ~ C_cod + F_cod, zlab = "Diâmetro", col = "lightblue")

##### ==========================================================
##### Modelo reduzido automático com hierarquia (via BIC)
##### ==========================================================
lin    <- c("C_cod","F_cod","D_cod","V_cod")
ints   <- combn(lin, 2, FUN=function(x) paste(x, collapse=":"))
quads  <- paste0("I(", lin, "^2)")
form_str <- paste("Diametro ~", paste(c(lin, ints, quads), collapse=" + "))
lm_full  <- lm(as.formula(form_str), data = codificados)  # “full” em lm (para comparações e intervalos)

# (b) seleção backward por BIC
lm_step <- step(lm_full, direction="backward", k=log(nrow(codificados)), trace=0)

# (c) impor hierarquia
get_terms <- function(obj) attr(terms(obj), "term.labels")
terms_sel <- get_terms(lm_step)
pais_interacoes  <- unique(unlist(strsplit(terms_sel[grepl(":", terms_sel)], ":")))
pais_quadraticos <- gsub("^I\\((.*)\\^2\\)$", "\\1", terms_sel[grepl("^I\\(.*\\^2\\)$", terms_sel)])
pais_necessarios <- unique(c(pais_interacoes, pais_quadraticos))
lineares_faltantes <- setdiff(pais_necessarios, terms_sel)
termos_finais <- unique(c(terms_sel, lineares_faltantes))

form_red_str <- paste("Diametro ~", paste(termos_finais, collapse=" + "))
modelo_red   <- lm(as.formula(form_red_str), data=codificados)

cat("\n=== SUMMARY (modelo_red) ===\n"); print(summary(modelo_red))
cat("\n=== ANOVA (modelo_red) ===\n");  print(anova(modelo_red))

# Diagnóstico e métricas do reduzido
op <- par(mfrow=c(2,2)); plot(modelo_red); par(op)
pred_red <- predict(modelo_red)
obs      <- dados$Diametro
err      <- pred_red - obs
mse_r    <- mean(err^2)
rmse_r   <- sqrt(mse_r)
mape_r   <- mean(abs(err/obs))*100
r2aj_r   <- summary(modelo_red)$adj.r.squared
cat("\n=== Métricas (modelo_red) ===\n")
cat(sprintf("R²aj: %.4f | MAPE (%%): %.3f | MSE: %.5f | RMSE: %.5f\n", r2aj_r, mape_r, mse_r, rmse_r))

##### ==========================================================
##### Predições em novos pontos
##### ==========================================================
novos <- data.frame(
  Concentracao = c(0.09, 0.10, 0.11),
  Fluxo        = c(0.30, 0.50, 0.70),
  Distancia    = rep(17.5, 3),
  Voltagem     = c(18, 20, 22)
)
novos_cod <- within(novos, {
  C_cod <- (Concentracao - centros$Concentracao)/intervalos$Concentracao
  F_cod <- (Fluxo        - centros$Fluxo)       /intervalos$Fluxo
  D_cod <- (Distancia    - centros$Distancia)   /intervalos$Distancia
  V_cod <- (Voltagem     - centros$Voltagem)    /intervalos$Voltagem
})

# FIX: intervalos só com modelos lm
pred_full_int <- predict(lm_full, newdata=novos_cod, interval="prediction")
pred_red_int  <- predict(modelo_red, newdata=novos_cod, interval="prediction")
cat("\n=== Predições (lm_full) ===\n"); print(cbind(novos, pred_full_int))
cat("\n=== Predições (modelo_red) ===\n");  print(cbind(novos, pred_red_int))

##### ==========================================================
##### Comparações coerentes
##### ==========================================================
# FIX: comparar lm vs lm
cat("\n=== anova(modelo_red, lm_full) ===\n"); print(anova(modelo_red, lm_full))
cat("\n=== AIC/BIC ===\n"); print(AIC(lm_full, modelo_red)); print(BIC(lm_full, modelo_red))
cat("\n=== R²aj ===\n"); print(summary(modelo_full)$adj.r.squared); print(summary(modelo_red)$adj.r.squared)

# Métricas novamente para tabela
pred_full <- predict(lm_full)
erro_full <- pred_full - obs
rmse_full <- sqrt(mean(erro_full^2))
mape_full <- mean(abs(erro_full/obs))*100

comparacao <- data.frame(
  Modelo = c("Completo (lm_full)","Reduzido"),
  R2_aj  = c(summary(lm_full)$adj.r.squared, summary(modelo_red)$adj.r.squared),
  RMSE   = c(rmse_full, rmse_r),
  MAPE   = c(mape_full, mape_r),
  AIC    = c(AIC(lm_full), AIC(modelo_red)),
  BIC    = c(BIC(lm_full), BIC(modelo_red))
)
print(comparacao)

# FIX: escolha final com teste coerente (linha 2 do anova)
p_lr <- anova(modelo_red, lm_full)$`Pr(>F)`[2]
if (!is.na(p_lr) && p_lr > 0.05 && BIC(modelo_red) <= BIC(lm_full)) {
  modelo_para_mapas <- modelo_red
  modelo_nome <- "Modelo reduzido (lm)"
} else {
  modelo_para_mapas <- modelo_full
  modelo_nome <- "Modelo completo (rsm)"
}
cat(sprintf("\n*** Usando %s para mapas/superfícies. ***\n", modelo_nome))

##### ==========================================================
##### Diagnóstico de resíduos
##### ==========================================================
# Garantir que temos um modelo escolhido
if (!exists("modelo_para_mapas")) {
  if (exists("modelo_red")) {
    modelo_para_mapas <- modelo_red
  } else if (exists("modelo_full")) {
    modelo_para_mapas <- modelo_full
  } else {
    stop("Nenhum modelo encontrado para diagnóstico de resíduos.")
  }
}

# Preditos, resíduos e observado
pred_best  <- as.numeric(predict(modelo_para_mapas))
resid_best <- as.numeric(residuals(modelo_para_mapas))
aux <- data.frame(pred = pred_best, resid = resid_best, obs = dados$Diametro)

# Gráficos
g_res1 <- ggplot(aux, aes(x = pred, y = resid)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "Resíduos vs Ajustado", x = "Ajustado", y = "Resíduo") +
  theme_minimal(base_size = 12)

g_res2 <- ggplot(aux, aes(sample = resid)) +
  stat_qq() +
  stat_qq_line() +
  labs(title = "QQ-Plot dos Resíduos") +
  theme_minimal(base_size = 12)

g_res3 <- ggplot(aux, aes(x = resid)) +
  geom_histogram(bins = 15, color = "black") +
  labs(title = "Histograma dos Resíduos", x = "Resíduo", y = "Frequência") +
  theme_minimal(base_size = 12)

g_res4 <- ggplot(aux, aes(x = pred, y = obs)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(title = "Observado vs Ajustado", x = "Ajustado", y = "Observado") +
  theme_minimal(base_size = 12)

# Painel 2x2
print((g_res1 + g_res2) / (g_res3 + g_res4))

##### ==========================================================
##### Mapa de erro da simplificação
##### ==========================================================
# FIX: gerar grade real e codificar
codificar <- function(x, centro, intervalo) (x - centro)/intervalo
fatores_reais <- c("Concentracao","Fluxo","Distancia","Voltagem")
# usaremos par (Concentracao, Fluxo) como exemplo de auditoria
seqC <- seq(centros$Concentracao - intervalos$Concentracao, centros$Concentracao + intervalos$Concentracao, length.out=40)
seqF <- seq(centros$Fluxo        - intervalos$Fluxo,        centros$Fluxo        + intervalos$Fluxo,        length.out=40)
grid  <- expand.grid(Concentracao=seqC, Fluxo=seqF)
grid$Distancia <- centros$Distancia
grid$Voltagem  <- centros$Voltagem
grid$C_cod <- codificar(grid$Concentracao, centros$Concentracao, intervalos$Concentracao)
grid$F_cod <- codificar(grid$Fluxo,        centros$Fluxo,        intervalos$Fluxo)
grid$D_cod <- codificar(grid$Distancia,    centros$Distancia,    intervalos$Distancia)
grid$V_cod <- codificar(grid$Voltagem,     centros$Voltagem,     intervalos$Voltagem)

grid$y_full <- predict(lm_full,   newdata = grid)
grid$y_red  <- predict(modelo_red, newdata = grid)
grid$delta  <- grid$y_full - grid$y_red
# (visualização opcional)
# ggplot(grid, aes(Concentracao, Fluxo, fill=delta)) + geom_tile() + scale_fill_viridis_c()

##### ==========================================================
##### Superfícies 3D e contornos 2D
##### ==========================================================
nomes_formatados <- list(
  Concentracao = "Concentração (g/mL)",
  Fluxo        = "Fluxo (mL/h)",
  Distancia    = "Distância (cm)",
  Voltagem     = "Voltagem (kV)"
)
fatores_reais <- c("Concentracao","Fluxo","Distancia","Voltagem")
pares <- combn(fatores_reais, 2, simplify = FALSE)

`%||%` <- function(a, b) if (is.null(a) || is.na(a)) b else a

superficie_3d <- function(f1, f2, modelo, centros, intervalos, nomes_formatados, n = 40) {
  seq1 <- seq(centros[[f1]] - intervalos[[f1]], centros[[f1]] + intervalos[[f1]], length.out = n)
  seq2 <- seq(centros[[f2]] - intervalos[[f2]], centros[[f2]] + intervalos[[f2]], length.out = n)
  grid <- expand.grid(f1 = seq1, f2 = seq2); names(grid) <- c(f1, f2)
  for (f in setdiff(names(centros), c(f1, f2))) grid[[f]] <- centros[[f]]
  grid$C_cod <- codificar(grid$Concentracao, centros$Concentracao, intervalos$Concentracao)
  grid$F_cod <- codificar(grid$Fluxo,        centros$Fluxo,        intervalos$Fluxo)
  grid$D_cod <- codificar(grid$Distancia,    centros$Distancia,    intervalos$Distancia)
  grid$V_cod <- codificar(grid$Voltagem,     centros$Voltagem,     intervalos$Voltagem)

  # usar o modelo selecionado para mapas
  yhat <- if (inherits(modelo, "rsm")) predict(modelo, newdata = grid) else as.numeric(predict(modelo, newdata = grid))
  grid <- grid[order(match(grid[[f1]], seq1), match(grid[[f2]], seq2)), ]
  zmat <- matrix(yhat, nrow = n, ncol = n, byrow = FALSE)

  t1 <- as.character(nomes_formatados[[f1]] %||% f1)
  t2 <- as.character(nomes_formatados[[f2]] %||% f2)

  plot_ly(x = ~seq1, y = ~seq2, z = ~zmat) |>
    add_surface(colorscale = "Viridis") |>
    layout(
      title = paste(t1, "vs", t2, sprintf("— %s", if (inherits(modelo, "rsm")) "rsm" else "lm")),
      scene = list(
        xaxis = list(title = t1),
        yaxis = list(title = t2),
        zaxis = list(title = "Diâmetro (µm)")
      )
    )
}

# gerar superfícies com o modelo escolhido
modelo_escolhido <- if (exists("modelo_para_mapas")) modelo_para_mapas else modelo_full
for (parf in pares) {
  f1 <- parf[1]; f2 <- parf[2]
  print(superficie_3d(f1, f2, modelo_escolhido, centros, intervalos, nomes_formatados, n = 40))
}

# Painel 2D de contornos preenchidos
graficos <- list()
for (parf in pares) {
  f1 <- parf[1]; f2 <- parf[2]
  seq1 <- seq(centros[[f1]] - intervalos[[f1]], centros[[f1]] + intervalos[[f1]], length.out = 60)
  seq2 <- seq(centros[[f2]] - intervalos[[f2]], centros[[f2]] + intervalos[[f2]], length.out = 60)
  grid <- expand.grid(f1=seq1, f2=seq2); names(grid) <- c(f1, f2)
  for (f in setdiff(names(centros), c(f1, f2))) grid[[f]] <- centros[[f]]
  grid$C_cod <- codificar(grid$Concentracao, centros$Concentracao, intervalos$Concentracao)
  grid$F_cod <- codificar(grid$Fluxo,        centros$Fluxo,        intervalos$Fluxo)
  grid$D_cod <- codificar(grid$Distancia,    centros$Distancia,    intervalos$Distancia)
  grid$V_cod <- codificar(grid$Voltagem,     centros$Voltagem,     intervalos$Voltagem)

  grid$yhat <- if (inherits(modelo_escolhido, "rsm")) predict(modelo_escolhido, newdata = grid) else as.numeric(predict(modelo_escolhido, newdata = grid))

  g <- ggplot(grid, aes_string(x = f1, y = f2, fill = "yhat")) +
    geom_tile() +
    geom_contour(aes(z = yhat), color = "black", alpha = 0.6) +
    scale_fill_viridis_c(option = "C") +
    labs(title = paste(nomes_formatados[[f1]], "vs", nomes_formatados[[f2]], sprintf("— %s", if (inherits(modelo_escolhido,"rsm") ) "rsm" else "lm")),
         x = nomes_formatados[[f1]], y = nomes_formatados[[f2]], fill = "Diâmetro (µm)") +
    theme_minimal(base_size = 11)
  graficos[[paste(f1,f2,sep="_")]] <- g
}
print(wrap_plots(graficos, ncol = 2))
