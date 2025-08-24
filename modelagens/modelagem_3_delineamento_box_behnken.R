##### ==========================================================
##### 0) Pacotes
##### ==========================================================
use_or_install <- function(pkgs) {
  miss <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
  if (length(miss)) install.packages(miss, quiet = TRUE)
  invisible(lapply(pkgs, library, character.only = TRUE))
}
use_or_install(c("rsm","ggplot2","viridis","patchwork"))
plotly_ok <- requireNamespace("plotly", quietly = TRUE)
if (plotly_ok) library(plotly)

theme_set(theme_minimal(base_size = 12))

##### ==========================================================
##### 1) Dados (BBD)
##### ==========================================================
dados <- data.frame(
  Concentracao = c(0.08, 0.08, 0.08, 0.08, 0.08, 0.08,
                   0.10, 0.10, 0.10, 0.10, 0.10, 0.10,
                   0.10, 0.10, 0.10, 0.10, 0.10, 0.10,
                   0.12, 0.12, 0.12, 0.12, 0.12, 0.12, 0.10),
  Fluxo = c(0.2, 0.5, 0.5, 0.5, 0.5, 0.8,
            0.2, 0.2, 0.2, 0.2, 0.5, 0.5,
            0.5, 0.5, 0.8, 0.8, 0.8, 0.8,
            0.2, 0.5, 0.5, 0.5, 0.5, 0.8, 0.5),
  Distancia = c(17.5, 15.0, 17.5, 17.5, 20.0, 17.5,
                15.0, 17.5, 17.5, 20.0, 15.0, 15.0,
                20.0, 20.0, 15.0, 17.5, 17.5, 20.0,
                17.5, 15.0, 17.5, 17.5, 20.0, 17.5, 17.5),
  Voltagem = c(20, 20, 15, 25, 20, 20,
               20, 15, 25, 20, 15, 25,
               15, 25, 20, 15, 25, 20,
               20, 20, 15, 25, 20, 20, 20),
  Diametro = c(2.622, 3.293, 3.363, 2.899, 2.829, 2.626,
               1.885, 3.667, 2.804, 3.120, 2.928, 2.993,
               3.121, 2.486, 1.957, 2.425, 2.933, 1.871,
               2.154, 3.101, 3.300, 4.007, 3.108, 3.078, 2.108)
)
stopifnot(all(complete.cases(dados)))

##### ==========================================================
##### 2) Codificação
##### ==========================================================
centros <- c(Concentracao = 0.10, Fluxo = 0.50, Distancia = 17.5, Voltagem = 20)
passos  <- c(Concentracao = 0.02, Fluxo = 0.30, Distancia = 2.5, Voltagem = 5)

dados_coded <- coded.data(
  dados,
  C_cod ~ (Concentracao - centros["Concentracao"]) / passos["Concentracao"],
  F_cod ~ (Fluxo        - centros["Fluxo"])        / passos["Fluxo"],
  D_cod ~ (Distancia    - centros["Distancia"])    / passos["Distancia"],
  V_cod ~ (Voltagem     - centros["Voltagem"])     / passos["Voltagem"]
)

decode_real <- function(C_cod=0, F_cod=0, D_cod=0, V_cod=0) {
  c(Concentracao = centros["Concentracao"] + passos["Concentracao"]*C_cod,
    Fluxo        = centros["Fluxo"]        + passos["Fluxo"]       *F_cod,
    Distancia    = centros["Distancia"]    + passos["Distancia"]   *D_cod,
    Voltagem     = centros["Voltagem"]     + passos["Voltagem"]    *V_cod)
}

##### ==========================================================
##### 3) Modelo quadrático completo (RSM)
##### ==========================================================
modelo_full <- rsm(Diametro ~ SO(C_cod, F_cod, D_cod, V_cod), data = dados_coded)

cat("\n=== SUMMARY (modelo_full) ===\n"); print(summary(modelo_full))
cat("\n=== ANOVA (modelo_full) ===\n"); print(anova(modelo_full))

# Ponto estacionário (codificado e real)
xs <- canonical(modelo_full)$xs
cat("\n--- Ponto estacionário (CODIFICADO) ---\n"); print(xs)
xs_real <- decode_real(C_cod = xs["C_cod"], F_cod = xs["F_cod"],
                       D_cod = xs["D_cod"], V_cod = xs["V_cod"])
cat("\n--- Ponto estacionário (REAIS) ---\n"); print(xs_real)

cat("\n--- Equação (variáveis REAIS) ---\n")
print(summary(modelo_full, decode = TRUE))

##### ==========================================================
##### 4) Métricas de erro (MAPE, MSE, RMSE)
##### ==========================================================
pred_full <- as.numeric(predict(modelo_full))
obs <- as.numeric(dados$Diametro)
erro <- pred_full - obs
MAPE_full <- mean(abs(erro/obs))*100
MSE_full  <- mean(erro^2)
RMSE_full <- sqrt(MSE_full)
cat(sprintf("\nFULL  -> MAPE=%.2f%%  MSE=%.4f  RMSE=%.4f\n", MAPE_full, MSE_full, RMSE_full))

##### ==========================================================
##### 5) Modelo reduzido (C e V)
##### ==========================================================
modelo_reduzido <- rsm(Diametro ~ FO(C_cod, V_cod) + PQ(C_cod, V_cod),
                       data = dados_coded)

cat("\n=== SUMMARY (modelo_reduzido) ===\n"); print(summary(modelo_reduzido))
cat("\n=== ANOVA (modelo_reduzido) ===\n"); print(anova(modelo_reduzido))

pred_red <- as.numeric(predict(modelo_reduzido))
erro_red <- pred_red - obs
MAPE_red <- mean(abs(erro_red/obs))*100
MSE_red  <- mean(erro_red^2)
RMSE_red <- sqrt(MSE_red)
cat(sprintf("\nRED   -> MAPE=%.2f%%  MSE=%.4f  RMSE=%.4f\n", MAPE_red, MSE_red, RMSE_red))

# Ponto estacionário do reduzido
xc <- canonical(modelo_reduzido)$xs
cat("\n--- Ponto estacionário (REDUZIDO, CODIFICADO) ---\n"); print(xc)
xc_real <- decode_real(C_cod = xc["C_cod"], F_cod = 0, D_cod = 0, V_cod = xc["V_cod"])
cat("\n--- Ponto estacionário (REDUZIDO, REAIS) ---\n"); print(xc_real)

cat("\n--- Equação (REDUZIDO em variáveis REAIS) ---\n")
print(summary(modelo_reduzido, decode = TRUE))

##### ==========================================================
##### 6) Comparação e escolha do modelo para mapas
##### ==========================================================
comp_tab <- data.frame(
  Modelo = c("Completo","Reduzido"),
  R2_aj  = c(summary(modelo_full)$adj.r.squared, summary(modelo_reduzido)$adj.r.squared),
  RMSE   = c(RMSE_full, RMSE_red),
  MAPE   = c(MAPE_full, MAPE_red),
  AIC    = c(AIC(modelo_full), AIC(modelo_reduzido)),
  BIC    = c(BIC(modelo_full), BIC(modelo_reduzido))
)
print(comp_tab)

teste_lr <- anova(modelo_reduzido, modelo_full)
cat("\n=== Teste de razão de verossimilhança (reduzido vs completo) ===\n"); print(teste_lr)

usa_reduzido <- !is.na(teste_lr$`Pr(>F)`[2]) && teste_lr$`Pr(>F)`[2] > 0.05 && BIC(modelo_reduzido) <= BIC(modelo_full)
modelo_para_mapas <- if (usa_reduzido) modelo_reduzido else modelo_full
modelo_nome <- if (usa_reduzido) "Modelo reduzido" else "Modelo completo"
cat(sprintf("\n*** Usando %s para mapas/superfícies. ***\n", modelo_nome))

##### ==========================================================
##### 7) Diagnósticos de resíduos (apenas melhor modelo)
##### ==========================================================
pred_best  <- as.numeric(predict(modelo_para_mapas))
resid_best <- as.numeric(residuals(modelo_para_mapas))
aux <- data.frame(pred = pred_best, resid = resid_best, obs = dados$Diametro)

g_res1 <- ggplot(aux, aes(x = pred, y = resid)) +
  geom_point() + geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title="Resíduos vs Ajustado", x="Ajustado", y="Resíduo")

g_res2 <- ggplot(aux, aes(sample = resid)) +
  stat_qq() + stat_qq_line() + labs(title="QQ-Plot dos Resíduos")

g_res3 <- ggplot(aux, aes(x = resid)) +
  geom_histogram(bins = 15, color="black") +
  labs(title="Histograma dos Resíduos", x="Resíduo", y="Frequência")

g_res4 <- ggplot(aux, aes(x = pred, y = obs)) +
  geom_point() + geom_abline(slope=1, intercept=0, linetype="dashed") +
  labs(title="Observado vs Ajustado", x="Ajustado", y="Observado")

print((g_res1 + g_res2) / (g_res3 + g_res4))

##### ==========================================================
##### 8) Funções auxiliares para contorno e superfície
##### ==========================================================
gera_grade_real <- function(f1, f2, n = 60, k = 1.0) {
  seq1 <- seq(centros[f1] - k*passos[f1]*2, centros[f1] + k*passos[f1]*2, length.out = n)
  seq2 <- seq(centros[f2] - k*passos[f2]*2, centros[f2] + k*passos[f2]*2, length.out = n)
  grid <- expand.grid(x1 = seq1, x2 = seq2, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  names(grid) <- c(f1, f2)
  outros <- setdiff(names(centros), c(f1, f2))
  for (nm in outros) grid[[nm]] <- centros[nm]
  grid
}

codificar_grade <- function(grid) {
  transform(grid,
    C_cod = (Concentracao - centros["Concentracao"])/passos["Concentracao"],
    F_cod = (Fluxo        - centros["Fluxo"]       )/passos["Fluxo"],
    D_cod = (Distancia    - centros["Distancia"]   )/passos["Distancia"],
    V_cod = (Voltagem     - centros["Voltagem"]    )/passos["Voltagem"]
  )
}

plot_contorno <- function(modelo, f1, f2, n = 80, k = 1.0) {
  stopifnot(f1 %in% names(centros), f2 %in% names(centros), f1 != f2)
  grid <- codificar_grade(gera_grade_real(f1, f2, n = n, k = k))
  grid$yhat <- as.numeric(predict(modelo, newdata = grid))
  ggplot(grid, aes_string(x = f1, y = f2, z = "yhat")) +
    geom_contour_filled(bins = 12) +
    scale_fill_viridis_d(name = "Diâmetro (µm)", option = "D") +
    labs(title = paste0("Mapa de contorno: ", f1, " × ", f2, " (", modelo_nome, ")"),
         x = c(Concentracao="Concentração (g/mL)", Fluxo="Fluxo (mL/h)",
               Distancia="Distância (cm)", Voltagem="Voltagem (kV)")[f1],
         y = c(Concentracao="Concentração (g/mL)", Fluxo="Fluxo (mL/h)",
               Distancia="Distância (cm)", Voltagem="Voltagem (kV)")[f2]) +
    theme(legend.position = "right")
}

plot_superficie <- function(modelo, f1, f2, n = 60, k = 1.0) {
  grid <- codificar_grade(gera_grade_real(f1, f2, n = n, k = k))
  grid$yhat <- as.numeric(predict(modelo, newdata = grid))
  seq1 <- sort(unique(grid[[f1]]))
  seq2 <- sort(unique(grid[[f2]]))
  Z <- matrix(grid$yhat, nrow = length(seq1), ncol = length(seq2))
  if (plotly_ok) {
    plot_ly(x = seq1, y = seq2, z = Z) |>
      add_surface() |>
      layout(
        title = paste0("Superfície: ", f1, " × ", f2, " (", modelo_nome, ")"),
        scene = list(
          xaxis = list(title = c(Concentracao="Concentração (g/mL)", Fluxo="Fluxo (mL/h)",
                                 Distancia="Distância (cm)", Voltagem="Voltagem (kV)")[f1]),
          yaxis = list(title = c(Concentracao="Concentração (g/mL)", Fluxo="Fluxo (mL/h)",
                                 Distancia="Distância (cm)", Voltagem="Voltagem (kV)")[f2]),
          zaxis = list(title = "Diâmetro (µm)")
        )
      )
  } else {
    persp(x = seq1, y = seq2, z = Z, theta = 45, phi = 25,
          xlab = f1, ylab = f2, zlab = "Diâmetro (µm)",
          ticktype = "detailed", col = "lightblue",
          main = paste0("Superfície: ", f1, " × ", f2, " (", modelo_nome, ")"))
  }
}

##### ==========================================================
##### 9) Mapas de contorno e superfícies (somente melhor modelo)
##### ==========================================================
fatores <- c("Concentracao","Fluxo","Distancia","Voltagem")
pares <- combn(fatores, 2, simplify = FALSE)

plots_contorno <- lapply(pares, function(p) plot_contorno(modelo_para_mapas, p[1], p[2], n=100, k=1.0))
if (length(plots_contorno) == 6) {
  suppressWarnings({
    print((plots_contorno[[1]] | plots_contorno[[2]]) /
          (plots_contorno[[3]] | plots_contorno[[4]]) /
          (plots_contorno[[5]] | plots_contorno[[6]]))
  })
} else {
  lapply(plots_contorno, print)
}

superficies <- lapply(pares, function(p) plot_superficie(modelo_para_mapas, p[1], p[2], n=70, k=1.0))
invisible(superficies)

##### ==========================================================
##### 10) Previsões de referência
##### ==========================================================
y_centro <- as.numeric(predict(modelo_para_mapas, newdata = data.frame(C_cod=0, F_cod=0, D_cod=0, V_cod=0)))
cat(sprintf("\nPrevisão no ponto central: %.3f µm\n", y_centro))

y_xs <- as.numeric(predict(modelo_para_mapas, newdata = data.frame(
  C_cod = xs["C_cod"], F_cod = xs["F_cod"], D_cod = xs["D_cod"], V_cod = xs["V_cod"]
)))
cat(sprintf("Previsão no ponto estacionário do modelo completo: %.3f µm\n", y_xs))
