# ==============================================================================
# TFG ECONOMÍA — Bitcoin y Oro como Activos Refugio
# Autor: Javier Lara
# Periodo de análisis: 01/01/2018 - 28/02/2026
# Lenguaje: R
# ==============================================================================
# ESTRUCTURA DEL SCRIPT
# 0. Paquetes y preparación
# 1. Funciones auxiliares
# 2. Descarga y preparación de datos
# 3. Tasa libre de riesgo
# 4. Tabla 1 — Estadísticos descriptivos y Ratio de Sharpe
# 5. Tabla 2 — Test ADF
# 6. Tabla 3 — Tests Ljung-Box y ARCH-LM
# 7. Tabla 4 — Correlaciones por subperiodo y correlación móvil
# 8. Tabla 5 — Modelo principal con dummies (periodo completo)
# 9. Tabla 6 — Robustez con VIX
# 10. Tablas 7a y 7b — Subperiodos
# 11. Figuras 1 a 4
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. PAQUETES Y PREPARACIÓN
# ------------------------------------------------------------------------------

# Instalar si faltan:
# install.packages(c(
#   "quantmod", "xts", "zoo", "dplyr", "tidyr", "lmtest", "sandwich",
#   "tseries", "ggplot2", "scales", "FinTS"
# ))

suppressPackageStartupMessages({
  library(quantmod)
  library(xts)
  library(zoo)
  library(dplyr)
  library(tidyr)
  library(lmtest)
  library(sandwich)
  library(tseries)
  library(ggplot2)
  library(scales)
  library(FinTS)
})

options(scipen = 999)

dir.create("salidas", showWarnings = FALSE)
dir.create("figuras", showWarnings = FALSE)

# ------------------------------------------------------------------------------
# 1. FUNCIONES AUXILIARES
# ------------------------------------------------------------------------------

sig_stars <- function(p) {
  dplyr::case_when(
    is.na(p)      ~ NA_character_,
    p < 0.01      ~ "***",
    p < 0.05      ~ "**",
    p < 0.10      ~ "*",
    TRUE          ~ "n.s."
  )
}

fmt_p <- function(p, digits = 4) {
  ifelse(
    is.na(p),
    NA_character_,
    ifelse(
      p < 0.0001,
      "<0.0001",
      format(round(p, digits), nsmall = digits, scientific = FALSE, trim = TRUE)
    )
  )
}

clasificar_activo <- function(beta_crisis, p_beta_crisis_bi, tol = 1e-12) {
  if (is.na(beta_crisis) || is.na(p_beta_crisis_bi)) return("Sin evidencia")

  if (beta_crisis < -tol && p_beta_crisis_bi < 0.05) {
    return("Refugio fuerte")
  } else if (beta_crisis <= tol && p_beta_crisis_bi >= 0.05) {
    return("Refugio débil")
  } else if (beta_crisis > tol && p_beta_crisis_bi < 0.05) {
    return("Procíclico")
  } else {
    return("Sin evidencia")
  }
}

estimar_dummy <- function(activo, indice, umbral_pct,
                          nombre_activo, nombre_indice,
                          datos_uso) {

  r_a <- as.numeric(datos_uso[, activo])
  r_m <- as.numeric(datos_uso[, indice])

  n <- length(r_a)
  if (n < 30) stop("Muestra demasiado pequeña para estimar el modelo.")

  umbral <- as.numeric(quantile(r_m, probs = umbral_pct, na.rm = TRUE))
  D_t <- as.numeric(r_m <= umbral)
  inter <- D_t * r_m

  df_mod <- data.frame(r_a = r_a, r_m = r_m, inter = inter)
  mod <- lm(r_a ~ r_m + inter, data = df_mod)

  k <- length(coef(mod))
  df_res <- n - k
  lag_nw <- max(floor(4 * (n / 100)^(2/9)), 1)

  V <- tryCatch(
    NeweyWest(mod, lag = lag_nw, prewhite = FALSE),
    error = function(e) vcov(mod)
  )

  b <- coef(mod)
  b0 <- unname(b[1])
  b1 <- unname(b[2])
  b2 <- unname(b[3])

  se0 <- sqrt(V[1, 1])
  se1 <- sqrt(V[2, 2])
  se2 <- sqrt(V[3, 3])

  t0 <- b0 / se0
  t1 <- b1 / se1
  t2 <- b2 / se2

  p0_bi <- 2 * pt(abs(t0), df = df_res, lower.tail = FALSE)
  p1_bi <- 2 * pt(abs(t1), df = df_res, lower.tail = FALSE)
  p2_bi <- 2 * pt(abs(t2), df = df_res, lower.tail = FALSE)

  beta_crisis <- b1 + b2

  var_beta_crisis <- V[2, 2] + V[3, 3] + 2 * V[2, 3]
  if (var_beta_crisis < 0 && abs(var_beta_crisis) < 1e-12) {
    var_beta_crisis <- 0
  }

  se_beta_crisis <- sqrt(var_beta_crisis)
  t_beta_crisis <- beta_crisis / se_beta_crisis

  p_beta_crisis_bi <- 2 * pt(abs(t_beta_crisis), df = df_res, lower.tail = FALSE)

  clasificacion <- clasificar_activo(beta_crisis, p_beta_crisis_bi)

  data.frame(
    Activo = nombre_activo,
    Indice = nombre_indice,
    Umbral = paste0(umbral_pct * 100, "%"),
    N_crisis = sum(D_t),

    Beta0 = b0,
    p_Beta0_bi = p0_bi,

    Beta1 = b1,
    p_Beta1_bi = p1_bi,
    Sig_Beta1 = sig_stars(p1_bi),

    Beta2 = b2,
    p_Beta2_bi = p2_bi,
    Sig_Beta2 = sig_stars(p2_bi),

    Beta_crisis = beta_crisis,
    SE_Beta_crisis = se_beta_crisis,
    p_Beta_crisis_bi = p_beta_crisis_bi,
    Sig_Beta_crisis = sig_stars(p_beta_crisis_bi),

    Clasificacion = clasificacion,
    stringsAsFactors = FALSE
  )
}

estimar_vix <- function(activo, nombre_activo, datos_uso) {

  r_a <- as.numeric(datos_uso[, activo])
  r_m <- as.numeric(datos_uso[, "SP500"])
  D_t <- as.numeric(datos_uso[, "DUMMY_VIX"])

  n <- length(r_a)
  inter <- D_t * r_m

  df_mod <- data.frame(r_a = r_a, r_m = r_m, inter = inter)
  mod <- lm(r_a ~ r_m + inter, data = df_mod)

  k <- length(coef(mod))
  df_res <- n - k
  lag_nw <- max(floor(4 * (n / 100)^(2/9)), 1)

  V <- tryCatch(
    NeweyWest(mod, lag = lag_nw, prewhite = FALSE),
    error = function(e) vcov(mod)
  )

  b <- coef(mod)
  b1 <- unname(b[2])
  b2 <- unname(b[3])

  se1 <- sqrt(V[2, 2])
  se2 <- sqrt(V[3, 3])

  t1 <- b1 / se1
  t2 <- b2 / se2

  p1_bi <- 2 * pt(abs(t1), df = df_res, lower.tail = FALSE)
  p2_bi <- 2 * pt(abs(t2), df = df_res, lower.tail = FALSE)

  beta_crisis <- b1 + b2
  var_beta_crisis <- V[2, 2] + V[3, 3] + 2 * V[2, 3]
  if (var_beta_crisis < 0 && abs(var_beta_crisis) < 1e-12) {
    var_beta_crisis <- 0
  }
  se_beta_crisis <- sqrt(var_beta_crisis)
  t_beta_crisis <- beta_crisis / se_beta_crisis
  p_beta_crisis_bi <- 2 * pt(abs(t_beta_crisis), df = df_res, lower.tail = FALSE)

  clasificacion <- clasificar_activo(beta_crisis, p_beta_crisis_bi)

  data.frame(
    Activo = nombre_activo,
    N_crisis = sum(D_t),

    Beta1 = b1,
    p_Beta1_bi = p1_bi,
    Sig_Beta1 = sig_stars(p1_bi),

    Beta2 = b2,
    p_Beta2_bi = p2_bi,
    Sig_Beta2 = sig_stars(p2_bi),

    Beta_crisis = beta_crisis,
    SE_Beta_crisis = se_beta_crisis,
    p_Beta_crisis_bi = p_beta_crisis_bi,
    Sig_Beta_crisis = sig_stars(p_beta_crisis_bi),

    Clasificacion = clasificacion,
    stringsAsFactors = FALSE
  )
}

calcular_cor_movil <- function(activo, indice, ventana, datos_uso) {
  rollapply(
    datos_uso[, c(activo, indice)],
    width = ventana,
    FUN = function(x) cor(x[, 1], x[, 2]),
    by.column = FALSE,
    align = "right"
  )
}

# ------------------------------------------------------------------------------
# 2. DESCARGA Y PREPARACIÓN DE DATOS
# ------------------------------------------------------------------------------

getSymbols(
  c("BTC-USD", "GLD", "^GSPC", "^STOXX50E", "^N225", "^IRX", "^VIX"),
  from = "2018-01-01",
  to   = "2026-02-28",
  src  = "yahoo",
  auto.assign = TRUE
)

# Precios de cierre ajustados y cierres
btc_raw  <- Ad(`BTC-USD`)
gld_raw  <- Ad(GLD)
sp_raw   <- Ad(GSPC)
euro_raw <- Ad(STOXX50E)
nik_raw  <- Ad(N225)
irx_raw  <- Cl(IRX)
vix_raw  <- Cl(VIX)

# Rendimientos logarítmicos diarios
btc_ret  <- diff(log(btc_raw), lag = 1)
gld_ret  <- diff(log(gld_raw), lag = 1)
sp_ret   <- diff(log(sp_raw), lag = 1)
euro_ret <- diff(log(euro_raw), lag = 1)
nik_ret  <- diff(log(nik_raw), lag = 1)

# Homogeneización estricta por intersección de fechas
lista_series <- list(
  btc_ret, gld_ret, sp_ret, euro_ret, nik_ret, irx_raw, vix_raw
)

datos_raw <- Reduce(function(x, y) merge(x, y, join = "inner"), lista_series)
colnames(datos_raw) <- c("BTC", "ORO", "SP500", "EURO", "NIKKEI", "IRX", "VIX")

datos <- na.omit(datos_raw)

cat("Observaciones totales:", nrow(datos), "\n")
cat("Periodo:", as.character(index(datos)[1]), "a", as.character(index(datos)[nrow(datos)]), "\n")

# ------------------------------------------------------------------------------
# 3. TASA LIBRE DE RIESGO
# ------------------------------------------------------------------------------

tasa_rf_anual  <- mean(as.numeric(datos[, "IRX"]), na.rm = TRUE) / 100
tasa_rf_diaria <- tasa_rf_anual / 252

cat("Tasa libre de riesgo anual media:", round(tasa_rf_anual * 100, 4), "%\n")
cat("Tasa libre de riesgo diaria:", round(tasa_rf_diaria * 100, 6), "%\n")

# ------------------------------------------------------------------------------
# 4. TABLA 1 — ESTADÍSTICOS DESCRIPTIVOS Y RATIO DE SHARPE
# ------------------------------------------------------------------------------

activos_desc <- c("BTC", "ORO", "SP500")
nombres_desc <- c("Bitcoin", "Oro", "S&P 500")

tabla1_completa <- data.frame()

for (i in seq_along(activos_desc)) {
  serie <- as.numeric(datos[, activos_desc[i]])
  media <- mean(serie, na.rm = TRUE)
  desv  <- sd(serie, na.rm = TRUE)
  asim  <- mean((serie - media)^3, na.rm = TRUE) / desv^3
  kurt_exc <- mean((serie - media)^4, na.rm = TRUE) / desv^4 - 3
  volat_anual <- desv * sqrt(252) * 100
  sharpe <- (media - tasa_rf_diaria) / desv * sqrt(252)

  tabla1_completa <- rbind(
    tabla1_completa,
    data.frame(
      Activo = nombres_desc[i],
      Media_diaria = media,
      Desv_Est = desv,
      Asimetria = asim,
      Curtosis = kurt_exc,
      Volat_anual_pct = volat_anual,
      Sharpe_anual = sharpe,
      stringsAsFactors = FALSE
    )
  )
}

tabla1_completa <- tabla1_completa %>%
  mutate(
    Media_diaria = round(Media_diaria, 6),
    Desv_Est = round(Desv_Est, 6),
    Asimetria = round(Asimetria, 4),
    Curtosis = round(Curtosis, 4),
    Volat_anual_pct = round(Volat_anual_pct, 2),
    Sharpe_anual = round(Sharpe_anual, 4)
  )

tabla1_cuerpo <- tabla1_completa %>%
  mutate(
    Media_diaria = ifelse(Activo == "S&P 500", "-", sprintf("%.6f", Media_diaria)),
    Desv_Est = ifelse(Activo == "S&P 500", "-", sprintf("%.6f", Desv_Est)),
    Asimetria = ifelse(Activo == "S&P 500", "-", sprintf("%.4f", Asimetria)),
    Curtosis = ifelse(Activo == "S&P 500", "-", sprintf("%.4f", Curtosis)),
    Volat_anual_pct = ifelse(Activo == "S&P 500", "-", paste0(sprintf("%.2f", Volat_anual_pct), "%")),
    Sharpe_anual = sprintf("%.4f", Sharpe_anual)
  )

cat("\n--- TABLA 1: ESTADÍSTICOS DESCRIPTIVOS Y RATIO DE SHARPE ---\n")
print(tabla1_cuerpo, row.names = FALSE)

write.csv(tabla1_completa, "salidas/tabla1_descriptivos_completa.csv", row.names = FALSE)
write.csv(tabla1_cuerpo, "salidas/tabla1_descriptivos_cuerpo.csv", row.names = FALSE)

# ------------------------------------------------------------------------------
# 5. TABLA 2 — TEST DE DICKEY-FULLER AUMENTADO
# ------------------------------------------------------------------------------

activos_adf <- c("BTC", "ORO", "SP500", "EURO", "NIKKEI")
nombres_adf <- c("Bitcoin", "Oro", "S&P 500", "Eurostoxx 50", "Nikkei 225")

tabla2_adf <- data.frame()

for (i in seq_along(activos_adf)) {
  serie <- as.numeric(datos[, activos_adf[i]])
  res_adf <- adf.test(serie, alternative = "stationary")

  tabla2_adf <- rbind(
    tabla2_adf,
    data.frame(
      Activo = nombres_adf[i],
      Estadistico_ADF = round(res_adf$statistic, 3),
      P_valor = ifelse(res_adf$p.value < 0.01, "<0.01", round(res_adf$p.value, 4)),
      Resultado = ifelse(res_adf$p.value < 0.05, "Estacionaria", "No estacionaria"),
      stringsAsFactors = FALSE
    )
  )
}

cat("\n--- TABLA 2: TEST ADF ---\n")
print(tabla2_adf, row.names = FALSE)

write.csv(tabla2_adf, "salidas/tabla2_adf.csv", row.names = FALSE)

# ------------------------------------------------------------------------------
# 6. TABLA 3 — TESTS LJUNG-BOX Y ARCH-LM
# ------------------------------------------------------------------------------

tabla3_arch <- data.frame()

for (i in seq_along(activos_adf)) {
  serie <- as.numeric(datos[, activos_adf[i]])

  lb_test  <- Box.test(serie^2, lag = 10, type = "Ljung-Box")
  arch_lm  <- FinTS::ArchTest(serie, lags = 5)

  tabla3_arch <- rbind(
    tabla3_arch,
    data.frame(
      Activo = nombres_adf[i],
      Ljung_Box = round(lb_test$statistic, 2),
      p_LB = ifelse(lb_test$p.value < 0.0001, "<0.0001", round(lb_test$p.value, 4)),
      ARCH_LM = round(arch_lm$statistic, 2),
      p_ARCH = ifelse(arch_lm$p.value < 0.0001, "<0.0001", round(arch_lm$p.value, 4)),
      Heterocedasticidad = ifelse(arch_lm$p.value < 0.05, "Sí", "No"),
      stringsAsFactors = FALSE
    )
  )
}

cat("\n--- TABLA 3: TESTS LJUNG-BOX Y ARCH-LM ---\n")
print(tabla3_arch, row.names = FALSE)

write.csv(tabla3_arch, "salidas/tabla3_arch_ljungbox.csv", row.names = FALSE)

# ------------------------------------------------------------------------------
# 7. TABLA 4 — CORRELACIONES POR SUBPERIODO Y CORRELACIÓN MÓVIL
# ------------------------------------------------------------------------------

fase1_inicio <- as.Date("2018-01-01")
fase1_fin    <- as.Date("2020-02-28")
fase2_inicio <- as.Date("2020-03-01")
fase2_fin    <- as.Date("2022-12-31")
fase3_inicio <- as.Date("2023-01-01")
fase3_fin    <- as.Date("2026-02-28")

fechas <- index(datos)
f1 <- fechas >= fase1_inicio & fechas <= fase1_fin
f2 <- fechas >= fase2_inicio & fechas <= fase2_fin
f3 <- fechas >= fase3_inicio & fechas <= fase3_fin

activos_cor <- c("BTC", "ORO")
indices_cor <- c("SP500", "EURO", "NIKKEI")
nombres_act <- c("Bitcoin", "Oro")
nombres_idx <- c("S&P 500", "Eurostoxx 50", "Nikkei 225")

tabla4_subperiodos <- data.frame()

for (a in seq_along(activos_cor)) {
  for (i in seq_along(indices_cor)) {
    serie_a <- as.numeric(datos[, activos_cor[a]])
    serie_i <- as.numeric(datos[, indices_cor[i]])

    tabla4_subperiodos <- rbind(
      tabla4_subperiodos,
      data.frame(
        Activo = nombres_act[a],
        Indice = nombres_idx[i],
        Fase1_PreCOVID = round(cor(serie_a[f1], serie_i[f1]), 3),
        Fase2_Pandemia = round(cor(serie_a[f2], serie_i[f2]), 3),
        Fase3_ETFs = round(cor(serie_a[f3], serie_i[f3]), 3),
        Total = round(cor(serie_a, serie_i), 3),
        stringsAsFactors = FALSE
      )
    )
  }
}

cat("\n--- TABLA 4: CORRELACIONES POR SUBPERIODO ---\n")
print(tabla4_subperiodos, row.names = FALSE)

write.csv(tabla4_subperiodos, "salidas/tabla4_correlaciones_subperiodo.csv", row.names = FALSE)

# Correlaciones móviles para las figuras
cor_btc_sp_30  <- calcular_cor_movil("BTC", "SP500", 30, datos)
cor_btc_sp_60  <- calcular_cor_movil("BTC", "SP500", 60, datos)
cor_btc_sp_120 <- calcular_cor_movil("BTC", "SP500", 120, datos)
cor_oro_sp_60  <- calcular_cor_movil("ORO", "SP500", 60, datos)

cor_btc_eu_30  <- calcular_cor_movil("BTC", "EURO", 30, datos)
cor_btc_eu_60  <- calcular_cor_movil("BTC", "EURO", 60, datos)
cor_btc_eu_120 <- calcular_cor_movil("BTC", "EURO", 120, datos)
cor_oro_eu_60  <- calcular_cor_movil("ORO", "EURO", 60, datos)

cor_btc_nk_30  <- calcular_cor_movil("BTC", "NIKKEI", 30, datos)
cor_btc_nk_60  <- calcular_cor_movil("BTC", "NIKKEI", 60, datos)
cor_btc_nk_120 <- calcular_cor_movil("BTC", "NIKKEI", 120, datos)
cor_oro_nk_60  <- calcular_cor_movil("ORO", "NIKKEI", 60, datos)

# ------------------------------------------------------------------------------
# 8. TABLA 5 — MODELO PRINCIPAL CON DUMMIES (PERIODO COMPLETO)
# ------------------------------------------------------------------------------

activos <- c("BTC", "ORO")
indices <- c("SP500", "EURO", "NIKKEI")
umbrales <- c(0.01, 0.05, 0.10)

nomb_act <- c("Bitcoin", "Oro")
nomb_idx <- c("S&P 500", "Eurostoxx 50", "Nikkei 225")

resultados_principal <- list()

for (a in 1:2) {
  for (i in 1:3) {
    for (u in 1:3) {
      clave <- paste0(activos[a], "_", indices[i], "_p", umbrales[u] * 100)
      resultados_principal[[clave]] <- estimar_dummy(
        activo = activos[a],
        indice = indices[i],
        umbral_pct = umbrales[u],
        nombre_activo = nomb_act[a],
        nombre_indice = nomb_idx[i],
        datos_uso = datos
      )
    }
  }
}

tabla5_tecnica <- do.call(rbind, resultados_principal)

tabla5_cuerpo <- tabla5_tecnica %>%
  transmute(
    Activo,
    Indice,
    Umbral,
    N_crisis,
    Beta_normal = round(Beta1, 4),
    p_Beta_normal = fmt_p(p_Beta1_bi),
    Beta_crisis = round(Beta_crisis, 4),
    p_Beta_crisis = fmt_p(p_Beta_crisis_bi),
    Clasificacion
  )

cat("\n--- TABLA 5: MODELO PRINCIPAL (PERIODO COMPLETO) ---\n")
print(tabla5_cuerpo, row.names = FALSE)

write.csv(tabla5_tecnica, "salidas/tabla5_modelo_principal_tecnica.csv", row.names = FALSE)
write.csv(tabla5_cuerpo, "salidas/tabla5_modelo_principal_cuerpo.csv", row.names = FALSE)

# ------------------------------------------------------------------------------
# 9. TABLA 6 — ROBUSTEZ CON DUMMY VIX
# ------------------------------------------------------------------------------

vix_serie <- as.numeric(datos[, "VIX"])
vix_p95 <- as.numeric(quantile(vix_serie, probs = 0.95, na.rm = TRUE))

datos_vix <- cbind(datos, DUMMY_VIX = as.numeric(vix_serie >= vix_p95))

cat("\nUmbral VIX p95:", round(vix_p95, 2), "\n")
cat("Días con VIX extremo:", sum(as.numeric(datos_vix[, "DUMMY_VIX"])), "\n")

tabla6_tecnica <- rbind(
  estimar_vix("BTC", "Bitcoin", datos_vix),
  estimar_vix("ORO", "Oro", datos_vix)
)

tabla6_cuerpo <- tabla6_tecnica %>%
  transmute(
    Activo,
    Beta_normal = round(Beta1, 4),
    p_Beta_normal = fmt_p(p_Beta1_bi),
    Beta_crisis = round(Beta_crisis, 4),
    p_Beta_crisis = fmt_p(p_Beta_crisis_bi),
    Clasificacion
  )

cat("\n--- TABLA 6: ROBUSTEZ CON DUMMY VIX (95%) ---\n")
print(tabla6_cuerpo, row.names = FALSE)

write.csv(tabla6_tecnica, "salidas/tabla6_robustez_vix_tecnica.csv", row.names = FALSE)
write.csv(tabla6_cuerpo, "salidas/tabla6_robustez_vix_cuerpo.csv", row.names = FALSE)

# ------------------------------------------------------------------------------
# 10. TABLAS 7a Y 7b — ANÁLISIS POR SUBPERIODOS
# ------------------------------------------------------------------------------

fases <- list(
  list(nombre = "1. Pre-COVID", inicio = fase1_inicio, fin = fase1_fin),
  list(nombre = "2. Pandemia",  inicio = fase2_inicio, fin = fase2_fin),
  list(nombre = "3. ETFs",      inicio = fase3_inicio, fin = fase3_fin)
)

resultados_sub <- list()

for (fase in fases) {
  rango <- paste0(as.character(fase$inicio), "/", as.character(fase$fin))
  datos_sub <- na.omit(datos[rango])

  for (a in 1:2) {
    for (i in 1:3) {
      res <- estimar_dummy(
        activo = activos[a],
        indice = indices[i],
        umbral_pct = 0.05,
        nombre_activo = nomb_act[a],
        nombre_indice = nomb_idx[i],
        datos_uso = datos_sub
      )

      res$Fase <- fase$nombre
      res$Obs_subperiodo <- nrow(datos_sub)

      clave <- paste0(fase$nombre, "_", activos[a], "_", indices[i])
      resultados_sub[[clave]] <- res
    }
  }

  cat("✓ Completada:", fase$nombre, "\n")
}

tabla7_tecnica <- do.call(rbind, resultados_sub)

tabla7a_bitcoin <- tabla7_tecnica %>%
  filter(Activo == "Bitcoin") %>%
  transmute(
    Fase,
    Indice,
    Beta_normal = round(Beta1, 4),
    p_Beta_normal = fmt_p(p_Beta1_bi),
    Beta_crisis = round(Beta_crisis, 4),
    p_Beta_crisis = fmt_p(p_Beta_crisis_bi),
    Clasificacion
  )

tabla7b_oro <- tabla7_tecnica %>%
  filter(Activo == "Oro") %>%
  transmute(
    Fase,
    Indice,
    Beta_normal = round(Beta1, 4),
    p_Beta_normal = fmt_p(p_Beta1_bi),
    Beta_crisis = round(Beta_crisis, 4),
    p_Beta_crisis = fmt_p(p_Beta_crisis_bi),
    Clasificacion
  )

cat("\n--- TABLA 7a: BITCOIN POR FASES (5%) ---\n")
print(tabla7a_bitcoin, row.names = FALSE)

cat("\n--- TABLA 7b: ORO POR FASES (5%) ---\n")
print(tabla7b_oro, row.names = FALSE)

write.csv(tabla7_tecnica, "salidas/tabla7_subperiodos_tecnica.csv", row.names = FALSE)
write.csv(tabla7a_bitcoin, "salidas/tabla7a_bitcoin_subperiodos.csv", row.names = FALSE)
write.csv(tabla7b_oro, "salidas/tabla7b_oro_subperiodos.csv", row.names = FALSE)

cat("\n✓ ANÁLISIS NUMÉRICO FINALIZADO\n")

# ==============================================================================
# 11. GENERACIÓN DE FIGURAS
# ==============================================================================

# ------------------------------------------------------------------------------
# FIGURA 1 — DISTRIBUCIÓN DE RENDIMIENTOS
# ------------------------------------------------------------------------------

btc_vec <- as.numeric(datos[, "BTC"])
oro_vec <- as.numeric(datos[, "ORO"])

df_dist <- data.frame(
  rendimiento = c(btc_vec, oro_vec),
  activo = rep(c("Bitcoin", "Oro"), c(length(btc_vec), length(oro_vec)))
)

params_dist <- df_dist %>%
  group_by(activo) %>%
  summarise(
    media = mean(rendimiento),
    desv = sd(rendimiento),
    .groups = "drop"
  )

fig1 <- ggplot(df_dist, aes(x = rendimiento)) +
  geom_histogram(
    aes(y = after_stat(density), fill = activo),
    bins = 100, alpha = 0.5, color = NA
  ) +
  geom_density(aes(color = activo), linewidth = 0.8) +
  stat_function(
    data = df_dist %>% filter(activo == "Bitcoin"),
    fun = dnorm,
    args = list(
      mean = params_dist$media[params_dist$activo == "Bitcoin"],
      sd   = params_dist$desv[params_dist$activo == "Bitcoin"]
    ),
    color = "#B71C1C", linewidth = 0.8, linetype = "dashed"
  ) +
  stat_function(
    data = df_dist %>% filter(activo == "Oro"),
    fun = dnorm,
    args = list(
      mean = params_dist$media[params_dist$activo == "Oro"],
      sd   = params_dist$desv[params_dist$activo == "Oro"]
    ),
    color = "#E65100", linewidth = 0.8, linetype = "dashed"
  ) +
  facet_wrap(~ activo, ncol = 2, scales = "free") +
  scale_fill_manual(values = c("Bitcoin" = "#2196F3", "Oro" = "#FF9800"), guide = "none") +
  scale_color_manual(values = c("Bitcoin" = "#1565C0", "Oro" = "#E65100"), guide = "none") +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Figura 1 — Distribución de rendimientos diarios: Bitcoin y Oro",
    subtitle = "Histograma de densidad empírica y distribución normal teórica superpuesta",
    x = "Rendimiento logarítmico diario",
    y = "Densidad",
    caption = "Fuente: elaboración propia. Datos Yahoo Finance 2018-2026."
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 8, color = "grey40"),
    strip.text = element_text(face = "bold", size = 10),
    panel.grid.minor = element_blank()
  )

ggsave("figuras/figura1_distribuciones.png", plot = fig1,
       width = 12, height = 5, dpi = 300, bg = "white")

# ------------------------------------------------------------------------------
# FUNCIÓN AUXILIAR — FIGURAS 2a, 2b y 2c
# ------------------------------------------------------------------------------

figura_correlacion <- function(cor_btc_30, cor_btc_60, cor_btc_120,
                               cor_oro_60, titulo, archivo) {

  df_cor_xts <- Reduce(
    function(x, y) merge(x, y, join = "inner"),
    list(cor_btc_30, cor_btc_60, cor_btc_120, cor_oro_60)
  )

  colnames(df_cor_xts) <- c("BTC_30", "BTC_60", "BTC_120", "ORO_60")
  df_cor <- data.frame(fecha = index(df_cor_xts), coredata(df_cor_xts)) %>% na.omit()

  fig <- ggplot(df_cor, aes(x = fecha)) +
    annotate("rect", xmin = fase1_inicio, xmax = fase1_fin,
             ymin = -Inf, ymax = Inf, fill = "#E3F2FD", alpha = 0.4) +
    annotate("rect", xmin = fase2_inicio, xmax = fase2_fin,
             ymin = -Inf, ymax = Inf, fill = "#FFF3E0", alpha = 0.4) +
    annotate("rect", xmin = fase3_inicio, xmax = fase3_fin,
             ymin = -Inf, ymax = Inf, fill = "#E8F5E9", alpha = 0.4) +
    geom_line(aes(y = BTC_30), color = "#90CAF9", linewidth = 0.4, alpha = 0.7) +
    geom_line(aes(y = BTC_60), color = "#1565C0", linewidth = 0.8) +
    geom_line(aes(y = BTC_120), color = "#0D47A1", linewidth = 0.4, linetype = "dashed") +
    geom_line(aes(y = ORO_60), color = "#FF9800", linewidth = 0.8) +
    geom_hline(yintercept = 0, linewidth = 0.5, linetype = "dashed", color = "black") +
    scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.25)) +
    scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
    labs(
      title = titulo,
      subtitle = "Bitcoin: 30d / 60d / 120d | Oro: 60d",
      x = NULL,
      y = "Coeficiente de correlación de Pearson",
      caption = "Fuente: elaboración propia. Datos Yahoo Finance 2018-2026."
    ) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", size = 11),
      plot.subtitle = element_text(size = 8, color = "grey40"),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )

  ggsave(archivo, plot = fig, width = 14, height = 5, dpi = 300, bg = "white")
  invisible(fig)
}

figura_correlacion(
  cor_btc_30 = cor_btc_sp_30,
  cor_btc_60 = cor_btc_sp_60,
  cor_btc_120 = cor_btc_sp_120,
  cor_oro_60 = cor_oro_sp_60,
  titulo = "Figura 2a — Correlación móvil con el S&P 500",
  archivo = "figuras/figura2a_correlacion_sp500.png"
)

figura_correlacion(
  cor_btc_30 = cor_btc_eu_30,
  cor_btc_60 = cor_btc_eu_60,
  cor_btc_120 = cor_btc_eu_120,
  cor_oro_60 = cor_oro_eu_60,
  titulo = "Figura 2b — Correlación móvil con el Eurostoxx 50",
  archivo = "figuras/figura2b_correlacion_eurostoxx.png"
)

figura_correlacion(
  cor_btc_30 = cor_btc_nk_30,
  cor_btc_60 = cor_btc_nk_60,
  cor_btc_120 = cor_btc_nk_120,
  cor_oro_60 = cor_oro_nk_60,
  titulo = "Figura 2c — Correlación móvil con el Nikkei 225",
  archivo = "figuras/figura2c_correlacion_nikkei.png"
)

# ------------------------------------------------------------------------------
# FUNCIÓN AUXILIAR — FIGURAS 3a, 3b y 3c
# ------------------------------------------------------------------------------

figura_zonas_crisis <- function(indice, nombre_display, archivo) {
  serie_i <- as.numeric(datos[, indice])
  fechas_i <- index(datos)

  u10 <- quantile(serie_i, 0.10, na.rm = TRUE)
  u5  <- quantile(serie_i, 0.05, na.rm = TRUE)
  u1  <- quantile(serie_i, 0.01, na.rm = TRUE)

  y_min <- min(serie_i) * 1.05
  df_plot <- data.frame(fecha = fechas_i, rendimiento = serie_i)

  fig <- ggplot(df_plot, aes(x = fecha, y = rendimiento)) +
    annotate("rect", xmin = min(fechas_i), xmax = max(fechas_i),
             ymin = u10, ymax = u5, fill = "#FFF9C4", alpha = 0.8) +
    annotate("rect", xmin = min(fechas_i), xmax = max(fechas_i),
             ymin = u5, ymax = u1, fill = "#FFE0B2", alpha = 0.8) +
    annotate("rect", xmin = min(fechas_i), xmax = max(fechas_i),
             ymin = y_min, ymax = u1, fill = "#FFCDD2", alpha = 0.8) +
    geom_hline(yintercept = u10, color = "#F9A825", linewidth = 0.7) +
    geom_hline(yintercept = u5,  color = "#E65100", linewidth = 0.7) +
    geom_hline(yintercept = u1,  color = "#B71C1C", linewidth = 0.7) +
    geom_line(color = "#37474F", linewidth = 0.35, alpha = 0.9) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.4, linetype = "dashed") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                       limits = c(y_min, max(serie_i) * 1.05)) +
    scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
    labs(
      title = paste0("Figura 3 — Zonas de crisis por intensidad: ", nombre_display),
      subtitle = "Zona amarilla: 10% | Zona naranja: 5% | Zona roja: 1%",
      x = NULL,
      y = "Rendimiento logarítmico diario",
      caption = "Fuente: elaboración propia. Datos Yahoo Finance 2018-2026."
    ) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", size = 11),
      plot.subtitle = element_text(size = 8, color = "grey40"),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )

  ggsave(archivo, plot = fig, width = 14, height = 5, dpi = 300, bg = "white")
  invisible(fig)
}

figura_zonas_crisis("SP500", "S&P 500", "figuras/figura3a_zonas_sp500.png")
figura_zonas_crisis("EURO",  "Eurostoxx 50", "figuras/figura3b_zonas_eurostoxx.png")
figura_zonas_crisis("NIKKEI", "Nikkei 225", "figuras/figura3c_zonas_nikkei.png")

# ------------------------------------------------------------------------------
# FIGURA 4 — COMPARACIÓN DE beta_crisis POR SUBPERIODOS (UMBRAL 5%)
# ------------------------------------------------------------------------------

df_fig4 <- tabla7_tecnica %>%
  mutate(
    Fase = factor(Fase, levels = c("1. Pre-COVID", "2. Pandemia", "3. ETFs")),
    Indice = factor(Indice, levels = c("S&P 500", "Eurostoxx 50", "Nikkei 225")),
    Activo = factor(Activo, levels = c("Bitcoin", "Oro")),
    etiqueta = sprintf("%.2f", Beta_crisis),
    y_texto = ifelse(Beta_crisis >= 0, Beta_crisis + 0.07, Beta_crisis - 0.07)
  )

lim_inf <- min(df_fig4$Beta_crisis, na.rm = TRUE) - 0.15
lim_sup <- max(df_fig4$Beta_crisis, na.rm = TRUE) + 0.15

fig4 <- ggplot(df_fig4, aes(x = Fase, y = Beta_crisis, fill = Activo)) +
  geom_col(position = position_dodge(width = 0.70), width = 0.55, alpha = 0.90) +
  geom_hline(yintercept = 0, linewidth = 0.5, linetype = "dashed", color = "black") +
  geom_text(
    aes(label = etiqueta, y = y_texto, group = Activo),
    position = position_dodge(width = 0.70),
    size = 3,
    fontface = "bold"
  ) +
  facet_wrap(~ Indice, ncol = 3) +
  scale_fill_manual(values = c("Bitcoin" = "#2196F3", "Oro" = "#FF9800"), name = "Activo") +
  scale_y_continuous(limits = c(lim_inf, lim_sup)) +
  labs(
    title = "Figura 4 — Comparación del efecto total en crisis por subperiodos (umbral 5%)",
    subtitle = "Se representa beta_crisis = beta1 + beta2",
    x = NULL,
    y = "beta_crisis",
    caption = "Fuente: elaboración propia a partir de las Tablas 7a y 7b."
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 8, color = "grey40"),
    strip.text = element_text(face = "bold", size = 10),
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 8)
  )

ggsave("figuras/figura4_beta_crisis_subperiodos.png", plot = fig4,
       width = 14, height = 6, dpi = 300, bg = "white")

cat("\n✓ TODAS LAS FIGURAS GENERADAS\n")
cat("Tablas guardadas en: /salidas\n")
cat("Figuras guardadas en: /figuras\n")
