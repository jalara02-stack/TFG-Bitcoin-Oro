# ==============================================================================
# TRABAJO DE FIN DE GRADO
# Bitcoin como Activo Refugio: Comparación con el Oro
# Autor: Javier Lara
# Periodo de análisis: 01/01/2018 - 28/02/2026
# Lenguaje: R
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. LIBRERÍAS
# ------------------------------------------------------------------------------
# Instalar si no están disponibles:
install.packages(c("quantmod","xts","zoo","dplyr","lmtest","sandwich",
                    "tseries","ggplot2","scales","patchwork","tidyr"))

library(quantmod)      # Descarga de datos financieros de Yahoo Finance
library(xts)           # Series temporales
library(zoo)           # Funciones de ventanas móviles
library(dplyr)         # Manipulación de datos
library(lmtest)        # Tests de hipótesis sobre modelos de regresión
library(sandwich)      # Estimador de Newey-West
library(tseries)       # Test ADF y Ljung-Box
library(ggplot2)       # Gráficos
library(scales)        # Formato de ejes en ggplot2
library(patchwork)     # Combinación de gráficos
library(tidyr)         # Transformación de datos

# ------------------------------------------------------------------------------
# 2. DESCARGA Y PREPARACIÓN DE DATOS
# ------------------------------------------------------------------------------
# Fuente: Yahoo Finance
# Periodo: 01/01/2018 - 28/02/2026
# Activos: Bitcoin (BTC-USD), Oro físico ETF (GLD), S&P 500 (^GSPC),
#          Eurostoxx 50 (^STOXX50E), Nikkei 225 (^N225),
#          T-Bill 3 meses (^IRX), VIX (^VIX)

getSymbols(c("BTC-USD", "GLD", "^GSPC", "^STOXX50E", "^N225", "^IRX", "^VIX"),
           from = "2018-01-01",
           to   = "2026-02-28",
           src  = "yahoo",
           auto.assign = TRUE)

# Extraemos precios de cierre ajustados
btc_raw  <- Ad(`BTC-USD`)
gld_raw  <- Ad(GLD)
sp_raw   <- Ad(GSPC)
euro_raw <- Ad(STOXX50E)
nik_raw  <- Ad(N225)
irx_raw  <- Cl(IRX)
vix_raw  <- Cl(VIX)

# Calculamos rendimientos logarítmicos diarios: r_t = ln(P_t / P_{t-1})
btc_ret  <- diff(log(btc_raw),  lag = 1)
gld_ret  <- diff(log(gld_raw),  lag = 1)
sp_ret   <- diff(log(sp_raw),   lag = 1)
euro_ret <- diff(log(euro_raw), lag = 1)
nik_ret  <- diff(log(nik_raw),  lag = 1)

# Fusionamos todas las series en un único objeto xts
# La fusión elimina automáticamente los días en que alguna serie no tiene dato
# (fines de semana, festivos occidentales y japoneses)
datos_raw <- merge(btc_ret, gld_ret, sp_ret, euro_ret, nik_ret,
                   irx_raw, vix_raw,
                   join = "inner")

colnames(datos_raw) <- c("BTC", "ORO", "SP500", "EURO", "NIKKEI",
                         "IRX", "VIX")

# Eliminamos NAs residuales
datos <- na.omit(datos_raw)

cat("Observaciones totales tras homogeneización:", nrow(datos), "\n")
cat("Periodo:", as.character(index(datos)[1]),
    "a", as.character(index(datos)[nrow(datos)]), "\n")

# ------------------------------------------------------------------------------
# 3. TASA LIBRE DE RIESGO
# ------------------------------------------------------------------------------
# Se emplea el T-Bill a 3 meses (^IRX) como tasa libre de riesgo anual.
# Se convierte a rendimiento diario para el cálculo del Ratio de Sharpe.

tasa_rf_anual <- mean(as.numeric(datos[, "IRX"]), na.rm = TRUE) / 100
tasa_rf_diaria <- tasa_rf_anual / 252

cat("Tasa libre de riesgo anual media (T-Bill 3m):",
    round(tasa_rf_anual * 100, 4), "%\n")
cat("Tasa libre de riesgo diaria equivalente:",
    round(tasa_rf_diaria * 100, 6), "%\n")

# ------------------------------------------------------------------------------
# 4. ESTADÍSTICOS DESCRIPTIVOS Y RATIO DE SHARPE — TABLA 1
# ------------------------------------------------------------------------------

activos_desc <- c("BTC", "ORO", "SP500")
nombres_desc <- c("Bitcoin", "Oro (GLD)", "S&P 500")

tabla_descriptivos <- data.frame()

for (i in 1:length(activos_desc)) {
  serie <- as.numeric(datos[, activos_desc[i]])
  n     <- length(serie)
  media <- mean(serie, na.rm = TRUE)
  desv  <- sd(serie, na.rm = TRUE)
  asim  <- mean((serie - media)^3, na.rm = TRUE) / desv^3
  kurt_exc <- mean((serie - media)^4, na.rm = TRUE) / desv^4 - 3
  volat_anual <- desv * sqrt(252) * 100
  sharpe <- (mean(serie, na.rm = TRUE) - tasa_rf_diaria) /
    sd(serie, na.rm = TRUE) * sqrt(252)
  
  tabla_descriptivos <- rbind(tabla_descriptivos, data.frame(
    Activo       = nombres_desc[i],
    Media_diaria = round(media, 6),
    Desv_Est     = round(desv, 6),
    Asimetria    = round(asim, 4),
    Curtosis_exc = round(kurt_exc, 4),
    Volat_anual_pct = round(volat_anual, 2),
    Sharpe_anual = round(sharpe, 4)
  ))
}

cat("\n--- TABLA 1: ESTADÍSTICOS DESCRIPTIVOS Y RATIO DE SHARPE ---\n")
print(tabla_descriptivos, row.names = FALSE)

# ------------------------------------------------------------------------------
# 5. TEST DE DICKEY-FULLER AUMENTADO — TABLA 2
# ------------------------------------------------------------------------------
# H0: la serie tiene raíz unitaria (no estacionaria)
# Rechazo de H0 confirma estacionariedad

cat("\n--- TABLA 2: TEST ADF ---\n")

activos_adf  <- c("BTC", "ORO", "SP500", "EURO", "NIKKEI")
nombres_adf  <- c("Bitcoin", "Oro (GLD)", "S&P 500", "Eurostoxx 50",
                  "Nikkei 225")

tabla_adf <- data.frame()
for (i in 1:length(activos_adf)) {
  serie  <- as.numeric(datos[, activos_adf[i]])
  result <- adf.test(serie, alternative = "stationary")
  tabla_adf <- rbind(tabla_adf, data.frame(
    Activo     = nombres_adf[i],
    Estadistico_ADF = round(result$statistic, 3),
    P_valor    = ifelse(result$p.value < 0.01, "<0.01",
                        round(result$p.value, 4)),
    Resultado  = ifelse(result$p.value < 0.05,
                        "Estacionaria", "No estacionaria")
  ))
}
print(tabla_adf, row.names = FALSE)

# ------------------------------------------------------------------------------
# 6. TESTS DE LJUNG-BOX Y ARCH-LM — TABLA 3
# ------------------------------------------------------------------------------
# Ljung-Box sobre rendimientos al cuadrado: detecta volatility clustering
# ARCH-LM: confirma heterocedasticidad condicional
# Ambos justifican el uso del estimador Newey-West

cat("\n--- TABLA 3: TESTS LJUNG-BOX Y ARCH-LM ---\n")

tabla_arch <- data.frame()
for (i in 1:length(activos_adf)) {
  serie   <- as.numeric(datos[, activos_adf[i]])
  lb_test <- Box.test(serie^2, lag = 10, type = "Ljung-Box")
  arch_lm <- FinTS::ArchTest(serie, lags = 5)
  tabla_arch <- rbind(tabla_arch, data.frame(
    Activo    = nombres_adf[i],
    LB_Q      = round(lb_test$statistic, 2),
    LB_pvalor = ifelse(lb_test$p.value < 0.0001, "<0.0001",
                       round(lb_test$p.value, 4)),
    ARCH_LM   = round(arch_lm$statistic, 2),
    ARCH_pvalor = ifelse(arch_lm$p.value < 0.0001, "<0.0001",
                         round(arch_lm$p.value, 4)),
    Heterocedasticidad = ifelse(arch_lm$p.value < 0.05, "Sí", "No")
  ))
}
print(tabla_arch, row.names = FALSE)

# ------------------------------------------------------------------------------
# 7. CORRELACIÓN MÓVIL — TABLA 4
# ------------------------------------------------------------------------------
# Coeficiente de Pearson en ventanas de 30, 60 y 120 días
# Calculado para cada activo (BTC, ORO) frente a cada índice

# Definición de fases temporales
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

cat("\n--- TABLA 4: CORRELACIONES MEDIAS POR SUBPERIODO ---\n")
tabla_subperiodos <- data.frame()

for (a in 1:length(activos_cor)) {
  for (i in 1:length(indices_cor)) {
    serie_a <- as.numeric(datos[, activos_cor[a]])
    serie_i <- as.numeric(datos[, indices_cor[i]])
    
    tabla_subperiodos <- rbind(tabla_subperiodos, data.frame(
      Activo          = nombres_act[a],
      Indice          = nombres_idx[i],
      Fase1_PreCOVID  = round(cor(serie_a[f1], serie_i[f1]), 3),
      Fase2_Pandemia  = round(cor(serie_a[f2], serie_i[f2]), 3),
      Fase3_ETFs      = round(cor(serie_a[f3], serie_i[f3]), 3),
      Total_periodo   = round(cor(serie_a, serie_i), 3)
    ))
  }
}
print(tabla_subperiodos, row.names = FALSE)

# Correlaciones móviles a 30, 60 y 120 días (para gráficos)
calcular_cor_movil <- function(activo, indice, ventana) {
  rollapply(datos[, c(activo, indice)],
            width      = ventana,
            FUN        = function(x) cor(x[, 1], x[, 2]),
            by.column  = FALSE,
            align      = "right")
}

cor_btc_sp_60  <- calcular_cor_movil("BTC", "SP500",  60)
cor_btc_sp_30  <- calcular_cor_movil("BTC", "SP500",  30)
cor_btc_sp_120 <- calcular_cor_movil("BTC", "SP500", 120)
cor_oro_sp_60  <- calcular_cor_movil("ORO", "SP500",  60)

cor_btc_eu_60  <- calcular_cor_movil("BTC", "EURO",   60)
cor_btc_eu_30  <- calcular_cor_movil("BTC", "EURO",   30)
cor_btc_eu_120 <- calcular_cor_movil("BTC", "EURO",  120)
cor_oro_eu_60  <- calcular_cor_movil("ORO", "EURO",   60)

cor_btc_nk_60  <- calcular_cor_movil("BTC", "NIKKEI", 60)
cor_btc_nk_30  <- calcular_cor_movil("BTC", "NIKKEI", 30)
cor_btc_nk_120 <- calcular_cor_movil("BTC", "NIKKEI",120)
cor_oro_nk_60  <- calcular_cor_movil("ORO", "NIKKEI", 60)

# ------------------------------------------------------------------------------
# 8. MODELO CON DUMMIES — TABLAS 5a, 5b Y 6
# ------------------------------------------------------------------------------
# Modelo: r_{i,t} = β0 + β1·r_{m,t} + β2·(Dt·r_{m,t}) + ε_t
# Dt = 1 si r_{m,t} <= percentil umbral (crisis), 0 en caso contrario
# Errores estándar Newey-West | p-valor β2 unilateral (H1: β2 > 0)

estimar_dummy <- function(activo, indice, umbral_pct,
                          nombre_activo, nombre_indice,
                          datos_uso = datos) {
  
  r_a <- as.numeric(datos_uso[, activo])
  r_i <- as.numeric(datos_uso[, indice])
  n   <- length(r_a)
  
  umbral <- quantile(r_i, probs = umbral_pct, na.rm = TRUE)
  Dt     <- as.numeric(r_i <= umbral)
  inter  <- Dt * r_i
  
  mod    <- lm(r_a ~ r_i + inter)
  lag_nw <- max(floor(4 * (n / 100)^(2/9)), 1)
  vcov_nw <- tryCatch(
    sandwich::NeweyWest(mod, lag = lag_nw, prewhite = FALSE),
    error = function(e) vcov(mod)
  )
  
  b0 <- coef(mod)[1]; b1 <- coef(mod)[2]; b2 <- coef(mod)[3]
  se0 <- sqrt(vcov_nw[1,1])
  se1 <- sqrt(vcov_nw[2,2])
  se2 <- sqrt(vcov_nw[3,3])
  t2  <- b2 / se2
  
  pval_b0 <- 2 * pt(abs(b0/se0), df = n-3, lower.tail = FALSE)
  pval_b1 <- 2 * pt(abs(b1/se1), df = n-3, lower.tail = FALSE)
  pval_b2_bil <- 2 * pt(abs(t2), df = n-3, lower.tail = FALSE)
  pval_b2_uni <- ifelse(b2 > 0, pval_b2_bil / 2, 1 - pval_b2_bil / 2)
  
  clasificacion <- dplyr::case_when(
    b2 < 0 & pval_b2_uni < 0.05  ~ "Refugio fuerte",
    b2 < 0 & pval_b2_uni >= 0.05 ~ "Refugio débil",
    b2 > 0 & pval_b2_uni < 0.05  ~ "Procíclico",
    b2 > 0 & pval_b2_uni >= 0.05 ~ "Sin evidencia",
    TRUE ~ "Refugio débil"
  )
  
  sig <- dplyr::case_when(
    pval_b2_uni < 0.01 ~ "***",
    pval_b2_uni < 0.05 ~ "**",
    pval_b2_uni < 0.10 ~ "*",
    TRUE               ~ "n.s."
  )
  
  data.frame(
    Activo        = nombre_activo,
    Indice        = nombre_indice,
    Umbral        = paste0(umbral_pct * 100, "%"),
    N_crisis      = sum(Dt),
    Beta0         = round(b0, 6),
    Beta1         = round(b1, 4),
    p_beta1       = round(pval_b0, 4),
    Beta2         = round(b2, 4),
    p_beta2_uni   = round(pval_b2_uni, 4),
    Sig           = sig,
    Clasificacion = clasificacion,
    stringsAsFactors = FALSE
  )
}

# Estimación de 18 modelos: 2 activos × 3 índices × 3 umbrales
activos  <- c("BTC", "ORO")
indices  <- c("SP500", "EURO", "NIKKEI")
umbrales <- c(0.01, 0.05, 0.10)
nomb_act <- c("Bitcoin", "Oro")
nomb_idx <- c("S&P 500", "Eurostoxx 50", "Nikkei 225")

resultados_principal <- list()
for (a in 1:2) {
  for (i in 1:3) {
    for (u in 1:3) {
      clave <- paste0(activos[a], "_", indices[i], "_p", umbrales[u]*100)
      resultados_principal[[clave]] <- estimar_dummy(
        activos[a], indices[i], umbrales[u], nomb_act[a], nomb_idx[i]
      )
    }
  }
}

tabla_principal <- do.call(rbind, resultados_principal)

cat("\n--- TABLAS 5a y 5b: MODELO CON DUMMIES PERÍODO COMPLETO ---\n")
cat("Bitcoin:\n")
print(tabla_principal[tabla_principal$Activo == "Bitcoin",
                      c("Indice","Umbral","N_crisis","Beta0","Beta1",
                        "Beta2","p_beta2_uni","Sig","Clasificacion")],
      row.names = FALSE)
cat("\nOro:\n")
print(tabla_principal[tabla_principal$Activo == "Oro",
                      c("Indice","Umbral","N_crisis","Beta0","Beta1",
                        "Beta2","p_beta2_uni","Sig","Clasificacion")],
      row.names = FALSE)

# Tabla 6 — síntesis por intensidad
cat("\n--- TABLA 6: SÍNTESIS β2 POR INTENSIDAD DE CRISIS ---\n")
tabla6 <- tabla_principal %>%
  select(Activo, Indice, Umbral, Beta2, Sig) %>%
  tidyr::pivot_wider(names_from  = Umbral,
                     values_from = c(Beta2, Sig))
print(tabla6, row.names = FALSE)

# ------------------------------------------------------------------------------
# 9. ANÁLISIS DE ROBUSTEZ CON VIX — TABLA 7
# ------------------------------------------------------------------------------
# Dt_VIX = 1 cuando VIX >= percentil 95 de su distribución histórica
# Solo se estima frente al S&P 500

cat("\n--- TABLA 7: ROBUSTEZ CON DUMMY VIX (percentil 95) ---\n")

vix_serie  <- as.numeric(datos[, "VIX"])
vix_p95    <- quantile(vix_serie, probs = 0.95, na.rm = TRUE)
Dt_VIX     <- as.numeric(vix_serie >= vix_p95)

cat("Umbral VIX p95:", round(vix_p95, 2), "\n")
cat("Días de pánico extremo:", sum(Dt_VIX), "\n")

estimar_vix <- function(activo, nombre_activo) {
  r_a   <- as.numeric(datos[, activo])
  r_sp  <- as.numeric(datos[, "SP500"])
  n     <- length(r_a)
  inter <- Dt_VIX * r_sp
  
  mod    <- lm(r_a ~ r_sp + inter)
  lag_nw <- floor(4 * (n / 100)^(2/9))
  vcov_nw <- sandwich::NeweyWest(mod, lag = lag_nw, prewhite = FALSE)
  
  b1  <- coef(mod)[2]; b2 <- coef(mod)[3]
  se1 <- sqrt(vcov_nw[2,2]); se2 <- sqrt(vcov_nw[3,3])
  t2  <- b2 / se2
  
  pval_b1     <- 2 * pt(abs(b1/se1), df = n-3, lower.tail = FALSE)
  pval_b2_bil <- 2 * pt(abs(t2), df = n-3, lower.tail = FALSE)
  pval_b2_uni <- ifelse(b2 > 0, pval_b2_bil / 2, 1 - pval_b2_bil / 2)
  
  sig <- dplyr::case_when(
    pval_b2_uni < 0.01 ~ "***",
    pval_b2_uni < 0.05 ~ "**",
    pval_b2_uni < 0.10 ~ "*",
    TRUE               ~ "n.s."
  )
  
  clasificacion <- dplyr::case_when(
    b2 < 0 & pval_b2_uni < 0.05  ~ "Refugio fuerte",
    b2 < 0 & pval_b2_uni >= 0.05 ~ "Refugio débil",
    b2 > 0 & pval_b2_uni < 0.05  ~ "Procíclico",
    TRUE                          ~ "Sin evidencia"
  )
  
  cat("\n---", nombre_activo, "---\n")
  cat(sprintf("  β1 = %.4f (p = %.4f)\n", b1, pval_b1))
  cat(sprintf("  β2 = %.4f (p_uni = %.4f) %s — %s\n",
              b2, pval_b2_uni, sig, clasificacion))
}

estimar_vix("BTC", "Bitcoin")
estimar_vix("ORO", "Oro")

# ------------------------------------------------------------------------------
# 10. ANÁLISIS POR SUBPERIODOS — TABLAS 8a Y 8b
# ------------------------------------------------------------------------------
# Replica el modelo con dummies para cada fase temporal
# Los percentiles se calculan DENTRO de cada subperiodo

fases <- list(
  list(nombre = "1. Pre-COVID",  inicio = fase1_inicio, fin = fase1_fin),
  list(nombre = "2. Pandemia",   inicio = fase2_inicio, fin = fase2_fin),
  list(nombre = "3. ETFs",       inicio = fase3_inicio, fin = fase3_fin)
)

resultados_sub <- list()

for (fase in fases) {
  rango     <- paste0(as.character(fase$inicio), "/",
                      as.character(fase$fin))
  datos_sub <- na.omit(datos[rango])
  
  for (a in 1:2) {
    for (i in 1:3) {
      res <- estimar_dummy(
        activo        = activos[a],
        indice        = indices[i],
        umbral_pct    = 0.05,
        nombre_activo = nomb_act[a],
        nombre_indice = nomb_idx[i],
        datos_uso     = datos_sub
      )
      res$Fase <- fase$nombre
      clave <- paste0(fase$nombre, "_", activos[a], "_", indices[i])
      resultados_sub[[clave]] <- res
    }
  }
  cat("✓ Completada:", fase$nombre, "\n")
}

tabla_sub <- do.call(rbind, resultados_sub)

cat("\n--- TABLA 8a: BITCOIN POR FASES (umbral 5%) ---\n")
tabla8a <- tabla_sub %>%
  filter(Activo == "Bitcoin") %>%
  select(Fase, Indice, Beta1, Beta2, p_beta2_uni, Sig, Clasificacion)
print(tabla8a, row.names = FALSE)

cat("\n--- TABLA 8b: ORO POR FASES (umbral 5%) ---\n")
tabla8b <- tabla_sub %>%
  filter(Activo == "Oro") %>%
  select(Fase, Indice, Beta1, Beta2, p_beta2_uni, Sig, Clasificacion)
print(tabla8b, row.names = FALSE)

cat("\n✓ ANÁLISIS FINALIZADO\n")





































# ==============================================================================
# GENERACIÓN DE FIGURAS
# Trabajo de Fin de Grado — Bitcoin como Activo Refugio
# Autor: Javier Lara
# Descripción: Genera y guarda todas las figuras del trabajo en alta resolución.
#              Requiere haber ejecutado previamente el Bloque 1 de análisis.
# Directorio de salida: carpeta de trabajo activa (getwd())
# ==============================================================================

# Verificar que los objetos necesarios están disponibles
stopifnot(exists("datos"), exists("tabla_sub"), exists("tabla_principal"))

colores_fase <- c("Fase 1" = "#E3F2FD",
                  "Fase 2" = "#FFF3E0",
                  "Fase 3" = "#E8F5E9")

# ------------------------------------------------------------------------------
# FIGURA 1 — Distribución de rendimientos: Bitcoin y Oro
# ------------------------------------------------------------------------------

btc_vec <- as.numeric(datos[, "BTC"])
oro_vec <- as.numeric(datos[, "ORO"])

df_dist <- data.frame(
  rendimiento = c(btc_vec, oro_vec),
  activo      = rep(c("Bitcoin", "Oro"),
                    c(length(btc_vec), length(oro_vec)))
)

params_dist <- df_dist %>%
  group_by(activo) %>%
  summarise(media = mean(rendimiento), desv = sd(rendimiento),
            .groups = "drop")

fig1 <- ggplot(df_dist, aes(x = rendimiento)) +
  geom_histogram(aes(y = after_stat(density), fill = activo),
                 bins = 100, alpha = 0.5, color = NA) +
  geom_density(aes(color = activo), linewidth = 0.8) +
  stat_function(
    data = df_dist %>% filter(activo == "Bitcoin"),
    fun  = dnorm,
    args = list(mean = params_dist$media[params_dist$activo == "Bitcoin"],
                sd   = params_dist$desv[params_dist$activo  == "Bitcoin"]),
    color = "#B71C1C", linewidth = 0.8, linetype = "dashed"
  ) +
  stat_function(
    data = df_dist %>% filter(activo == "Oro"),
    fun  = dnorm,
    args = list(mean = params_dist$media[params_dist$activo == "Oro"],
                sd   = params_dist$desv[params_dist$activo  == "Oro"]),
    color = "#E65100", linewidth = 0.8, linetype = "dashed"
  ) +
  facet_wrap(~ activo, ncol = 2, scales = "free") +
  scale_fill_manual(values  = c("Bitcoin" = "#2196F3", "Oro" = "#FF9800"),
                    guide   = "none") +
  scale_color_manual(values = c("Bitcoin" = "#1565C0", "Oro" = "#E65100"),
                     guide  = "none") +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title    = "Figura 1 — Distribución de rendimientos diarios: Bitcoin y Oro",
    subtitle = paste0("Histograma de densidad empírica (sólido) vs ",
                      "distribución normal teórica (discontinuo)\n",
                      "Las colas pesadas son visibles en ambos activos, ",
                      "siendo más pronunciadas en el Bitcoin"),
    x        = "Rendimiento logarítmico diario",
    y        = "Densidad",
    caption  = "Fuente: elaboración propia. Datos Yahoo Finance 2018-2026."
  ) +
  theme_minimal(base_size = 10) +
  theme(plot.title    = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 8, color = "grey40"),
        strip.text    = element_text(face = "bold", size = 10),
        panel.grid.minor = element_blank())

ggsave("figura1_distribuciones.png", plot = fig1,
       width = 12, height = 5, dpi = 300, bg = "white")
cat("✓ Figura 1 guardada\n")


# ------------------------------------------------------------------------------
# FUNCIÓN AUXILIAR — Figura de correlación móvil
# Usada para las Figuras 2a, 2b y 2c
# ------------------------------------------------------------------------------

figura_correlacion <- function(cor_btc_30, cor_btc_60, cor_btc_120,
                               cor_oro_60, titulo, archivo) {
  fechas_f <- index(cor_btc_60)
  
  df_cor <- data.frame(
    fecha    = fechas_f,
    BTC_30   = as.numeric(cor_btc_30)[
      (length(as.numeric(cor_btc_30)) -
         length(fechas_f) + 1):length(as.numeric(cor_btc_30))],
    BTC_60   = as.numeric(cor_btc_60),
    BTC_120  = as.numeric(cor_btc_120)[
      (length(as.numeric(cor_btc_120)) -
         length(fechas_f) + 1):length(as.numeric(cor_btc_120))],
    ORO_60   = as.numeric(cor_oro_60)
  ) %>% na.omit()
  
  fig <- ggplot(df_cor, aes(x = fecha)) +
    annotate("rect", xmin = fase1_inicio, xmax = fase1_fin,
             ymin = -Inf, ymax = Inf, fill = "#E3F2FD", alpha = 0.4) +
    annotate("rect", xmin = fase2_inicio, xmax = fase2_fin,
             ymin = -Inf, ymax = Inf, fill = "#FFF3E0", alpha = 0.4) +
    annotate("rect", xmin = fase3_inicio, xmax = fase3_fin,
             ymin = -Inf, ymax = Inf, fill = "#E8F5E9", alpha = 0.4) +
    annotate("text", x = as.Date("2019-01-01"), y = 0.90,
             label = "Fase 1", size = 2.8, color = "#1565C0",
             fontface = "bold") +
    annotate("text", x = as.Date("2021-06-01"), y = 0.90,
             label = "Fase 2", size = 2.8, color = "#E65100",
             fontface = "bold") +
    annotate("text", x = as.Date("2024-07-01"), y = 0.90,
             label = "Fase 3", size = 2.8, color = "#2E7D32",
             fontface = "bold") +
    geom_line(aes(y = BTC_30),  color = "#90CAF9",
              linewidth = 0.4, alpha = 0.7) +
    geom_line(aes(y = BTC_60),  color = "#1565C0", linewidth = 0.8) +
    geom_line(aes(y = BTC_120), color = "#0D47A1",
              linewidth = 0.4, linetype = "dashed") +
    geom_line(aes(y = ORO_60),  color = "#FF9800", linewidth = 0.8) +
    geom_hline(yintercept = 0, linewidth = 0.5,
               linetype = "dashed", color = "black") +
    scale_y_continuous(limits = c(-1, 1),
                       breaks = seq(-1, 1, 0.25)) +
    scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
    labs(
      title    = titulo,
      subtitle = paste0(
        "Bitcoin: azul claro = 30d | azul = 60d | azul oscuro = 120d\n",
        "Oro: naranja = 60d | Áreas: Fase 1 (azul), ",
        "Fase 2 (naranja), Fase 3 (verde)"),
      x        = NULL,
      y        = "Coeficiente de correlación de Pearson",
      caption  = "Fuente: elaboración propia. Datos Yahoo Finance 2018-2026."
    ) +
    theme_minimal(base_size = 10) +
    theme(plot.title       = element_text(face = "bold", size = 11),
          plot.subtitle    = element_text(size = 7.5, color = "grey40"),
          panel.grid.minor = element_blank(),
          axis.text.x      = element_text(angle = 45, hjust = 1))
  
  ggsave(archivo, plot = fig,
         width = 14, height = 5, dpi = 300, bg = "white")
  cat("✓ Guardada:", archivo, "\n")
  invisible(fig)
}

# Figura 2a — S&P 500
figura_correlacion(
  cor_btc_30  = cor_btc_sp_30,
  cor_btc_60  = cor_btc_sp_60,
  cor_btc_120 = cor_btc_sp_120,
  cor_oro_60  = cor_oro_sp_60,
  titulo  = "Figura 2a — Correlación móvil con el S&P 500",
  archivo = "figura2a_correlacion_sp500.png"
)

# Figura 2b — Eurostoxx 50
figura_correlacion(
  cor_btc_30  = cor_btc_eu_30,
  cor_btc_60  = cor_btc_eu_60,
  cor_btc_120 = cor_btc_eu_120,
  cor_oro_60  = cor_oro_eu_60,
  titulo  = "Figura 2b — Correlación móvil con el Eurostoxx 50",
  archivo = "figura2b_correlacion_eurostoxx.png"
)

# Figura 2c — Nikkei 225
figura_correlacion(
  cor_btc_30  = cor_btc_nk_30,
  cor_btc_60  = cor_btc_nk_60,
  cor_btc_120 = cor_btc_nk_120,
  cor_oro_60  = cor_oro_nk_60,
  titulo  = "Figura 2c — Correlación móvil con el Nikkei 225",
  archivo = "figura2c_correlacion_nikkei.png"
)

# ------------------------------------------------------------------------------
# FUNCIÓN AUXILIAR — Zonas de crisis por intensidad
# Usada para las Figuras 3a, 3b y 3c
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
    annotate("rect",
             xmin = min(fechas_i), xmax = max(fechas_i),
             ymin = u10, ymax = u5,
             fill = "#FFF9C4", alpha = 0.8) +
    annotate("rect",
             xmin = min(fechas_i), xmax = max(fechas_i),
             ymin = u5, ymax = u1,
             fill = "#FFE0B2", alpha = 0.8) +
    annotate("rect",
             xmin = min(fechas_i), xmax = max(fechas_i),
             ymin = y_min, ymax = u1,
             fill = "#FFCDD2", alpha = 0.8) +
    geom_hline(yintercept = u10, color = "#F9A825",
               linewidth = 0.7) +
    geom_hline(yintercept = u5,  color = "#E65100",
               linewidth = 0.7) +
    geom_hline(yintercept = u1,  color = "#B71C1C",
               linewidth = 0.7) +
    annotate("label", x = max(fechas_i), y = (u10 + u5) / 2,
             label = "10%\nmoderado", color = "#F57F17",
             fill = "white", size = 2.5, hjust = 1, label.size = 0) +
    annotate("label", x = max(fechas_i), y = (u5 + u1) / 2,
             label = "5%\ngrave", color = "#E65100",
             fill = "white", size = 2.5, hjust = 1, label.size = 0) +
    annotate("label", x = max(fechas_i), y = u1 * 1.3,
             label = "1%\nextremo", color = "#B71C1C",
             fill = "white", size = 2.5, hjust = 1, label.size = 0) +
    annotate("rect", xmin = fase1_inicio, xmax = fase1_fin,
             ymin = 0, ymax = 0.005, fill = "#1565C0", alpha = 0.5) +
    annotate("rect", xmin = fase2_inicio, xmax = fase2_fin,
             ymin = 0, ymax = 0.005, fill = "#E65100", alpha = 0.5) +
    annotate("rect", xmin = fase3_inicio, xmax = fase3_fin,
             ymin = 0, ymax = 0.005, fill = "#2E7D32", alpha = 0.5) +
    annotate("text", x = as.Date("2019-01-01"), y = 0.007,
             label = "Fase 1", size = 2.2, color = "#1565C0",
             fontface = "bold") +
    annotate("text", x = as.Date("2021-06-01"), y = 0.007,
             label = "Fase 2", size = 2.2, color = "#E65100",
             fontface = "bold") +
    annotate("text", x = as.Date("2024-07-01"), y = 0.007,
             label = "Fase 3", size = 2.2, color = "#2E7D32",
             fontface = "bold") +
    geom_line(color = "#37474F", linewidth = 0.35, alpha = 0.9) +
    geom_hline(yintercept = 0, color = "black",
               linewidth = 0.4, linetype = "dashed") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                       limits = c(y_min, max(serie_i) * 1.05)) +
    scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
    labs(
      title   = paste0("Zonas de crisis por intensidad — ", nombre_display),
      subtitle = paste0("Zona amarilla: 10% moderado | ",
                        "Zona naranja: 5% grave | ",
                        "Zona roja: 1% extremo"),
      x        = NULL,
      y        = "Rendimiento logarítmico diario",
      caption  = paste0("Fuente: elaboración propia. ",
                        "Datos Yahoo Finance 2018-2026.")
    ) +
    theme_minimal(base_size = 10) +
    theme(plot.title       = element_text(face = "bold", size = 11),
          plot.subtitle    = element_text(size = 7.5, color = "grey40"),
          panel.grid.minor = element_blank(),
          axis.text.x      = element_text(angle = 45, hjust = 1))
  
  ggsave(archivo, plot = fig,
         width = 14, height = 5, dpi = 300, bg = "white")
  cat("✓ Guardada:", archivo, "\n")
  invisible(fig)
}

figura_zonas_crisis("SP500", "S&P 500",      "figura3a_zonas_sp500.png")
figura_zonas_crisis("EURO",  "Eurostoxx 50", "figura3b_zonas_eurostoxx.png")
figura_zonas_crisis("NIKKEI","Nikkei 225",   "figura3c_zonas_nikkei.png")

# ------------------------------------------------------------------------------
# FIGURA 4 — Comparación β₂ Bitcoin vs Oro por subperiodos (umbral 5%)
# ------------------------------------------------------------------------------

df_fig4 <- tabla_sub %>%
  filter(Umbral == "5%") %>%
  mutate(
    Fase_label = dplyr::case_when(
      grepl("Pre-COVID", Fase) ~ "F1\nPre-COVID",
      grepl("Pandemia",  Fase) ~ "F2\nPandemia",
      grepl("ETFs",      Fase) ~ "F3\nETFs"
    ),
    Fase_label = factor(Fase_label,
                        levels = c("F1\nPre-COVID", "F2\nPandemia", "F3\nETFs"))
  )

fig4 <- ggplot(df_fig4,
               aes(x = Fase_label, y = Beta2,
                   fill = Activo, group = Activo)) +
  geom_col(position = position_dodge(width = 0.65),
           width = 0.55, alpha = 0.85) +
  geom_hline(yintercept = 0, linewidth = 0.5,
             linetype = "dashed", color = "black") +
  geom_text(
    aes(label = sprintf("%.2f", Beta2),
        vjust = ifelse(Beta2 >= 0, -0.3, 1.3)),
    position = position_dodge(width = 0.65),
    size = 2.8, fontface = "bold"
  ) +
  scale_fill_manual(
    values = c("Bitcoin" = "#2196F3", "Oro" = "#FF9800"),
    name   = "Activo"
  ) +
  facet_wrap(~ Indice, ncol = 3) +
  scale_y_continuous(breaks = seq(-1, 2.5, 0.5),
                     limits = c(-1.0, 1.5)) +
  labs(
    title    = expression(
      paste("Figura 4 — Comparación ",
            beta[2], " Bitcoin vs Oro por subperiodos (Umbral 5%)")),
    subtitle = paste0(
      "Bitcoin (azul) vs Oro (naranja) | ",
      "B2 positivo = comportamiento prociclico en crisis | ",
      "B2 negativo = comportamiento refugio en crisis"),
    x        = NULL,
    y        = expression(paste(beta[2],
                                " — coeficiente de interaccion en crisis")),
    caption  = paste0(
      "Fuente: elaboracion propia. Datos Yahoo Finance 2018-2026.\n",
      "Significatividad: ***p<0.01 **p<0.05 *p<0.10 n.s.")
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title       = element_text(face = "bold", size = 11),
    plot.subtitle    = element_text(size = 8, color = "grey40"),
    strip.text       = element_text(face = "bold", size = 10),
    legend.position  = "bottom",
    panel.grid.minor = element_blank(),
    axis.text.x      = element_text(size = 8)
  )

ggsave("figura4_comparacion_beta2.png", plot = fig4,
       width = 14, height = 6, dpi = 300, bg = "white")
cat("✓ Figura 4 guardada\n")

cat("\n✓ TODAS LAS FIGURAS GENERADAS Y GUARDADAS\n")
cat("Directorio:", getwd(), "\n")

