# ==============================================================================
# SCRIPT TFG: BITCOIN VS ORO COMO ACTIVO REFUGIO
# ==============================================================================
# Autor: [Javier Lara]
# Fecha de extracción de datos: Enero 2018 - Febrero 2026
# Descripción: Análisis econométrico para evaluar la condición de activo 
# refugio del Bitcoin y el Oro frente a caídas extremas del mercado.
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. LIBRERÍAS Y DESCARGA DE DATOS
# ------------------------------------------------------------------------------
library(quantmod)
library(PerformanceAnalytics)
library(rugarch)
library(lmtest)
library(sandwich)
library(tseries)
library(zoo)

# Descarga de series temporales (Yahoo Finance)
getSymbols(c("BTC-USD", "GLD", "^GSPC", "^VIX"), from = "2018-01-01", to = "2026-02-24")

# Alineación de fechas (elimina fines de semana del BTC para cuadrar con bolsa)
datos_precios <- na.omit(merge(Ad(`BTC-USD`), Ad(GLD), Ad(GSPC), Cl(VIX)))
colnames(datos_precios) <- c("Precio_BTC", "Precio_ORO", "Precio_SP500", "Nivel_VIX")

# Cálculo de rendimientos logarítmicos
rendimientos_temp <- na.omit(CalculateReturns(datos_precios[, 1:3], method = "log"))
base_datos_final <- na.omit(merge(rendimientos_temp, datos_precios$Nivel_VIX))
colnames(base_datos_final) <- c("BTC", "ORO", "SP500", "VIX")

cat("\n--- Días hábiles totales analizados:", nrow(base_datos_final), "---\n")

# ------------------------------------------------------------------------------
# 2. ESTADÍSTICOS DESCRIPTIVOS Y TESTS DE ESTACIONARIEDAD
# ------------------------------------------------------------------------------
cat("\n=== ESTADÍSTICOS DESCRIPTIVOS ===\n")
# Bitcoin
cat("BTC -> Media:", mean(base_datos_final$BTC), "| Desv.Est:", sd(base_datos_final$BTC), 
    "| Asim:", skewness(base_datos_final$BTC, method="moment"), 
    "| Curtosis:", kurtosis(base_datos_final$BTC, method="moment") + 3, "\n")
# Oro
cat("ORO -> Media:", mean(base_datos_final$ORO), "| Desv.Est:", sd(base_datos_final$ORO), 
    "| Asim:", skewness(base_datos_final$ORO, method="moment"), 
    "| Curtosis:", kurtosis(base_datos_final$ORO, method="moment") + 3, "\n")

cat("\n=== TEST DE DICKEY-FULLER AUMENTADO (Estacionariedad) ===\n")
print(adf.test(base_datos_final$BTC))
print(adf.test(base_datos_final$ORO))

cat("\n=== TEST DE LJUNG-BOX (Autocorrelación / Efectos ARCH) ===\n")
print(Box.test(base_datos_final$BTC^2, type = "Ljung-Box", lag = 10))
print(Box.test(base_datos_final$ORO^2, type = "Ljung-Box", lag = 10))

# ------------------------------------------------------------------------------
# 3. RIESGO: MAXIMUM DRAWDOWN Y RATIO DE SHARPE
# ------------------------------------------------------------------------------
cat("\n=== RIESGO Y DESEMPEÑO ===\n")
cat("Maximum Drawdown (SP500, BTC, ORO):", 
    round(maxDrawdown(base_datos_final$SP500)*100, 2), "%,",
    round(maxDrawdown(base_datos_final$BTC)*100, 2), "%,",
    round(maxDrawdown(base_datos_final$ORO)*100, 2), "%\n")

cat("Sharpe Ratio Anualizado (SP500, BTC, ORO):", 
    round(SharpeRatio.annualized(base_datos_final$SP500, Rf=0), 4), ",",
    round(SharpeRatio.annualized(base_datos_final$BTC, Rf=0), 4), ",",
    round(SharpeRatio.annualized(base_datos_final$ORO, Rf=0), 4), "\n")

# ------------------------------------------------------------------------------
# 4. VOLATILIDAD CONDICIONAL: MODELO GARCH(1,1)
# ------------------------------------------------------------------------------
especificacion_garch <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
                                   mean.model = list(armaOrder = c(0, 0), include.mean = TRUE),
                                   distribution.model = "norm")

ajuste_garch_btc <- ugarchfit(spec = especificacion_garch, data = base_datos_final$BTC)
ajuste_garch_oro <- ugarchfit(spec = especificacion_garch, data = base_datos_final$ORO)

cat("\n=== COEFICIENTES GARCH (1,1) ===\n")
cat("BTC:\n"); print(coef(ajuste_garch_btc))
cat("ORO:\n"); print(coef(ajuste_garch_oro))

# ------------------------------------------------------------------------------
# 5. CONTRASTE OLS Y ERRORES ROBUSTOS HAC (S&P 500 y VIX)
# ------------------------------------------------------------------------------
# 5.1 Umbrales y Dummies (S&P 500 al 5% y VIX al 95%)
q5 <- quantile(base_datos_final$SP500, 0.05)
D5 <- ifelse(base_datos_final$SP500 <= q5, 1, 0)
interaccion_5_pura <- D5 * base_datos_final$SP500

q95_vix <- quantile(base_datos_final$VIX, 0.95)
D_vix <- ifelse(base_datos_final$VIX >= q95_vix, 1, 0)
interaccion_vix_pura <- D_vix * base_datos_final$SP500

# 5.2 Modelos de Regresión
mod_btc_5 <- lm(base_datos_final$BTC ~ base_datos_final$SP500 + interaccion_5_pura)
mod_oro_5 <- lm(base_datos_final$ORO ~ base_datos_final$SP500 + interaccion_5_pura)

mod_btc_vix <- lm(base_datos_final$BTC ~ base_datos_final$SP500 + interaccion_vix_pura)
mod_oro_vix <- lm(base_datos_final$ORO ~ base_datos_final$SP500 + interaccion_vix_pura)

# 5.3 Test de Heterocedasticidad (Breusch-Pagan)
cat("\n=== TEST DE BREUSCH-PAGAN ===\n")
print(bptest(mod_btc_5)); print(bptest(mod_oro_5))

# 5.4 Resultados con Errores Robustos Newey-West (HAC) y Test de Wald
cat("\n=== MODELOS OLS (S&P 500 al 5%) CON HAC ===\n")
print(coeftest(mod_btc_5, vcov = vcovHAC)); print(waldtest(mod_btc_5, vcov = vcovHAC))
print(coeftest(mod_oro_5, vcov = vcovHAC)); print(waldtest(mod_oro_5, vcov = vcovHAC))

cat("\n=== MODELOS OLS (VIX al 95%) CON HAC ===\n")
print(coeftest(mod_btc_vix, vcov = vcovHAC))
print(coeftest(mod_oro_vix, vcov = vcovHAC))

# ------------------------------------------------------------------------------
# 6. ANÁLISIS ESTRUCTURAL POR SUBPERIODOS
# ------------------------------------------------------------------------------
fechas <- index(base_datos_final)
fase1 <- fechas >= as.Date("2018-01-01") & fechas <= as.Date("2020-02-28")
fase2 <- fechas >= as.Date("2020-03-01") & fechas <= as.Date("2022-12-31")
fase3 <- fechas >= as.Date("2023-01-01") & fechas <= as.Date("2026-12-31")

cat("\n=== CORRELACIONES S&P 500 POR SUBPERIODOS ===\n")
cat("BTC Fase 1 (Pre-COVID):", cor(base_datos_final$BTC[fase1], base_datos_final$SP500[fase1]), "\n")
cat("BTC Fase 2 (Pandemia): ", cor(base_datos_final$BTC[fase2], base_datos_final$SP500[fase2]), "\n")
cat("BTC Fase 3 (ETFs):     ", cor(base_datos_final$BTC[fase3], base_datos_final$SP500[fase3]), "\n")

cat("\nORO Fase 1 (Pre-COVID):", cor(base_datos_final$ORO[fase1], base_datos_final$SP500[fase1]), "\n")
cat("ORO Fase 2 (Pandemia): ", cor(base_datos_final$ORO[fase2], base_datos_final$SP500[fase2]), "\n")
cat("ORO Fase 3 (ETFs):     ", cor(base_datos_final$ORO[fase3], base_datos_final$SP500[fase3]), "\n")