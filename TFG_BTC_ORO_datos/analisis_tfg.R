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

# ------------------------------------------------------------------------------
# 7. GENERACIÓN DE GRÁFICOS.
# ------------------------------------------------------------------------------
# Configuración inicial para guardar o visualizar
# (Si deseas guardarlos automáticamente, puedes usar png("nombre.png") antes de cada plot)
# FIGURA 1: Distribución modelos vs Normal (Campanas de Gauss)
par(mfrow = c(1, 2)) # Dividir pantalla en 2 columnas
# Bitcoin
hist(base_datos_final$BTC, breaks = 50, prob = TRUE, col = "lightblue",
 main = "A) Bitcoin: Normal vs Empírica", xlab = "Rendimientos", ylab = "Densidad")
lines(density(base_datos_final$BTC), col = "darkblue", lwd = 2)
curve(dnorm(x, mean = mean(base_datos_final$BTC), sd = sd(base_datos_final$BTC)),
 col = "red", lwd = 2, lty = 2, add = TRUE)
# Oro
hist(base_datos_final$ORO, breaks = 50, prob = TRUE, col = "navajowhite",
 main = "B) Oro: Normal vs Empírica", xlab = "Rendimientos", ylab = "Densidad")
lines(density(base_datos_final$ORO), col = "darkorange", lwd = 2)
curve(dnorm(x, mean = mean(base_datos_final$ORO), sd = sd(base_datos_final$ORO)),
 col = "red", lwd = 2, lty = 2, add = TRUE)
par(mfrow = c(1, 1)) # Restaurar pantalla
# FIGURA 2: Volatilidad Condicional Estimada GARCH (1,1)
vol_btc <- sigma(ajuste_garch_btc)
vol_oro <- sigma(ajuste_garch_oro)
plot(vol_btc, type = "l", col = "blue", ylab = "Volatilidad (Riesgo Diario)",
 xlab = "Fecha", main = "Volatilidad Condicional Estimada (GARCH 1,1)",
 ylim = c(0, max(vol_btc)))
lines(vol_oro, col = "darkorange", lwd = 2)
legend("topright", legend = c("Bitcoin", "Oro"), col = c("blue", "darkorange"), lwd = 2)
# FIGURA 3: Correlación Móvil (Ventana de 60 días)
cor_movil_btc <- rollapply(base_datos_final, width = 60,
 FUN = function(x) cor(x[, "BTC"], x[, "SP500"]),
by.column = FALSE, align = "right")
cor_movil_oro <- rollapply(base_datos_final, width = 60,
 FUN = function(x) cor(x[, "ORO"], x[, "SP500"]),
by.column = FALSE, align = "right")
plot(cor_movil_btc, type = "l", col = "blue", ylim = c(-1, 1),
 main = "Correlación Móvil (60 días): BTC y Oro vs S&P 500",
 ylab = "Nivel de Correlación", xlab = "Fecha")
lines(cor_movil_oro, col = "darkorange", lwd = 2)
abline(h = 0, col = "black", lty = 2)
legend("bottomleft", legend = c("BTC", "Oro"), col = c("blue", "darkorange"), lwd = 2)
# FIGURA 4: Estudio de Robustez (Ventanas de 30, 60 y 120 días)
# (Ejemplo simplificado para BTC)
cor_movil_btc_30 <- rollapply(base_datos_final, width = 30, FUN = function(x) cor(x[, "BTC"],
x[, "SP500"]), by.column = FALSE, align = "right")
cor_movil_btc_120 <- rollapply(base_datos_final, width = 120, FUN = function(x) cor(x[, "BTC"],
x[, "SP500"]), by.column = FALSE, align = "right")
plot(cor_movil_btc_30, type = "l", col = "lightgray", ylim = c(-0.5, 1),
 main = "Robustez Bitcoin: Correlación a 30, 60 y 120 días", ylab = "Correlación BTCS&P500")
lines(cor_movil_btc, col = "red", lwd = 2)
lines(cor_movil_btc_120, col = "darkgreen", lwd = 2)
abline(h = 0, lty = 2)
legend("topleft", legend = c("30 días", "60 días", "120 días"), col = c("lightgray", "red",
"darkgreen"), lwd = 2)
# FIGURA 5: Correlación SIN el 5% de peores caídas
datos_sin_crisis <- base_datos_final[base_datos_final$SP500 > q5, ]
cor_sin_crisis_btc <- rollapply(datos_sin_crisis, width = 60, FUN = function(x) cor(x[, "BTC"],
x[, "SP500"]), by.column = FALSE, align = "right")
cor_sin_crisis_oro <- rollapply(datos_sin_crisis, width = 60, FUN = function(x) cor(x[, "ORO"],
x[, "SP500"]), by.column = FALSE, align = "right")
plot(cor_sin_crisis_btc, type = "l", col = "blue", ylim = c(-1, 1),
 main = "Correlación Móvil (60 días): SIN EL 5% DE PEORES CAÍDAS", ylab = "Correlación")
lines(cor_sin_crisis_oro, col = "darkorange", lwd = 2)
abline(h = 0, lty = 2)
# FIGURA 6: Gráfico de Dispersión frente al Pánico (VIX)
par(mfrow = c(1, 2))
# Bitcoin
plot(as.numeric(base_datos_final$VIX), as.numeric(base_datos_final$BTC), col = "darkblue", pch
= 20,
 xlab = "Nivel del Índice VIX (Miedo)", ylab = "Rendimientos Diarios BTC", main = "A)
Bitcoin frente al Pánico (VIX)")
abline(v = q95_vix, lty = 2, lwd = 2)
abline(lm(base_datos_final$BTC ~ base_datos_final$VIX), col = "red", lwd = 2)
# Oro
plot(as.numeric(base_datos_final$VIX), as.numeric(base_datos_final$ORO), col = "darkorange",
pch = 20,
 xlab = "Nivel del Índice VIX (Miedo)", ylab = "Rendimientos Diarios Oro", main = "B) Oro
frente al Pánico (VIX)")
abline(v = q95_vix, lty = 2, lwd = 2)
abline(lm(base_datos_final$ORO ~ base_datos_final$VIX), col = "red", lwd = 2)
par(mfrow = c(1, 1))
# FIGURA 7: Evolución Estructural (Gráfico de Barras)
cor_fases_btc <- c(0.01, 0.43, 0.32) # Valores extraídos del análisis por subperiodos
cor_fases_oro <- c(-0.13, 0.18, 0.08)
datos_barras <- rbind(cor_fases_btc, cor_fases_oro)
colnames(datos_barras) <- c("Fase 1: Pre-COVID", "Fase 2: Pandemia", "Fase 3: ETFs Wall St")
barplot(datos_barras, beside = TRUE, col = c("royalblue", "darkorange"),
 main = "Evolución Estructural de la Correlación frente al S&P 500",
 ylab = "Coeficiente de Correlación de Pearson", ylim = c(-0.2, 0.5),
 legend.text = c("Bitcoin (BTC)", "Oro Físico (GLD)"),
 args.legend = list(x = "topleft"))
abline(h = 0)
cat("\n=== GRÁFICOS GENERADOS CON ÉXITO ===\n")
