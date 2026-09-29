# 0. PAQUETES
library(lavaan)
library(MASS)
library(dplyr)
library(ggplot2)

set.seed(12345)


# NOTA: la simulación no utiliza los datos reales. Solo genera datos
# simulados con los mismos nombres de variables del modelo.
#
# Ejecutar con la raíz del repositorio como directorio de trabajo.
# Los resultados se guardan en Analisis/MonteCarlo/.

dir_montecarlo <- file.path("Analisis", "MonteCarlo")
dir.create(dir_montecarlo, showWarnings = FALSE, recursive = TRUE)


# 3. MODELO DE MEDIACIÓN----

modelo <- '

  # Modelo del mediador
  PSS_DIRt ~ a1*MAIA_Autorregulacion_DIRd + a2*Slope__eeg_global

  # Modelo del outcome
  GHQ12_DIRt ~ b*PSS_DIRt + c1*MAIA_Autorregulacion_DIRd + c2*Slope__eeg_global

  # Efectos indirectos
  ind_MAIA := a1*b
  ind_Slope := a2*b

  # Efectos totales
  total_MAIA := c1 + (a1*b)
  total_Slope := c2 + (a2*b)

'


# 4. DEFINIR LOS TRES ESCENARIOS----

#
# Todos los coeficientes se interpretan como
# coeficientes estandarizados poblacionales hipotéticos.
#
# IMPORTANTE:
# Estos valores NO provienen de tus resultados.
# Son escenarios de sensibilidad a priori.
#
# a1 = MAIA -> PSS
# a2 = Slope -> PSS
# b  = PSS -> GHQ12
# c1 = MAIA -> GHQ12 directo
# c2 = Slope -> GHQ12 directo
#
# r_X1X2 = correlación poblacional hipotética
#          entre MAIA y Slope.

escenarios <- data.frame(
  
  escenario = c(
    "Pequeno",
    "Moderado",
    "Grande"
  ),
  
  # MAIA -> PSS
  a1 = c(
    -0.15,
    -0.25,
    -0.35
  ),
  
  # Slope -> PSS
  a2 = c(
    -0.15,
    -0.25,
    -0.35
  ),
  
  # PSS -> GHQ12
  b = c(
    0.40,
    0.50, 
    0.55
  ),
  
  # MAIA -> GHQ12
  # efecto directo hipotético
  c1 = c(
    -0.10,
    -0.20,
    -0.30
  ),
  
  # Slope -> GHQ12
  # efecto directo hipotético
  c2 = c(
    -0.10,
    -0.20,
    -0.30
  ),
  
  # Correlación MAIA-Slope
  r_X1X2 = c(
    0.10,
    0.15,
    0.20
  )
)


# 5. CALCULAR LOS EFECTOS INDIRECTOS POBLACIONALES ----

escenarios <- escenarios %>%
  
  mutate(
    
    ind_MAIA_poblacional =
      a1 * b,
    
    ind_Slope_poblacional =
      a2 * b
    
  )


cat("\nEscenarios simulados:\n")

print(escenarios)


# 6. FUNCIÓN PARA GENERAR UNA BASE SIMULADA----
#
# Las cuatro variables tendrán aproximadamente:
#
# media = 0
# varianza = 1
#
# Por tanto, los coeficientes introducidos pueden
# interpretarse como coeficientes estandarizados.
#
# ============================================================

simular_datos <- function(
    N,
    a1,
    a2,
    b,
    c1,
    c2,
    r_X1X2
) {
  
  # ==========================================================
  # 6.1. Generar los dos predictores
  # ==========================================================
  
  Sigma_X <- matrix(
    c(
      1,       r_X1X2,
      r_X1X2,  1
    ),
    nrow = 2,
    byrow = TRUE
  )
  
  
  X <- MASS::mvrnorm(
    n = N,
    mu = c(0, 0),
    Sigma = Sigma_X
  )
  
  
  MAIA_Autorregulacion_DIRd <- X[, 1]
  
  Slope__eeg_global <- X[, 2]
  
  
  # ==========================================================
  # 6.2. Generar PSS
  # ==========================================================
  #
  # PSS = a1*MAIA + a2*Slope + error
  #
  # Elegimos la varianza del error para que:
  #
  # Var(PSS) ≈ 1
  #
  # ==========================================================
  
  R2_PSS <-
    a1^2 +
    a2^2 +
    2 * a1 * a2 * r_X1X2
  
  
  var_error_PSS <- 1 - R2_PSS
  
  
  if(var_error_PSS <= 0) {
    
    stop(
      "La varianza residual de PSS es <= 0.
       Revisa los parámetros especificados."
    )
    
  }
  
  
  PSS_DIRt <-
    
    a1 * MAIA_Autorregulacion_DIRd +
    
    a2 * Slope__eeg_global +
    
    rnorm(
      N,
      mean = 0,
      sd = sqrt(var_error_PSS)
    )
  
  
  # ==========================================================
  # 6.3. Covarianzas necesarias para generar GHQ12
  # ==========================================================
  
  cov_MAIA_PSS <-
    a1 +
    a2 * r_X1X2
  
  
  cov_Slope_PSS <-
    a2 +
    a1 * r_X1X2
  
  
  # ==========================================================
  # 6.4. Varianza explicada de GHQ12
  # ==========================================================
  #
  # GHQ12 =
  # b*PSS +
  # c1*MAIA +
  # c2*Slope +
  # error
  #
  # ==========================================================
  
  var_pred_GHQ <-
    
    b^2 +
    
    c1^2 +
    
    c2^2 +
    
    2 * b * c1 * cov_MAIA_PSS +
    
    2 * b * c2 * cov_Slope_PSS +
    
    2 * c1 * c2 * r_X1X2
  
  
  var_error_GHQ <- 1 - var_pred_GHQ
  
  
  if(var_error_GHQ <= 0) {
    
    stop(
      "La varianza residual de GHQ12 es <= 0.
       Revisa los parámetros especificados."
    )
    
  }
  
  
  # ==========================================================
  # 6.5. Generar GHQ12
  # ==========================================================
  
  GHQ12_DIRt <-
    
    b * PSS_DIRt +
    
    c1 * MAIA_Autorregulacion_DIRd +
    
    c2 * Slope__eeg_global +
    
    rnorm(
      N,
      mean = 0,
      sd = sqrt(var_error_GHQ)
    )
  
  
  # ==========================================================
  # 6.6. Crear dataframe
  # ==========================================================
  
  datos_simulados <- data.frame(
    
    MAIA_Autorregulacion_DIRd =
      MAIA_Autorregulacion_DIRd,
    
    Slope__eeg_global =
      Slope__eeg_global,
    
    PSS_DIRt =
      PSS_DIRt,
    
    GHQ12_DIRt =
      GHQ12_DIRt
    
  )
  
  
  return(datos_simulados)
}


# ============================================================
# 7. FUNCIÓN PARA CALCULAR POTENCIA EN UN N
# ============================================================
#
# La potencia es la proporción de simulaciones en las
# que p < .05.
#
# Calculamos:
#
# 1. Potencia indirecto MAIA
# 2. Potencia indirecto Slope
# 3. Potencia para ambos indirectos simultáneamente
#
# ============================================================

calcular_potencia <- function(
    N,
    parametros,
    nsim = 5000,
    alpha = 0.05
) {
  
  sig_MAIA <- rep(NA, nsim)
  
  sig_Slope <- rep(NA, nsim)
  
  sig_ambos <- rep(NA, nsim)
  
  convergio <- rep(FALSE, nsim)
  
  
  for(i in seq_len(nsim)) {
    
    
    # ========================================================
    # 7.1. Simular datos
    # ========================================================
    
    datos_sim <- simular_datos(
      
      N = N,
      
      a1 = parametros$a1,
      
      a2 = parametros$a2,
      
      b = parametros$b,
      
      c1 = parametros$c1,
      
      c2 = parametros$c2,
      
      r_X1X2 = parametros$r_X1X2
    )
    
    
    # ========================================================
    # 7.2. Ajustar modelo
    # ========================================================
    
    fit <- try(
      
      lavaan::sem(
        
        model = modelo,
        
        data = datos_sim,
        
        estimator = "ML",
        
        fixed.x = FALSE
        
      ),
      
      silent = TRUE
      
    )
    
    
    # ========================================================
    # 7.3. Comprobar convergencia
    # ========================================================
    
    if(!inherits(fit, "try-error")) {
      
      
      conv <- try(
        lavInspect(fit, "converged"),
        silent = TRUE
      )
      
      
      if(
        !inherits(conv, "try-error") &&
        isTRUE(conv)
      ) {
        
        convergio[i] <- TRUE
        
        
        # ====================================================
        # 7.4. Extraer parámetros
        # ====================================================
        
        parametros_estimados <-
          parameterEstimates(
            fit,
            standardized = TRUE
          )
        
        
        # ====================================================
        # 7.5. p-value indirecto MAIA
        # ====================================================
        
        p_MAIA <-
          
          parametros_estimados$pvalue[
            
            parametros_estimados$lhs ==
              "ind_MAIA" &
              
              parametros_estimados$op ==
              ":="
            
          ]
        
        
        # ====================================================
        # 7.6. p-value indirecto Slope
        # ====================================================
        
        p_Slope <-
          
          parametros_estimados$pvalue[
            
            parametros_estimados$lhs ==
              "ind_Slope" &
              
              parametros_estimados$op ==
              ":="
            
          ]
        
        
        # ====================================================
        # 7.7. Guardar significación
        # ====================================================
        
        if(
          length(p_MAIA) == 1 &&
          length(p_Slope) == 1
        ) {
          
          sig_MAIA[i] <-
            p_MAIA < alpha
          
          
          sig_Slope[i] <-
            p_Slope < alpha
          
          
          sig_ambos[i] <-
            (p_MAIA < alpha) &&
            (p_Slope < alpha)
          
        }
        
      }
      
    }
    
  }
  
  
  # ==========================================================
  # 7.8. Resultados
  # ==========================================================
  
  validos <- convergio &
    !is.na(sig_MAIA) &
    !is.na(sig_Slope)
  
  
  if(sum(validos) == 0) {
    
    return(
      
      data.frame(
        
        N = N,
        
        power_ind_MAIA = NA,
        
        power_ind_Slope = NA,
        
        power_ambos = NA,
        
        tasa_convergencia = 0
        
      )
      
    )
    
  }
  
  
  data.frame(
    
    N = N,
    
    power_ind_MAIA =
      mean(sig_MAIA[validos]),
    
    power_ind_Slope =
      mean(sig_Slope[validos]),
    
    power_ambos =
      mean(sig_ambos[validos]),
    
    tasa_convergencia =
      mean(convergio)
    
  )
  
}


# ============================================================
# 8. TAMAÑOS MUESTRALES A PROBAR
# ============================================================
#
# Primero hacemos una búsqueda gruesa cada 10 sujetos.
#
# ============================================================

tamanos_muestrales <- seq(
  from = 80,
  to = 120,
  by = 10
)


# ============================================================
# 9. NÚMERO DE SIMULACIONES
# ============================================================
#
# Recomendación:
#
# 1000 = exploración inicial
# 5000 = resultado final más estable
# 10000 = alta precisión
#
# Empieza con 1000.
#
# ============================================================

NSIM <- 5000


# ============================================================
# 10. EJECUTAR MONTE CARLO
# ============================================================

lista_resultados <- list()

contador <- 1


for(s in seq_len(nrow(escenarios))) {
  
  
  cat(
    "\n\n==========================================\n"
  )
  
  cat(
    "ESCENARIO:",
    escenarios$escenario[s],
    "\n"
  )
  
  cat(
    "a1 =", escenarios$a1[s],
    "| a2 =", escenarios$a2[s],
    "| b =", escenarios$b[s],
    "\n"
  )
  
  cat(
    "Indirecto MAIA =",
    escenarios$ind_MAIA_poblacional[s],
    "\n"
  )
  
  cat(
    "Indirecto Slope =",
    escenarios$ind_Slope_poblacional[s],
    "\n"
  )
  
  cat(
    "==========================================\n"
  )
  
  
  parametros_actuales <-
    escenarios[s, ]
  
  
  for(N_actual in tamanos_muestrales) {
    
    
    cat(
      "Simulando N =",
      N_actual,
      "\n"
    )
    
    
    resultado_actual <-
      calcular_potencia(
        
        N = N_actual,
        
        parametros =
          parametros_actuales,
        
        nsim = NSIM,
        
        alpha = .05
        
      )
    
    
    resultado_actual$escenario <-
      escenarios$escenario[s]
    
    
    resultado_actual$a1 <-
      escenarios$a1[s]
    
    
    resultado_actual$a2 <-
      escenarios$a2[s]
    
    
    resultado_actual$b <-
      escenarios$b[s]
    
    
    resultado_actual$ind_MAIA_poblacional <-
      escenarios$ind_MAIA_poblacional[s]
    
    
    resultado_actual$ind_Slope_poblacional <-
      escenarios$ind_Slope_poblacional[s]
    
    
    lista_resultados[[contador]] <-
      resultado_actual
    
    
    contador <- contador + 1
    
  }
  
}


# ============================================================
# 11. UNIR TODOS LOS RESULTADOS
# ============================================================

resultados_finales <-
  
  bind_rows(
    lista_resultados
  )


# Ordenar
resultados_finales <-
  
  resultados_finales %>%
  
  arrange(
    escenario,
    N
  )


cat(
  "\n\nRESULTADOS COMPLETOS:\n"
)

print(resultados_finales)

# ============================================================
# 12. N MÍNIMO:
# POTENCIA >= .80 PARA CADA INDIRECTO
# ============================================================
#
# Este criterio exige:
#
# power_ind_MAIA >= .80
# Y
# power_ind_Slope >= .80
#
# ============================================================

N_minimo_indirectos <-
  
  resultados_finales %>%
  
  filter(
    
    power_ind_MAIA >= .80,
    
    power_ind_Slope >= .80
    
  ) %>%
  
  group_by(
    escenario
  ) %>%
  
  slice_min(
    order_by = N,
    n = 1,
    with_ties = FALSE
  ) %>%
  
  ungroup()


cat(
  "\n\n==========================================\n"
)

cat(
  "N MÍNIMO PARA QUE AMBOS INDIRECTOS\n"
)

cat(
  "TENGAN POTENCIA INDIVIDUAL >= .80\n"
)

cat(
  "==========================================\n"
)

print(
  N_minimo_indirectos
)


# ============================================================
# 13. N MÍNIMO:
# PROBABILIDAD >= .80 DE DETECTAR LOS DOS
# EN EL MISMO ESTUDIO
# ============================================================
#
# Este es un criterio más exigente.
#
# ============================================================

N_minimo_ambos <-
  
  resultados_finales %>%
  
  filter(
    power_ambos >= .80
  ) %>%
  
  group_by(
    escenario
  ) %>%
  
  slice_min(
    order_by = N,
    n = 1,
    with_ties = FALSE
  ) %>%
  
  ungroup()


cat(
  "\n\n==========================================\n"
)

cat(
  "N MÍNIMO PARA TENER >= .80 DE PROBABILIDAD\n"
)

cat(
  "DE DETECTAR AMBOS INDIRECTOS SIMULTÁNEAMENTE\n"
)

cat(
  "==========================================\n"
)

print(
  N_minimo_ambos
)


# ============================================================
# 14. GRÁFICO:
# POTENCIA DE CADA EFECTO INDIRECTO
# ============================================================

resultados_largos <-
  
  resultados_finales %>%
  
  select(
    escenario,
    N,
    power_ind_MAIA,
    power_ind_Slope
  ) %>%
  
  tidyr::pivot_longer(
    
    cols = c(
      power_ind_MAIA,
      power_ind_Slope
    ),
    
    names_to = "efecto",
    
    values_to = "potencia"
    
  )


grafico_indirectos <-
  
  ggplot(
    resultados_largos,
    aes(
      x = N,
      y = potencia,
      linetype = efecto
    )
  ) +
  
  geom_line(
    linewidth = 0.9
  ) +
  
  geom_hline(
    yintercept = .80,
    linetype = "dashed"
  ) +
  
  facet_wrap(
    ~ escenario
  ) +
  
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, .10)
  ) +
  
  labs(
    
    title =
      "Potencia Monte Carlo de los efectos indirectos",
    
    subtitle =
      "Modelo de mediación con MAIA y exponente aperiódico como predictores",
    
    x =
      "Tamaño muestral (N)",
    
    y =
      "Potencia estadística",
    
    linetype =
      "Efecto indirecto"
    
  ) +
  
  theme_classic(
    base_size = 12
  )


print(
  grafico_indirectos
)


# ============================================================
# 15. GRÁFICO:
# PROBABILIDAD DE DETECTAR AMBOS INDIRECTOS
# ============================================================

grafico_ambos <-
  
  ggplot(
    resultados_finales,
    aes(
      x = N,
      y = power_ambos,
      linetype = escenario
    )
  ) +
  
  geom_line(
    linewidth = 1
  ) +
  
  geom_hline(
    yintercept = .80,
    linetype = "dashed"
  ) +
  
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, .10)
  ) +
  
  labs(
    
    title =
      "Potencia para detectar ambos efectos indirectos",
    
    x =
      "Tamaño muestral (N)",
    
    y =
      "Potencia estadística",
    
    linetype =
      "Escenario"
    
  ) +
  
  theme_classic(
    base_size = 12
  )


print(
  grafico_ambos
)


# ============================================================
# 16. GUARDAR RESULTADOS
# ============================================================

write.csv(
  
  resultados_finales,
  
  file.path(dir_montecarlo, "MonteCarlo_resultados_completos.csv"),
  
  row.names = FALSE
  
)


write.csv(
  
  N_minimo_indirectos,
  
  file.path(dir_montecarlo, "MonteCarlo_N_minimo_indirectos.csv"),
  
  row.names = FALSE
  
)


write.csv(
  
  N_minimo_ambos,
  
  file.path(dir_montecarlo, "MonteCarlo_N_minimo_ambos.csv"),
  
  row.names = FALSE
  
)


# ============================================================
# FIN
# ============================================================
