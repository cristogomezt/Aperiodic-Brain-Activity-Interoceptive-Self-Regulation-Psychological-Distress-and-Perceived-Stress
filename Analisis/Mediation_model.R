library(lavaan)
library(tidyverse)


# Ejecutar con la raíz del repositorio como directorio de trabajo
# (en RStudio: Session > Set Working Directory, o abrir el proyecto en la raíz).

# 1. Cargar base de datos ----
# Generada por fooof/processing.py
data <- read.csv("Analisis/df.csv")

# Carpeta de salida para las tablas
dir_tablas <- file.path("Analisis", "Tablas")
dir.create(dir_tablas, showWarnings = FALSE, recursive = TRUE)

# 2. Seleccionar variables ----

df <- na.omit(data[, c(
  "MAIA_Autorregulacion_DIRd",
  "Exponent__Global",
  "PSS_DIRt",
  "GHQ_DIRt"
)])     

# 3. Modelo ----
modelo <- '
  # Modelo del mediador
  PSS_DIRt ~ a1*MAIA_Autorregulacion_DIRd + a2*Exponent__Global

  # Modelo del outcome
  GHQ_DIRt ~ b*PSS_DIRt + c1*MAIA_Autorregulacion_DIRd + c2*Exponent__Global

  # Efectos indirectos
  ind_MAIA := a1*b
  ind_Exponent := a2*b

  # Efectos totales
  total_MAIA := c1 + (a1*b)
  total_Exponent := c2 + (a2*b)
'

set.seed(1234)

fit <- sem(
  modelo,
  data = df,
  se = "bootstrap",
  bootstrap = 5000,
  fixed.x = FALSE
)

res <- parameterEstimates(
  fit,
  standardized = TRUE,
  ci = TRUE,
  boot.ci.type = "bca.simple"
)
### SEGUNDA PARTE: GENERAR TABLA DE RESULTADOS----
library(flextable)
library(officer)
library(pagedown)

tabla_modelo <- res %>%
  filter(op %in% c("~", ":=")) %>%
  mutate(
    Effect = case_when(
      label == "a1" ~ "MAIA Self-regulation → Perceived stress",
      label == "a2" ~ "1/f exponent → Perceived stress",
      label == "b" ~ "Perceived stress → Psychological distress",
      label == "c1" ~ "MAIA Self-regulation → Psychological distress",
      label == "c2" ~ "1/f exponent → Psychological distress",
      label == "ind_MAIA" ~ "Indirect effect: MAIA Self-regulation → PSS → GHQ-12",
      label == "ind_Exponent" ~ "Indirect effect: 1/f exponent → PSS → GHQ-12",
      label == "total_MAIA" ~ "Total effect: MAIA Self-regulation → GHQ-12",
      label == "total_Exponent" ~ "Total effect: 1/f exponent → GHQ-12",
      TRUE ~ paste(lhs, op, rhs)
    ),
    B = round(est, 3),
    SE = round(se, 3),
    z = round(z, 3),
    sig = case_when(
      is.na(pvalue) ~ "",
      pvalue < .001 ~ "***",
      pvalue < .01  ~ "**",
      pvalue < .05  ~ "*",
      TRUE ~ ""
    ),
    
    p = case_when(
      is.na(pvalue) ~ "",
      pvalue < .001 ~ "< .001***",
      TRUE ~ paste0(sprintf("%.3f", pvalue), sig)
    ),
    `95% BCa CI` = paste0("[", round(ci.lower, 3), ", ", round(ci.upper, 3), "]"),
    beta = round(std.all, 3)
  ) %>%
  select(Effect, B, SE, z, p, `95% BCa CI`, beta)

# 4. Crear flextable ----

ft <- flextable(tabla_modelo) %>%
  set_header_labels(
    Effect = "Effect",
    B = "B",
    SE = "SE",
    z = "z",
    p = "p",
    `95% BCa CI` = "95% BCa CI",
    beta = "β"
  ) %>%
  theme_booktabs() %>%
  autofit() %>%
  align(align = "center", part = "all") %>%
  align(j = "Effect", align = "left", part = "all") %>%
  fontsize(size = 9, part = "all") %>%
  bold(part = "header")

save_as_docx(
  ft,
  path = file.path(dir_tablas, "tabla_modelo_mediacion.docx")
)

# 5. Exportar a PDF ----

save_as_html(
  ft,
  path = file.path(dir_tablas, "tabla_modelo_mediacion.html")
)

# Requiere Google Chrome o Chromium instalado
pagedown::chrome_print(
  input = file.path(dir_tablas, "tabla_modelo_mediacion.html"),
  output = file.path(dir_tablas, "tabla_modelo_mediacion.pdf")
)

# Matriz de correlacion---- 
library(Hmisc) #Para los Pvalue
library(corrplot) #Graficos
vars_matrix <- as.matrix(df)
# Correlaciones Pearson con p-values
res <- rcorr(vars_matrix, type = "pearson")

# Matriz de correlaciones 
res$r

# Matriz de p-values
res$P

# R^2 ----
r2 <- lavInspect(fit, "rsquare")

r2

summary(
  fit,
  standardized = TRUE,
  fit.measures = TRUE,
  rsquare = TRUE
)

# Estadisticos descriptivos ----
psych::describe(df)

