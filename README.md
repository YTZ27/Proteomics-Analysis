# Proteomics-Analysis

#Contexto del estudio
#Los datos corresponden a una cohorte con mediciones en tres timepoints (Baseline, Week.6, Week.12) bajo dos condiciones (Treated / Untreated), #reclutada en múltiples sitios. Las proteínas se midieron en dos paneles: Cardiometabolic e Inflammation.
#Dado que los datos provienen de dos corridas distintas, existe un batch effect sistemático entre proyectos que se corrige mediante bridge #normalization, usando muestras medidas en ambas corridas como referencia.


#Estructura del pipeline
olink_pipeline
│
├── 1.  Configuración del entorno
├── 2.  Carga de datos
├── 3.  Exploración pre-normalización       ← PCA, distribución NPX, correlación
├── 4.  Bridge normalization                ← olink_normalization_bridge()
├── 5.  Exploración post-normalización      ← PCA, scree plot, CV antes/después
├── 6.  Filtrado de proteínas               ← LOD ≤ 30% + SD ≥ 1.0
├── 7.  Detección de outliers               ← mediana ± 3×IQR por proteína/timepoint
├── 8.  Exploración adicional               ← UMAP por Treatment y Time
├── 9.  Expresión diferencial               ← LMM con olink_lmer()
└── 10. Exportación                         ← CSV normalizados y resultados

#Decisiones metodológicas clave
Bridge normalization
Se usó olink_normalization_bridge() tomando data1 como proyecto de referencia. El factor de ajuste (Adj_factor) se calcula proteína a proteína a partir de las bridge samples compartidas. Las bridge samples de data2 se excluyen del dataset final para evitar pseudoreplicación.

#Por qué scale. = TRUE en el PCA
El PCA post-bridge se corre con scale. = TRUE para ser consistente con el comportamiento interno de olink_pca_plot. Con datos Olink de múltiples paneles mezclados, la varianza se distribuye entre ~180 proteínas y los PCs individuales capturan poca varianza cada uno (~2-3%). Esto es matemáticamente esperado y no es un error — por eso el scree plot se genera por panel separado.

#Filtrado de proteínas
LOD ≤ 30%: proteínas con más del 30% de muestras bajo el límite de detección tienen señal insuficiente para análisis estadístico.
SD ≥ 1.0: proteínas con poca variabilidad entre muestras no aportarán poder estadístico en los modelos.

#Detección de outliers
El criterio IQR se aplica por proteína y por timepoint (group = c("OlinkID", "Time")). Agrupar solo por proteína mezclaría mediciones de Baseline con Week.12, distorsionando el rango normal. Una muestra se considera outlier global si >20% de sus proteínas caen fuera del rango mediana ± 3×IQR.


#De los modelos de expresion diferencial, El Modelo 3 es el apropiado para este diseño porque la pregunta biológica central es si el efecto del tratamiento difiere según el momento de la visita (Treatment:Time). El efecto aleatorio por sujeto (1|Subject) controla la correlación entre medidas repetidas del mismo individuo.

#Clonar el repositorio
git clone https://github.com/YTZ27/olink-pipeline.git
cd olink-pipeline

#Colocar los archivos de entrada en el directorio de trabajo
#Ajustar setwd() en la línea 26 del script

#Ejecutar el pipeline completo
source("olink_pipeline.R")
