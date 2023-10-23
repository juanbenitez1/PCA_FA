---
title: "FINAL AEM 2023 | UTDT"
author:
- "Juan Pablo Benitez"
date: "23/10/2023"
output: pdf_document
fontsize: 9pt
editor_options: 
  markdown: 
    wrap: 72
---

```{r setup}
rm(list = ls())
knitr::opts_knit$set(
  root.dir = "C:/Users/Acer Nitro 5/iCloudDrive/Di Tella/Analisis Multivariado/Examen")
knitr::opts_chunk$set(
  echo = TRUE,
  warnings = FALSE,
  message = FALSE)
options(scipen = 999)

suppressWarnings({
  library("dplyr")
  library("grid")
  library("gridExtra")
  library("ggpubr")
  library("tidyverse")
  library("kableExtra")
  library("knitr")
  library("ggplot2")
  library("factoextra")
  library("cluster")
  library("MVN")
  library("corrplot")
})
```

*Definimos las funciones a utilizar*

Función de summary:
```{r}
summ <- function(data) {
  data %>%
    select_if(is.numeric) %>%
    summarize_all(list(Media = ~ round(mean(.), 3),
                       Mediana = ~ round(median(.), 3),
                       Desvio = ~ round(sd(.), 3),
                       CV = ~ round((sd(.) / mean(.)), 3),
                       Min = ~ round(min(.), 3),
                       P25 = ~ round(quantile(., probs = 0.25), 3),
                       P50 = ~ round(median(.), 3),
                       P75 = ~ round(quantile(., probs = 0.75), 3),
                       Max = ~ round(max(.), 3))) %>%
    pivot_longer(everything(), names_sep = "_", names_to = c("variable", ".value"))
}
```

Función de summary para PCA:
```{r}
summ2 <- function(data) {
  data %>%
    select_if(is.numeric) %>%
    summarize_all(list(Media = ~ round(mean(.), 3),
                       Mediana = ~ round(median(.), 3),
                       Desvio = ~ round(sd(.), 3),
                       Min = ~ round(min(.), 3),
                       P25 = ~ round(quantile(., probs = 0.25), 3),
                       P50 = ~ round(median(.), 3),
                       P75 = ~ round(quantile(., probs = 0.75), 3),
                       Max = ~ round(max(.), 3))) %>%
    pivot_longer(everything(), names_sep = "_", names_to = c("variable", ".value"))
}
```

Función de ftest:
```{r}
fa_ftest <- function(dat,m,alpha){
  modelo <- psych::fa(scale(dat), nfactors=m,rotate='varimax', fm='ml', max.iter=1000)
  modelo$PVAL
  
  # Escalares
  p <- ncol(dat)
  n <- nrow(dat)
  
  # Matriz de correlacion Sigma estimada
  Co <- var(scale(dat)) 
  
  # Uniqueness
  e2 <- 1-modelo$communality
  
  # Matriz de cargas
  L  <- as.matrix(unclass(modelo$loadings))
  LL_tr <- L %*% t(L) # Lambda * Lambda' 
  
  # Varianzas especificas
  uu  <- diag(e2)
  Vo  <- LL_tr+uu  # Lambda*Lambda' + Psi
  
  # Armamos el estadístico lambda
  detLLPsi <- det(Vo)
  detS <- det(Co)
  
  # Corrección de Bartlett
  num <- ((n-1)-(2*p+4*m+5)/6)
  gl  <-  ((p-m)^2 - (p + m))/2
  
  #est <-  num*(ln0 - ln1)
  est <- n * log(detLLPsi/detS)
  estBartlett <- num * log(detLLPsi/detS) / n 
  
  # Computamos criterios de info
  aic <-  est - 2*gl
  bic <-  est - gl*log(n)
  
  # Para test
  v_c <- qchisq(1-alpha,gl)
  pval <- pchisq(est,gl,lower.tail = F)
  return(data.frame("Factors"=m,
                    "Stat"=est,
                    "Crit"=v_c,
                    "DoF"=gl,
                    "P-value"=pval,
                    "AIC" = aic,
                    "BIC" = bic))
}
```

\
\
\
Funcion gráfico matriz varianzas y covarianzas:
```{r}
cov_plot <- function(matrix_data,main_size,numbers_size,titulo){
  # Definimos la matriz de varianzas y covarianzas
  cov_matrix_rounded <- matrix_data
  
  # Establecer el tamaño del gráfico y la fuente
  par(mar = c(10, 10, 3, 3))  # Ajustar los márgenes
  cex.main <- num(main_size)  # Tamaño del título
  cex.axis <- num(numbers_size)  # Tamaño de la fuente de los números en las celdas
  
  # Crear un gráfico en blanco sin ejes
  plot(1, type = "n", 
       xlim = c(0.5, ncol(cov_matrix_rounded) + 0.5), 
       ylim = c(0.5, nrow(cov_matrix_rounded) + 0.5),
       xaxt = "n", 
       yaxt = "n", 
       xlab = "", 
       ylab = "", 
       xaxs = "i", 
       yaxs = "i")
  
  # Dibujar líneas divisorias de la matriz
  abline(h = 1:(nrow(cov_matrix_rounded)) - 0.5, col = "black", lwd = 2)
  abline(v = 1:(ncol(cov_matrix_rounded)) - 0.5, col = "black", lwd = 2)
  
  # Agregar valores numéricos en las celdas sin torcerlos y ajustar el tamaño de la fuente
  for (i in 1:ncol(cov_matrix_rounded)) {
    for (j in 1:nrow(cov_matrix_rounded)) {
      text(i, nrow(cov_matrix_rounded) - j + 1, cov_matrix_rounded[j, i], cex = cex.axis)
    }
  }
  
  # Agregar etiquetas de filas y columnas y ajustar el tamaño de la fuente
  axis(1, at = 1:ncol(cov_matrix_rounded), 
       labels = colnames(cov_matrix_rounded), 
       las = 2, cex.axis)
  
  axis(2, at = nrow(cov_matrix_rounded):1, 
       labels = rownames(cov_matrix_rounded), 
       las = 1, cex.axis)
  
  # Colorear solo la diagonal principal de celeste
  for (i in 1:ncol(cov_matrix_rounded)) {
    rect(i - 0.5, 
         nrow(cov_matrix_rounded) - i + 1.5, 
         i + 0.5, 
         nrow(cov_matrix_rounded) - i + 0.5, 
         col = rgb(0, 0, 1, alpha = 0.3), 
         border = "black")
  }
  # Agregar el título
  title(titulo, cex.main=cex.main)
}
```


(1) *A partir de la base de datos disponible elige un subconjunto de 5 a 20 indicadores y un año de referencia para efectuar tu análisis. Esta será la matriz de datos que utilizarás durante todo el examen. Cualquier subconjunto de información es válido, no hay elecciones incorrectas; sin embargo, verifica que la cantidad de observaciones disponibles supere a la cantidad de variables elegidas.*

Cargamos los datos
```{r}
data <- read.csv("./Data_Extract_from_WBWDI_Latam/WBWDI_Data.csv")
```
Generamos nuestro dataset filtrando por:

- Países del MERCOSUR (Argentina, Brasil, Paraguay, Uruguay y Venezuela) + 
  Asociados (Chile, Colombia, Ecuador, Perú y Bolivia) - Año 2008
  
- Variables:

    - EG.CFT.ACCS.ZS: Access to clean fuels and technologies for cooking (% of population)

    - AG.LND.AGRI.ZS: Agricultural land (% of land area)

    - EN.ATM.METH.AG.ZS: Agricultural methane emissions (% of total)

    - EN.ATM.NOXE.AG.ZS: Agricultural nitrous oxide emissions (% of total)

    - EG.USE.COMM.CL.ZS: Alternative and nuclear energy (% of total energy use)

    - EN.CO2.ETOT.ZS: CO2 emissions from electricity and heat production, total (% of total fuel combustion)

    - EN.CO2.MANF.ZS: CO2 emissions from manufacturing industries and construction (% of total fuel combustion)

    - EN.CO2.TRAN.ZS: CO2 emissions from transport (% of total fuel combustion)


```{r}
series_codes <- c("EG.CFT.ACCS.ZS", "AG.LND.AGRI.ZS", "EN.ATM.METH.AG.ZS",
                  "EN.ATM.NOXE.AG.ZS", "EG.USE.COMM.CL.ZS", "EN.CO2.ETOT.ZS", 
                  "EN.CO2.MANF.ZS", "EN.CO2.TRAN.ZS")

country_codes <- c("ARG","BRA","CHL","COL","ECU","PRY","PER","URY", "VEN", "BOL")

filtered_dataset <- data %>%
  filter(Series.Code %in% series_codes, Country.Code %in% country_codes) %>%
  select("Country.Name", "Country.Code", "Series.Name", "Series.Code", "X2008..YR2008.")

filtered_dataset <- filtered_dataset %>%
  rename("2008" = "X2008..YR2008.")

filtered_dataset$`2008` <- as.numeric(filtered_dataset$`2008`)
summary(filtered_dataset)
```
```{r}
unpivoted_dataset <- filtered_dataset %>%
  select(-"Country.Code",-"Series.Name")
  
unpivoted_dataset <- unpivoted_dataset %>%
  pivot_wider(
    names_from = Series.Code,
    values_from = '2008'
  )

segment1 <- unpivoted_dataset[, 1:5]
segment2 <- unpivoted_dataset[, c(1, 6:9)]

suppressWarnings(knitr::kable(segment1, caption = "Dataset (Part 1)", format ="latex", 
                              longtable=TRUE)%>%
  kableExtra::kable_styling(latex_options = "scale_down"))

suppressWarnings(knitr::kable(segment2, caption = "Dataset (Part 2)", format ="latex", 
                              longtable=TRUE)%>%
  kableExtra::kable_styling(latex_options = "scale_down"))
```

(2) *Realiza un breve análisis descriptivo de la información disponible para el conjunto de variables y año elegidos. Calcula y describe distintas medidas de variabilidad global de los datos.*

Describimos la base

Nombres de las columnas:
```{r}
names(unpivoted_dataset)
```
A continuación se presenta la tabla con las estadísticas descriptivas las variables numéricas del dataset, presentando medidas de dispersión como el desvío estandar y el coeficiente de variación, medidas de tendencia central como la Media y la Mediana, valores máximos y mínimos, y medidas de ubicación como los percentiles 25, 50 y 75.
```{r}
suppressWarnings(knitr::kable(summ(unpivoted_dataset), caption = "Summary", 
                              format ="latex", longtable=TRUE)%>%
  kableExtra::kable_styling(latex_options = "scale_down"))
```
Podemos observar que las variables EG.CFT.ACCS.ZS y EN.ATM.NOXE.AG.ZS tienen una media significativamente superiores a las demás, lo cual puede tener impacto en el análisis de componentes principales. 

Otra observación relevante es que, a excepción de las variables EG.CFT.ACCS.ZS y EN.CO2.TRAN.ZS, las medias de las variables están muy cercanas a la mediana, lo cual nos brinda indicios de que puede no haber tanta dispersión. Esto lo podemos ver también cuando observamos el Coeficiente de Variación de las variables, en este caso casi todas las variables presentan un CV inferior a 0.5, a excepeción de AG.LND.AGRI.ZS que presenta un CV levemente superior (0.515).

Además, podemos ver que los valores máximos y mínimos de las variables EG.USE.COMM.CL.ZS, EN.CO2.ETOT.ZS y EN.CO2.MANF.ZS son significativamente diferentes a los de las demás variables. Esto nos brinda indicios de que pueden existir outliers en diversas variables y es un aspecto que deberá ser considerado en el momento de aplicar el Análisis de Componentes Principales, dado que el PCA es sensible a la escala de los datos y a la varianza en las variables, por lo que la presencia de valores extremadamente altos en algunas variables puede influir en los resultados de manera desproporcionada.


(3) *Calcula la matriz de variancias y covariancias y la matriz de correlaciones muestrales e interpreta brevemente sus resultados.*

__Matriz de varianzas y covarianzas__

La matriz de varianzas y covarianzas, también llamada matriz S, presenta en la diagonal principal la varianza de cada variable, mientras que los elementos fuera de la diagonal son las covarianzas entre cada par de variables.

Cuando observamos las covarianzas, analizamos si existe una relación lineal entre cada par de variables. Para ello, hay ciertos aspectos que debemos tener en cuenta:

- Signo: Si la covarianza es positiva, indica que las dos variables tienden a aumentar juntas. Si la covarianza es negativa, indica que cuando una variable aumenta, la otra tiende a disminuir.

- Magnitud: La magnitud de la covarianza indica cuán fuerte es la relación lineal entre las variables. Cuanto mayor sea el valor absoluto de la covarianza, más fuerte será la relación lineal. En este sentido, un valor cercano a cero indica que no hay una relación lineal fuerte entre las variables.

- Escala de las variables: Para analizar la covarianza entre variables es importante que las mismas tengan la misma escala, de lo contrario puede ser muy difícil realizar una comparación entre datos que se encuentran medidos en diferentes ecalas.

Con respecto al último punto, en nuestro conjunto de datos todas las variables se encuentran expresadas en la misma unidad de medida, por lo que no tendremos inconvenientes por la escala al realizar las comparaciones.

```{r}
cov_matrix1 <- cov(unpivoted_dataset[, 2:ncol(unpivoted_dataset)])
matrix_data <- round(cov_matrix1, 3)
cov_plot(matrix_data, 1.2, 0.7, "Matriz de Varianzas y Covarianzas Muestral")
```
A partir de los valores de la diagonal principal podemos observar que las variables con mayor varianza, y por lo tanto mayor dispersión, son EG.CFT.ACCS.ZS, AG.LND.AGRI.ZS, EN.ATM.METH.AG.ZS, EN.CO2.ETOT.ZS y EN.CO2.TRAN.ZS.

Observando las covarianzas podemos ver que existen relaciones lineales fuertes entre algunas variables:

- EN.CO2.ETOT.ZS y EG.CFT.ACCS.ZS: en este caso la covarianza es de 164.838, lo cual puede ser considerado como una magnitud elevada teniendo en cuenta el resto de los valores de covarianzas. Con respecto al signo, la misma es positiva, por lo que se espera que estas dos variables tiendan a aumentar juntas.

- EN.CO2.TRAN.ZS y EG.CFT.ACCS.ZS: la covarianzas para este par de variables es una magnitud elevada y de signo negativo (-244.107), por lo que podemos esperar que cuando una variable aumente, la otra tienda disminuir.

Del mismo modo, podemos observar pares de variables donde la magnitud de la covarianza es muy baja:

- EG.USE.COMM.CL.ZS y AG.LND.AGRI.ZS: 0.544

- EN.CO2.TRAN.ZS y EN.ATM.NOXE.AG.ZS: 0.163

En ambos casos la covarianza es muy cercana a cero, por lo que tenemos indicios para pensar que no existe una relación lineal fuerte entre ambos pares de variables.

__Matriz de correlaciones muestrales__

A continuación se presenta la matriz de correlaciones muestrales, también llamada matriz R, que tal como su nombre lo indica nos muestra el grado de correlación entre las variables.

La interpretación de la matriz R es similar a la que realizamos para la matriz S en términos de signos y magnitudes. Sin embargo, hay dos diferencias importantes a señalar:

1) Mientras en la matriz S analizamos si existe una relación lineal entre dos variables, en la matriz R analizamos el grado de correlación existente.

2) En la matriz R los datos están escalados, por lo que el análisis no se vería afectado en caso de que las variables estuviesen expresadas en unidades de medida diferentes.

```{r}
cor_matrix <- cor(unpivoted_dataset[, 2:ncol(unpivoted_dataset)])
cor_matrix_rounded <- round(cor_matrix, 3)
corrplot(cor_matrix_rounded,
         title = "Matriz de correlaciones (R)",
         type = "lower",
         method = "color",
         diag = FALSE,
         tl.col = "black",
         tl.srt = 45,
         addCoef.col = "black",
         col = c("red", "green"),
         mar = c(0, 0, 1, 0),
         addgrid.col = "black",
         tl.cex = 0.7,
         number.cex = 0.8)
```
En el gráfico de la matriz se ha excluído la diagonal principal ya que la correlación de una variable con respecto a sí misma es igual a 1 y no aporta información extra. De manera similar, se ha excluído la diagonal superior de la matriz dado que es un espejo de la diagonal inferior.

Analizando los signos de las correlaciones, se presentan en color rojo las correlaciones negativas y en verde las correlaciones positivas.

Viendo las magnitudes, podemos ver que hay tres pares de variables cuya correlación puede ser considerada "fuerte", y son aquellas cuyas magnitud es mayor a 0.7:

- EN.CO2.TRAN.ZS y EG.CFT.ACCS.ZS: -0.77
- EN.CO2.TRAN.ZS y EN.CO2.ETOT.ZS: -0.87
- EN.CO2.ETOT.ZS y EG.CFT.ACCS.ZS: 0.75

En los dos primeros casos la correlación es negativa, lo que indica que cuando una de las variables aumenta, la otra tiende a disminuir. En el tercer caso, la correlación es positiva, por lo que ambas variables tienden a aumentar o disminuir juntas.

Luego, podemos observar casos de correlación "moderada" (valores de correlación entre 0.40 y 0.69), algunos de ellos son:

- EN.ATM.METH.AG.ZS y AG.LND.AGRI.ZS: 0.69
- EN.ATM.NOXE.AG.ZS y AG.LND.AGRI.ZS: 0.64
- EN.CO2.MANF.ZS y EN.CO2.ETOT.ZS: 0.31

También es importante resaltar que hay casos donde la correlación es muy débil o nula (valores entre 0 y 0.19), por ejemplo:

- EN.CO2.TRAN.ZS y EN.ATM.NOXE.AG.ZS: 0
- EN.CO2.TRAN.ZS y AG.LND.AGRI.ZS: 0.4
- EN.CO2.MANF.ZS y EG.USE.COMM.CL.ZS: 0.16

Es importante prestar atención a los pares de variables donde la correlación es muy débil, ya que la existencia de estos puede generar complicaciones al momento de realizar la aplicacion de modelos como PCA y Análisis Factorial.

(4) *Realiza un análisis de componentes principales utilizando tanto la matriz de variancias y covariancias como la matriz de correlaciones. Describe de manera general los resultados obtenidos en cada caso y compara los resultados destacando diferencias y/o similitudes entre ambos casos. En particular, describe la proporción de variancia explicada por cada componente principal. Elije los resultados de uno de estos análisis para continuar con la próxima consigna.*

Preparamos nuestro dataset para implementar análisis de componentes principales (PCA)
```{r}
dataset <- unpivoted_dataset[,2:ncol(unpivoted_dataset)]
dim(dataset)
```
Dado que la cantidad de observaciones es mayor que la cantidad de variables utilizadas, podemos implementar PCA.

Análisis de componentes principales utilizando matriz de correlaciones muestrales:
```{r}
pc <- princomp(dataset, cor=T)
summary(pc)

pc$loadings
```

```{r}
fviz_eig(pc, choice = c("variance"), addlabels = T)
```

```{r}
resumen <- summary(pc)

resumen_manual <- data.frame(
  Comp.1 = c(pc$sdev[1], pc$sdev[1]^2 / sum(pc$sdev^2), 
             cumsum(pc$sdev^2 / sum(pc$sdev^2))[1]),
  
  Comp.2 = c(pc$sdev[2], pc$sdev[2]^2 / sum(pc$sdev^2), 
             cumsum(pc$sdev^2 / sum(pc$sdev^2))[2])
)

rownames(resumen_manual) <- c("Standard Deviation", 
                              "Proportion of Variance", 
                              "Cumulative Proportion")

suppressWarnings(knitr::kable(summ(resumen_manual), caption = "Summary (PCA 1)", 
                              format ="latex", longtable=TRUE)%>%
  kableExtra::kable_styling(latex_options = "scale_down"))
```
\
\
Resultados relevantes:

- Cuando se observa la desviación estándar, el componente 1 siempre explicará la mayor variabilidad del dataset, seguido por el componente 2 y así sucesivamente. En este caso, el Componente 1 tiene una desviación estándar de 1.79 mientras que el Componente 2 tiene una desviación estándar de 1.57.

- Con respecto a la proporción de Varianza Explicada, se observa que el Componente 1 explica el 40.25% de la varianza total, mientras que el Componente 2 explica el 30.85% de la varianza total.

- Si contemplamos solamente los primeros dos componentes principales obtenidos a partir del modelo, se observa que entre ambos explican aproximadamente el 71% de la varianza de nuestro conjunto de datos. 

Cuando observamos el Scree Plot vemos que, si bien es un buen enfoque considerar solamente los dos primeros componentes principales, no estaría mal si consideráramos 3 componentes en lugar de 2. De hecho, si se agregamos un componente más a nuestro análisis (Componente 3), la varianza acumulada entre los tres sería de aproximadamente el 85%. Sin embargo, a fines prácticos y de acuerdo a lo solicitado en la consigna, se ha procedido conservando únicamente los dos primeros componentes.


Los loadings indican cómo cada variable original se relaciona con cada componente principal. Los valores en esta tabla muestran la fuerza y la dirección de la relación.

Analizando los resultados del Componente 1 podemos observar que se encuentra fuertemente relacionado con las siguientes variables:

- EN.CO2.TRAN.ZS: CO2 emissions from transport (% of total fuel combustion) (Relacion negativa)

- EN.CO2.ETOT.ZS: CO2 emissions from electricity and heat production, total (% of total fuel combustion) (Relacion positiva)

- EG.CFT.ACCS.ZS: Access to clean fuels and technologies for cooking (% of population) (Relacion positiva)

Por otro lado, las variables que menor peso tienen en el Componente 1 son:

- AG.LND.AGRI.ZS: Agricultural land (% of land area) (Relación negativa)

- EN.ATM.NOXE.AG.ZS: Agricultural nitrous oxide emissions (% of total) (Relación negativa)

- EN.ATM.METH.AG.ZS: Agricultural methane emissions (% of total) (Relación negativa)

- EG.USE.COMM.CL.ZS: Alternative and nuclear energy (% of total energy use) (Relación negativa)

Por lo tanto, podríamos intuir que el Componente 1 representa las emisiones de gases de efecto invernadero relacionadas al consumo de energía en aspectos del día a día, como ser transporte, cocina, calefacción.

En el caso del Componente 2, las variables con mayor peso en el mismo son:

- AG.LND.AGRI.ZS: Agricultural land (% of land area) (Relación positiva)

- EN.ATM.NOXE.AG.ZS: Agricultural nitrous oxide emissions (% of total) (Relación positiva)

- EN.ATM.METH.AG.ZS: Agricultural methane emissions (% of total) (Relación positiva)

Y las de menor peso son:

- EN.CO2.TRAN.ZS: CO2 emissions from transport (% of total fuel combustion) (Relación negativa)

- EN.CO2.ETOT.ZS: CO2 emissions from electricity and heat production, total (% of total fuel combustion) (Relación positiva)

- EG.CFT.ACCS.ZS: Access to clean fuels and technologies for cooking (% of population) (Relación positiva)

- EG.USE.COMM.CL.ZS: Alternative and nuclear energy (% of total energy use) (Relación negativa)

En este caso, podríamos pensar que el Componente 2 representa las emisiones de CO2 derivadas de las actividades  agrícolas.

Un aspecto importante a resaltar es que, en este caso, las variables que mayor peso tienen en el Componente 1 son las que menor peso tienen en el Componente 2, y viceversa.

Por otro lado, la variable EG.USE.COMM.CL.ZS: Alternative and nuclear energy (% of total energy use) parece no tener mucha relación con los Componentes 1 y 2, pero viendo los resultados se observa que si estaría fuertemente relacionada con el Componente 3 en caso de contemplarlo en nuestro análisis.
  

Análisis de componentes principales utilizando matriz de varianzas y covarianzas:
```{r}
pc2 <- princomp(dataset)
summary(pc2)
pc2$loadings
```

```{r}
fviz_eig(pc2, choice = c("variance"), addlabels = T)
```

```{r}
resumen <- summary(pc2)

resumen_manual <- data.frame(
  Comp.1 = c(pc2$sdev[1], pc2$sdev[1]^2 / sum(pc2$sdev^2), 
             cumsum(pc2$sdev^2 / sum(pc2$sdev^2))[1]),
  
  Comp.2 = c(pc2$sdev[2], pc2$sdev[2]^2 / sum(pc2$sdev^2), 
             cumsum(pc2$sdev^2 / sum(pc2$sdev^2))[2])
)

rownames(resumen_manual) <- c("Standard Deviation", 
                              "Proportion of Variance", 
                              "Cumulative Proportion")

suppressWarnings(knitr::kable(summ(resumen_manual), caption = "Summary (PCA 2)", 
                              format ="latex", longtable=TRUE)%>%
  kableExtra::kable_styling(latex_options = "scale_down"))
```
\
\
Resultados relevantes:

Al momento de observar la desviación estándar, varianza explicada y varianza acumulada, las interpretaciones para el caso de la matriz de varianzas y covarianzas serán similares a las que obtuvimos en el caso del PCA calculado a partir de la matriz de correlaciones. Sin embargo, existen algunas diferencias:

- Cuando se observa la desviación estándar, se mantiene la idea de que el componente 1 siempre explicará la mayor variabilidad del dataset, seguido por el componente 2 y así sucesivamente. Pero en este caso se puede observar que, en general, el desvío estándar de los distintos componentes es significativamente mayor a los obtenidos en el caso de la matriz de correlación. Esto ocurre porque la matriz de varianzas y covarianzas conserva la escala y la varianza original de las variables. Por lo tanto, las desviaciones estándar son más grandes en este contexto. Entonces, vemos que el componente 1 tiene un desvío estándar de aproximadamente 25.90 y el componente 2 tiene un desvío estándar de aproximadamente 22.70

- La interpretación de la varianza explicada será similar pero también veremos como afecta el hecho de que la matriz de varianzas y covarianzas no estandarice las variables. En este caso, el Componente 1 explica el 47.81% de la varianza total mientras el Componente 2 explica el 36.72% de la varianza total. Ambos valores son mayores en comparación con el análisis basado en la matriz de correlaciones. Esto ocurre porque, al no estandarizar las variables, la varianza original de cada variable tiene un mayor impacto en la explicación de la varianza total. 

- La varianza acumulada explicada por los componentes 1 y 2 es aún mayor en este caso, explicando en conjunto un 84.53% de la varianza total de los datos. Esto nso brinda indicios de que estos dos componentes principales retienen una mayor cantidad de información sobre la varianza y la covarianza original en los datos. 
A diferencia del análisis realizado con la matriz de correlación, en este caso tenemos casi un 85% de varianza explicada por los primeros dos componentes y, si bien el Scree Plot brinda indicios de que podríamos considerar un tercer componente, se puede considerar como suficiente para avanzar con el análisis. 

Ahora que hemos utilizado la matriz de varianzas y covarianzas, los loadings aún indican cómo cada variable original se relaciona con cada componente principal pero ahora lo hacen en términos de covarianzas. Es decir, los loadings en el PCA basado en la matriz de varianzas y covarianzas representan las covarianzas entre las variables originales y los componentes.

Analizando los resultados del Componente 1 podemos observar que se encuentra fuertemente relacionado con las siguientes variables:

- EN.CO2.TRAN.ZS: CO2 emissions from transport (% of total fuel combustion) (Covarianza negativa)

- EG.CFT.ACCS.ZS: Access to clean fuels and technologies for cooking (% of population) (Covarianza positiva)

- EN.CO2.ETOT.ZS: CO2 emissions from electricity and heat production, total (% of total fuel combustion) (Covarianza positiva)

Se observa que son las mismas tres variables que se encontraban fuertemente asociadas con el componente 1 cuando corrimos el PCA basado en la matriz de correlaciones, con la diferencia de que ahora el acceso a energías limpias y tecnologías para cocinar tienen un mayor peso que las emisiones de CO2 derivadas de la electricidad y generación de calor.

Por otro lado, las variables que menor peso tienen en el Componente 1 son:

- AAG.LND.AGRI.ZS: Agricultural land (% of land area) (Covarianza nula)

- EN.ATM.NOXE.AG.ZS: Agricultural nitrous oxide emissions (% of total) (Covarianza nula)

- EG.USE.COMM.CL.ZS: Alternative and nuclear energy (% of total energy use) (Covarianza nula)

- EN.ATM.METH.AG.ZS: Agricultural methane emissions (% of total) (Covarianza negativa)

Nuevamente, se observa que las variables con menor fuerza de asociación o menor covarianza con el componente 1 son las mismas que presentaban la menor relación en el caso del PCA basado en la matriz de correlaciones.

Una diferencia importante es que en este caso las tres primeras variables mencionadas tienen una covarianza nula y la cuarta tiene una covarianza cercana a cero, por lo que podemos concluir que estas 4 variables prácticamente no influyen significativamente en la definición del componente.

Por lo tanto, al igual que en el PCA basado en la matriz de correlaciones, podríamos intuir que el Componente 1 representa las emisiones de CO2 relacionadas al consumo de energía en aspectos del día a día, como ser transporte, cocina, calefacción.

En el caso del Componente 2, las variables con mayor peso en el mismo son:

- AG.LND.AGRI.ZS: Agricultural land (% of land area) (Covarianza positiva)

- EN.ATM.METH.AG.ZS: Agricultural methane emissions (% of total) (Covarianza positiva)

Podemos observar dos de las tres variables que se encontraban fuertemente asociadas al componente 2 cuando corrimos el PCA basado en la matriz de correlaciones. En este caso, las variables que presentan una mayor covarianza con el componente 2 son “Agricultural land” y “Agricultural methane emissions”.

Y las de menor peso son:

- EG.USE.COMM.CL.ZS: Alternative and nuclear energy (% of total energy use) (Covarianza nula)

- EN.CO2.ETOT.ZS: CO2 emissions from electricity and heat production, total (% of total fuel combustion) (Covarianza nula)

- EN.CO2.TRAN.ZS: CO2 emissions from transport (% of total fuel combustion) (Covarianza nula)

Similar al caso del componente 1, podemos ver que las tres variables con menor grado de asociación con el componente 2 presentan una covarianza nula o cercana a cero, por lo que prácticamente no influyen en la definición de este componente.

En este caso, podríamos pensar que el Componente 2 representa las emisiones de CO2 derivadas de las actividades agrícolas.


Para poder definir si utilizamos la matriz de varianzas y covarianzas o la matriz de correlaciones muestrales para el PCA debemos tener en cuenta principalmente dos cuestiones:

1) Unidad de medida de las variables.

2) Presencia de outliers.

Estos aspectos son importantes a considerar si queremos utilizar la matriz de varianzas y covarianzas, ya que para ello necesitamos que todas las variables se encuentren en la misma unidad de medida y que no existan outliers. Esto último se debe a que la presencia de outliers genera varianzas muy elevadas con respecto a las varianzas de las demás variables, ocasionando que el primer componente obtenido del PCA se vea fuertemente representado solo por la/las variables con outliers y dejando con poco peso a las demás variables.

Sabemos que todas las variables seleccionadas en el presente trabajo se encuentran presentadas como proporciones, por lo que se cumple el primer item, es decir, tienen todas las misma unidad de medida.

A continuación, procedemos a elaborar un boxplot para cada variable a fin de observar si existen outliers:
```{r}
par(mfrow = c(3, 3))
par(mar = c(1, 4, 1, 1))

column_names <- names(dataset)

for (i in 1:8) {
  bp <- boxplot(dataset[, i], 
                main = column_names[i], 
                ylab = "Valor", 
                outline = TRUE)
  for (outlier in bp$out) {
    points(outlier, col = "red", pch = 20, cex = 2)
  }
}

par(mfrow = c(1, 1))
```
De acuerdo a lo presentado en los boxplots podemos ver que las siguientes tres variables presentan outliers:

- AG.LND.AGRL.ZS: Agricultural land (% of land area)

- EN.ATM.NOXE.AG.ZS: Agricultural nitrous oxide emissions (% of total)

- EN.CO2.TRAN.ZS: CO2 emissions from transport (% of total fuel combustion)

Según lo evaluado en los análisis de componentes principales, las tres variables son significativas en la determinación del componente 1 o del componente 2, por lo que no es posible eliminarlas de nuestro modelo.

Por lo tanto, se utilizarán los resultados obtenidos en el primer Análisis de Componentes Principales realizado a partir de la Matriz de Correlaciones Muestrales.

(5) *Tomando en cuenta las primeras dos componentes principales del análisis elegido en el punto anterior calcula los coeficientes de correlación de cada componente principal con respecto a cada variable original. Interpreta la representación de las primeras 2 componentes principales en términos de su correlación con las variables que las definen.*

\
\
Tomamos los dos primeros componentes:
```{r}
dos_comp <- pc$scores[, 1:2]
suppressWarnings(knitr::kable(dos_comp, caption = "Scores Comp.1 y Comp.2", 
                              format ="latex", longtable=TRUE)%>%
  kableExtra::kable_styling(latex_options = "scale_down"))
```

```{r}
pc$loadings
```
Calculamos la matriz de correlación de ambos componentes con las variables originales:

```{r}
# Crea la matriz de correlación
cor_matrix <- cor(dos_comp, dataset)

# Redondea los valores de la matriz de correlación
cor_matrix_rounded <- round(cor_matrix, 3)

# Crea la tabla de correlación
# Crea la tabla de correlación
corrplot(cor_matrix_rounded,
         title = "Matriz de correlaciones (R)",
         type = "full",
         method = "color",
         diag = TRUE,
         tl.col = "black",
         tl.srt = 45,
         addCoef.col = "black",
         col = c("red", "green"),
         mar = c(0, 0, 1, 0),
         addgrid.col = "black",
         tl.cex = 0.7,
         number.cex = 0.8,
         number.digits = 3)
```
Podemos ver que el Componente 1 se encuentra fuertemente correlacionado con tres variables:

1) EN.CO2.TRAN.ZS: CO2 emissions from transport (% of total fuel combustion)

2) EG.CFT.ACCS.ZS: Access to clean fuels and technologies for cooking (% of population)

3) EN.CO2.ETOT.ZS: CO2 emissions from electricity and heat production, total (% of total fuel combustion)

Siendo negativa la correlación en el primer caso y positiva en los dos casos restantes.

Por otro lado, tenemos tres variables que presentan un grado de correlación que podríamos considerar moderado:

1) EN.CO2.MANF.ZS: CO2 emissions from manufacturing industries and construction (% of total fuel combustion)

2) EN.ATM.METH.AG.ZS: Agricultural methane emissions (% of total)

3) EG.USE.COMM.CL.ZS: Alternative and nuclear energy (% of total energy use)

En el primer caso la correlación es positiva y en los dos casos restantes es negativa.

Finalmente, el Componente 1 presenta un grado de correlación "bajo" con las siguientes variables:

1) AG.LND.AGRI.ZS: Agricultural land (% of land area)

2) EN.ATM.NOXE.AG.ZS: Agricultural nitrous oxide emissions (% of total)

Realizando un análisis similar para la correlación entre el Componente 2 y las variables originales podemos ver que este componente se encuentra fuertemente correlacionado con dos variables:

1) AG.LND.AGRI.ZS: Agricultural land (% of land area)

2) EN.ATM.NOXE.AG.ZS: Agricultural nitrous oxide emissions (% of total)

En ambos casos la correlación es negativa.

Luego, se pueden ver dos variables con un grado de correlación moderado:

1) EN.ATM.METH.AG.ZS: Agricultural methane emissions (% of total)

2) EG.CFT.ACCS.ZS: Access to clean fuels and technologies for cooking (% of population)

Presentando correlación positiva en ambos casos.

Finalmente, las 4 variables restantes presentan un grado de correlación bajo:

1) EN.CO2.TRAN.ZS: CO2 emissions from transport (% of total fuel combustion)

2) EG.USE.COMM.CL.ZS: Alternative and nuclear energy (% of total energy use)

3) EN.CO2.MANF.ZS: CO2 emissions from manufacturing industries and construction (% of total fuel combustion)

4) EN.CO2.ETOT.ZS: CO2 emissions from electricity and heat production, total (% of total fuel combustion)

Donde las tres primeras variables presentan correlación negativa y la restante es positiva.

(6)*Calcula el valor de las primeras 2 componentes principales para cada una de las observaciones (países). Realiza un breve análisis de estadística descriptiva de cada una de ellas.*

```{r}
two_components <- pc$scores[, 1:2]
two_components_df <- data.frame("Country" = unpivoted_dataset$Country.Name, 
                                two_components)

two_components_df$Comp.1 <- as.numeric(two_components_df$Comp.1)
two_components_df$Comp.2 <- as.numeric(two_components_df$Comp.2)

suppressWarnings(knitr::kable(two_components_df, caption = "Comp.1 y Comp.2 por país", 
                              format ="latex", longtable=TRUE)%>%
  kableExtra::kable_styling(latex_options = "scale_down"))
```
Estadísticas descriptivas de los dos componentes:
```{r}
# Utilizamos la función summ2 definida al principio
suppressWarnings(knitr::kable(summ2(two_components_df), caption = "Summary", 
                              format ="latex", longtable=TRUE)%>%
  kableExtra::kable_styling(latex_options = "scale_down"))
```
Sabemos que en una combinación lineal la suma de los componentes es 0. Por lo tanto, dado que los componentes principales son combinaciones lineales, es de esperar que la media de los mismos sea igual a 0 también. En este caso, no tiene sentido calcular el coeficiente de variación ya que el cálculo del mismo implica utilizar a la media.

Por otro lado, viendo el desvío estandar se observa que el componente 1 presenta una mayor dispersión que el componente 2. Esto también se ver representado en la diferencia que existe entre los valores mínimos y máximos. Cabe señalar que este resultado es lo esperado, ya que en el análisis de componentes principales el componente 1 es el que captura la mayor varianza del conjunto de datos, seguido por el componente 2 y los demás.

(7)*Realiza un gráfico biplot para la representación conjunta de filas y columnas (equivalente a observaciones y variables) de la matriz de datos y describe los resultados obtenidos en términos de la distribución de las observaciones (países) en el subespacio de las primeras dos componentes.*

```{r}
dataset2 <- unpivoted_dataset
dataset2 <- column_to_rownames(dataset2, "Country.Name")
pc1 <- princomp(dataset2, cor=T)

fviz_pca_biplot(pc1, 
                axes = c(1,2),
                xlim = c(-4.5, 4),
                ylim = c(-2, 4))
```
A partir del Biplot podemos realizar distintos análisis:

En primer lugar, podemos ver la relación existente entre las variables del dataset y los componentes principales 1 y 2. Según el ángulo que presente el vector de una variable con respecto a uno de los componentes, podemos señalar que dicha variable está asociada más a un componente que al otro. Además, la longitud del vector nos brinda indicios de cuán fuerte es la relación que tiene la variable con el componente asociado.

En este caso, podemos ver que las siguientes variables presentan un ángulo cercano al eje de las abscisas, por lo que presentan una mayor asociación con el Componente Principal 1:

- EN.CO2.TRAN.ZS: CO2 emissions from transport

- EG.CFT.ACCS.ZS: Access to clean fuels and technologies for cooking

- EN.CO2.ETOT.ZS: CO2 emissions from electricity and heat production

Teniendo en cuenta el ángulo de los vectores, podemos señalar que las emisiones de CO2 derivadas del transporte son las que presentan una relación mas fuerte con el Componente 1. Esto se encuentra alineado con el análisis de los loadings realizado previamente.

Cuando observamos las variables mas cercanas al eje de las ordenadas, y por lo tanto mas asociadas al Componente Principal 2, son:

- AG.LND.AGRI.ZS: Agricultural land

- EN.ATM.NOXE.AG.ZS: Agricultural nitrous oxide emissions

- EN.ATM.METH.AG.ZS: Agricultural methane emissions

En ese orden según el peso que tienen en el Componente Principal 2 según el ángulo de los vectores. Al igual que en el caso analizado del Componente Principal 1, este análisis se encuentra alineado al análisis realizado de los loadings.

Por otro lado, la proximidad de los vectores entre sí nos brinda indicios de que las variables cercanas pueden estar fuertemente relacionadas entre sí. Por ello podemos suponer que existe una fuerte relación entre las siguientes variables:

- EN.CO2.TRAN.ZS y EG.USE.COMM.CL.ZS

- EG.CFT.ACCS.ZS y EN.CO2.ETOT.ZS

- EN.ATM.METH.AG.ZS, EN.ATM.NOXE.AG.ZS y AG.LND.AGRI.ZS

En los casos en que los vectores forman un ángulo próximo a 90° podemos suponer que las variables representadas por dichos vectores son independientes, por ejemplo:

- EN.ATM.NOXE.AG.ZS y EG.USE.COMM.CL.ZS

- EN.ATM.NOXE.AG.ZS y EG.CFT.ACCS.ZS

Esto tiene sentido cuando observamos la matriz de correlaciones de las variables originales y vemos que para ambos casos la correlación es cercana a cero, siendo -0.091 en el primer caso y 0.07 en el segundo.


Cuando observamos la distribución de los países en el Biplot podemos señalar algunos resultados importantes:

- Los países Brasil, Colombia, Ecuador, Bolivia y Perú están cercanos al origen, lo que nos brinda indicios de que pueden no estar caracterizados por ninguna de las variables en cuestión.

- Cuando dos vectores tienen un ángulo próximo a 180° decimos que dichas variables presentan una correlación negativa. En este sentido podemos señalar que Argentina se caracteriza tener mayor acceso a energías limpias y tecnologías para cocinar (EG.CFT.ACCS.ZS) y también por emitir CO2 derivado de la generación de eletricidad y calor (EN.CO2.ETOT.ZS) pero, a su vez, se caracteriza por emitir poco CO2 derivado de los transportes con respecto a Paraguay (EN.CO2.TRAN.ZS).

- En sentido inverso podemos señalar que Paraguay se caracteriza por tener un mayor grado de emisión de CO2 derivado del transporte (EN.CO2.TRAN.ZS), pero también se caracteriza por producir un grado menor de CO2 derivado de la producción de electricidad y calor (EN.CO2.ETOT.ZS) con respecto a Argentina.

- En el caso de Uruguay, podemos ver que se caracteriza por tener un mayor grado de tierras dedicadas a la agricultura (AG.LND.AGRI.ZS) y también un mayor grado de emisión de gases de efecto invernadero como el Óxido Nitroso (EN.ATM.NOXE.AG.ZS) y Metano (EN.ATM.METH.AG.ZS)

- Observando a Chile, se puede ver que se encuentra caracterizado principalmente por la emisión de CO2 derivada de los procesos de manufactura (EN.CO2.MANF.ZS)


(8)*Utilizando la representación biplot adecuada calcula una aproximación de dimensión 2 para la matriz de variancias y covariancias (o matriz de correlaciones, según corresponda) original e interpreta los resultados obtenidos.*

En este caso se ha decidido avanzar con el análisis de componentes principales efectuado con la matriz de correlaciones. Por lo tanto, a continuación se realizará la aproximación de la matriz de correlaciones (R).

Una limitación importante a resaltar es que en R no podemos extraer los datos directamente del biplot. Sin embargo, al extraer directamente de pc1, obtenemos las cargas y puntuaciones de los dos primeros componentes principales, lo cual es esencialmente lo mismo que se obtendría del gráfico biplot. Esto se debe a que el biplot combina esta información en un gráfico que visualiza las relaciones entre las variables y las observaciones en función de los componentes principales. 

Por lo tanto, no deberíamos perder información relevante al calcular la aproximación de la matriz de correlaciones original utilizando los resultados extraídos directamente de pc1.

```{r}
# Extraemos los desvíos estandar y los autovalores de pc1
SD <- diag(pc1$sdev, nrow = length(pc1$sdev), ncol = length(pc1$sdev))
Autovalor <- pc1$sdev

# Definimos la matriz D: matriz diagonal compuesta por los dos primeros autovalores
D <- diag(Autovalor[1:2])

# Calculamos la matriz A a partir de los loadings
A <- pc1$loadings[, 1:2]

# Calculamos la matriz C = D*A'
C <- D %*% t(A)

# Calculamos la matriz R aproximada
R_aprox <- t(C) %*% C

#R_aprox
corrplot(R_aprox,
         title = "Matriz de correlaciones aproximada (R aprox)",
         type = "lower",
         method = "color",
         diag = TRUE,
         tl.col = "black",
         tl.srt = 45,
         addCoef.col = "black",
         col = c("red", "green"),
         mar = c(0, 0, 1, 0),
         addgrid.col = "black",
         tl.cex = 0.7,
         number.cex = 0.8)
```
A partir de los datos extraídos del PCA podemos aproximar la matriz de correlación (R tilde) considerando la siguiente fórmula: *A x D x A^T* donde la matriz A representa los loadings del PCA (autovectores) y D es una matriz diagonal de 2x2 cuyos elementos diagonales son los autovalores de los respectivos componentes principales 1 y 2.

Analizando las correlaciones de las variables, podemos observar que las aproximaciones de las variables relacionadas al componente 1 y al componente 2 son muy buenas:

- Componente 1:
  
  - EN.CO2.TRAN.ZS: 0.94
  - EN.CO2.ETOT.ZS: 0.88
  - EG.CFT.ACCS.ZS: 0.75
  
- Componente 2:
  - AG.LND.AGRI.ZS: 0.84
  - EN.ATM.NOXE.AG.ZS: 0.67
  - EN.ATM.METH.AG.ZS: 0.74

**Tomando en cuenta las primeras 2 (dos) componentes principales que representaban el mayor porcentaje de variancia explicada y que calculaste en los ítems anteriores, resuelve los siguientes puntos:**

(9)*Realiza un análisis de clusters para hallar la jerarquía de agrupación de los países considerados en tu matriz de datos, en base a la información de las dos primeras componentes principales. Utiliza dos métodos distintos (puedes elegir entre los encadenamientos simple, completo, promedio, etc.).*

De acuerdo a lo estudiado en la bibliografía, en general los métodos con mejor desempeño suelen ser el método de encadenamiento promedio y el método de Ward.

**Método: Encadenamiento simple**

Este método calcula las distancias promedio entre todos los pares de puntos de los grupos que queremos agrupar y es útil cuando se buscan agrupamientos compactos en el espacio de los dos primeros componentes de un análisis PCA.

```{r}
returns <- column_to_rownames(two_components_df,var="Country")

d <- dist(returns)
fit <- hclust(d, method = "average")

print(fit)
```
**Método: Ward**

El método de Ward minimiza la varianza dentro de los grupos que se están agrupando. 

```{r}
returns <- column_to_rownames(two_components_df,var="Country")

d <- dist(returns)
fit1 <- hclust(d, method = "ward.D")

print(fit1)
```
(10)*En cada caso realiza una representación gráfica de los resultados obtenidos utilizando un dendrograma.*

**Dendrograma: Encadenamiento promedio**
```{r}
plot(as.dendrogram(fit), type="rectangle")
rect.hclust(fit, k = 5, border = 2:7)
abline(h = 2.1, col = 'red')
```

**Dendrograma: Método de Ward**
```{r}
plot(as.dendrogram(fit1), type="rectangle")
rect.hclust(fit1, k = 5, border = 2:7)
abline(h = 2.4, col = 'red')
```
Se han elaborado dos dendrogramas a partir de un análisis de clusters jerárquico donde se ha implementado el método de encadenamiento simple en primer lugar y el método de Ward en el segundo caso.

La primer observación importante que se puede realizar es que ambos métodos brindan agrupamientos idénticos, de modo que los clústers de países para los Componentes 1 y 2 del análisis PCA realizado previamente se conforman de la siguiente manera:

- Chile y Venezuela

- Perú, Brasil, Colombia y Ecuador

- Paraguay

- Argentina y Bolivia

Cuando comparamos estos agrupamientos con el Biplot elaborado anteriormente en el punto 7 vemos que los agrupamientos obtenidos por estos métodos tienen sentido. El único país del que podríamos dudar con respecto al agrupamiento es Bolivia, ya que en el Biplot parece estar mas próximo al grupo de Perú, Brasil, Colombia y Ecuador que a Argentina.

Sin embargo, podemos considerar que el análisis de clústers se asemeja a los agrupamientos que podemos observar en el biplot realizado a partir del análisis PCA.

**Por último**

(13)*Efectúa un análisis factorial para describir la variabilidad común entre las variables de tu matriz de datos. El punto de partida será un modelo con un solo factor.*

Para llevar a cabo un análisis factorial, en primer lugar debemos testear si nuestros datos cumplen con el supuesto de normalidad multivariada, para ello aplicaremos dos pruebas para testear la siguiente hipótesis:

\
\
Ho: los datos siguen una distribución normal multivariada.

Ha: los datos no siguen una distribución normal multivariada.

__Prueba de Mardia__

La prueba de Mardia se basa en estadísticas de momentos multivariados: Skewness y Kurtosis.

Para esta prueba, la función calcula los coeficientes de Skewness y Kurtosis de Mardia, así como su significación estadística correspondiente.
```{r}
result1 <- mvn(data = unpivoted_dataset[,2:ncol(unpivoted_dataset)], mvnTest = "mardia")

suppressWarnings(knitr::kable(result1$multivariateNormality, caption = "Prueba de Mardia", 
                              format ="latex", longtable=TRUE)%>%
  kableExtra::kable_styling(latex_options = "scale_down"))
```
Si bien observamos que en la última columna se indica que el dataset sigue una distribución normal multivariada, también lo podemos observar a partir de los p-valores de la Skewness y la Kurtosis, donde en ambos casos el p-valor es mayor al nivel de significancia 0.05, lo cual no permite rechazar la hipótesis nula de que los datos siguen una distribución normal multivariada.

__Prueba de Henze-Zirkler__

La prueba de Henze-Zirkler se basa en una distancia multivariada entre la distribución de los datos y una distribución normal multivariada. Evalúa si los datos se distribuyen de manera similar a una distribución multivariada normal estándar.
```{r}
result2 <- mvn(data = unpivoted_dataset[,2:ncol(unpivoted_dataset)], mvnTest = "hz")

suppressWarnings(knitr::kable(result2$multivariateNormality, caption = "Prueba de Henze-Zirkler", 
                              format ="latex", longtable=TRUE)%>%
  kableExtra::kable_styling(latex_options = "scale_down"))
```
Del mismo modo que hemos analizado la prueba de Mardia, vemos que la columna MVN nos indica que los datos siguen una distribución normal multivariada, lo cual podemos validar con el p-valor = 0.08 > alpha = 0.05.

Dado que tenemos suficiente evidencia para considerar que nuestros datos tienen distribución normal multivariada, podemos avanzar con el desarrollo del Análisis Factorial utilizando el método de máxima verosimilitud para el cálculo de los factores.

En primer lugar, se realizará la Prueba de Esfericidad de Barlett, la cual contrasta la hipótesis nula de que la matriz de correlaciones es una matriz identidad, en cuyo caso no existirían correlaciones significativas entre las variables y el modelo factorial no sería pertinente.

```{r}
library(psych)
mat_cor <- cor(dataset)
cortest.bartlett(unpivoted_dataset[,2:ncol(unpivoted_dataset)], n = 10)
```
Observamos que el p-valor de la prueba de esfericidad de Bartlett es igual a 0.07, el cual se encuentra ligeramente por encima del nivel de significancia de 0.05. Esto nos indica que el conjunto de datos no es adecuado para aplicar Análisis Factorial, esto se debe a que no podemos asumir que los datos cumplan estrictamente con la esfericidad necesaria para llevar a cabo un análisis factorial. Lo que significa que la correlación entre las variables podría no ser constante a lo largo de todas las dimensiones.

```{r}
KMO(mat_cor)
```

A partir del KMO, se refuerza el resultado observado en el test de Bartlett de que el conjunto de datos no es adecuado para aplicar Análisis Factorial.

Por lo tanto, de ahora en adelante es importante tener prudencia al momento de interpretar y utilizar los resultados.

*Análisis Factorial*:

***Modelo con un solo factor:***
```{r}
fa <- psych::fa(scale(dataset), nfactors=1, rotate ='varimax', fm ='ml', max.iter=1000)
fa
```
(14)*Estima la matriz de variancias y covariancias que surge del modelo factorial, la cual se descompone en variancia común (debido a los factores) y variancia específica. Compara estos resultados con respecto a la matriz S de variancias y covariancias muestral.*
```{r}
# Calcular la matriz de covarianzas
cov_matrix1 <- cov(dataset)

# Calcular la matriz de correlaciones
cor_matrix1 <- cor(dataset)

# Crear la matriz diagonal SD
SD <- diag(sqrt(diag(cov_matrix1)))

# Matriz L
L <- fa$loadings

# Calcular (L * L')
LL <- L %*% t(L)

# Calcular la matriz Psi (varianzas únicas)
Psi <- diag(diag(cor_matrix1) - diag(LL))

# Matriz de correlaciones aproximadas
R_Aprox1 <- LL + Psi

# Matriz de varianzas aproximadas
S_Aprox1 <- SD %*% R_Aprox1 %*% t(SD)

# Crear un vector de nombres de variables
nombres_var <- colnames(dataset)

# Asignar los nombres de las variables a la matriz
rownames(S_Aprox1) <- nombres_var
colnames(S_Aprox1) <- nombres_var

# Mostrar la matriz S_Aprox1
#print(S_Aprox1)
matrix_data <- round(S_Aprox1, 3)
cov_plot(matrix_data, 1.2, 0.7, "Matriz de Varianzas y Covarianzas Aproximada")

# Mostrar matriz S muestral
cov_matrix1 <- cov(dataset)
matrix_data <- round(cov_matrix1, 3)
cov_plot(matrix_data, 1.2, 0.7, "Matriz de Varianzas y Covarianzas Muestral")
#print(cov(dataset))
```
Si observamos la diagonal principal de la matriz de varianzas y covarianzas aproximada podremos notar que es igual a la diagonal de la matriz de varianzas y covarianzas muestrales, esto es algo nomral y esperado.

El análisis factorial no afecta a las varianzas individuales de las variables originales, es por ello que tiene sentido que las diagonales de ambas matrices sean similares.

A continuación se presenta la matriz cuyos elementos son los resultados de restar:

Matriz de Variancias y Covariancias muestrales - Matriz de Variancias y Covariancias aproximadas
```{r}
residual_cov <- cov(dataset)-S_Aprox1
matrix_data <- round(residual_cov, 3)
cov_plot(matrix_data, 1.2, 0.7, "Diferencia entre Matriz de Covarianzas Muestral y Aprox.")
```
De acuerdo a lo mencionado anteriormente, es de esperar que la diagonal principal sea 0 para todos los valores, ya que ambas matrices comparten los mismos valores.

Observando las covarianzas, podemos ver que el modelo ha logrado capturar muy bien algunas covarianzas, lo cual podemos comprobar en los casos donde la diferencia es cercana a cero, por ejemplo:

- EN.CO2.TRAN.ZS y EG.CFT.ACCS.ZS: 0.098

- EG.USE.COMM.CL.ZS y AG.LND.AGRI.ZS: -0.382

- EN.CO2.TRAN.ZS y EG.USE.COMM.CL.ZS: 0.057

Por otro lado, hay covarianzas muy elevadas en otros casos, lo que nos da un indicio de que el modelo no las ha capturado muy bien:

- AG.LND.AGRI.ZS y EN.ATM.METH.AG.ZS: 215.765

- EN.CO2.MANF.ZS y AG.LND.AGRI.ZS: -79.354

(15)*Repite los dos puntos anteriores agregando un nuevo factor al modelo. Se producen cambios en las comunalidades de las variables con respecto al modelo estimado en el punto 13? Determina cuál de las dos especificaciones del modelo factorial (m=1 ó m =2) resulta más adecuada para representar la estructura de asociación entre las variables. Es posible que en función de las variables y/o países elegidos en tu matriz de datos la especificación del modelo factorial con dos factores no sea la óptima. En ese caso agrega los comentarios que creas conveniente. *

***Modelo con dos factores:***
```{r}
fa2 <- psych::fa(scale(dataset), nfactors=2, rotate ='varimax', fm ='ml', max.iter=1000)
fa2
```
**Estimación de matriz de varianzas y covarianzas del modelo factorial**
```{r}
# Matriz L
L <- fa2$loadings

# Calcular (L * L')
LL <- L %*% t(L)

# Calcular la matriz Psi (varianzas únicas)
Psi <- diag(diag(cor_matrix) - diag(LL))

# Matriz de correlaciones aproximadas
R_Aprox2 <- LL + Psi

# Matriz de varianzas aproximadas
S_Aprox2 <- SD %*% R_Aprox2 %*% SD

# Crear un vector de nombres de variables
nombres_var <- colnames(dataset)

# Asignar los nombres de las variables a la matriz
rownames(S_Aprox2) <- nombres_var
colnames(S_Aprox2) <- nombres_var

# Mostrar la matriz S_Aprox2
#print(S_Aprox2)
matrix_data <- round(S_Aprox2, 3)
cov_plot(matrix_data, 1.2, 0.7, "Matriz de Varianzas y Covarianzas Aproximada")

# Mostrar matriz S muestral
cov_matrix1 <- cov(dataset)
matrix_data <- round(cov_matrix1, 3)
cov_plot(matrix_data, 1.2, 0.7, "Matriz de Varianzas y Covarianzas Muestral")
```

Restamos: Matriz de Variancias y Covariancias muestrales - Matriz de Variancias y Covariancias aproximadas
```{r}
residual_cov <- cov(dataset)-S_Aprox2
matrix_data <- round(residual_cov, 3)
cov_plot(matrix_data, 1.2, 0.7, 
         "Diferencia entre Matriz de Covarianzas Muestral y Aprox.")
```
Podemos analizar la matriz de diferencias para comparar que ocurre con los mismos pares de variables que analizamos en el caso del modelo de un solo factor:

- EN.CO2.TRAN.ZS y EG.CFT.ACCS.ZS: 0.18 (0.098 en modelo de 1 factor)

- EN.CO2.TRAN.ZS y EG.USE.COMM.CL.ZS: 0.119 (0.057 en modelo de 1 factor)

- EG.USE.COMM.CL.ZS y AG.LND.AGRI.ZS: -0.154 (-0.382 en modelo de 1 factor)

Si bien es cierto que en los dos primeros casos la covarianza capturada por el analisis factorial de 2 factores es menor que la capturada por el de 1 solo factor, la diferencia es muy leve. Además, podemos ver que ha mejorado en lo que respecta al tercer par de variables.

Sin embargo, la diferencia más notoria aparece cuando analizamos aquellas covarianzas que en el modelo de 1 factor arrojaban una diferencia muy grande entre lo muestral y lo aproximado:

- AG.LND.AGRI.ZS y EN.ATM.METH.AG.ZS: 0.08 (215.765 en modelo de 1 factor)

- EN.CO2.MANF.ZS y AG.LND.AGRI.ZS: -0.039 (-79.354 en modelo de 1 factor)

Observando estos resultados (y los demás que pueden verse en la comparación de matrices), tenemos evidencia para suponer que el modelo de dos factores es superior al modelo de un solo factor ya que captura mejor la covarianza de los datos.

**Comparación de comunalidades**
```{r}
sort(fa$communality,decreasing = T)->c1
sort(fa2$communality,decreasing = T)->c2
result_matrix <- cbind(c1, c2, "c2 - c1" = c2 - c1)
#result_matrix
suppressWarnings(knitr::kable(result_matrix, caption = "Comunalidades", format ="latex", 
             longtable=TRUE)%>%
  kableExtra::kable_styling(latex_options = "scale_down"))
```
Las comunalidades representan la proporción de varianza en una variable que puede atribuirse a los factores subyacentes en un modelo de análisis factorial. En general, al comparar dos modelos basados en las comunalidades, se busca determinar cuál de los modelos tiene comunalidades más altas. Comunalidades más altas indican que el modelo es más eficaz en capturar la estructura subyacente de los datos.

En nuestro análisis de comparación, hemos calculado la diferencia entre las comunalidades del segundo modelo y las comunalidades del primer modelo y la hemos mostrado en la tercera columna. En todos los casos, esta diferencia es positiva, lo que sugiere que el segundo modelo es más efectivo para capturar la estructura subyacente de los datos.

Por lo tanto, existen indicios para considerar que el modelo con dos factores es superior al modelo con un solo factor.

**Comparación de las unicidades**
```{r}
sort(fa$uniquenesses,decreasing = T)->u1
sort(fa2$uniquenesses,decreasing = T)->u2
result_matrix <- cbind(u1, u2, "u2 - u1" = u2 - u1)
#result_matrix
suppressWarnings(knitr::kable(result_matrix, caption = "Unicidades", format ="latex", 
             longtable=TRUE)%>%
  kableExtra::kable_styling(latex_options = "scale_down"))
```
Las unicidades se pueden considerar como el complemento de las comunalidades, es decir, la varianza que no es explicada por los factores subyacentes. En este contexto, es deseable que las unicidades de las variables sean lo más bajas posible.

Si observamos la tercera columna, que representa la diferencia entre las unicidades del segundo modelo y las unicidades del primer modelo, podemos notar que todos los valores son negativos. Esto tiene sentido, ya que si las comunalidades aumentaron en el modelo 2, las unicidades correspondientes deberían disminuir.

Esta revisión mantiene el contenido original pero busca mejorar la claridad y la coherencia en la exposición de ideas.

**Computamos criterios de información: Criterio de Schwarz (BIC) y Criterio de Akaike (AIC)**

Los criterios de información BIC y AIC son indicadores de ajuste del modelo, en ambos criterios siempre preferiremos el modelo que tenga un valor más bajo en ambos estadísticos.
```{r}
#Comparación de modelos
m_max <- ncol(dataset)

lista_resultados <- list()
for(j in 1:2){
  lista_resultados[[j]] <- fa_ftest(dataset,j,0.05)
}
tabla_ftest <- data.table::rbindlist(lista_resultados)

#tabla_ftest
suppressWarnings(knitr::kable(tabla_ftest, caption = "Resultados FTest", format ="latex", 
             longtable=TRUE)%>%
  kableExtra::kable_styling(latex_options = "scale_down"))
```
Podemos ver que:

- BIC_Modelo2 = -7.98 < BIC_Modelo1 = -4.77

- AIC_Modelo2 = -4.04 < AIC_Modelo1 = 1.27

Por lo tanto, siguiendo estos criterios, el Modelo 2 es preferible por sobre el Modelo 1, ya que tiene valores más bajos tanto de AIC como de BIC, lo que sugiere que se ajusta mejor a los datos y es potencialmente un modelo más parsimonioso.

A raíz de todas las pruebas implementadas, podemos considerar que el modelo de Análisis Factorial que incluye 2 factores es superior al modelo con un solo factor.

(16)*Interpreta brevemente los resultados obtenidos, tomando en cuenta las eventuales diferencias y semejanzas con los resultados del cálculo de componentes principales.*

Cuando trabajamos con modelos de Análisis Factorial, normalmente las matrices que obtenemos no son fáciles de interpretar. Por lo tanto, una práctica usual es transformar la matriz para simplificar la estructura de los factores, principalmente la rotación de los factores.

Cuando aplicamos una rotación a la matriz las comunalidades no varían, pero si se ve afectada la varianza explicada por cada factor.

Es por ello que al momento de correr el modelo se ha indicado que se debe hacer aplicando una rotación "Varimax", la cual implica una rotación ortogonal que minimiza el número de variables que tienen saturaciones altas en cada factor y simplifica la interpretación de los factores.

A continuación se presentan los diagramas de árbol que facilitan la interpretación del modelo de análisis factorial. Estos diagramas muestran la relación existente entre los factores y las variables observadas, también conocida como carga factorial, la cual se representa por flechas que salen del factor y llegan a las distintas variables.

Es común visualizar distintos formatos de flechas que permiten denotar el tipo de relación. De este modo, en los diagramas que se muestran a continuacíon se observan dos tipos de flechas:

- Flechas sólidas de color negro: indican una carga factorial positiva. Es decir, implican que cuando el factor aumenta, también lo hace la variable.

- Flechas punteadas de color rojo: al contrario del caso anterior, indican una carga factorial negativa. Por lo tanto, cuando el factor aumenta, la variable disminuye.

Como criterio para definir si una carga factorial es fuerte o no podemos considerar como "pesadas" a las cargas factoriales mayores a 0.5 y como "leves" a las que se encuentran por debajo de este valor.
```{r}
# Configuración del gráfico
layout(matrix(1:2, nrow = 1))
aspect_ratio <- par('pin')[1] / par('pin')[2]
scale_factor <- aspect_ratio/1.2
par(cex = scale_factor)

# Dibujamos los diagramos
fa.diagram(fa, rsize = 0.4, l.cex=1.2, main = "")
text(x = -2.5, 
     y = 9, 
     "Analisis Factorial - 1 Factor", 
     adj = c(0, 0.5), 
     font = 2,
     cex = 1.2)

fa.diagram(fa2, rsize = 0.4, l.cex=1.2, main = "")
text(x = -2.8,
     y = 9, 
     "Analisis Factorial - 2 Factores", 
     adj = c(0, 0.5), 
     font = 2,
     cex = 1.2)

# Restaura el tamaño de la fuente a su valor original
par(cex = 1)
```
En el diagrama de árbol del análisis factorial ajustado con un solo factor podemos ver que el único factor generado tiene cargas factoriales pesadas con 5 variables:

- EN.CO2.TRAN.ZS: -1 (Carga factorial negativa)
- EN.CO2.ETOT.ZS: 0.9 (Carga factorial)
- EG.CFT.ACCS.ZS: 0.8 (Carga factorial positiva)
- EN.CO2.MANF.ZS: 0.6 (Carga factorial positiva)
- EG.USE.COMM.CL.ZS: -0.5 (Carga factorial negativa)

En el caso del análisis factorial ajustado con dos factores podemos encontrar que el factor 1 presenta cargas factoriales fuertes con 4 variables:

- EN.CO2.TRAN.ZS: -1 (Carga factorial negativa)
- EN.CO2.ETOT.ZS: 0.9 (Carga factorial)
- EG.CFT.ACCS.ZS: 0.8 (Carga factorial positiva)
- EG.USE.COMM.CL.ZS: -0.5 (Carga factorial negativa)

Podemos observar que la carga factorial del factor 1 con estas 4 variables es exactamente la misma que presenta en el modelo de un solo factor, la diferencia es que ahora la variable EN.CO2.MANF.ZS se encuentra vinculada con el factor 2, siendo también 4 las variables asociadas con el factor 2:
- AG.LND.AGRI.ZS: 1 (Carga factorial positiva)
- EN.ATM.METH.AG.ZS: 0.7 (Carga factorial positiva)
- EN.ATM.NOXE.AG.ZS: 0.6 (Carga factorial positiva)
- EN.CO2.MANF.ZS: -0.6 (Carga factorial negativa)

Si comparamos los resultados del modelo de análisis factorial con dos factores con los obtenidos en el Análisis de Componentes Principales, podemos ver que el Factor 1 presenta asociación con las mismas tres variables que se encontraban relacionadas con el Componente 1 del PCA: EN.CO2.TRAN.ZS, EN.CO2.ETOT.ZS y EG.CFT.ACCS.ZS.

Del mismo modo, podemos ver que el Factor 2 se encuentra asociado con las mismas tres variables que se relacionan con el Componente 2 del PCA: AG.LND.AGRI.ZS, EN.ATM.METH.AG.ZS y EN.ATM.NOXE.AG.ZS.

A diferencia del PCA, el Análisis Factorial ha incluído una cuarta variable asociada a cada factor, EG.USE.COMM.CL.ZS en el caso del Factor 1 y EN.CO2.MANF.ZS en el caso del Factor 2.

Por lo tanto, podríamos intuir que el Componente 1 representa las emisiones de gases de efecto invernadero relacionadas al consumo de energía en aspectos del día a día, como ser transporte, cocina, calefacción.

Por lo tanto, podríamos intuir que los factores latentes son:

- ML1: Emisiones de gases de efecto invernadero derivados del consumo de energía diaria (transporte, comida, calefacción)

- ML2: Emisiones de gases de efecto invernadero derivados de la actividad agrícola.

Si bien es cierto que los resultados obtenidos por el análisis factorial son similares a los obtenidos por el PCA, es importante resaltar que estos resultados deben ser considerados con cautela debido a que el conjunto de datos no pasó satisfactoriamente la prueba de Bartlett ni el KMO, lo cual nos indica que el análisis factorial puede no ser la técnica más adecuada para estos datos.