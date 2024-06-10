if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GEOquery")
BiocManager::install("affy")
BiocManager::install("Biobase")
BiocManager::install("genefilter")
BiocManager::install("hgu133plus2.db") # Affymetrix HG-U133_Plus_2 Array annotation data (chip hgu133plus2) https://www.bioconductor.org/packages/release/BiocViews.html#___AnnotationData
BiocManager::install("hgu133plus2cdf") # Affymetrix HG-U133_Plus_2 Array annotation data (chip hgu133plus2)

install.packages("factoextra")
install.packages("compareGroups")
install.packages("janitor")
install.packages("caret")
install.packages("compareGroups")
install.packages("scales")
install.packages("ISLR", dependencies = TRUE)
install.packages("e1071", dependencies = TRUE)

## Carga de paquetes al entorno de trabajo:
### Básicos
library(GEOquery) # Para descargar datos de GEO
library(affy) # Para trabajar con Arrays de Affymetrix
library(ggplot2)
library(reshape2)
library(factoextra)
library(preprocessCore)
library(genefilter)
library(annotate)
library(caret)
library(compareGroups)
library(tidyverse)
library(compareGroups)
library(glmnet)
library(scales)
library(ISLR)
library(e1071)


###Anotación
# Anotación [HG-U133_Plus_2] Affymetrix Human Genome U133 Plus 2.0 Array; GEO Platform: GPL570 	
# Bioconductor Annotación
library(hgu133plus2.db)
library(hgu133plus2cdf)


## Descarga matriz de expresión 

Ya pre-procesada.
Filas: genes
Columnas: Muestras
Valores: expresión génica

rel_muestras <- GEOquery::getGEO("GSE42568", GSEMatrix =TRUE, getGPL=TRUE) 
if (length(rel_muestras) > 1) idx <- grep("GPL570", attr(rel_muestras, "names")) else idx <- 1
eset_muestras <- rel_muestras[[idx]]

## Acceso información biologica
# Información muestral (fenotipos), anotada en el objeto ExpressionSet:

pheno_muestras <- pData(eset_muestras)

## Información tecnica de las sondas 

### Correspondiencia sondas y genes

SondasMetadata <- fData(eset_muestras)

## Descarga datos crudos (.CEL)

**Descarga de los archivos**
options(timeout = max(300, getOption("timeout")))
setwd("...")
gcel=getGEOSuppFiles("GSE42568")


##Cambio directorio
setwd("GSE42568")
##Descomprimir archivos
system(" tar xvf GSE42568_RAW.tar")

**Lectura de los ficheros y almacenamiento en un objeto AffyBatch**
GSE42568  = ReadAffy()

# Pre-procesado de los datos
# Normalización método MAS5 (de la casa comercial de Affymetrix)
## Lectura de intensidades normalizadas y asignación de Ausencia/Presencia

eset_mas5 <- mas5(GSE42568)
call_mas5 <- mas5calls(GSE42568)

## Filtrado de genes

eset.filt1 = nsFilter(eset_mas5,var.func=IQR,var.cutoff=0.5,require.GOBP=TRUE)
eset.filt2 = nsFilter(eset_mas5,var.func=median,var.cutoff=0.5,require.GOBP=TRUE)
sel = intersect(featureNames(eset.filt1),featureNames(eset.filt2))
eset_mas5_filt = eset_mas5[sel,]
call_mas5_filt = call_mas5[sel,]

RawValues <- as.data.frame(list(exprs(eset_mas5_filt)))
RAw_presence <- as.data.frame(list(exprs(call_mas5_filt)))

RawTable <- cbind(RawValues,RAw_presence)



##Gráfico mosaico de Presencias/Ausencias**
  
#Gráfico mosaico
grade <- pheno_muestras$characteristics_ch1.4
#Trasponemos los datos
presence_t = as.data.frame(t(RAw_presence))
#Añadimos el grado a los datos
presence_plus_grade = cbind(presence_t, grade)
#Separamos los datos por clase
split_presence = split(presence_plus_grade, presence_plus_grade$grade)

#Tablas de contingencia detección sonda/grado en frecuencia absoluta
presence_table <- data.frame(rbind(table(unlist(split_presence[[1]]))[c(1,3,4)],
                                   table(unlist(split_presence[[2]]))[c(1,3,4)],
                                   table(unlist(split_presence[[3]]))[c(1,3,4)],
                                   table(unlist(split_presence[[4]]))[c(1,3,4)]))
rownames(presence_table) <- c("grade: 1", "grade: 2", "grade: 3", "grade: NA")

mosaicplot(presence_table, main=NULL, color=c("#E1F7E6","#2ACCCB","#02547d"),
           xlab = "Grade", ylab="Detección de sonda",
           border = c("#E1F7E6","#2ACCCB","#02547d"), cex.axis = 1,
           las=par(cex.lab=1.2))

## Eliminacion de sondas presentes de expresión ausente

# Criterio: presencia > 10 muestras
presence_filter = RawValues[which(rowSums(RawTable=="P")>10),]

## Normalización por quantil de los datos
raw_matrix <- as.matrix(presence_filter)
normalized_data<-normalize.quantiles(raw_matrix)
normalized_data<-as.data.frame(normalized_data)
names(normalized_data)<-names(presence_filter)
rownames(normalized_data)<-rownames(presence_filter)

##Anotación del expression set

ID = rownames(normalized_data)
Genes <- getSYMBOL(ID,"hgu133plus2.db")
normalized_data <- cbind(normalized_data, Genes)
normalized_data <-aggregate(normalized_data, by = list(normalized_data$Genes), FUN = median)
rownames(normalized_data) <- normalized_data$Group.1
normalized_data <- normalized_data[,-1]
normalized_data <- normalized_data[,-ncol(normalized_data)]


##Transformación logaritmica de los datos normalizados
norm_log_data<-log2(normalized_data)


#Análisis descriptivo de los datos

## Boxplot distribución de los datos normalizados vs no normalizados**
par(mfrow=c(1,2))
plotDensity(raw_matrix, main="Datos sin normalizar"); grid()
plotDensity(norm_log_data, main="Datos normalizados"); grid()

set.seed(123)
x<-sample(ncol(RawValues),15)

ggplot(data=melt(log2(RawValues[,x])), aes(x=variable, y=value)) +
  geom_boxplot(fill="#2ACCCB", alpha=0.7) +
  theme_bw() + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=10, angle=90, hjust=1),
        plot.title = element_text(hjust = 0.5, size=15, margin = ggplot2::margin(b=12)))+
  ggtitle("Datos sin normalizar"); grid()

ggplot(data=melt(log2(norm_log_data[,x])), aes(x=variable, y=value)) +
  geom_boxplot(fill="#2ACCCB", alpha=0.7) +
  theme_bw() + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=10, angle=90, hjust=1),
        plot.title = element_text(hjust = 0.5, size=15, margin = ggplot2::margin(b=12)))+
  ggtitle("Datos normalizados"); grid()


# Análisis de componentes principales

#Anotamos informacion muestral
#PCA datos normalizados y transformados logaritmicamente

data_t <- as.data.frame(t(norm_log_data))
data_plus_grade <- cbind(data_t, grade)
pca_analysis<-prcomp(data_plus_grade[,-ncol(data_plus_grade)])


## Porcentajes de varianza explicada por cada PC

screeplot_1 <- fviz_eig(pca_analysis, addlabels = TRUE, ylim = c(0, 25), barfill = '#2AcccB', barcolor="black")

ggpubr::ggpar(screeplot_1, title = "", ggtheme = theme_bw(), 
              ylab = "Porcentaje de la varianza explicado",
              xlab = "Componente principal",
              font.x = 15,
              font.y = 15,
              font.tickslab = 12)

pca_plot1 <- fviz_pca_ind(pca_analysis, axes = c(1, 2), geom.ind = "point", col.ind = "black",
                          palette = c("#8B1C42","#0A4A75","#36D6DA","#F7B409","#F03B0E"),
                          pointshape = 23,
                          pointsize = 2.5,
                          fill.ind = data_plus_grade$grade,
                          addEllipses = TRUE,
                          mean.point = FALSE,
                          alpha.ind = 0.7)

ggpubr::ggpar(pca_plot1, title = "", legend.title = "Grade", legend = c(0.90, 0.3),
              ggtheme = theme_bw(), font.tickslab = 10, xlab = "PC1", ylab = "PC2", 
              font.x = 16, font.y = 16, font.legend = 14)

#Contribucion de cada gen a la PC1
fviz_contrib(pca_analysis, choice="var", axes = 1, top = 20) 
#Contribucion de cada gen a la PC2
fviz_contrib(pca_analysis, choice="var", axes = 2, top = 20)
```

#Repetimos el ejercicio con las componentes 2 y 3:
  

pca_plot2 <- fviz_pca_ind(pca_analysis, axes= c(2,3), geom.ind = "point", col.ind="black",
                          palette = c("#8B1C42","#0A4A75","#36D6DA","#F7B409","#F03B0E"),
                          pointshape = 23,
                          pointsize = 2.5,
                          fill.ind = data_plus_grade$grade,
                          addEllipses = TRUE,
                          mean.point = FALSE,
                          alpha.ind = 0.7)

ggpubr::ggpar(pca_plot2, title = "", legend.title = "Grade", legend = c(0.90, 0.3),
              ggtheme = theme_bw(), font.tickslab = 10, xlab = "PC2", ylab = "PC3", 
              font.x = 16, font.y = 16, font.legend = 14)

pca_plot2 <- fviz_pca_ind(pca_analysis, axes= c(1,3), geom.ind = "point", col.ind="black",
                          palette = c("#8B1C42","#0A4A75","#36D6DA","#F7B409","#F03B0E"),
                          pointshape = 23,
                          pointsize = 2.5,
                          fill.ind = data_plus_grade$grade,
                          addEllipses = TRUE,
                          mean.point = FALSE,
                          alpha.ind = 0.7)

ggpubr::ggpar(pca_plot2, title = "", legend.title = "Grade", legend = c(0.90, 0.3),
              ggtheme = theme_bw(), font.tickslab = 10, xlab = "PC1", ylab = "PC3", 
              font.x = 16, font.y = 16, font.legend = 14)

pca_plot2 <- fviz_pca_ind(pca_analysis, axes= c(1,4), geom.ind = "point", col.ind="black",
                          palette = c("#8B1C42","#0A4A75","#36D6DA","#F7B409","#F03B0E"),
                          pointshape = 23,
                          pointsize = 2.5,
                          fill.ind = data_plus_grade$grade,
                          addEllipses = TRUE,
                          mean.point = FALSE,
                          alpha.ind = 0.7)

ggpubr::ggpar(pca_plot2, title = "", legend.title = "Grade", legend = c(0.90, 0.3),
              ggtheme = theme_bw(), font.tickslab = 10, xlab = "PC1", ylab = "PC4", 
              font.x = 16, font.y = 16, font.legend = 14)

#Clasificador multiclase basado en SVMs

data_plus_grade_new <- data_plus_grade

# Convertir la variable dependiente a un factor
data_plus_grade_new$grade <- factor(data_plus_grade_new$grade, 
                                    levels = c('grade: NA', 'grade: 1', 'grade: 2', 'grade: 3'))

# Calcular la matriz de correlación
# Excluir la variable dependiente
correlation_matrix <- cor(data_plus_grade_new[, -ncol(data_plus_grade_new)])  

# Encontrar las variables altamente correlacionadas
highlyCorrelated <- findCorrelation(correlation_matrix, cutoff = 0.9)

# Mostrar las variables que se eliminarán
print(names(data_plus_grade_new)[highlyCorrelated])  

# Eliminar las variables para eliminar la correlación por pares
df_reduced <- data_plus_grade_new[, -highlyCorrelated]

# Configurar el índice de partición
set.seed(123)
trainIndex <- createDataPartition(df_reduced$grade, p = .7, 
                                  list = FALSE, 
                                  times = 1)
data_train <- df_reduced[trainIndex,]
data_test <- df_reduced[-trainIndex,]

svm_cv <- tune("svm", grade ~ ., data = data_train, kernel="linear", 
               ranges = list(cost = c(0.0001, 0.0005, 0.001, 0.01, 0.1, 1)))


ggplot(data = svm_cv$performances, aes(x = cost, y = error)) +
  geom_line() +
  geom_point() +
  labs(title = "Error de clasificación vs hiperparámetro C") +
  theme_bw()


svm_cv$best.parameters

modelo_svm <- svm_cv$best.model

# Aciertos del modelo con los datos de entrenamiento
paste("Error de entrenamiento:", 100*mean(data_train$grade != modelo_svm$fitted), "%")

table(prediccion = modelo_svm$fitted, clase_real = data_train$grade)

#Evaluación con datos de entrenamiento
predicciones <- predict(object = modelo_svm, newdata = data_test)

paste("Error de test:", 100 * mean(data_test$grade != predicciones), "%")

table(prediccion = predicciones, clase_real = data_test$grade)

#Evaluacion sobre todo el conjunto
predicciones <- predict(object = modelo_svm, newdata = df_reduced)

paste("Error de test:", 100 * mean(df_reduced$grade != predicciones), "%")

table(prediccion = predicciones, clase_real = df_reduced$grade)

# Probamos a reducir a dis clases (si grado = 3 G3==1 si no G3==0)

data_plus_grade_new$G3<-ifelse(data_plus_grade_new$grade == "grade: 3", 1, 0)
data_plus_grade_new <- data_plus_grade_new[-(ncol(data_plus_grade_new)-1)]
# Convertir la variable dependiente a un factor
data_plus_grade_new$G3 <- factor(data_plus_grade_new$G3, 
                                 levels = c(0,1))


# Calcular la matriz de correlación
# Excluir la variable dependiente
correlation_matrix <- cor(data_plus_grade_new[, -ncol(data_plus_grade_new)])  

# Encontrar las variables altamente correlacionadas
highlyCorrelated <- findCorrelation(correlation_matrix, cutoff = 0.9)

# Mostrar las variables que se eliminarán
print(names(data_plus_grade_new)[highlyCorrelated])  

# Eliminar las variables altamente correlacionadas
df_reduced <- data_plus_grade_new[, -highlyCorrelated]

# Configurar el índice de partición
set.seed(123)
trainIndex <- createDataPartition(df_reduced$G3, p = .7, 
                                  list = FALSE, 
                                  times = 1)
data_train <- df_reduced[trainIndex,]
data_test <- df_reduced[-trainIndex,]

svm_cv2 <- tune("svm", G3 ~ ., data = data_train, kernel="linear", 
                ranges = list(cost = c(0.0001, 0.0005, 0.001, 0.01, 0.1, 1)))

ggplot(data = svm_cv2$performances, aes(x = cost, y = error)) +
  geom_line() +
  geom_point() +
  labs(title = "Error de clasificación vs hiperparámetro C") +
  theme_bw()


svm_cv2$best.parameters

modelo_svm2 <- svm_cv2$best.model

# Aciertos del modelo con los datos de entrenamiento
paste("Error de entrenamiento:", 100*mean(data_train$G3 != modelo_svm2$fitted), "%")

table(prediccion = modelo_svm2$fitted, clase_real = data_train$G3)

#Evaluación con datos de entrenamiento
predicciones <- predict(object = modelo_svm2, newdata = data_test)

paste("Error de test:", 100 * mean(data_test$G3 != predicciones), "%")

table(prediccion = predicciones, clase_real = data_test$G3)



# Pruebas del Modelo de regresion logistica
#Para trabajar con la funcion glmnet hay que convertir la variable objetivo a datos de tipo numerico.

data_train <- data_train %>%
  mutate(G3 = ifelse(G3 == "0", 0, 1))

data_test <- data_test %>%
  mutate(G3 = ifelse(G3 == "0", 0, 1))

modelo.glmnet <- glmnet(x = data_train[,-ncol(data_train)], y = data_train$G3, family = "binomial", alpha = 0, nlambda = 100, standardize = TRUE)

# Evolución de los coeficientes en función de lambda
# ==============================================================================
regularizacion <- modelo.glmnet$beta %>% 
  as.matrix() %>%
  t() %>% 
  as_tibble() %>%
  mutate(lambda = modelo.glmnet$lambda)

regularizacion <- regularizacion %>%
  pivot_longer(
    cols = !lambda, 
    names_to = "predictor",
    values_to = "coeficientes"
  )

regularizacion %>%
  ggplot(aes(x = lambda, y = coeficientes, color = predictor)) +
  geom_line() +
  scale_x_log10(
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x))
  ) +
  labs(title = "Coeficientes del modelo en función de la regularización") +
  theme_bw() +
  theme(legend.position = "none")

# Evolución del error en función de lambda
# ==============================================================================
set.seed(123)
cv_error <- cv.glmnet(
  x = as.matrix(data_train[,-ncol(data_train)]), 
  y = data_train$G3,
  alpha  = 0,
  nfolds = 10,
  type.measure = "mse",
  standardize  = TRUE
)

plot(cv_error)

# Mejor valor lambda encontrado
# ==============================================================================
paste("Mejor valor de lambda encontrado:", cv_error$lambda.min)

# Mejor valor lambda encontrado + 1sd
# ==============================================================================
# Mayor valor de lambda con el que el test-error no se aleja más de 1sd del mínimo.
paste("Mejor valor de lambda encontrado + 1 desviación estándar:", cv_error$lambda.1se)

# Mejor modelo lambda óptimo
# ==============================================================================
modelo <- glmnet(
  x = as.matrix(data_train[,-ncol(data_train)]), 
  y = data_train$G3,
  alpha       = 0,
  lambda      = cv_error$lambda.1se,
  standardize = TRUE
)

# Coeficientes del modelo
# ==============================================================================
df_coeficientes <- coef(modelo) %>%
  as.matrix() %>%
  as_tibble(rownames = "predictor") %>%
  rename(coeficiente = s0)

df_coeficientes %>%
  filter(predictor != "(Intercept)") %>%
  ggplot(aes(x = predictor, y = coeficiente)) +
  geom_col() +
  labs(title = "Coeficientes del modelo Ridge") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 6, angle = 45))

## Variables con mayor peso en el modelo Ridge
summary(abs(df_coeficientes$coeficiente))
df_coeficientes %>%
  mutate(abs.coef = abs(coeficiente)) %>% 
  filter(abs.coef >= 0.00045)


# Predicciones de entrenamiento
# ==============================================================================
predicciones_train <- predict(modelo, newx = as.matrix(data_train[-ncol(data_train)]))

# MSE de entrenamiento
# ==============================================================================
training_mse <- mean((predicciones_train - data_train$G3)^2)
paste("Error (mse) de entrenamiento:", training_mse)

# Predicciones de test
# ==============================================================================
predicciones_test <- predict(modelo, newx = as.matrix(data_test[,-ncol(data_test)]))

# MSE de test
# ==============================================================================
test_mse <- mean((predicciones_test - data_test$G3)^2)
paste("Error (mse) de test:", test_mse)

# Convertir probabilidades a clases (0.5 es el umbral)
pred_class <- ifelse(predicciones_test > 0.5, "S", "N")

# Crear una tabla de confusión
y_ref <- ifelse(data_test$G3 == 0, "N", "S")
confusionMatrix(data = factor(pred_class, levels = c("N", "S")), 
                reference = as.factor(y_ref))

