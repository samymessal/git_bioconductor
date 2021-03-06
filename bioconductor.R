if (!require("BiocManager"))
  +     install.packages("BiocManager")
BiocManager::install("Biobase")

#Creacion del expressionSet a partir de los ficheros del ordenador

#AssayData
dataDirectory <- system.file("extdata", package="Biobase") #obtenemos la ruta del fichero exdata
exprsFile <- file.path(dataDirectory, "exprsData.txt")  #Obtenemos el documento exprsData.txt del fichero exdata
exprs <- as.matrix(read.table(exprsFile, header = TRUE, 
                              sep="\t", row.names = 1, as.is=TRUE)) #Introducimos el contenido del documento
                                                                    #en una matriz
exprs
head(exprs)
class(exprs)
dim(exprs)

minimalSet <- ExpressionSet(assayData = exprs)  #Introducimos la matriz creada en la clase ExpressionSet
class(minimalSet)

#Phenotypic data

pDataFile <- file.path(dataDirectory, "pData.txt") #Obtencion del documento pData
pData <- read.table(pDataFile, header=TRUE, row.names = 1, sep ="\t") #Introducimos en un dataframe el documento

class(pData)
names(pData)
dim(pData)  #Notese que el numero de lineas de pData corresponde al numero de columnas del assayData

summary(pData)
pData
sapply(pData, class)

pData[15, 1:2]
pData[20, 1:2]

pData[pData$score>0.8,]

# A veces el nombre de las columnas puede ser un codigo dificil de entender para el usuario,
# por ello, podemos encontrar mas informacion en un fichero metadata o crear uno a partir de la informacion 
# encontrada en internet.

metadata <- data.frame(labelDescription = 
                        c("Patient gender", "Case/Control", "Tumor progress on XYZ scale"),
                      row.names = c("gender", "type", "score"))
metadata

#AnnotatedDataFrame Annotations

phenoData <- new("AnnotatedDataFrame", data=pData, varMetadata=metadata)
phenoData
phenoData[,1]

#Estas son tres de las funciones que vienen con el objeto de ExpressionSet "AnnotatedDataFrame"
pData(phenoData)
sampleNames(phenoData)
varMetadata(phenoData)

#Anadimos la informacion sobre el chib que estamos utilizando

annotation <- "hgu95av2"
#Experiment description MIAME

#Creamos un objeto de ExpressionSet "MIAME"
experimentData <- new("MIAME",
                      name="Samy Messal",
                      lab = "The zwin lab",
                      contact = "samy.messal.sm@gmail.com",
                      title = "Smoking-Cancer Experiment",
                      abstract = "An example ExpressionSet",
                      url = "www.lab.zwin.exist",
                      other = list(notes="Created from text files"))

#Ensamblaje del ExpressionSet

exampleSet <- ExpressionSet(assayData = exprs, phenoData = phenoData, 
                            experimentData=experimentData, annotation = annotation)
exampleSet

help("ExpressionSet-class")  #Obtendremos ayuda 

featureNames(exampleSet)
varLabels(exampleSet)
annotation(exampleSet)

#subset of the exampleSet

mok <- exampleSet[1:8, 1:5]
mok
featureNames(mok)
dim(mok)
