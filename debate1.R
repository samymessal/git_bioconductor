#Creamos la primera version del documento R debate1

#Ejercicio 1

#Create an ExpressionSet from scratch using, for example, information from a study downladed from GEO.
#Utilizaré el estudio que he presentado en el formulario.

#Para obtener el ExpressionSet a partir de los datos proporcionados en GEO entry: 
#GSE65480 en el formato .CEL
#Utilizaré otros paquetes de Bioconductor: affy, lima, hgu95a.db, annotate.
#Siguiendo este tutorial http://rstudio-pubs-static.s3.amazonaws.com/16793_fba4b435fe4e4c17bf6a13b5d8d05eec.html

BiocManager::install(c("affy","lima", "hgu95a.db", "annotate"))

require(affy)
require(lima)
require(hgu95a.db)
require(annotate)
require(Biobase)
cel_file <- ReadAffy(celfile.path = "./GSE65480_RAW", compress = TRUE)

cel_file
colnames(cel_file)
class(cel_file)
dim(cel_file)

str(cel_file)

annotation(cel_file)

pData(cel_file) #Podemos ver que no tenemos informacion sobre los pacientes

varMetadata(cel_file) 

featureData(cel_file)

#assay data
assayData_1 <- assayData(cel_file)
#Pheno data

pheno_1 <- phenoData(cel_file)


#crearé un dataframe para poder crear otra clase "AnnotatedDataFrame".

pheno2 <- data.frame(sampleNames(assayData_1))
colnames(pheno2) <- "sample"

meta <- varMetadata(cel_file)

class(meta)
dim(meta)
phenoData_creada <- new("AnnotatedDataFrame", data=pheno2, varMetadata = meta)

#Annotation

annotation_1 <- annotation(cel_file)

#Experiment data

experimentData(cel_file) #El fichero no contiene informacion sore el experimento, buscaremos en internet

experimentData_1 <- new("MIAME",
                      name="Takaaki Kobayashi",
                      lab = "Kyorin university",
                      contact = "	ck9t-kbys@asahi-net.or.jp",
                      title = "	Expression data at each site in colon cancer",
                      )

#Por ultimo ensamblamos todos los objetos creados anteriormente en la clase ExpressionSet.

ExpressionSet_colon <- ExpressionSet(assayData=assayData_1, phenoData = pheno_1, 
                                     annotation = annotation_1, experimentData = experimentData_1)

#Cuando intentamos crear el ExpressionSet con el la instancia de "AnottedDataFrame" creada
#manualmente, nos sale un mensaje de error debido a que los nombres de las muestras difieren

sampleNames(assayData_1) == sampleNames(phenoData_creada)
sampleNames(assayData_1)[1]

#Find out how you can access and change ...
# the expression values or the covariates in the phenoData

dim(ExpressionSet_colon)
dim(cel_file)
class(cel_file)
class(ExpressionSet_colon)
