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

cel_file <- ReadAffy(celfile.path = "./GSE65480_RAW", compress = TRUE)

cel_file

class(cel_file)
dim(cel_file)

str(cel_file)

annotation(cel_file)

pData(cel_file) #Podemos ver que no tenemos informacion sobre los pacientes

varMetadata(cel_file) 

featureData(cel_file)

#assay data
assayData <- assayData(cel_file)
#Pheno data

pheno1 <- phenoData(cel_file)
pheno1

#crearé un dataframe para poder crear otra clase "AnnotatedDataFrame".

pheno2 <- data.frame(sampleNames(assayData(cel_file)))
colnames(pheno2) <- "sampleID"
pheno2
meta <- varMetadata(cel_file)
meta
class(meta)
dim(meta)
sampleNames(assayData(cel_file)) == sampleNames(phenoData(cel_file))
phenoData <- new("AnnotatedDataFrame", data=pheno2, varMetadata = meta)

#Annotation

annotation <- annotation(cel_file)

#Experiment data

experimentData(cel_file) #El fichero no contiene informacion sore el experimento, buscaremos en internet

experimentData <- new("MIAME",
                      name="Takaaki Kobayashi",
                      lab = "Kyorin university",
                      contact = "	ck9t-kbys@asahi-net.or.jp",
                      title = "	Expression data at each site in colon cancer",
                      )

#Por ultimo ensamblamos todos los objetos creados anteriormente en la clase ExpressionSet.

ExpressionSet_colon <- ExpressionSet(assayData=assayData, phenoData = phenoData, 
                                     annotation = annotation, experimentData = experimentData)





