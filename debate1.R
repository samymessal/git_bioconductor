#Creamos la primera version del documento R debate1

#Ejercicio 1

#Create an ExpressionSet from scratch using, for example, information from a study downladed from GEO.
#Utilizaré el estudio que he presentado en el formulario.

#Para obtener el ExpressionSet a partir de los datos proporcionados en GEO entry, descargaremos
#el series_matrix.txt. El identificador del estudio es GSE65480.

datadir <- "."  #El proyecto esta creado en este directorio
info <- readLines(file.path(datadir, "GSE65480_series_matrix.txt"), n=62)
# A partir de comandos del terminal, utilizando grep, sabemos que el !series_matrix_begin
#se encuentra en la linea 63.

rows2read <- 33362 - 63 - 2 #33362 son el numero de lineas totales, 63 las de la info y
                            #2 del titulo y final.

matriz_GSE <- read.table(file.path(datadir, "GSE65480_series_matrix.txt"), skip=63, 
                         header = TRUE, sep="\t", row.names=1, nrows=rows2read)

dim(matriz_GSE)

#Para facilitar la visualizacion de los datos, cambiaremos el nombre de las columnas
#Las separaremos en central y frontal, que son los sitios de donde se han extraido las muestras.
#A partir de la informacion obtenida en GEO, podemos ver que las muestras centrales corresponden
# a identificadores pares.
copia_matriz <- matriz_GSE
count <- 1
count_name <- 1

for (i in colnames(copia_matriz)){
  if (as.integer(substr(i, nchar(i)-1, nchar(i))) %% 2 == 0){
    colnames(copia_matriz)[count] <- paste("Central",count_name,sep="_")
  }
  else {
    colnames(copia_matriz)[count]<- paste("Frontal", count_name, sep="_")
    count_name <- count_name + 1
  
  }
  count <- count + 1
}
colnames(copia_matriz)

maximos <- apply(copia_matriz, 2, max)
minimos <- apply(copia_matriz, 2, min)

par(mfrow=c(2,1))
plot(maximos, type = "l", xaxt="n", xlab="", col="blue")
axis(1, at = seq(1, 40), las=2, labels = colnames(copia_matriz) )

plot(minimos, type= "l", col="red", xaxt="n", xlab="")
axis(1, at = seq(1,40), las=2, labels=colnames(copia_matriz))

par(mfrow=c(1,1))
boxplot(copia_matriz, las=2)

cluster <- hclust(dist(t(copia_matriz)), method = "ward.D2")
plot(cluster, hang=-1) #Los grupos central y frontal parecen agrupados entre si
                        #No parece haber clara diferencia en la expresion
                        #entre las muestras centrales y frontales


#Hemos analizado la matriz de expresion, ahora pasamos a obtener el phenoData
#La informacion normalmente contenida en targets se encuentra en el series matrix
#despues de la linea !Sample_. A partir del terminal buscaremos la posicion de esta.

rows2read_sup <- 1
sample_suplementary_file <- read.table(file.path(datadir, "GSE65480_series_matrix.txt"), skip = 60, 
                                       sep="\t", row.names = 1, nrows = rows2read_sup)

dim(sample_suplementary_file) #tenemos que convertirlo a 40 filas
sample_suplementary_file <- t(sample_suplementary_file)
dim(sample_suplementary_file)

targets <- data.frame(sampleNames = colnames(matriz_GSE), group = colnames(copia_matriz),
                      sup_file = sample_suplementary_file, row.names = 1)

info_columns <- data.frame(labelDescription = c("Type of Sample", "More info link"))

phenoData_colon <- new("AnnotatedDataFrame", data=targets, varMetadata= info_columns)


info <- new("MIAME",name = "Takaaki Kobayashi", lab = "Kyorin university", 
             contact = "ck9t-kbys@asahi-net.or.jp", title= "Expression data at each site in colon cancer")

#Para utilizar la funcion ExpressionSet necesitamos convertir la matriz_GSE (es un dataframe),
# en matriz

matriz_GSE <- data.matrix(matriz_GSE)
myEset <- ExpressionSet(assayData = matriz_GSE, phenoData = phenoData_colon,
                        experimentData = info, )

class(myEset)

dim(exprs(myEset))
head(pData(myEset))

#Rompemos el expresionset en grupo central y frontal.
#Para ello creamos una funcion para obtener los indices de cada muestra.

get_index <- function(x, pair=0){
  count <- 1
  pair_index <- c()
  odd_index <- c()
  for (i in x){
    if(as.integer(substr(i, nchar(i)-1, nchar(i))) %% 2 == 0){
      pair_index <- c(pair_index, count)
    }
    else{
      odd_index <- c(odd_index, count)
    }
    count<- count + 1
  }
  if (pair==0){
    return(pair_index)
  }
  else{
    return(odd_index)
  }
}
samples <- sampleNames(myEset)

subset_central <- myEset[, get_index(samples)]

subset_frontal <- myEset[, get_index(samples,1)]

head(exprs(subset_central))
head(exprs(subset_frontal))

head(pData(subset_central)) #podemos ver que el phenoData también ha sido dividido

#Podemos cambiar los datos de expresion del subset_central

exprs(subset_central)[1,1] <- 0.0098
head(exprs(subset_central)) #El primer valor de expresion ha sido cambiado a 0.0098

pData(subset_central)$group[1] <- "lo que me apetezca"
head(pData(subset_central))     #Podemos cambiar los grupos en el pData tambien

dif <- exprs(subset_central)/exprs(subset_frontal)
p_dif <- data.frame(sampleNames = paste(sampleNames(subset_central), " - ", 
                                            sampleNames(subset_frontal)), row.names = 1,
                    group = paste(pData(subset_central)$group," - ", pData(subset_frontal)$group))

var_dif <- data.frame(labelDescription = c("Samples used"), row.names = 1)

pheno_dif <- new("AnnotatedDataFrame", data = p_dif, varMetadata = var_dif )
dif_eset <- ExpressionSet(assayData = dif)
phenoData(dif_eset) <- pheno_dif
head(pData(dif_eset))
min(exprs(dif_eset))
dif_sorted <- sort(dif)



#Ejericicio 2 diapositivas
#GEO query

BiocManager::install("GEOquery")
require(GEOquery)

browseVignettes("GEOquery")
gse <- getGEO("GSE65480")
class(gse)
gse[[1]]

ExpressionSet_colonGEO <- gse[[1]]
dim(ExpressionSet_colonGEO)


