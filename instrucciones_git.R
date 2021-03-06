#To use the save of versions with git in a local environment 

# Seleccionamos el fichero en el cuadrante alto derecho, le damos a commit para guardar.

#Para encontrar antiguas versionas, seleccionamos diff y luego history. 
#Copiamos el SHA de la version buscada y utilizamos el siguiente comando de terminal.
#git reset --hard SEGUIDO DEL SHA

#Para crear un repositorio en github, necesitamos un access token que obtenemos 
#a partir de nuestra cuenta github

library(usethis)

#Una vez obtenido el tokken, lo guardamos en nuestro environment con la funcion

edit_r_environ() #Bajo el nombre GITHUB_PAT = "tokken"

#Para crear el repo

use_github(protocol = "https", auth_token = Sys.getenv("GITHUB_PAT"))

#Modificamos ahora para ver que pasa