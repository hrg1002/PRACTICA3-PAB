Ejecutar en la terminal una vez estamos dentro de la carpeta con los archivos siguiendo esta estructura:
python Pipeline.py -i <archivo.fasta> -m <matriz de alineamiento> -g <penalización de gap> -o <nombre archivo salida> -v

Ejemplo:
python Pipeline.py -i data/gene.fasta -m BLOSUM62 -g -10.0 -o salida.txt -v


Para ver el accesion number se crea el archivo txt en la carpeta 'data'