import argparse
# Definimos el parser
parser = argparse.ArgumentParser(description='Alineamiento de secuencias con el algoritmo de Smith-Waterman')

# Mis argumentos
# Definimos los argumentos que vamos a recibir por la línea de comandos
'''
Los argumentos serán los siguientes:
• --input / -i Fichero con la secuencia de entrada. Por no complicar mucho m´as el comando, supongamos
que siempre recibiremos el fichero en formato fasta.
• --matrix / -m Nombre de la matriz a utilizar en el alineamiento.
• --gap / -g Valor de penalización de gap en el alineamiento.
• --output / -o Fichero donde se almacenará la salida.

'''
parser.add_argument('--input', '-i', type=str, help='Fichero con la secuencia de entrada')
parser.add_argument('--matrix', '-m', type=str, help='Nombre de la matriz a utilizar en el alineamiento')
parser.add_argument('--gap', '-g', type=float, default=-10.0, help='Valor de penalización de gap en el alineamiento')
parser.add_argument('--output', '-o', type=str, help='Fichero donde se almacenará la salida')
parser.add_argument('--verbose', '-v', action='store_true', help='Modo detallado')


args = parser.parse_args ()
if args.verbose :
    print ("Modo detallado activado")
    print (f" Archivo de entrada : { args.input }")
    print (f" Matriz de alineamiento : { args.matrix }")
    print (f" Penalización de gap : { args.gap }")
    print (f" Archivo de salida : { args.output }")

# ejemplo de uso
# python argumentos.py -i data/gene.fasta -m BLOSUM62 -g -10.0 -o salida.txt -v 