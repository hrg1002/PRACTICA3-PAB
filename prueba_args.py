import argparse
parser = argparse . ArgumentParser ( description =" Ejemplo con varios argumentos ")
parser . add_argument (" --input ", "-i", type =str , help =" Ruta al archivo de entrada ")

parser . add_argument (" -- output ", "-o", type =str , default =" salida . txt",
help =" Ruta al archivo de salida ( por defecto : salida . txt )")

# Argumento opcional como bandera
parser . add_argument (" -- verbose ", "-v", action =" store_true ",
help =" Activa la salida detallada ")

args = parser . parse_args ()
if args . verbose :
    print (" Modo detallado activado ")
    print (f" Archivo de entrada : { args . input_file }")
    print (f" Archivo de salida : { args . output }")
else :
    print (f" Procesando { args . input_file }... ")