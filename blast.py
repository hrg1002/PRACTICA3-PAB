'''
A traves del uso de BLAST, se obtendria el accesion number de la secuencia proporcionada en el fichero gene.fasta.
Podeis utilizar cualquiera de las alternativas de uso de BLAST estudiadas en el Notebook 2B de la asignatura
(en local o contra el servidor de NCBI).
Almacenad el resultado obtenido con BLAST en el fichero que se ha indicado en el parametro --output.
Como seleccionar el accesion number del hit adecuado:
• Deberia tener un e-value menor de 0.01
• Entre los que cumplan la condicion anterior seleccionaremos el que mayor score tenga.
'''
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SeqIO

#leemos la secuencia del fichero gene.fasta
registro= SeqIO.read("data/gene.fasta", "fasta")

# ejecutamos el BLAST
blast_result = NCBIWWW.qblast("blastp", "nr", registro.seq, expect=0.01, hitlist_size=10)

# Guardamos el resultado en un fichero 
with open("data/salida.xml", "w") as archivo_blast:
    archivo_blast.write(blast_result.read())


# Cerramos el resultado del BLAST
blast_result.close()

# Leemos el resultado del BLAST
# queremos el accession number del hit adecuado y que mayor score tenga

with open("data/salida.xml") as archivo_blast:
    blast_record = NCBIXML.read(archivo_blast)

    # Inicializamos el mejor hit
    mejor_hit = None
    mejor_score = 0

    # Iteramos sobre los hits
    for hit in blast_record.alignments:
        for hsp in hit.hsps:
            if hsp.expect < 0.01 and hsp.score > mejor_score:
                mejor_hit = hit
                mejor_score = hsp.score


    # Guardar el accession number del mejor hit en un archivo de salida
    with open("data/salida.txt", "w") as salida:
        if mejor_hit:
            salida.write(f"Accession number: {mejor_hit.accession}\n")
        else:
            salida.write("No se encontraton\n")
