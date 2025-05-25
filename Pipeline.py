import argparse
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SeqIO 
from itertools import combinations
import numpy as np
from Bio import pairwise2
from Bio.Align.substitution_matrices import load
from utils import *
from itertools import combinations
from Bio import Entrez

CORREO = "raul.ortegare@gmail.com"
TU_API = "90fba815d4e3b0c6891f4558feedd9ad5e09"

class Pipeline:
    def __init__(self):
        """
        Inicializa el pipeline.
        """
        self.sequences = []

    def load_sequences(self, file_path, file_format, verbose=False):
        """
        Carga secuencias desde un archivo.

        Args:
            file_path (str): Ruta del archivo.
            file_format (str): Formato del archivo.
            verbose (bool): Si se imprime un resumen de las secuencias.
        """
        self.sequences = load_sequences(file_path, file_format, verbose)

    @staticmethod
    @contar_tiempo # Decorador para medir el tiempo de ejecución de smith_waterman
    def smith_waterman(seq1, seq2, matrix_name, gap_penalty):
        """
        Implementación del algoritmo Smith-Waterman.

        Args:
            seq1 (str): Primera secuencia.
            seq2 (str): Segunda secuencia.
            matrix_name (str): Nombre de la matriz de sustitución.
            gap_penalty (int): Penalización por huecos.

        Returns:
            tuple: (alineamiento seq1, alineamiento seq2, puntuación)
        """
        matrix = load(matrix_name)
        n = len(seq1) + 1
        m = len(seq2) + 1
        score_matrix = np.zeros((n, m), dtype=int)
        traceback = np.zeros((n, m), dtype=str)

        max_score = 0
        max_pos = (0, 0)

        for i in range(1, n):
            for j in range(1, m):
                match = matrix[seq1[i - 1], seq2[j - 1]]
                diag = score_matrix[i-1][j-1] + match
                up = score_matrix[i-1][j] + gap_penalty
                left = score_matrix[i][j-1] + gap_penalty
                score_matrix[i][j] = max(0, diag, up, left)

                if score_matrix[i][j] == diag:
                    traceback[i][j] = 'D'
                elif score_matrix[i][j] == up:
                    traceback[i][j] = 'U'
                elif score_matrix[i][j] == left:
                    traceback[i][j] = 'L'
                else:
                    traceback[i][j] = '0'

                if score_matrix[i][j] >= max_score:
                    max_score = score_matrix[i][j]
                    max_pos = (i, j)
        # Backtracking
        i, j = max_pos
        aligned1, aligned2 = [], []

        while score_matrix[i][j] != 0:
            if traceback[i][j] == 'D':
                aligned1.append(seq1[i - 1])
                aligned2.append(seq2[j - 1])
                i -= 1
                j -= 1
            elif traceback[i][j] == 'U':
                aligned1.append(seq1[i - 1])
                aligned2.append('-')
                i -= 1
            elif traceback[i][j] == 'L':
                aligned1.append('-')
                aligned2.append(seq2[j - 1])
                j -= 1

        aligned1.reverse()
        aligned2.reverse()

        return ''.join(aligned1), ''.join(aligned2), max_score
    
    @staticmethod
    @contar_tiempo
    def pairwise2_alignment(seq1, seq2, matrix_name, gap_open=-10, gap_extend=-0.5):
        """
        Realiza un alineamiento utilizando Biopython.

        Args:
            seq1 (str): Primera secuencia.
            seq2 (str): Segunda secuencia.
            matrix_name (str): Nombre de la matriz de sustitución.
            gap_open (int): Penalización por abrir un hueco.
            gap_extend (int): Penalización por extender un hueco.

        Returns:
            list: Lista de alineamientos.
        """
        alignments = pairwise2.align.globalds(seq1, seq2, load(matrix_name), gap_open, gap_extend)
        return alignments
    
    @staticmethod
    def blast(seq, output_xml="data/salida.xml"):
        '''
        Ejecuta BLAST para obtener el accession number de la secuencia proporcionada.
        Args:
            seq (str): Secuencia a buscar.
            output_xml (str): Ruta del archivo XML de salida.
        Returns:
            str: Accession number del mejor hit encontrado.  
        '''
        # Ejecutar BLAST
        blast_result = NCBIWWW.qblast("blastp", "nr", seq, expect=0.01, hitlist_size=10)
        with open(output_xml, "w") as archivo_blast:
            archivo_blast.write(blast_result.read())
        blast_result.close()

        # Leer el resultado BLAST
        with open(output_xml) as archivo_blast:
            blast_record = NCBIXML.read(archivo_blast)
            mejor_hit = None
            mejor_score = 0
            for hit in blast_record.alignments:
                for hsp in hit.hsps:
                    if hsp.expect < 0.01 and hsp.score > mejor_score:
                        mejor_hit = hit
                        mejor_score = hsp.score
            if mejor_hit:
                return mejor_hit.accession
            else:
                return None
@staticmethod
def busqueda(identificador):
    """
    Busca un identificador en la base de datos de NCBI y devuelve la secuencia asociada.
    
    Args:
        identificador (str): Identificador a buscar.
    
    Returns:
        
    """
    Entrez.email = CORREO
    Entrez.api_key = TU_API
    handle = Entrez.efetch(db = "protein", id = identificador, rettype = "gb", retmode = "text")
    record = SeqIO.read(handle, "genbank")
    handle.close()
    print(f"Accession number: {record.name}")
    organismo = record.annotations.get("organism")
    print(f"Organismo: {organismo}")

    return record.seq


if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description='Pipeline de alineamiento de secuencias')
    parser.add_argument('--input', '-i', type=str, required=True, help='Fichero con la secuencia de entrada (FASTA)')
    parser.add_argument('--matrix', '-m', type=str, required=True, help='Nombre de la matriz de alineamiento')
    parser.add_argument('--gap', '-g', type=float, default=-4.0, help='Penalización de gap')
    parser.add_argument('--output', '-o', type=str, required=True, help='Fichero de salida')
    parser.add_argument('--verbose', '-v', action='store_true', help='Modo detallado')

    args = parser.parse_args()

    pipeline = Pipeline()
    pipeline.load_sequences(args.input, file_format="fasta", verbose=args.verbose)

    print(f"\nCargando archivo: {args.input}")
    print(f"Usando matriz de alineamiento: {args.matrix}")

    sequences = pipeline.sequences
    matrix_name = args.matrix
    output_file = args.output

    # Ejecutar BLAST sobre la primera secuencia y guardar el accession number
    first_seq = sequences[0].seq
    accession = Pipeline.blast(str(first_seq))

    # Guardar el accession number en el archivo de salida antes de los alineamientos
    # Guardar el accession number del mejor hit en un archivo de salida
    with open("data/salida.txt", "w") as salida:
        if accession:
            salida.write(f"Accession number: {accession}\n")
        else:
            salida.write("No se encontró ningún hit significativo en BLAST\n")
    # Buscar la secuencia y obtener la secuencia asociada
    secuencia_buscada = busqueda(accession)
    aligned1, aligned2, score_sw = Pipeline.smith_waterman(first_seq, secuencia_buscada, matrix_name, args.gap)
    print("Alineamiento Smith-Waterman:")
    print(aligned1)
    print(aligned2)
    print("Puntuación:", score_sw)




