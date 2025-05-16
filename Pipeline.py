import sys
from itertools import combinations
import numpy as np
from Bio import pairwise2
from Bio.Align.substitution_matrices import load
from utils import *
from itertools import combinations

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

    def extract_subsequences(self, start_motif, end_motif):
        """
        Extrae subsecuencias entre dos motivos.

        Args:
            start_motif (str): Motivo de inicio.
            end_motif (str): Motivo de fin.

        Returns:
            dict: Diccionario con los identificadores de secuencia como claves y las subsecuencias extraídas como valores.
        """
        return extract_subsequences(self.sequences, start_motif, end_motif)

    def gc_content(self, pos_start=None, pos_end=None):
        """
        Calcula el contenido de GC para cada secuencia.

        Args:
            pos_start (int): Posición de inicio del cálculo.
            pos_end (int): Posición de fin del cálculo.

        Returns:
            list: Lista de tuplas con los IDs de las secuencias y su contenido de GC.
        """
        gc_contents = []
        for seq in self.sequences:
            gc_percent = gc_content(str(seq.seq), pos_start, pos_end)
            gc_contents.append((seq.id, gc_percent))
        return gc_contents
    
    def count_motifs(self, motifs, overlapping=False):
        """
        Cuenta la frecuencia de patrones de ADN específicos en las secuencias.

        Args:
            motifs (list): Lista de motivos a buscar.
            overlapping (bool): Si se permite solapamiento en la búsqueda.

        Returns:
            dict: Diccionario con los IDs de las secuencias como claves y un diccionario de motivos como valores.
        """
        motif_counts = {}
        for seq in self.sequences:
            seq_id = seq.id
            motif_counts[seq_id] = count_motifs(str(seq.seq), motifs, overlapping)
        return motif_counts
    
    def identify_codons(self, table=11):
        """
        Identifica los codones de inicio y parada en las secuencias de mRNA.

        Args:
            table (int): Código de la tabla de codones.

        Returns:
            dict: Diccionario con los IDs de las secuencias como claves y una lista de codones como valores.
        """
        codons = []
        for seq in self.sequences:
            codons_in_seq = identify_codons(str(seq.seq), table)
            codons.append((seq.id, codons_in_seq))
        return codons
    
    def check_translations(self, peptide):
        """
        Comprueba si las secuencias de ADN se traducen correctamente a un péptido.
        
        Args:
            motifs (list): Lista de motivos a buscar.
            overlapping (bool): Si se permite solapamiento en la búsqueda.
        
        Returns:
            dict: Diccionario con los IDs de las secuencias como claves y un diccionario de motivos como valores.
        """
        translations = []
        for seq in self.sequences:
            iguales, discrepancies = check_translation(str(seq.seq), peptide)
            translations.append((seq.id, iguales, discrepancies))
        return translations
    
    def save_sequences(self, output_file):
        """
        Guarda secuencias en un archivo.
        """
        return save_sequences(self.sequences, output_file)

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
    
    def save_alignment_results(self, results, output_file):
        """
        Guarda los resultados de los alineamientos en un archivo.
        Args:
            results (list): lista con los resultados.
            output_file (str): nombre del fichero de salida
        """
        with open(output_file, "w") as f:
            for result in results:
                id1, id2, aligned1, aligned2, score = result
                f.write(f"Comparacion: {id1} vs {id2}\n")
                f.write(f"Alineamiento 1: {aligned1}\n")
                f.write(f"Alineamiento 2: {aligned2}\n")
                f.write(f"Puntuacio4n: {score}\n")
                f.write("-" * 80 + "\n")
        print(f"Resultados guardados en {output_file}") 

if __name__ == "__main__":
    if len(sys.argv) < 4: # Comprobar argumentos
        print("Uso: python Pipeline.py <archivo.fasta> <matriz_alineamiento> <output_file.txt>")
        sys.exit(1) 
    
    fasta_file = sys.argv[1] # Nombre del archivo fasta, el primer argumento
    matrix_name = sys.argv[2] # Nombre de la matriz de alineamiento que queramos usar
    output_file = sys.argv[3] # Nombre del fichero de salida
    pipeline = Pipeline()
    sequences=load_sequences(fasta_file, file_format="fasta", verbose=True)

    print(f"\nCargando archivo: {fasta_file}")
    print(f"Usando matriz de alineamiento: {matrix_name}")

    # Comparación de tiempos por combinación
    results = []
    for i, (seq1, seq2) in enumerate(combinations(sequences, 2), 1):
        print(f"\n🔹 Comparando par {i}: {seq1.id} vs {seq2.id}")
        # Nuestra implementación
        aligned1, aligned2, score_sw = Pipeline.smith_waterman(str(seq1.seq), str(seq2.seq), matrix_name, -4)
        print(f"Smith-Waterman → Puntuación: {score_sw}")
        results.append((seq1.id, seq2.id, aligned1, aligned2, score_sw))
        # Biopython
        alignments = Pipeline.pairwise2_alignment(str(seq1.seq), str(seq2.seq), matrix_name)
    
    # Comparción de tiempos totales
        # Nuestra implementación
    start_time = time.time()
    for seq1, seq2 in combinations(sequences, 2):
        aligned1, aligned2, score_sw = Pipeline.smith_waterman.__wrapped__(str(seq1.seq), str(seq2.seq), matrix_name, -4)
    end_time = time.time()
    total_time = end_time - start_time
    print(f"\n Tiempo total de ejecución para todas las combinaciones con nuestra implementación: {total_time:.4f} segundos")
        # Usando Biopython
    start_time_2 = time.time()
    for seq1, seq2 in combinations(sequences, 2):
        alignments = Pipeline.pairwise2_alignment.__wrapped__(str(seq1.seq), str(seq2.seq), matrix_name)
    end_time_2 = time.time()
    total_time_2 = end_time_2 - start_time_2
    print(f"\n Tiempo total de ejecución para todas las combinaciones con la dunción de Biopython: {total_time_2:.4f} segundos")

    # Guardado de los resultados
    pipeline.save_alignment_results(results,output_file)

