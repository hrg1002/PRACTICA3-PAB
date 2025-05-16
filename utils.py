from Bio.Seq import Seq
from Bio.Align.substitution_matrices import load
import numpy as np
import time
from Bio import SeqIO
from functools import wraps

def get_seq_length(dna_seq):
    """
    Funcion que devuelve la longitud de una secuencia de DNA
    Args :
    dna_seq ( str ): Secuencia de DNA
    Returns :
    int : Longitud de la secuencia de DNA
    """
    return len( dna_seq )

# Ejercicio 1
def count_motifs(dna_seq, motifs, overlapping):
    """
    Cuenta la frecuencia de patrones de ADN específicos en la secuencia.
    Args:
        dna_seq (str): Secuencia de DNA.
        motifs (str): Patrón de ADN a buscar.
        overlapping (bool): Si se permite solapamiento en la búsqueda.
    Returns:
        int: Frecuencia del motivo en la secuencia.
    """
    contador = 0
    salto = 1 if overlapping else len(motifs)
    for i in range(0, len(dna_seq) - len(motifs) + 1, salto):
        if dna_seq[i:i+len(motifs)] == motifs:
            contador += 1
    return contador

# Ejercicio 2
def gc_content(dna_seq, pos_start, pos_end):
    """
    Calcula el porcentaje de GC en una secuencia de ADN.

    Args:
        dna_seq (str): Secuencia de ADN.
        pos_start (int): Posición inicial (opcional).
        pos_end (int): Posición final (opcional).

    Returns:
        float: Porcentaje de GC.
    """
    if pos_start is not None and (pos_start < 0 or pos_start >= len(dna_seq)):
        return print(f"pos_start={pos_start} está fuera del rango de la secuencia.")
    if pos_end is not None and (pos_end < 0 or pos_end > len(dna_seq)):
        return print(f"pos_end={pos_end} está fuera del rango de la secuencia.")
    if pos_start is not None and pos_end is not None:
        dna_seq = dna_seq[pos_start:pos_end]
    gc_count = dna_seq.count('G') + dna_seq.count('C')
    return (gc_count / len(dna_seq)) * 100

# Ejercicio 3
def get_coding_strand(mrna_seq):
    """
    Obtiene la hebra codificante de ADN a partir de una secuencia de mRNA.

    Args:
        mrna_seq (str): Secuencia de mRNA.

    Returns:
        str: Hebra codificante de ADN.
    """
    return str(Seq(mrna_seq).reverse_complement())

# Ejercicio 4
def identify_codons(mrna_seq, table):
    """
    Identifica los codones de inicio y parada de una secuencia de mRNA.

    Args:
        mrna_seq (str): Secuencia de mRNA.
        table (dict): Tabla de codones.

    Returns:
        list: Lista de tuplas con posiciones de inicio y parada, y la secuencia codificante.
    """
    codones = []
    for start_codon in table['start']:
        for stop_codon in table['stop']:
            start_pos = mrna_seq.find(start_codon)
            stop_pos = mrna_seq.find(stop_codon, start_pos)
            if start_pos != -1 and stop_pos != -1:
                codones.append((start_pos, stop_pos, mrna_seq[start_pos:stop_pos]))
    return codones

# Ejercicio 5
def check_translation(dna, peptide):
    """
    Comprueba si una secuencia de ADN se traduce correctamente a un péptido.

    Args:
        dna (str): Secuencia de ADN.
        peptide (str): Péptido.

    Returns:
        bool: (bool, list) Indicando si coinciden y las posiciones de discrepancia.
    """
    translation_seq = Seq(dna).translate()
    discrepancias = []
    for i in range(min(len(translation_seq),len(peptide))):
        if translation_seq[i] != peptide[i]:
            discrepancias.append(i+1)

    return (translation_seq == peptide, discrepancias)

# Ejercicio 6
def load_sequences(file_path, file_format, verbose=False):
    """
    Carga secuencias de un archivo.

    Args:
        file_path (str): Ruta del archivo.
        file_format (str): Formato del archivo.
        verbose (bool): Si se imprime un resumen de las secuencias.

    Returns:
        list: Lista de secuencias.
    """
    sequences = list(SeqIO.parse(file_path, file_format))
    if verbose:
        print(f"Secuencias totales: {len(sequences)}")
        for seq_record in sequences:
            print(f"ID: {seq_record.id}")
            print(f'Descripción del registro: {seq_record.description}')
            print(f"Secuencia: {repr(seq_record.seq)}")
            print(f"Longitud: {len(seq_record)}")
            
    return sequences

# Ejercicio 7
def filter_sequences_by_length(min_length, max_length, sequences):
    """
    Filtra secuencias por longitud.

    Args:
        min_length (int): Longitud mínima.
        max_length (int): Longitud máxima.
        sequences (list): Lista de secuencias.

    Returns:
        list: Lista de secuencias filtradas.
    """
    return [seq for seq in sequences if min_length <= len(seq) <= max_length]

# Ejercicio 8
def extract_subsequences(sequences, start_motif, end_motif):
    """
    Extrae subsecuencias entre dos motivos.
    
    Args:
        sequences (list): Lista de secuencias.
        start_motif (str): Motivo de inicio.
        end_motif (str): Motivo de fin.
        
    Returns:
        dic: Diccionario con los identificadores de secuencia como claves y las subsecuencias extraídas como valores.
    """
    subsecuencias = {}
    for seq in sequences:
        start_pos = seq.seq.find(start_motif)
        end_pos = seq.seq.find(end_motif, start_pos+len(start_motif))
        if start_pos != -1 and end_pos != -1:
            subsecuencias[seq.id] = seq.seq[start_pos:end_pos+len(end_motif)]

    return subsecuencias

# Ejercicio 9
def index_sequences(file_path, file_format):
    """
    Indexa las secuencias basado en sus identificadores.

    Args:
        file_path (str): Ruta del archivo.
        file_format (str): Formato del archivo.

    Returns:
        dict: Diccionario con los identificadores de secuencia como claves y las secuencias como valores
    """
    return SeqIO.index(file_path, file_format)

def get_sequence_by_id(seq_id, seqs_index):
    """
    Obtiene una secuencia basada en su identificador.

    Args:
        seq_id (str): Identificador de la secuencia.
        seqs_index (dict): Diccionario con los identificadores de secuencia como claves y las secuencias como valores.

    Returns:
        SeqRecord: Secuencia correspondiente al identificador.
    """
    return seqs_index[seq_id]

# Ejercicio 10
def save_sequences(output_file, sequences):
    """
    Guarda secuencias en un archivo.
    
    Args:
        output_file (str): Ruta del archivo.
        sequences (list): Lista de secuencias.

    Returns:
        None
    """
    SeqIO.write(sequences, output_file, "fasta")
    print(f"Número total de secuencias almacenadas: {len(sequences)}")
    print(f"Longitud media: {sum(len(seq) for seq in sequences) / len(sequences)}")

def contar_tiempo(funcion):
        """
        Decorador para medir el tiempo de ejecución de una función.
        """
        @wraps(funcion)
        def wrapper(*args, **kwargs):
            start_time = time.time()
            result = funcion(*args, **kwargs)
            end_time = time.time()
            # Calcula el tiempo de ejecución
            print(f"Tiempo de ejecución de {funcion.__name__}: {end_time - start_time:.2f} segundos")
            return result
        return wrapper