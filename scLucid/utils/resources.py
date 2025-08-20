import logging
from typing import Dict, List

log = logging.getLogger(__name__)

__all__ = ["load_gmt_file"]

def load_gmt_file(file_path: str) -> Dict[str, List[str]]:
    """
    Parses a .gmt file into a dictionary of gene sets.

    Args:
        file_path: Path to the .gmt file.

    Returns:
        A dictionary where keys are gene set names and values are lists of genes.
    """
    log.info(f"Loading gene sets from GMT file: {file_path}")
    gene_sets = {}
    try:
        with open(file_path, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                set_name = parts[0]
                genes = parts[2:]
                gene_sets[set_name] = genes
        log.info(f"Successfully loaded {len(gene_sets)} gene sets.")
        return gene_sets
    except FileNotFoundError:
        log.error(f"GMT file not found at: {file_path}")
        raise
    except Exception as e:
        log.error(f"Failed to parse GMT file: {e}")
        raise