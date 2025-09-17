import subprocess
from typing import Dict, List

def check_r_environment(packages: List[str]) -> Dict[str, bool]:
    """
    Checks for the presence of R, rpy2, and a list of R packages.

    Returns:
        A dictionary mapping each checked item to a boolean (True=found, False=missing).
    """
    status = {}
    # Check for R itself
    try:
        subprocess.check_output(["R", "--version"])
        status["R"] = True
    except FileNotFoundError:
        status["R"] = False

    # Check for rpy2
    try:
        import rpy2
        status["rpy2"] = True
    except ImportError:
        status["rpy2"] = False

    # Check for R packages
    for pkg in packages:
        cmd = f'Rscript -e "if (!requireNamespace(\'{pkg}\', quietly=TRUE)) quit(status=1)"'
        try:
            # Use subprocess.run for better control
            subprocess.run(cmd, shell=True, check=True, capture_output=True)
            status[pkg] = True
        except subprocess.CalledProcessError:
            status[pkg] = False
            
    return status

def print_r_environment_status():
    """Prints a formatted report of the R environment status."""
    pkgs_to_check = ["Seurat", "SingleCellExperiment", "monocle3", "CellChat", "CopyKAT", "BayesPrism", "DWLS", "DESeq2", "slingshot"]
    status = check_r_environment(pkgs_to_check)
    
    print("--- scLucid R Environment Check ---")
    print(f"R installed: {'✅' if status.get('R') else '❌'}")
    print(f"rpy2 installed: {'✅' if status.get('rpy2') else '❌'}")
    print("\nRequired R Packages:")
    for pkg in pkgs_to_check:
        print(f"  - {pkg}: {'✅' if status.get(pkg) else '❌'}")
    
    if not all(status.values()):
        print("\n[Warning] Your R environment is not fully configured for all tools.")
        print("Please install missing components. For R packages, use BiocManager::install('...') or install.packages('...').")