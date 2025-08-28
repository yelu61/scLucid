import subprocess

def check_r():
    try:
        output = subprocess.check_output(["R", "--version"])
        print("R found:", output.decode().split('\n')[0])
        return True
    except FileNotFoundError:
        print("R is not installed or not in PATH.")
        return False

def check_rpy2():
    try:
        import rpy2
        print("rpy2 found:", rpy2.__version__)
        return True
    except ImportError:
        print("rpy2 is not installed.")
        return False

def check_r_package(pkg):
    cmd = f'Rscript -e "if (!requireNamespace(\'{pkg}\', quietly=TRUE)) quit(status=1)"'
    try:
        subprocess.check_call(cmd, shell=True)
        print(f"R package {pkg}: OK")
        return True
    except subprocess.CalledProcessError:
        print(f"R package {pkg} not found. Install with: R -e 'install.packages(\"{pkg}\")' or BiocManager.")
        return False

if __name__ == "__main__":
    check_r()
    check_rpy2()
    for pkg in ["Seurat", "SingleCellExperiment", "monocle3", "CellChat", "CopyKAT", "BayesPrism", "DWLS", "DESeq2"]:
        check_r_package(pkg)