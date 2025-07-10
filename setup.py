from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="scRNA-toolkit",
    version="0.1.0",
    author="Ye Lu",
    author_email="foliageyelu@gmail.com",
    description="A toolkit for single-cell RNA-seq data analysis",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/yelu61/scRNA-toolkit.git",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    python_requires=">=3.7",
    install_requires=[
        "scanpy>=1.9.0",
        "anndata>=0.8.0",
        "numpy>=1.20.0",
        "pandas>=1.3.0",
        "matplotlib>=3.4.0",
        "seaborn>=0.11.0",
        "scikit-learn>=1.0.0",
        "scrublet>=0.2.3",
        "infercnvpy>=0.4.0",
        "scvi-tools>=0.16.0",
    ],
    entry_points={
        'console_scripts': [
            'scrna-toolkit=scRNA.cli:main', 
        ],
    },
)