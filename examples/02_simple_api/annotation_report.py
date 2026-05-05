"""Generate a publication-style annotation review report with scLucid."""

from pathlib import Path

import scanpy as sc

from scLucid.analysis import AnnotationConfig, run_annotation
from scLucid.plotting import export_annotation_report


def main() -> None:
    data_path = Path("data/annotated_example.h5ad")
    output_dir = Path("results/annotation_report")
    output_dir.mkdir(parents=True, exist_ok=True)

    adata = sc.read_h5ad(data_path)

    config = AnnotationConfig(
        cluster_key="leiden_clusters",
        marker_species="human",
        run_scoring=True,
        run_celltypist=True,
        final_method="hybrid",
        key_added="cell_type_hybrid",
        save_dir=str(output_dir),
        plot=False,
        report=True,
    )
    adata = run_annotation(adata, config=config)

    export_annotation_report(
        adata,
        annotation_key="cell_type_hybrid",
        cluster_key="leiden_clusters",
        save=str(output_dir / "cell_type_hybrid_annotation_report.png"),
        export_formats=("png", "pdf"),
        write_sidecars=True,
        show=False,
    )


if __name__ == "__main__":
    main()
