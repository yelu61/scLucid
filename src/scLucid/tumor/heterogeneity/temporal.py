"""
Temporal dynamics and longitudinal heterogeneity analysis.

This module provides tools for tracking tumor heterogeneity
dynamics over time and treatment response.
"""

import numpy as np
import pandas as pd
from typing import Dict, List, Optional, Tuple, Union
from anndata import AnnData
import logging

log = logging.getLogger(__name__)


class TemporalAnalyzer:
    """
    Analyze temporal dynamics of tumor heterogeneity.

    Parameters
    ----------
    time_key : str
        Column containing time points

    Attributes
    ----------
    dynamics_ : pd.DataFrame
        Temporal dynamics metrics
    """

    def __init__(
        self,
        time_key: str = "time_point",
    ):
        self.time_key = time_key
        self.dynamics_: Optional[pd.DataFrame] = None

    def track_temporal_dynamics(
        self,
        adata: AnnData,
        clone_key: str = "clone_id",
        patient_key: Optional[str] = None,
    ) -> pd.DataFrame:
        """
        Track clone dynamics over time.

        Parameters
        ----------
        adata : AnnData
            Expression data with multiple time points
        clone_key : str
            Column containing clone IDs
        patient_key : str, optional
            Column containing patient IDs

        Returns
        -------
        pd.DataFrame
            Temporal dynamics metrics
        """
        results = []

        if patient_key is None:
            patients = ["patient_1"]
        else:
            patients = adata.obs[patient_key].unique()

        for patient in patients:
            if patient_key:
                patient_mask = adata.obs[patient_key] == patient
                patient_adata = adata[patient_mask]
            else:
                patient_adata = adata

            time_points = sorted(patient_adata.obs[self.time_key].unique())

            for time_point in time_points:
                mask = patient_adata.obs[self.time_key] == time_point
                time_adata = patient_adata[mask]

                # Calculate clonal composition
                clone_counts = time_adata.obs[clone_key].value_counts()
                clone_props = clone_counts / clone_counts.sum()

                # Calculate diversity
                shannon = -np.sum(clone_props * np.log(clone_props + 1e-10))
                simpson = 1 - np.sum(clone_props ** 2)

                result = {
                    "patient": patient,
                    "time_point": time_point,
                    "n_cells": time_adata.n_obs,
                    "n_clones": len(clone_counts),
                    "shannon_diversity": shannon,
                    "simpson_diversity": simpson,
                    "dominant_clone": clone_counts.index[0],
                    "dominant_clone_freq": clone_props.iloc[0],
                }

                # Add clone frequencies
                for clone, freq in clone_props.items():
                    result[f"clone_{clone}_freq"] = freq

                results.append(result)

        self.dynamics_ = pd.DataFrame(results)
        return self.dynamics_

    def analyze_treatment_response_trajectory(
        self,
        adata: AnnData,
        treatment_key: str = "treatment",
        response_key: str = "response",
        clone_key: str = "clone_id",
    ) -> pd.DataFrame:
        """
        Analyze treatment response trajectories.

        Parameters
        ----------
        adata : AnnData
            Expression data with treatment information
        treatment_key : str
            Column containing treatment information
        response_key : str
            Column containing response status
        clone_key : str
            Column containing clone IDs

        Returns
        -------
        pd.DataFrame
            Treatment response trajectories
        """
        results = []

        time_points = sorted(adata.obs[self.time_key].unique())

        for time_point in time_points:
            mask = adata.obs[self.time_key] == time_point
            time_adata = adata[mask]

            result = {
                "time_point": time_point,
                "n_cells": time_adata.n_obs,
            }

            # Treatment status
            if treatment_key in time_adata.obs.columns:
                result["treatment"] = time_adata.obs[treatment_key].iloc[0]

            # Response status
            if response_key in time_adata.obs.columns:
                result["response"] = time_adata.obs[response_key].iloc[0]

            # Clonal composition
            clone_counts = time_adata.obs[clone_key].value_counts()
            result["n_clones"] = len(clone_counts)
            result["shannon_diversity"] = -np.sum(
                (clone_counts / clone_counts.sum()) * np.log(clone_counts / clone_counts.sum() + 1e-10)
            )

            results.append(result)

        return pd.DataFrame(results)

    def calculate_clone_turnover(
        self,
        adata: AnnData,
        clone_key: str = "clone_id",
        patient_key: Optional[str] = None,
    ) -> pd.DataFrame:
        """
        Calculate clone turnover between time points.

        Parameters
        ----------
        adata : AnnData
            Expression data with multiple time points
        clone_key : str
            Column containing clone IDs
        patient_key : str, optional
            Column containing patient IDs

        Returns
        -------
        pd.DataFrame
            Clone turnover metrics
        """
        results = []

        if patient_key is None:
            patients = ["patient_1"]
        else:
            patients = adata.obs[patient_key].unique()

        for patient in patients:
            if patient_key:
                patient_mask = adata.obs[patient_key] == patient
                patient_adata = adata[patient_mask]
            else:
                patient_adata = adata

            time_points = sorted(patient_adata.obs[self.time_key].unique())

            for i in range(len(time_points) - 1):
                t1, t2 = time_points[i], time_points[i + 1]

                mask1 = patient_adata.obs[self.time_key] == t1
                mask2 = patient_adata.obs[self.time_key] == t2

                clones1 = set(patient_adata[mask1].obs[clone_key].unique())
                clones2 = set(patient_adata[mask2].obs[clone_key].unique())

                # Calculate turnover metrics
                shared = len(clones1 & clones2)
                lost = len(clones1 - clones2)
                gained = len(clones2 - clones1)

                jaccard = shared / len(clones1 | clones2) if len(clones1 | clones2) > 0 else 0

                results.append({
                    "patient": patient,
                    "time_from": t1,
                    "time_to": t2,
                    "shared_clones": shared,
                    "lost_clones": lost,
                    "gained_clones": gained,
                    "jaccard_index": jaccard,
                    "turnover_rate": (lost + gained) / len(clones1) if len(clones1) > 0 else 0,
                })

        return pd.DataFrame(results)


def track_temporal_dynamics(
    adata: AnnData,
    time_key: str = "time_point",
    clone_key: str = "clone_id",
    patient_key: Optional[str] = None,
) -> pd.DataFrame:
    """
    Track clone dynamics over time.

    Parameters
    ----------
    adata : AnnData
        Expression data with multiple time points
    time_key : str
        Column containing time points
    clone_key : str
        Column containing clone IDs
    patient_key : str, optional
        Column containing patient IDs

    Returns
    -------
    pd.DataFrame
        Temporal dynamics metrics
    """
    analyzer = TemporalAnalyzer(time_key=time_key)
    return analyzer.track_temporal_dynamics(adata, clone_key, patient_key)


def analyze_treatment_response_trajectory(
    adata: AnnData,
    time_key: str = "time_point",
    treatment_key: str = "treatment",
    response_key: str = "response",
    clone_key: str = "clone_id",
) -> pd.DataFrame:
    """
    Analyze treatment response trajectories.

    Parameters
    ----------
    adata : AnnData
        Expression data with treatment information
    time_key : str
        Column containing time points
    treatment_key : str
        Column containing treatment information
    response_key : str
        Column containing response status
    clone_key : str
        Column containing clone IDs

    Returns
    -------
    pd.DataFrame
        Treatment response trajectories
    """
    analyzer = TemporalAnalyzer(time_key=time_key)
    return analyzer.analyze_treatment_response_trajectory(
        adata, treatment_key, response_key, clone_key
    )


def detect_clonal_sweep(
    adata: AnnData,
    time_key: str = "time_point",
    clone_key: str = "clone_id",
    frequency_threshold: float = 0.5,
) -> pd.DataFrame:
    """
    Detect clonal sweeps (expansion of a single clone).

    Parameters
    ----------
    adata : AnnData
        Expression data with multiple time points
    time_key : str
        Column containing time points
    clone_key : str
        Column containing clone IDs
    frequency_threshold : float
        Threshold for detecting sweep

    Returns
    -------
    pd.DataFrame
        Detected clonal sweeps
    """
    results = []

    time_points = sorted(adata.obs[time_key].unique())

    for clone in adata.obs[clone_key].unique():
        frequencies = []

        for time_point in time_points:
            mask = adata.obs[time_key] == time_point
            time_adata = adata[mask]

            clone_freq = (time_adata.obs[clone_key] == clone).mean()
            frequencies.append(clone_freq)

        frequencies = np.array(frequencies)

        # Detect sweep: initial low frequency -> final high frequency
        if frequencies[0] < 0.2 and frequencies[-1] > frequency_threshold:
            results.append({
                "clone": clone,
                "initial_freq": frequencies[0],
                "final_freq": frequencies[-1],
                "fold_change": frequencies[-1] / (frequencies[0] + 1e-6),
                "is_sweep": True,
            })

    return pd.DataFrame(results)
