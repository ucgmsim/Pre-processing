#!/usr/bin/env python3
"""
This script allows the user to visualise a type5 realisation, including
fault geometry, jump points and planned rupture propagation.

Usage
-----
```
$ python visualise_type5_realisation.py path/to/realisation.yaml output_directory
```
"""

from pathlib import Path
from typing import Annotated

import numpy as np
import typer
from pygmt_helper import plotting
from source_modelling import sources
from workflow import realisations

app = typer.Typer()


def generate_figure(realisation_name: str, faults: list[sources.Fault]) -> None:
    corners = np.vstack([fault.corners for fault in faults])
    region = (
        corners[:, 1].min() - 0.5,
        corners[:, 1].max() + 0.5,
        corners[:, 0].min() - 0.25,
        corners[:, 0].max() + 0.25,
    )
    return plotting.gen_region_fig(realisation_name, region=region)


def plot_faults(fig, faults: list[sources.Fault]) -> None:
    for fault in faults:
        for plane in fault.planes:
            corners = plane.corners
            fig.plot(
                x=corners[:, 1].tolist() + [corners[0, 1]],
                y=corners[:, 0].tolist() + [corners[0, 0]],
                pen="1p",
                fill="red",
            )


def plot_on_fault(
    fig, fault: sources.Fault, local_coordinates: np.ndarray, **plot_args
) -> None:
    global_coordinates = fault.fault_coordinates_to_wgs_depth_coordinates(
        local_coordinates
    )
    fig.plot(x=global_coordinates[1], y=global_coordinates[0], **plot_args)


def plot_arc_between(fig, fault_a: sources.Fault, fault_b: sources.Fault) -> None:
    fault_a_centroid = fault_a.centroid[:2][::-1]
    fault_b_centroid = fault_b.centroid[:2][::-1]
    direction_vector = np.append(fault_a_centroid, fault_b_centroid).reshape((1, -1))
    fig.plot(
        data=direction_vector,
        style="=0.5c+s+e+a30+gblack+h0.5+p1p,black",
        pen="1p",
    )


@app.command()
def main(
    realisation_path: Annotated[
        Path,
        typer.Argument(
            help="Path to realisation file", exists=True, dir_okay=False, readable=True
        ),
    ],
    output_path: Annotated[
        Path, typer.Argument(help="Output plot path", dir_okay=False, writable=True)
    ],
):
    "Plot a realisation with pygmt."
    rupture_propagation_config: realisations.RupturePropagationConfig = (
        realisations.read_config_from_realisation(
            realisations.RupturePropagationConfig, realisation_path
        )
    )
    source_config: realisations.SourceConfig = (
        realisations.read_config_from_realisation(
            realisations.SourceConfig, realisation_path
        )
    )
    metadata: realisations.RealisationMetadata = (
        realisations.read_config_from_realisation(
            realisations.RealisationMetadata, realisation_path
        )
    )
    initial_fault_name = next(
        fault_name
        for fault_name, parent in rupture_propagation_config.rupture_causality_tree.items()
        if parent is None
    )
    faults = source_config.source_geometries
    fig = generate_figure(metadata.name, list(faults.values()))
    plot_faults(fig, list(faults.values()))
    initial_fault = faults[initial_fault_name]
    plot_on_fault(
        fig,
        initial_fault,
        rupture_propagation_config.hypocentre,
        style="a0.3c",
        fill="white",
        pen="black",
    )

    for fault, parent in rupture_propagation_config.rupture_causality_tree.items():
        if parent is None:
            continue
        plot_arc_between(fig, faults[parent], faults[fault])
        jump_point = rupture_propagation_config.jump_points[fault]
        plot_on_fault(fig, faults[parent], jump_point.from_point, style='c0.1c', fill='white', pen='black')
        plot_on_fault(fig, faults[fault], jump_point.to_point, style='c0.1c', fill='white', pen='black')
    fig.savefig(output_path)


if __name__ == "__main__":
    app()
