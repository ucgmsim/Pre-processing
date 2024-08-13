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

from srf_generation import realisation
from srf_generation.realisation import Realisation

app = typer.Typer()


def plot_realisation(realisation_obj: Realisation, output_path: Path):
    """Plot faults involved in a rupture scenario using pygmt.

    Parameters
    ----------
    realisation : Realisation
        The realisation to plot
    output_path : Path
        The path to output the realisation plot to.
    """

    corners = np.vstack([fault.corners() for fault in realisation_obj.faults.values()])
    region = (
        corners[:, 1].min() - 0.5,
        corners[:, 1].max() + 0.5,
        corners[:, 0].min() - 0.25,
        corners[:, 0].max() + 0.25,
    )
    fig = plotting.gen_region_fig(realisation_obj.name, region=region)

    for rupture_fault in realisation.topologically_sorted_faults(
        list(realisation_obj.faults.values())
    ):
        for plane in rupture_fault.planes:
            corners = plane.corners
            fig.plot(
                x=corners[:, 1].tolist() + [corners[0, 1]],
                y=corners[:, 0].tolist() + [corners[0, 0]],
                pen="1p",
                fill="red",
            )
        hypocentre = rupture_fault.fault_coordinates_to_wgsdepth_coordinates(
            np.array([rupture_fault.shyp, rupture_fault.dhyp])
        )
        if rupture_fault.parent:
            parent_centroid = (
                rupture_fault.parent.fault_coordinates_to_wgsdepth_coordinates(
                    np.array([0, 0])
                )
            )[:2][::-1]
            fault_centroid = rupture_fault.fault_coordinates_to_wgsdepth_coordinates(
                np.array([0, 0])
            )[:2][::-1]
            direction_vector = np.append(parent_centroid, fault_centroid).reshape(
                (1, -1)
            )
            fig.plot(
                data=direction_vector,
                style="=0.5c+s+e+a30+gblack+h0.5+p1p,black",
                pen="1p",
            )
            fig.plot(
                x=hypocentre[1],
                y=hypocentre[0],
                style="c0.1c",
                fill="white",
                pen="black",
            )
        else:
            fig.plot(
                x=hypocentre[1],
                y=hypocentre[0],
                style="a0.3c",
                fill="white",
                label="Hypocentre",
                pen="black",
            )
    fig.legend()
    fig.savefig(output_path)


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
    realisation_obj = realisation.read_realisation(realisation_path)
    plot_realisation(realisation_obj, output_path)


if __name__ == "__main__":
    app()
