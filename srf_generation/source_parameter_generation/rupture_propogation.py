"""
Rupture Propogation Module

This module provides functions for computing likely rupture paths from
information about the distances between faults.

Functions:
    - shaw_dieterich_distance_model:
      Compute fault jump probabilities using the Shaw-Dieterich distance model.
    - prune_distance_graph: Prune the distance graph based on a cutoff value.
    - probability_graph: 
      Convert a graph of distances between faults into a graph of jump
      probabilities using the Shaw-Dieterich model.
    - probabilistic_minimum_spanning_tree: Generate a probabilistic minimum spanning tree.
    - probabilistic_shortest_path: Generate the shortest path tree from the initial fault.

Typing Aliases:
    - DistanceGraph: A graph representing distances between faults.
    - ProbabilityGraph: A graph representing jump probabilities between faults.
    - RuptureCausalityTree: A tree representing the causality of ruptures between faults.
"""

from collections import defaultdict
from typing import Optional

import numpy as np

DistanceGraph = dict[str, dict[str, int]]
ProbabilityGraph = defaultdict[str, dict[str, float]]

RuptureCausalityTree = dict[str, Optional[str]]


def shaw_dieterich_distance_model(distance: float, d0, delta) -> float:
    """
    Compute fault jump probabilities using the Shaw-Dieterich distance model.

    Parameters
    ----------
    distance : float
        The distance between two faults.
    d0 : float
        The characteristic distance parameter.
    delta : float
        The characteristic slip distance parameter.

    Returns
    -------
    float
        The calculated probability.
    """
    return min(1, np.exp(-(distance - delta) / d0))


def prune_distance_graph(distances: DistanceGraph, cutoff: int) -> DistanceGraph:
    """
    Prune the distance graph based on a cutoff value.

    Parameters
    ----------
    distances : DistanceGraph
        The graph of distances between faults.
    cutoff : int
        The cutoff distance in metres.

    Returns
    -------
    DistanceGraph
        A copy of the input distance graph, keeping only edges that are less
        than the cutoff.
    """
    return {
        fault_u: {
            fault_v: distance
            for fault_v, distance in neighbours_fault_u.items()
            if distance < cutoff
        }
        for fault_u, neighbours_fault_u in distances.items()
    }


def probability_graph(
    distances: DistanceGraph, d0: float = 3, delta: float = 1
) -> ProbabilityGraph:
    """
    Convert a graph of distances between faults into a graph of jump
    probablities using the Shaw-Dieterich model.

    Parameters
    ----------
    distances : DistanceGraph
        The distance graph between faults.
    d0 : float, optional
        The d0 parameter for the Shaw_Dieterich model. See `shaw_dieterich_distance_model`.
    delta : float, optional
        The delta parameter for the Shaw_Dieterich model. See `shaw_dieterich_distance_model`.

    Returns
    -------
    ProbabilityGraph
        The graph with faults as vertices. Each edge (fault_u, fault_v)
        has a log-probability -p as a weight. The log-probability -p here
        is the negative of the log-probability a rupture propogates from
        fault_u to fault_v, relative to the probability it propogates to
        any of the other neighbours of fault_u.
    """
    probabilities_raw = {
        fault_u: {
            fault_v: shaw_dieterich_distance_model(distance, d0, delta)
            for fault_v, distance in neighbours_fault_u.items()
        }
        for fault_u, neighbours_fault_u in distances.items()
    }

    probabilities_log = defaultdict(dict)
    for fault_u, neighbours_fault_u in probabilities_raw.items():
        normalising_constant = sum(neighbours_fault_u.values())
        for fault_v, prob in neighbours_fault_u.items():
            probabilities_log[fault_u][fault_v] = normalising_constant - np.log(prob)
    return probabilities_log


def probabilistic_minimum_spanning_tree(
    probability_graph: ProbabilityGraph, initial_fault: str
) -> RuptureCausalityTree:
    """
    Generate a probabilistic minimum spanning tree.

    The minimum spanning tree of the probability graph represents rupture
    causality tree with the highest likelihood of occuring assuming that
    rupture jumps are independent.

    NOTE: While the overall probability of the minimum spanning tree is high,
    the paths from the initial fault may not be the most likely.

    Parameters
    ----------
    probability_graph : ProbabilityGraph
        The probability graph.
    initial_fault : str
        The initial fault.

    Returns
    -------
    RuptureCausalityTree
        The probabilistic minimum spanning tree.
    """
    rupture_causality_tree = {initial_fault: None}
    path_probabilities = defaultdict(lambda: np.inf)
    path_probabilities[initial_fault] = 0
    processed_faults = set()
    for _ in range(len(probability_graph)):
        current_fault = min(
            (fault for fault in probability_graph if fault not in processed_faults),
            key=lambda fault: path_probabilities[fault],
        )
        processed_faults.add(current_fault)
        for fault_neighbour, log_probability in probability_graph[
            current_fault
        ].items():
            if (
                fault_neighbour not in processed_faults
                and log_probability < path_probabilities[fault_neighbour]
            ):
                path_probabilities[fault_neighbour] = log_probability
                rupture_causality_tree[fault_neighbour] = current_fault
    return rupture_causality_tree


def probabilistic_shortest_path(
    probability_graph: ProbabilityGraph, initial_fault: str
) -> RuptureCausalityTree:
    """
    Given a probability graph, return the shortest path tree from the initial fault.

    The interpretation of this function in a probability theory context is:
    each path from the initial fault to all subsequent faults is the most
    likely.

    NOTE: While the paths from the initial fault to all subsequent faults are
    the most likely, the overall rupture causality may not be the most likely.

    Parameters
    ----------
    probability_graph : ProbabilityGraph
        The probability graph.
    initial_fault : str
        The initial fault.

    Returns
    -------
    RuptureCausalityTree
        The probabilistic shortest path tree.
    """
    rupture_causality_tree = {initial_fault: None}
    path_probabilities = defaultdict(lambda: np.inf)
    path_probabilities[initial_fault] = 0
    processed_faults = set()
    while len(processed_faults) < len(probability_graph):
        current_fault = min(
            (fault for fault in probability_graph if fault not in processed_faults),
            key=lambda fault: path_probabilities[fault],
        )
        for fault_neighbour, log_probability in probability_graph[
            current_fault
        ].items():
            candidate_path_probability = (
                path_probabilities[current_fault] + log_probability
            )
            if candidate_path_probability < path_probabilities[fault_neighbour]:
                path_probabilities[fault_neighbour] = candidate_path_probability
                rupture_causality_tree[fault_neighbour] = current_fault

        processed_faults.add(current_fault)
    return rupture_causality_tree
