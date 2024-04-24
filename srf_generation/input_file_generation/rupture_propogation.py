import numpy as np
from collections import defaultdict
import random


DistanceGraph = dict[str, dict[str, int]]
ProbabilityGraph = defaultdict[str, dict[str, float]]

RuptureCausalityTree = dict[str, str | None]


def shaw_dieterich_distance_model(distance: float, d0, delta) -> float:
    return min(1, np.exp(-(distance - delta) / d0))


def prune_distance_graph(distances: DistanceGraph, cutoff: int) -> ProbabilityGraph:
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
    probabilities_raw = {
        fault_u: {
            fault_v: shaw_dieterich_distance_model(distance, d0, delta)
            for fault_v, distance in neighbours_fault_u.items()
        }
        for fault_u, neighbours_fault_u in distances.items()
    }
    print(probabilities_raw)

    probabilities_log = defaultdict(lambda: dict())
    for fault_u, neighbours_fault_u in probabilities_raw.items():
        normalising_constant = sum(neighbours_fault_u.values())
        for fault_v, prob in neighbours_fault_u.items():
            probabilities_log[fault_u][fault_v] = normalising_constant - np.log(prob)
    return probabilities_log


def probabilistic_minimum_spanning_tree(probability_graph, initial_fault: str):
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
