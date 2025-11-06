'''
fin.py

翅片参数与节点初始化
Fin parameters and node initialization

Author: Yuang
Email: chris14658@naver.com
'''

from dataclasses import dataclass
from math import pi
from typing import List


@dataclass
class Node:
    """
    Node information: index, position, temperature
    """
    idx: int
    x: float
    T: float


@dataclass
class FinParams:
    length: float = 0.1         # Fin length (m)
    total_node: int = 50        # Number of nodes
    T0: float = 300.0           # Base temperature (K)
    Ta: float = 20.0            # Ambient temperature (K)
    hc: float = 100.0           # Convection coefficient (W/m²·K)
    D: float = 0.005            # Fin diameter (m)
    error: float = 1e-2         # Convergence threshold for energy balance
    delta: float = 3e-2         # Convergence factor
    max_step: int = 100_000     # Maximum iteration steps
    print_step: int = 10000     # Print interval
    kt_slope: float = (112.0 - 10.3) / (300.0 - 20.0)   # Slope dk/dT for thermal conductivity
    kt_bias: float = 10.3 - kt_slope * 20.0             # Bias for thermal conductivity

    def __post_init__(self):
        self.total_node += 1

    def geometry(self):
        """
        Calculate dx, cross-sectional perimeter, cross-sectional area
        """
        dx = self.length / (self.total_node - 1)
        P  = pi * self.D
        A  = pi * (self.D**2) / 4.0
        return dx, P, A


def init_nodes(p: FinParams) -> List[Node]:
    """
    init nodes with boundary condition x=0, T=T0; others T=Ta
    """
    dx, _, _ = p.geometry()
    nodes: List[Node] = []
    
    # Initialize all nodes at Ta
    for i in range(p.total_node):
        x = i * dx
        T = p.Ta    
        nodes.append(Node(idx=i, x=x, T=T))
    
    # First node with temperature T0
    nodes[0].T = p.T0
    return nodes