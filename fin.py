from dataclasses import dataclass
from math import pi
from typing import List

@dataclass
class Node:
    """一维鳍片的节点：序号、位置、温度"""
    idx: int
    x: float
    T: float

@dataclass
class FinParams:
    """鳍片物理与数值参数"""
    length: float = 0.1         # m, fin长度
    total_node: int = 50        # 节点数
    T0: float = 300.0           # K, 基部温度 (x=0) Dirichlet
    Ta: float = 20.0            # K, 环境温度
    hc: float = 100.0           # W/m^2-K, 对流换热系数
    D: float = 0.005            # m, 直径
    error: float = 1e-2         # 能量残差容差
    delta: float = 1e-2         # 残差下降步长
    max_step: int = 100_000
    print_step: int = 2000

    def __post_init__(self):
        self.total_node += 1

    def geometry(self):
        dx = self.length / (self.total_node - 1)
        perimeter  = pi * self.D
        area  = pi * (self.D**2) / 4.0
        return dx, perimeter, area

def thermal_conductivity(T: float) -> float:
    return (1017.0/2800.0)*T + (85.0/28.0)

def init_nodes(p: FinParams) -> List[Node]:
    """init nodes with boundary condition at x=0"""
    dx, _, _ = p.geometry()
    nodes: List[Node] = []
    
    for i in range(p.total_node):
        x = i * dx
        T = p.Ta     # 初始温度全部设为Ta
        nodes.append(Node(idx=i, x=x, T=T))
    
    nodes[0].T = p.T0   # 保留边界条件：基部温度固定为T0
    return nodes