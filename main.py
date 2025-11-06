"""
main.py

启动仿真并绘制温度分布
Run the simulation and generate the temperature distribution plot

Author: Yuang
Email: chris14658@naver.com
"""

import numpy as np
import matplotlib.pyplot as plt
from fin import FinParams, init_nodes, thermal_conductivity
from solver import energy_balance


def main():
    p = FinParams(length=0.1, total_node=50, T0=300.0, Ta=20.0, hc=100.0, D=0.005,
                  error=1e-2, delta=1e-2, max_step=100000, print_step=10000)
    nodes = init_nodes(p)
    energy_balance(nodes, p, k_of_T=thermal_conductivity)

    x = np.array([node.x for node in nodes])
    T = np.array([node.T for node in nodes])

    step = p.total_node // 10
    for i in range(0, p.total_node, step):
        print(f"node {i}, x={x[i]:.6f} m, T={T[i]:.2f} K")

    plt.plot(x, T, color='blue', linewidth=0.5, marker='o', 
                markersize=0.5, markerfacecolor='black', markeredgecolor='black')
    plt.text(0.98, 0.95, f"{p.total_node-1:.0f}nodes, energy error={p.error:.2e}",
         ha='right', va='top',
         transform=plt.gca().transAxes)
    plt.xlabel("Position along the fin (m)")
    plt.ylabel("Temperature (K)")
    plt.title("Temperature Distribution")
    plt.grid(True)
    plt.show()


if __name__ == "__main__":
    main()