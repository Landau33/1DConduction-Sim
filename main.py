"""
main.py

启动仿真并绘制温度分布
Run the simulation and generate the temperature distribution plot

Author: Yuang
Email: chris14658@naver.com
"""

import numpy as np
import matplotlib.pyplot as plt
import argparse
from fin import FinParams, init_nodes
from solver import residual_solver, picard_thomas_solver, newton_solver


def parse_args():
    parser = argparse.ArgumentParser(description="1D Fin Heat Conduction Simulation")

    parser.add_argument(
        "--solver",
        type=str,
        choices=["residual", "picard", "newton"],
        default="residual",
        help="Select solver: residual | picard | newton (default: residual)"
    )

    return parser.parse_args()


def main():
    args = parse_args()
    params = FinParams(length=0.1, total_node=50, T0=300.0, Ta=20.0, hc=100.0, D=0.005,
                  error=1e-2, delta=3e-2, max_step=100000, print_step=10000)
    nodes = init_nodes(params)

    if args.solver == "residual":
        print("Using Residual Solver...")
        residual_solver(nodes, params)
    elif args.solver == "picard":
        print("Using Picard-Thomas Solver...")
        picard_thomas_solver(nodes, params)
    elif args.solver == "newton":
        print("Using Newton Solver...")
        newton_solver(nodes, params)
    
    x = np.array([node.x for node in nodes])
    T = np.array([node.T for node in nodes])

    if params.total_node > 10:
        step = params.total_node // 10
        for i in range(0, params.total_node, step):
            print(f"node {i}, x={x[i]:.6f} m, T={T[i]:.2f} K")

    plt.plot(x, T, color='blue', linewidth=0.5, marker='o', 
                markersize=0.5, markerfacecolor='black', markeredgecolor='black')
    plt.text(0.98, 0.95, f"{params.total_node-1:.0f}nodes, energy thershold={params.error:.2e}",
         ha='right', va='top',
         transform=plt.gca().transAxes)
    plt.xlabel("Position along the fin (m)")
    plt.ylabel("Temperature (K)")
    plt.title("Temperature Distribution")
    plt.grid(True)
    plt.show()


if __name__ == "__main__":
    main()