"""
solver.py

能量平衡迭代求解器
Iterative energy-balance solver for the fin temperature field

Author: Yuang
Email: chris14658@naver.com
"""

from typing import List, Callable
from fin import Node, FinParams


def energy_balance(nodes: List[Node],
            params: FinParams,
            k_of_T: Callable[[float], float]) -> int:

    dx, P, A = params.geometry()
    current_step = 0

    while True:
        current_step += 1
        tol_count = 0
        max_energy = 0.0

        # calculate energy balance for each node
        for i in range(1, params.total_node):
            Ti  = nodes[i].T
            Tm1 = nodes[i-1].T

            # heat conduction via left control surface
            k_left = k_of_T(0.5*(Tm1 + Ti))
            Q_left = k_left * A/dx * (Tm1 - Ti)

            # convection
            conv = params.hc * P * dx * (params.Ta - Ti)

            # for nodes i from 2 to N-1
            if i <= params.total_node - 2:
                # heat conduction via right control surface
                Tp1 = nodes[i+1].T
                k_right = k_of_T(0.5*(Ti + Tp1))
                Q_right = k_right * A/dx * (Tp1 - Ti)
                # energy balance for node i
                total_energy = Q_left + Q_right + conv
            else:
                # energy balance for node N
                total_energy = Q_left + conv

            # update temperature of node i
            Ti_new = Ti + params.delta * total_energy
            if Ti_new > params.T0: Ti_new = params.T0
            if Ti_new < params.Ta: Ti_new = params.Ta
            nodes[i].T = Ti_new

            if abs(total_energy) <= params.error:
                tol_count += 1
            if abs(total_energy) > max_energy:
                max_energy = abs(total_energy)

        if params.print_step and current_step % params.print_step == 0:
            print(f"step={current_step}, max_energy={max_energy:.6f}, pass={tol_count}/{params.total_node-1}")
            if tol_count == (params.total_node - 1) or current_step >= params.max_step:
                return current_step