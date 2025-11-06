Version: v0.1.1
# 1D Fin Steady-State Heat Conduction Simulation

***This project was developed as an MA6801 course assignment at NTU for Fall 2025 (Semester 1).***

This project simulates the *steady-state temperature distribution on a 1D circular fin* with convection at its surface.  
The fin has temperature-dependent thermal conductivity `k(T)`, using two-point linear interpolation.

![Temperature Distribution](question.png)

---

## Mathematical Model

The fin is modeled by the energy balance at each control volume:

$$
Q_{\text{left}} + Q_{\text{right}} + Q_{\text{conv}} = 0
$$

where:

- Conduction:

$$
Q_{\text{left}} = k \frac{A}{dx}(T_{i-1} - T_i)
$$

$$
Q_{\text{right}} = k \frac{A}{dx}(T_{i+1} - T_i)
$$

- Convection:

$$
Q_{\text{conv}} = h P dx (T_a - T_i)
$$

Thermal conductivity is interpolated from two user-specified points:

$$
k(T) = k_1 + (k_2 - k_1)\frac{T - T_1}{T_2 - T_1}
$$

---

## Installation

```bash
pip install -r requirements.txt
```

---

## Project Structure

The simulation parameters are defined inside `main.py` through the `FinParams` class.  
Below is the default configuration used in this project:

```python
FinParams(
    length=0.1,          # Fin length (m)
    total_node=50,       # Number of grid nodes
    T0=300.0,            # Base temperature (K)
    Ta=20.0,             # Ambient temperature (K)
    hc=100.0,            # Convection coefficient (W/m²·K)
    D=0.005,             # Fin diameter (m)
    error=1e-2,          # Convergence threshold for energy balance
    delta=3e-2,          # Convergence factor
    max_step=100000,     # Maximum iteration steps
    print_step=10000     # Print interval
```
The temperature-dependent thermal conductivity is defined in `FinParams` class:

```python
    kt_slope=(112.0-10.3)/(300.0-20.0)  # Slope dk/dT for thermal conductivity
    kt_bias=10.3-kt_slope*20.0          # Bias for thermal conductivity
)           
```

---

## Usage

The program supports three solvers: `residual`, `picard`, and `newton` (default: residual).

Run the simulation:

```bash
python main.py --solver <solver>
```

Examples:

```bash
python main.py --solver residual
python main.py --solver picard
python main.py --solver newton
```

---
## Expected Output

```
Using Residual Solver...
[Residual] step 10000, max_energy=0.176659
[Residual] step 20000, max_energy=0.039235
[Residual] step 30000, max_energy=0.011394
[Residual] finished at step 31139
node 0, x=0.00 m, T=300.00 K
node 5, x=0.01 m, T=241.05 K
node 10, x=0.02 m, T=189.37 K
node 15, x=0.03 m, T=144.92 K
node 20, x=0.04 m, T=107.66 K
node 25, x=0.05 m, T=77.51 K
node 30, x=0.06 m, T=54.36 K
node 35, x=0.07 m, T=38.02 K
node 40, x=0.08 m, T=27.99 K
node 45, x=0.09 m, T=23.07 K
node 50, x=0.10 m, T=21.57 K
```