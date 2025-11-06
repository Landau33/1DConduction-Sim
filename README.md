# 1D Fin Steady-State Heat Conduction Simulation

This project simulates the **steady-state temperature distribution on a 1D circular fin** with convection at its surface.  
The fin has temperature-dependent thermal conductivity \(k(T)\), using **two-point linear interpolation**.

This is a simple and clear numerical model suitable for thermal engineering coursework, numerical methods practice, and fin performance analysis.
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
    error=1e-2,          # Convergence tolerance for residual
    delta=3e-2,          # Residual descent step
    max_step=100000,     # Maximum iteration steps
    print_step=1000      # Print residual every N steps
)
```
The temperature-dependent thermal conductivity is defined inside the function:

```python
def conductivity(T: float) -> float:
    return (1017.0/2800.0) * T + (85.0/28.0)
```

---

## Usage

Run the simulation:

```bash
python main.py
```