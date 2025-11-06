"""
solver.py

能量平衡迭代求解器
Iterative energy-balance solvers for the fin temperature field

Author: Yuang
Email: chris14658@naver.com
"""

from typing import List, Tuple
from fin import Node, FinParams


def residual_solver(nodes: List[Node],
                    params: FinParams) -> int:
    """
    Residual-balancing point-wise updater.
    Math model (steady 1D fin with convection along the surface):

        d/dx( k(T) * A * dT/dx ) - h * P * (T - Ta) = 0

    Control-volume balance at node i (i = 1..N-2):
        Q_left  + Q_right + Q_conv = 0

    where
        Q_left  = k_{i-1/2} * A/dx * (T_{i-1} - T_i)
        Q_right = k_{i+1/2} * A/dx * (T_{i+1} - T_i)
        Q_conv  = h * P * dx * (Ta - T_i)

    This routine updates T_i in-place by:
        T_i^{new} = T_i + delta * (Q_left + Q_right + Q_conv)
    """
    dx, P, A = params.geometry()
    step = 0

    while True:
        step += 1
        tol_count = 0
        max_energy = 0.0

        # calculate energy balance for each node
        for i in range(1, params.total_node):
            Ti  = nodes[i].T
            Tm1 = nodes[i-1].T

            # heat conduction via left control surface
            k_left = params.kt_slope * (0.5 * (Tm1 + Ti)) + params.kt_bias
            Q_left = k_left * A/dx * (Tm1 - Ti)

            # convection
            conv = params.hc * P * dx * (params.Ta - Ti)

            # for nodes i from 2 to N-1
            if i <= params.total_node - 2:
                # heat conduction via right control surface
                Tp1 = nodes[i+1].T
                k_right = params.kt_slope * (0.5 * (Ti + Tp1)) + params.kt_bias
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

        if params.print_step and step % params.print_step == 0:
            print(f"[Residual] step {step}, max_energy={max_energy:.6f}, pass={tol_count}/{params.total_node-1}")

        # Stop either when all nodes pass or reaching max iterations
        if tol_count == (params.total_node - 1) or step >= params.max_step:
            print(f"[Residual] finish at step {step}")
            return step
            

def thomas_tridiagonal(a, b, c, d):
    """
    Solve a tridiagonal linear system A x = d using the Thomas algorithm.
    
    The matrix A has the form:
    
        [ b0  c0   0   0   ...        ]
        [ a1  b1  c1   0   ...        ]
        [  0  a2  b2  c2   ...        ]
        [          ...                ]
        [              a_{n-1}  b_{n-1} ]
    
    Parameters
    ----------
    a : Sub-diagonal coefficients (a[0] unused; a[i] multiplies x[i-1])
    b : Main diagonal coefficients (must be non-zero)
    c : Super-diagonal coefficients (c[n-1] unused; c[i] multiplies x[i+1])
    d : Right-hand side vector
    
    Returns
    -------
    x : Solution vector satisfying A x = d
    
    Algorithm
    ---------
    Forward elimination:
        Modify coefficients to eliminate sub-diagonal entries,
        storing:
            cp[i] = modified super-diagonal (c')
            dp[i] = modified RHS (d')
        so the system becomes upper triangular.

    Back substitution:
        Solve for x[n-1], then x[n-2] ... x[0].

    Numerical Notes
    ---------------
    - This algorithm runs in O(n) time.
    - Assumes A is diagonally dominant or nonsingular.
    """
    n = len(b)
    cp = [0.0]*n
    dp = [0.0]*n

    beta = b[0]
    cp[0] = c[0] / beta if n > 1 else 0.0
    dp[0] = d[0] / beta
    for i in range(1, n):
        beta = b[i] - a[i] * cp[i-1]
        cp[i] = (c[i] / beta) if (i < n-1) else 0.0
        dp[i] = (d[i] - a[i]*dp[i-1]) / beta

    x = [0.0]*n
    x[-1] = dp[-1]
    for i in range(n-2, -1, -1):
        x[i] = dp[i] - cp[i]*x[i+1]
    return x


def picard_thomas_solver(nodes: List[Node],
                         params: FinParams,) -> Tuple[int, float]:
    """
    Picard (fixed-point) linearization + Thomas tridiagonal solve.
    Idea:
        - Freeze k(T) at current iterate using interface temperatures,
          i.e., k_{i±1/2}^{(n)} = k( (T_i^{(n)} + T_{i±1}^{(n)}) / 2 ).
        - With k frozen, the discrete equations are linear tridiagonal:
              a_i T_{i-1}^{(n+1)} + b_i T_i^{(n+1)} + c_i T_{i+1}^{(n+1)} = d_i

    Discrete coefficients (for i = 1..N-2):
        a_i =  k_{i-1/2} * A / dx
        c_i =  k_{i+1/2} * A / dx
        b_i = -(a_i + c_i) - h * P * dx
        d_i = - h * P * dx * Ta

    Boundary:
        i = 0: Dirichlet T_0 = T0  ->  b_0 = 1, d_0 = T0
        i = N-1 (tip): no right conduction term in this approximation:
              a_{N-1} = k_{N-3/2} * A / dx
              b_{N-1} = -a_{N-1} - h * P * dx
              d_{N-1} = - h * P * dx * Ta

    Convergence check:
        || T^{(n+1)} - T^{(n)} ||_inf <= error
    """
    N = params.total_node
    dx, P, A = params.geometry()

    T = [nd.T for nd in nodes]

    for step in range(1, params.max_step + 1):
        # 1) Freeze k at interfaces using current T
        k_left  = [0.0] * N    # used for i>=1   (left faces)
        k_right = [0.0] * N    # used for i<=N-2 (right faces)
        for i in range(1, N):
            k_left[i] = params.kt_slope * (0.5 * (T[i-1] + T[i])) + params.kt_bias
        for i in range(0, N-1):
            k_right[i] = params.kt_slope * (0.5 * (T[i] + T[i+1])) + params.kt_bias

        # 2) Assemble tridiagonal system A*T_new = d
        a = [0.0] * N
        b = [0.0] * N
        c = [0.0] * N
        d = [0.0] * N

        # Dirichlet at i=0
        a[0] = 0.0
        b[0] = 1.0
        c[0] = 0.0
        d[0] = params.T0

        hPdx = params.hc * P * dx

        # Interior nodes
        for i in range(1, N-1):
            a[i] = k_left[i]  * A / dx
            c[i] = k_right[i] * A / dx
            b[i] = -(a[i] + c[i]) - hPdx
            d[i] = -hPdx * params.Ta

        # Tip control volume (no right conduction)
        a[N-1] = k_left[N-1] * A / dx
        b[N-1] = -a[N-1] - hPdx
        c[N-1] = 0.0
        d[N-1] = -hPdx * params.Ta

        # 3) Solve the linear tridiagonal system
        T_new = thomas_tridiagonal(a, b, c, d)

        # 4) Convergence check: infinity norm of temperature change
        max_energy = max(abs(T_new[i] - T[i]) for i in range(N))

        for i in range(1, N):  # i=0 is fixed
            if T_new[i] > params.T0: T_new[i] = params.T0
            if T_new[i] < params.Ta: T_new[i] = params.Ta

        T = T_new

        if params.print_step and step % params.print_step == 0:
            print(f"[Picard] step {step}, max_energy={max_energy:.6e}")

        if max_energy <= params.error:
            for i in range(N):
                nodes[i].T = T[i]
            print(f"[Picard-Thomas] finish at step {step}")
            return step, max_energy

    # Not converged within max_step: still write back last iterate
    for i in range(N):
        nodes[i].T = T[i]
    return params.max_step, max_energy


def newton_solver(nodes: List[Node],
                 params: FinParams,) -> Tuple[int, float]:
    """
    Newton–Raphson nonlinear solve with a tridiagonal Jacobian.
    We solve R(T) = 0, where for node i (i=1..N-2):

        R_i(T) = Q_left + Q_right + Q_conv
               = k_{i-1/2}(T) * A/dx * (T_{i-1} - T_i)
               + k_{i+1/2}(T) * A/dx * (T_{i+1} - T_i)
               + h * P * dx * (Ta - T_i)

    with interface conductivities evaluated at the face-average temperature:
        k_{i-1/2} = k( (T_{i-1} + T_i)/2 ),  k_{i+1/2} = k( (T_i + T_{i+1})/2 )

    Newton step:
        J(T^{(n)}) * dT = - R(T^{(n)})
        T^{(n+1)} = T^{(n)} + dT

    Jacobian entries (1D three-diagonal) for interior nodes:
      Let m = dk/dT be constant (two-point linear k(T) makes m constant).
      Define L-face: TL = (T_{i-1}+T_i)/2, k_left = k(TL)
              R-face: TR = (T_i+T_{i+1})/2, k_right = k(TR)

      dR_i/dT_{i-1}  =  (A/dx) * [ k_left + (m/2)*(T_{i-1} - T_i) ]

      dR_i/dT_i      =  (A/dx) * [ -k_left + (m/2)*(T_i - T_{i-1})
                                   -k_right - (m/2)*(T_{i+1} - T_i) ]
                        - h * P * dx

      dR_i/dT_{i+1}  =  (A/dx) * [ k_right + (m/2)*(T_{i+1} - T_i) ]

    Boundary:
        i = 0: Dirichlet  => R_0 = T_0 - T0,  J_00 = 1
        i = N-1: tip uses the same insulated approximation (no right conduction term).
    """

    N = params.total_node
    dx, P, A = params.geometry()
    hPdx = params.hc * P * dx

    T = [nd.T for nd in nodes]

    for step in range(1, params.max_step + 1):
        # Residual vector and tridiagonal Jacobian (a,b,c)
        R = [0.0] * N
        a = [0.0] * N   # J_{i,i-1}
        b = [0.0] * N   # J_{i,i}
        c = [0.0] * N   # J_{i,i+1}

        # Dirichlet at i=0: R0 = T0 - T0_target -> here write as (T[0]-T0)=0, J00=1
        R[0] = T[0] - params.T0
        a[0], b[0], c[0] = 0.0, 1.0, 0.0

        # Interior and tip rows
        for i in range(1, N):
            Ti  = T[i]
            Tm1 = T[i-1]

            # Left face: k_left
            k_left = params.kt_slope * 0.5 * (Tm1 + Ti) + params.kt_bias

            # Left contribution to residual
            Q_left = (k_left * A / dx) * (Tm1 - Ti)

            # Left derivatives (from chain rule on k(TL) with dTL/dT = 1/2)
            dRi_im1 = (A/dx) * (k_left + 0.5 * params.kt_slope * (Tm1 - Ti))
            dRi_i   = (A/dx) * (-k_left - 0.5 * params.kt_slope * (Tm1 - T_i)) if False else (A/dx) * (-k_left + 0.5 * params.kt_slope * (Ti - Tm1))
            # (The latter form is algebraically identical; we keep the '+' variant for readability.)
            dRi_ip1 = 0.0

            # Right face only if not at the tip
            Q_right = 0.0
            if i <= N - 2:
                Tp1 = T[i+1]
                k_right  = params.kt_slope * 0.5 * (Ti + Tp1) + params.kt_bias
                Q_right  = (k_right * A / dx) * (Tp1 - Ti)

                # Right derivatives
                dRi_ip1 += (A/dx) * (k_right + 0.5 * params.kt_slope * (Tp1 - Ti))
                dRi_i   += (A/dx) * (-k_right - 0.5 * params.kt_slope * (Tp1 - Ti))

            # Convection at node i: hPdx*(Ta - Ti)
            Qconv = hPdx * (params.Ta - Ti)
            R[i]  = Q_left + Q_right + Qconv

            # Derivative of convection wrt Ti: -hPdx
            dRi_i += -hPdx

            # Assemble Jacobian tri-diagonals
            a[i] = dRi_im1
            b[i] = dRi_i
            c[i] = dRi_ip1

        # Solve J dT = -R with Thomas algorithm
        rhs = [-ri for ri in R]
        dT  = thomas_tridiagonal(a, b, c, rhs)

        # Update and clamp
        T_new = [0.0] * N
        T_new[0] = params.T0
        for i in range(1, N):
            T_new[i] = T[i] + dT[i]
            if T_new[i] > params.T0: T_new[i] = params.T0
            if T_new[i] < params.Ta: T_new[i] = params.Ta

        max_energy = max(abs(T_new[i] - T[i]) for i in range(N))
        T = T_new

        if params.print_step and step % params.print_step == 0:
            print(f"[Newton] step {step}, max_energy={max_energy:.6e}")

        if max_energy <= params.error:
            for i in range(N):
                nodes[i].T = T[i]
            print(f"[Newton] finish at step {step}")
            return step, max_energy

    # Not converged within max_step
    for i in range(N):
        nodes[i].T = T[i]
    return params.max_step, max_energy