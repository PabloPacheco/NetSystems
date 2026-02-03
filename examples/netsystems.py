import time
import numpy as np
from scipy import sparse
import scipy.sparse.linalg as splinalg
from scipy.optimize import minimize, least_squares, newton_krylov

import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
import pandas as pd
import seaborn as sns
from copy import deepcopy
from CoolProp.CoolProp import PropsSI


class Network:
    def __init__(self):
        # -------------------------
        # Topología y geometría
        # -------------------------
        self.n_nodes = 0
        self.n_elements = 0

        self.connectivity = None          # (n_elements, 2)
        self.node_coordinates = None      # (n_nodes, dim)
        self.element_lengths = None

        # -------------------------
        # Sistemas físicos
        # -------------------------
        self.systems = {}
        self.iteration = 0

    def calculate_element_lengths(self):
        self.element_lengths = np.zeros(self.n_elements)

        for e in range(self.n_elements):
            i, j = self.connectivity[e].astype(int)
            xi = self.node_coordinates[i]
            xj = self.node_coordinates[j]
            self.element_lengths[e] = np.linalg.norm(xj - xi)

    def set_connectivity(self, c):
        self.connectivity = c
        self.n_nodes = int(np.max(c) + 1)
        self.n_elements = c.shape[0]

    def add_system(self, system):
        system.network = self
        self.systems[system.name] = system

    def system(self, name):
        return self.systems[name]


    def get(self, var_name):
        matches = []
        for system in self.systems.values():
            if var_name in (system.x_name, system.b_name, system.element_variable_name):
                matches.append(system)

        if len(matches) == 1:
            return matches[0].get(var_name)
        elif len(matches) > 1:
            raise ValueError(
                f"Variable '{var_name}' es ambigua "
                f"({[s.name for s in matches]})"
            )
        else:
            raise KeyError(f"Variable '{var_name}' no encontrada en la red.")

    def __getitem__(self, var_name):
        return self.get(var_name)


class System:
    def __init__(self, name, system_type = "diffusive"):
        self.name = name
        self.type = system_type  # "diffusive" | "advective"

        self.network = None

        # -------------------------
        # Variables matemáticas
        # -------------------------
        self.x = None    # incógnita nodal
        self.b = None    # fuente / demanda nodal
        self.A = None
        self.initial_x = None

        # -------------------------
        # Metadatos físicos
        # -------------------------
        self.x_name = None
        self.b_name = None
        self.element_variable_name = None

        # -------------------------
        # Condiciones de borde
        # -------------------------
        self.known_x_nodes = np.array([], dtype=int)
        self.known_x_values = np.array([])

        self.known_b_nodes = np.array([], dtype=int)
        self.known_b_values = np.array([])

        # -------------------------
        # Propiedades de elementos
        # -------------------------
        self.element_conductance = None
        self.element_conductance_function = None

        self.element_flux = None
        self.element_flux_function = None

        self.nodal_gamma_function = None
        self.nodal_gamma = None

        # Variables de línea (ej: caudal)
        self.element_variable = None
        self.linear_solver = "cg"
        self.linear_solver_tolerance = 1e-6

    def get(self, var_name):
        if var_name == self.x_name:
            return self.x
        elif var_name == self.b_name:
            return self.b
        elif var_name == self.element_variable_name:
            return self.element_variable
        else:
            raise KeyError(
                f"Variable '{var_name}' no existe en el sistema '{self.name}'."
            )

    def __getitem__(self, var_name):
        return self.get(var_name)


def solve_linear_system(network, system_name):
    """
    Resuelve un sistema lineal nodal usando partición.
    """
    system = network.systems[system_name]

    n_nodes = network.n_nodes

    # 1. Inicializar vectores
    nodal_b = np.zeros(n_nodes)
    nodal_x = np.zeros(n_nodes)

    # 2. Aplicar fuentes conocidas
    nodal_b[system.known_b_nodes] = system.known_b_values

    # 3. Máscara de nodos con x conocido
    known_x_mask = np.zeros(n_nodes, dtype=bool)
    known_x_mask[system.known_x_nodes] = True

    # 4. Ensamblar matriz global
    global_matrix = assemble_global_matrix(network, system_name)

    # 5. Extraer particiones
    (K_uu, K_uk, K_kk,
     unknown_indices, known_indices) = extract_partitioned_matrices(
        global_matrix, known_x_mask)

    # 6. Vectores particionados
    b_u = nodal_b[unknown_indices]
    x_k = system.known_x_values.reshape(-1, 1)

    # 7. Sistema reducido
    rhs = b_u.reshape(-1, 1) - K_uk @ x_k

    # 8. Solver
    solver_fn = get_linear_solver_function(system.linear_solver)

    if system.linear_solver == "spsolve":
        x_u = solver_fn(K_uu, rhs.ravel())
        solver_info = 0
    else:
        K_uu = K_uu.tocsr() if hasattr(K_uu, "tocsr") else K_uu
        x_u, solver_info = solver_fn(
            K_uu, rhs, atol=system.linear_solver_tolerance
        )

    x_u = x_u.reshape(-1)

    # 9. Reconstruir solución completa
    nodal_x[system.known_x_nodes] = system.known_x_values
    nodal_x[unknown_indices] = x_u

    # 10. Calcular b en nodos con x conocido
    b_k = (K_kk @ x_k + K_uk.T @ x_u.reshape(-1, 1)).ravel()
    nodal_b[known_indices] = b_k

    # 11. Postproceso elemental (usa conectividad de la red)
    i_start = network.connectivity[:, 0]
    i_end = network.connectivity[:, 1]

    x_diff = nodal_x[i_start] - nodal_x[i_end]
    element_variable = system.element_conductance * x_diff

    # 12. Almacenar resultados
    system.x = nodal_x
    system.b = nodal_b
    system.x_difference = x_diff
    system.element_variable = element_variable

    return solver_info



def assemble_global_matrix(network, system_name):
    """
    Ensambla la matriz global nodal del sistema.
    """
    system = network.systems[system_name]

    n_nodes = network.n_nodes
    connectivity = network.connectivity

    global_matrix = sparse.lil_matrix((n_nodes, n_nodes))

    for i in range(network.n_elements):
        k = float(system.element_conductance[i])
        n1, n2 = connectivity[i]

        global_matrix[n1, n1] += k
        global_matrix[n1, n2] -= k
        global_matrix[n2, n1] -= k
        global_matrix[n2, n2] += k

    return global_matrix.tocsr()


def extract_partitioned_matrices(global_matrix, known_x_mask):
    """Extrae submatrices particionadas."""
    known = np.where(known_x_mask)[0]
    unknown = np.where(~known_x_mask)[0]

    K_uu = global_matrix[unknown, :][:, unknown]
    K_uk = global_matrix[unknown, :][:, known]
    K_kk = global_matrix[known, :][:, known]

    return K_uu, K_uk, K_kk, unknown, known


def get_linear_solver_function(solver_name):
    solver_map = {
        'spsolve': splinalg.spsolve,
        'cg': splinalg.cg,
        'bicg': splinalg.bicg,
        'bicgstab': splinalg.bicgstab,
        'cgs': splinalg.cgs,
        'gmres': splinalg.gmres,
        'minres': splinalg.minres,
        'lgmres': splinalg.lgmres,
        'qmr': splinalg.qmr,
        'gcrotmk': splinalg.gcrotmk,
        'tfqmr': splinalg.tfqmr
    }
    return solver_map.get(solver_name, splinalg.cg)


def friction_factor(Re, D, eps):
    """
    Vectorized Darcy friction factor.

    Parameters
    ----------
    Re : float or array_like
        Reynolds number.
    D : float or array_like
        Pipe diameter [m].
    eps : float or array_like
        Absolute roughness [m].

    Returns
    -------
    f : ndarray
        Darcy friction factor.
    """

    # Convert to arrays (enables broadcasting)
    Re  = np.asarray(Re, dtype=float)
    D   = np.asarray(D, dtype=float)
    eps = np.asarray(eps, dtype=float)

    f = np.empty_like(Re, dtype=float)

    # Regime masks
    mask_lam = Re < 2000.0
    mask_tur = Re > 4000.0
    mask_trn = (~mask_lam) & (~mask_tur)

    # -------------------
    # Laminar
    # -------------------
    if np.any(mask_lam):
        f[mask_lam] = 64.0 / Re[mask_lam]

    # -------------------
    # Turbulent (Swamee–Jain, log10)
    # -------------------
    if np.any(mask_tur):
        aa = eps[mask_tur] / (3.7 * D[mask_tur]) \
           + 5.74 / (Re[mask_tur]**0.9)

        f[mask_tur] = 0.25 / (np.log10(aa)**2)

    # -------------------
    # Transitional
    # -------------------
    if np.any(mask_trn):
        Re_t  = Re[mask_trn]
        D_t   = D[mask_trn]
        eps_t = eps[mask_trn]

        R = Re_t / 2000.0

        Y2 = eps_t / (3.7 * D_t) + 5.74 / (Re_t**0.9)
        Y3 = -0.86859 * np.log(
            eps_t / (3.7 * D_t) + 5.74 / (4000.0**0.9)
        )

        FA = Y3**(-2)
        FB = FA * (2.0 - 0.00514215 / (Y2 * Y3))

        X1 = 7.0 * FA - FB
        X2 = 0.128 - 17.0 * FA + 2.5 * FB
        X3 = -0.128 + 13.0 * FA - 2.0 * FB
        X4 = R * (0.032 - 3.0 * FA + 0.5 * FB)

        f[mask_trn] = X1 + R * (X2 + R * (X3 + X4))

    return f


def calculate_nonlinear_fluid_property(network):
    system = network.systems["fluid"]

    heads = system.x  # ya resuelto o estimado

    i = network.connectivity[:, 0].astype(int)
    j = network.connectivity[:, 1].astype(int)

    delta_h = heads[i] - heads[j]
    abs_delta_h = np.abs(delta_h) + 1e-12

    diameters = system.element_diameter[:network.n_elements]
    C = getattr(network, "chw", 110.0)

    if not hasattr(network, "element_lengths"):
        network.calculate_element_lengths_from_coordinates()

    k = (
        0.2784
        * C
        * diameters ** 2.63
        / network.element_lengths ** 0.54
        / abs_delta_h ** 0.46
    )

    return k


def get_initial_guess(network, system_name):
    system = network.systems[system_name]

    n_nodes = network.n_nodes

    # 1. Conductividad inicial
    system.element_conductance = np.ones(network.n_elements)

    # 2. Ensamblaje
    global_matrix = assemble_global_matrix(network, system_name)

    # 3. Vector b
    b = np.zeros(n_nodes)
    b[system.known_b_nodes] = system.known_b_values

    # 4. Máscara
    known_mask = np.zeros(n_nodes, dtype=bool)
    known_mask[system.known_x_nodes] = True

    u = np.where(~known_mask)[0]
    k = np.where(known_mask)[0]

    Auu = global_matrix[u][:, u]
    Auk = global_matrix[u][:, k]

    rhs = b[u] - Auk @ system.known_x_values

    Hu = splinalg.spsolve(Auu, rhs)

    # 5. Reconstrucción
    x = np.zeros(n_nodes)
    x[k] = system.known_x_values
    x[u] = Hu

    system.x = x

    # 6. Flujos iniciales
    i, j = network.connectivity[:, 0], network.connectivity[:, 1]
    system.element_variable = system.element_conductance * (x[i] - x[j])

    return Hu


def calculate_residual(network, system_name, unknown_x):
    system = network.systems[system_name]

    n_nodes = network.n_nodes

    known_mask = np.zeros(n_nodes, dtype=bool)
    known_mask[system.known_x_nodes] = True

    u = np.where(~known_mask)[0]
    k = np.where(known_mask)[0]

    x = np.zeros(n_nodes)
    x[k] = system.known_x_values
    x[u] = unknown_x

    system.x = x

    # 1. Propiedad no lineal
    system.element_conductance = system.element_conductance_function(network)

    # 2. Ensamblaje
    A = assemble_global_matrix(network, system_name)

    Auu = A[u][:, u]
    Auk = A[u][:, k]

    b = np.zeros(n_nodes)
    b[system.known_b_nodes] = system.known_b_values

    r = Auu @ unknown_x + Auk @ system.known_x_values - b[u]

    return np.asarray(r).ravel()


def finalize_solution(network, system_name, unknown_x):
    system = network.systems[system_name]

    n_nodes = network.n_nodes

    known_mask = np.zeros(n_nodes, dtype=bool)
    known_mask[system.known_x_nodes] = True

    u = np.where(~known_mask)[0]
    k = np.where(known_mask)[0]

    x = np.zeros(n_nodes)
    x[k] = system.known_x_values
    x[u] = unknown_x
    system.x = x

    system.element_conductance = system.element_conductance_function(network)

    A = assemble_global_matrix(network, system_name)

    b = np.zeros(n_nodes)
    b[system.known_b_nodes] = system.known_b_values

    bk = A[k][:, k] @ system.known_x_values + A[k][:, u] @ unknown_x
    b[k] = bk

    system.b = b

    i, j = network.connectivity[:, 0], network.connectivity[:, 1]

    system.element_variable = system.element_conductance * (x[i] - x[j])
    


def solve_least_squares_nonlinear(
    network,
    system_name,
    initial_guess,
    solver_kwargs=None,
    finalize=True
):
    """
    Resuelve el sistema no lineal usando scipy.optimize.least_squares.
    """

    if solver_kwargs is None:
        solver_kwargs = dict(
            method="trf",
            xtol=1e-8,
            ftol=1e-8,
            gtol=1e-8,
            max_nfev=1000
        )

    result = least_squares(
        fun=lambda x: calculate_residual(network, system_name, x),
        x0=np.asarray(initial_guess, dtype=float),
        **solver_kwargs
    )

    if finalize:
        finalize_solution(network, system_name, result.x)

    return result


def approximate_jacobian(function, vector, step_size=1e-6):
    """
    Aproxima el jacobiano usando diferencias finitas.
    """
    vector = np.asarray(vector, dtype=float)
    f0 = function(vector)

    n = len(vector)
    m = len(f0)

    J = np.zeros((m, n))

    for i in range(n):
        v = vector.copy()
        v[i] += step_size
        fi = function(v)
        J[:, i] = (fi - f0) / step_size

    return J


def solve_newton_raphson(
    network,
    system_name,
    initial_guess,
    tolerance=1e-8,
    max_iterations=50,
    verbose=True,
    finalize=True
):
    """
    Resuelve el sistema no lineal con Newton-Raphson
    usando jacobiano numérico.
    """

    x = np.asarray(initial_guess, dtype=float)

    for it in range(max_iterations):
        r = calculate_residual(network, system_name, x)
        norm_r = np.linalg.norm(r)

        if verbose:
            print(f"[Newton] Iter {it:02d} |R| = {norm_r:.3e}")

        if norm_r < tolerance:
            break

        J = approximate_jacobian(
            lambda y: calculate_residual(network, system_name, y),
            x
        )

        try:
            dx = np.linalg.solve(J, -r)
        except np.linalg.LinAlgError:
            dx, *_ = np.linalg.lstsq(J, -r, rcond=None)

        x += dx

    if finalize:
        finalize_solution(network, system_name, x)

    return x


def solve_newton_krylov_nonlinear(
    network,
    system_name,
    initial_guess,
    method="lgmres",
    tolerance=1e-8,
    max_iterations=200,
    finalize=True
):
    """
    Resuelve el sistema no lineal usando Newton-Krylov.
    """

    sol = newton_krylov(
        lambda x: calculate_residual(network, system_name, x),
        xin=np.asarray(initial_guess, dtype=float),
        method=method,
        maxiter=max_iterations,
        f_tol=tolerance
    )

    if finalize:
        finalize_solution(network, system_name, sol)

    return sol


def solve_minimize_nonlinear(
    network,
    system_name,
    initial_guess,
    method="BFGS",
    tol=1e-6,
    maxiter=1000,
    bounds=None,
    bounds_factor=6.0,
    min_abs_value=1.0,
    fallback_method="L-BFGS-B",
    verbose=False,
    finalize=True
):
    """
    Minimiza ||R(x)||_2 usando scipy.optimize.minimize
    """

    x0 = np.asarray(initial_guess, dtype=float)

    # -------------------------------------------------
    # Bounds automáticos ± factor · |x|
    # -------------------------------------------------
    if bounds is None:
        bounds = []
        for xi in x0:
            ref = max(abs(xi), min_abs_value)
            bounds.append((-bounds_factor * ref, bounds_factor * ref))

    methods_with_bounds = {
        "L-BFGS-B", "TNC", "SLSQP", "trust-constr"
    }

    if bounds is not None and method not in methods_with_bounds:
        if verbose:
            print(f"[INFO] '{method}' no soporta bounds → usando '{fallback_method}'")
        method = fallback_method

    # -------------------------------------------------
    # Función objetivo
    # -------------------------------------------------
    def objective(x):
        r = calculate_residual(network, system_name, x)
        return np.linalg.norm(r)

    result = minimize(
        fun=objective,
        x0=x0,
        method=method,
        bounds=bounds,
        tol=tol,
        options=dict(maxiter=maxiter, disp=verbose)
    )

    if finalize and result.success:
        finalize_solution(network, system_name, result.x)

    return result


def assemble_advective_system(network, system_name, flow):
    """
    Ensambla la matriz y el vector de un sistema advectivo nodal
    usando esquema upwind, con propiedades definidas en nodos.

    Args:
        network : Network
        system_name : str
        flow : array (n_elements,)  # caudal volumétrico por elemento

    Returns:
        A : csr_matrix
        b : ndarray
    """
    system = network.systems[system_name]

    n_nodes = network.n_nodes
    n_elem  = network.n_elements

    A = sparse.lil_matrix((n_nodes, n_nodes))
    b = np.zeros(n_nodes)

    # -------------------------------------------------
    # Element heat
    # -------------------------------------------------
    if system.element_flux_function is not None:
        system.element_flux = system.element_flux_function(net)

    # -------------------------------------------------
    # γ nodal
    # -------------------------------------------------
    if system.nodal_gamma_function is not None:
        gamma_nodal = system.nodal_gamma_function(network)
    else:
        gamma_nodal = system.nodal_gamma

    gamma_nodal = np.asarray(gamma_nodal, dtype=float)

    # -------------------------------------------------
    # Loop sobre elementos
    # -------------------------------------------------
    for e in range(n_elem):
        i, j = network.connectivity[e].astype(int)

        # γ elemental = promedio nodal
        gamma_e = 0.5 * (gamma_nodal[i] + gamma_nodal[j])

        Q = flow[e]
        adv_coeff = abs(Q) * gamma_e

        # -------------------------------
        # Upwind
        # -------------------------------
        if Q >= 0.0:
            up, down = i, j
        else:
            up, down = j, i

        A[down, down] += adv_coeff
        A[down, up]   -= adv_coeff

        # -------------------------------
        # Fuente elemental (opcional)
        # -------------------------------
        if system.element_flux is not None:
            b[down] += system.element_flux[e]

    system.A = A.tocsr()
    system.b = b

    return system.A, system.b




def solve_advective_system(network, system_name, flow):
    """
    Resuelve un sistema advectivo nodal lineal.

    Args:
        network : Network
        system_name : str
        flow : array (n_elements,)
    """
    system = network.systems[system_name]

    # -------------------------
    # Ensamblaje
    # -------------------------
    A, b = assemble_advective_system(network, system_name, flow)

    n_nodes = network.n_nodes

    known_mask = np.zeros(n_nodes, dtype=bool)
    known_mask[system.known_x_nodes] = True

    u = np.where(~known_mask)[0]
    k = np.where(known_mask)[0]

    Auu = A[u][:, u]
    Auk = A[u][:, k]

    rhs = b[u] - Auk @ system.known_x_values

    # -------------------------
    # Resolución lineal
    # -------------------------
    xu = splinalg.spsolve(Auu, rhs)

    # -------------------------
    # Reconstrucción
    # -------------------------
    x = np.zeros(n_nodes)
    x[k] = system.known_x_values
    x[u] = xu

    system.x = x

    return x


def solve_coupled_sequential(
    net,
    diffusive_system_name,
    advective_system_name,
    max_iter=20,
    tol=1e-4,
    report=True,
    hydraulic_tol=1e-6,
    hydraulic_max_iter=100,
):
    """
    Sequentially coupled diffusive–advective solver.

    Order:
        1) Diffusive system (e.g. hydraulics)
        2) Advective system (e.g. heat)

    Parameters
    ----------
    net : Network
        Network with registered systems.
    diffusive_system_name : str
        Name of the diffusive system (e.g. "fluid").
    advective_system_name : str
        Name of the advective system (e.g. "heat").
    """

    diff = net.systems[diffusive_system_name]
    adv  = net.systems[advective_system_name]

    # -------------------------------------------------
    # Initialization
    # -------------------------------------------------
    net.iteration = 0

    if adv.initial_x is None:
        adv.x = np.ones(net.n_nodes) * 300.
    if adv.x is None:
        adv.x = adv.initial_x
    if diff.initial_x is None:
        diff.x = np.ones(net.n_nodes)
    if diff.x is None:
        diff.x = diff.initial_x
    

    # Aplicar condiciones de contorno térmicas
    T = np.zeros(net.n_nodes)
    T[adv.known_x_nodes] = adv.known_x_values

    for it in range(max_iter):
        T_old = T.copy()

        # =================================================
        # 1) DIFFUSIVE SYSTEM (e.g. fluid)
        # =================================================
        guess = get_initial_guess(net, diffusive_system_name)

        res = solve_newton_raphson(
            net,
            system_name = diffusive_system_name, 
            initial_guess = guess,
            tolerance=hydraulic_tol,
            max_iterations=hydraulic_max_iter,
            verbose=False
        )

        finalize_solution(net, diffusive_system_name, res)

        # =================================================
        # 2) ADVECTIVE SYSTEM (e.g. heat)
        # =================================================
        solve_advective_system(net, advective_system_name, diff.element_variable)

        # =================================================
        # Convergence check (thermal)
        # =================================================
        T = adv.x
        err = np.linalg.norm(T - T_old)

        if report:
            print(f"[Coupled] Iter {it:02d}: ΔT = {err:.3e}")

        if err < tol:
            break

        net.iteration += 1

    return 


def plot_network_geometry(
    network,
    show_nodes = True,
    show_elements = True,
    show_node_numbers = True,
    show_element_numbers = True, 
    node_size = 120,
    node_color = "black",
    line_color = "black",
    line_width = 2.0,
    figure_size = (6, 6),
    title = "Network geometry",
    axis_off = True,
    # -------- Saving options --------
    save_figure = False,
    figure_name = "network_geometry.png",
    dpi = 300,
    pad_inches = 0.05,
):
    """
    Simple visualization of network geometry (nodes and connectivity),
    without physical results or color maps.
    """

    coords = network.node_coordinates
    connectivity = network.connectivity

    n_nodes = network.n_nodes
    n_elements = network.n_elements

    fig, ax = plt.subplots(figsize=figure_size)
    ax.set_aspect("equal")

    if axis_off:
        ax.axis("off")

    # -----------------------------
    # Elements (lines)
    # -----------------------------
    if show_elements:
        for e in range(n_elements):
            i, j = connectivity[e]
            x0, y0 = coords[i]
            x1, y1 = coords[j]

            ax.plot(
                [x0, x1],
                [y0, y1],
                color=line_color,
                linewidth=line_width,
                zorder=1
            )

            if show_element_numbers:
                xm = 0.5 * (x0 + x1)
                ym = 0.5 * (y0 + y1)
                ax.text(
                    xm,
                    ym,
                    f"{e}",
                    fontsize=9,
                    ha="center",
                    va="center",
                    color="black",
                    zorder=3
                )

    # -----------------------------
    # Nodes
    # -----------------------------
    if show_nodes:
        ax.scatter(
            coords[:, 0],
            coords[:, 1],
            s=node_size,
            color=node_color,
            edgecolors="black",
            zorder=2
        )

        if show_node_numbers:
            for i in range(n_nodes):
                ax.text(
                    coords[i, 0],
                    coords[i, 1],
                    f"{i}",
                    fontsize=9,
                    ha="center",
                    va="center",
                    color="white",
                    zorder=4
                )

    if title is not None:
        ax.set_title(title)

    # -----------------------------
    # Save figure (optional)
    # -----------------------------
    if save_figure:
        plt.savefig(
            figure_name,
            dpi=dpi,
            bbox_inches="tight",
            pad_inches=pad_inches
        )

    plt.show()



def plot_network_results(
    network,
    element_flow_name,
    node_flow_name,
    node_values=None,
    element_values=None,
    node_cmap="coolwarm",
    element_cmap="viridis",

    show_node_numbers=True,
    show_element_numbers=True,
    show_flow_arrows=True,

    # --------- Element value labels ----------
    show_element_values=False,
    element_value_format="{:.3f}",
    element_value_fontsize=9,
    element_value_color="black",   # "black" or "element"
    element_value_bbox=True,

    # --------- Node demand labels ----------
    show_node_demands=False,
    node_demand_format="{:.3f}",
    node_demand_fontsize=9,
    node_demand_color="black",     # "black" or "node"
    node_demand_offset=(0.03, 0.03),
    node_demand_bbox=True,

    # --------- Geometry ----------
    node_size=120,
    line_width=3.0,
    arrow_scale=1.0,
    arrow_width=0.006,
    figure_size=(8, 6),

    # --------- Labels ----------
    node_colorbar_label="Node value",
    element_colorbar_label="Element value",
    title=None,
    axis_off=True,

    # --------- Saving options ----------
    save_figure=False,
    figure_name="figure.png",
    dpi=300,
    pad_inches=0.05,
):
    """
    Generic visualization of nodal and element results for NetworkSimulation,
    including element values and nodal demands (node_flows).
    """

    # -----------------------------
    # Basic data
    # -----------------------------
    coords = network.node_coordinates
    connectivity = network.connectivity

    n_nodes = network.n_nodes
    n_elements = network.n_elements

    flows = network.get(element_flow_name)
    node_flows = network.get(node_flow_name)

    # -----------------------------
    # Figure
    # -----------------------------
    fig, ax = plt.subplots(figsize=figure_size)
    ax.set_aspect("equal")

    if axis_off:
        ax.axis("off")

    # -----------------------------
    # Nodes
    # -----------------------------
    node_norm = Normalize(
        vmin=np.nanmin(node_values),
        vmax=np.nanmax(node_values)
    )

    node_scatter = ax.scatter(
        coords[:, 0],
        coords[:, 1],
        c=node_values,
        cmap=node_cmap,
        norm=node_norm,
        s=node_size,
        edgecolors="black",
        zorder=3
    )

    cbar_nodes = fig.colorbar(node_scatter, ax=ax, shrink=0.8)
    cbar_nodes.set_label(node_colorbar_label)

    if show_node_numbers:
        for i in range(n_nodes):
            ax.text(
                coords[i, 0],
                coords[i, 1],
                f"{i}",
                fontsize=9,
                ha="center",
                va="center",
                color="black",
                zorder=4
            )

    # -----------------------------
    # Node demand labels (node_flows)
    # -----------------------------
    if show_node_demands and node_flows is not None:
        for i in range(n_nodes):
            if abs(node_flows[i]) > 0.0:
                dx, dy = node_demand_offset

                txt_color = "black"
                if node_demand_color == "node":
                    txt_color = node_scatter.cmap(
                        node_norm(node_values[i])
                    )

                bbox = dict(
                    boxstyle="round,pad=0.25",
                    facecolor="white",
                    edgecolor="black",   # <<< borde negro para nodos
                    linewidth=0.8
                ) if node_demand_bbox else None

                ax.text(
                    coords[i, 0] + dx,
                    coords[i, 1] + dy,
                    node_demand_format.format(node_flows[i]),
                    fontsize=node_demand_fontsize,
                    ha="left",
                    va="center",
                    color=txt_color,
                    bbox=bbox,
                    zorder=5
                )

    # -----------------------------
    # Elements
    # -----------------------------
    elem_norm = Normalize(
        vmin=np.nanmin(element_values),
        vmax=np.nanmax(element_values)
    )
    elem_mapper = ScalarMappable(norm=elem_norm, cmap=element_cmap)

    for e in range(n_elements):
        i, j = connectivity[e]
        x0, y0 = coords[i]
        x1, y1 = coords[j]

        color = elem_mapper.to_rgba(element_values[e])

        # ---- Draw element
        if show_flow_arrows and flows is not None:
            if flows[e] >= 0:
                xs, ys = x0, y0
                dx, dy = x1 - x0, y1 - y0
            else:
                xs, ys = x1, y1
                dx, dy = x0 - x1, y0 - y1

            ax.quiver(
                xs, ys,
                dx, dy,
                angles="xy",
                scale_units="xy",
                scale=arrow_scale,
                width=arrow_width,
                color=color,
                zorder=2
            )
        else:
            ax.plot(
                [x0, x1],
                [y0, y1],
                color=color,
                linewidth=line_width,
                zorder=1
            )

        # ---- Element index
        if show_element_numbers:
            xm = 0.5 * (x0 + x1)
            ym = 0.5 * (y0 + y1)
            ax.text(
                xm,
                ym,
                f"{e}",
                fontsize=8,
                ha="center",
                va="center",
                color="black",
                zorder=4
            )

        # ---- Element value label
        if show_element_values:
            xm = 0.5 * (x0 + x1)
            ym = 0.5 * (y0 + y1)

            txt_color = "black"
            if element_value_color == "element":
                txt_color = color

            bbox = dict(
                boxstyle="round,pad=0.25",
                facecolor="white",
                edgecolor="white",   # <<< borde negro para nodos
                linewidth=0.1
            ) if element_value_bbox else None

            ax.text(
                xm,
                ym,
                element_value_format.format(element_values[e]),
                fontsize=element_value_fontsize,
                ha="center",
                va="center",
                color=txt_color,
                bbox=bbox,
                zorder=5
            )

    # -----------------------------
    # Element colorbar
    # -----------------------------
    cbar_elements = fig.colorbar(
        elem_mapper,
        ax=ax,
        orientation="horizontal",
        fraction=0.05,
        pad=0.05
    )
    cbar_elements.set_label(element_colorbar_label)

    if title is not None:
        ax.set_title(title)

    # -----------------------------
    # Save figure
    # -----------------------------
    if save_figure:
        plt.savefig(
            figure_name,
            dpi=dpi,
            bbox_inches="tight",
            pad_inches=pad_inches
        )

    plt.show()

def get_global_matrix_csr(network, system_name, flow=None):
    """
    Obtiene la matriz global ensamblada para un sistema dado.
    
    Args:
        network: Objeto Network con sistemas cargados
        system_name: Nombre del sistema
        flow: Campo de flujo (solo requerido para sistemas advectivos)
    
    Returns:
        Matriz global en formato CSR
    """
    if system_name not in network.systems:
        raise ValueError(f"Sistema '{system_name}' no encontrado en la red")
    
    system = network.systems[system_name]
    
    if system.type == "diffusive":
        # Sistema difusivo: ensamblar directamente
        return assemble_global_matrix(network, system_name)
    elif system.type == "advective":
        # Sistema advectivo: necesita campo de flujo
        if flow is None:
            raise ValueError(f"Sistema advectivo '{system_name}' requiere parámetro 'flow'")
        A, _ = assemble_advective_system(network, system_name, flow)
        return A
    else:
        raise ValueError(f"Tipo de sistema desconocido: {system.type}")


def spy_matrix(A, title="Matriz Global", figsize=(8, 6), markersize=8):
    """
    Visualiza una matriz usando spy() para mostrar el patrón de dispersión.
    
    Args:
        A: Matriz en formato CSR o densa
        title: Título del gráfico
        figsize: Tamaño de la figura
        markersize: Tamaño de los marcadores en spy
    """
    # Convertir a matriz densa si es dispersa
    if hasattr(A, 'toarray'):
        A_dense = A.toarray()
    else:
        A_dense = np.array(A)
    
    fig, ax = plt.subplots(figsize=figsize)
    
    # Usar spy para mostrar el patrón
    ax.spy(A_dense, markersize=markersize)
    ax.set_title(title, fontsize=14)
    #ax.set_xlabel('Índice de columna', fontsize=12)
    #ax.set_ylabel('Índice de fila', fontsize=12)
    
    # Añadir cuadrícula
    ax.grid(True, which='both', color='gray', linestyle='-', linewidth=0.5, alpha=0.3)
    
    plt.tight_layout()
    plt.show()
    
    return fig, ax


def plot_matrix_with_colors(A, title="Matriz Global", figsize=(8, 6), cmap='viridis'):
    """
    Visualiza una matriz con colores según los valores.
    
    Args:
        A: Matriz en formato CSR o densa
        title: Título del gráfico
        figsize: Tamaño de la figura
        cmap: Mapa de colores
    """
    # Convertir a matriz densa si es dispersa
    if hasattr(A, 'toarray'):
        A_dense = A.toarray()
    else:
        A_dense = np.array(A)
    
    fig, ax = plt.subplots(figsize=figsize)
    
    # Crear mapa de calor
    im = ax.imshow(np.abs(A_dense), cmap=cmap, aspect='auto')
    
    # Añadir barra de colores
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('Valor absoluto', fontsize=11)
    
    # Configurar título y etiquetas
    ax.set_title(title, fontsize=14)
    #ax.set_xlabel('Índice de columna', fontsize=12)
    #ax.set_ylabel('Índice de fila', fontsize=12)
    
    # Mostrar valores si la matriz es pequeña
    n_rows, n_cols = A_dense.shape
    if n_rows <= 10 and n_cols <= 10:
        # Mostrar todos los índices
        ax.set_xticks(range(n_cols))
        ax.set_yticks(range(n_rows))
        ax.set_xticklabels(range(n_cols))
        ax.set_yticklabels(range(n_rows))
        
        # Añadir valores en las celdas
        for i in range(n_rows):
            for j in range(n_cols):
                val = A_dense[i, j]
                if abs(val) > 1e-10:  # Solo mostrar valores significativos
                    color = 'white' if abs(val) > 0.7 * np.abs(A_dense).max() else 'black'
                    ax.text(j, i, f'{val:.2f}', ha='center', va='center', color=color, fontsize=9)
    
    plt.tight_layout()
    plt.show()
    
    return fig, ax