## NetSystems - Fundamentals

Python program based on the code from the book:
- Brebbia, C. A., & Ferrante, A. J. (2013). *Computational hydraulics*. Butterworth-Heinemann.

<img src="BrebbiaBook.jpg" width="333">

NetSystems is a computational framework for solving systems of interconnected one-dimensional elements, as illustrated below:

<img src="netExample.svg" width="400">

Each one-dimensional element is mathematically represented as follows:

<img src="pipeElement.svg" width="300">

The fundamental relationship for each element is expressed through its local stiffness matrix:

$$
\begin{bmatrix}
b_k^i \\
b_j^i
\end{bmatrix}
= k^i
\begin{bmatrix}
1 & -1 \\
-1 & 1
\end{bmatrix}
\begin{bmatrix}
x_k^i \\
x_j^i
\end{bmatrix}
$$

where:
- $x_k^i$ and $x_j^i$ represent the nodal unknowns for element $i$
- $b_k^i$ and $b_j^i$ denote the corresponding nodal fluxes or source terms
- $k^i$ is the element conductivity (or stiffness) coefficient

A network of such one-dimensional elements consists of multiple nodes and connecting elements. By assembling individual element matrices according to their connectivity pattern, we obtain the global system matrix for the entire network. This global system takes the form:

$$ \mathbf{A} \mathbf{x} = \mathbf{b} $$

In many practical applications, the matrix $\mathbf{A}$ depends on the solution vector $\mathbf{x}$, leading to nonlinear problems. Net Simulation provides specialized solvers for such cases, which are represented as:

$$ \mathbf{A}(\mathbf{x}) \, \mathbf{x} = \mathbf{b} $$

where the global matrix $\mathbf{A}(\mathbf{x})$ varies as a function of the unknown vector $\mathbf{x}$.

### Solution Methodology

#### Linear Systems

The solution of linear systems employs the partition method. Depending on the physical nature of the problem and the boundary conditions, unknowns may reside either in vector $\mathbf{x}$ or vector $\mathbf{b}$. The global equation system is therefore partitioned accordingly:

$$ \mathbf{A} \mathbf{x} = \mathbf{b} \quad \rightarrow \quad 
\begin{bmatrix}
\mathbf{A}_E & \mathbf{A}_{EF} \\
\mathbf{A}_{EF}^T & \mathbf{A}_F
\end{bmatrix}
\begin{bmatrix}
\mathbf{x}_E \\
\mathbf{x}_F
\end{bmatrix}
= \begin{bmatrix}
\mathbf{b}_E \\
\mathbf{b}_F
\end{bmatrix}
$$

where:
- $\mathbf{x}_F$: nodes with unknown primary variables (to be solved)
- $\mathbf{x}_E$: nodes with prescribed primary variables (Dirichlet conditions)
- $\mathbf{b}_F$: nodes with known source terms (Neumann conditions)
- $\mathbf{b}_E$: nodes with unknown secondary variables (to be computed)

From this partitioned system, we derive the following relationship:

$$ \mathbf{A}_{EF}^T \mathbf{x}_E + \mathbf{A}_F \mathbf{x}_F = \mathbf{b}_F $$

Rearranging for solution yields:

$$\mathbf{A}_F \mathbf{x}_F = \mathbf{b}_F - \mathbf{A}_{EF}^T \mathbf{x}_E $$

This reduced system can be solved using various linear algebra techniques to obtain $\mathbf{x}_F$. Subsequently, the complete solution can be recovered by computing:

$$ \mathbf{b}_E = \mathbf{A}_E \mathbf{x}_E + \mathbf{A}_{EF} \mathbf{x}_F $$

#### Nonlinear Systems

For nonlinear problems where $\mathbf{A}$ depends on $\mathbf{x}$, Net Simulation implements several iterative solution strategies:
- **Newton-Raphson method** with numerical Jacobian approximation
- **Least-squares optimization** approaches
- **Newton-Krylov methods** for large-scale problems
- **Sequential coupling algorithms** for multi-physics applications

These methods handle the nonlinear dependence through successive linearization and iteration until convergence criteria are satisfied.

## Applications

The framework is particularly suited for:
- Pipe network hydraulics (water distribution systems)
- Thermal networks (district heating/cooling)
- Coupled fluid-thermal simulations
- Transport phenomena in network structures

# NetSystems: Modular Framework for Coupled Physical Systems Simulation

## Overview

NetSystems is a flexible and extensible Python framework for simulating complex, coupled physical systems. Built on a modular architecture, it enables modeling, solving, and analyzing multiple physical phenomena (fluid flow, heat transfer, diffusion, etc.) within arbitrary networks while maintaining clear separation between network topology and the governing physical systems.

## Core Philosophy

NetSystems emerges from the need to unify simulation of diverse physical phenomena under a common architecture. Instead of creating specific solvers for each problem, it provides a framework where:

1. **The network is independent** of the physical phenomena occurring within it
2. **Each system defines its own variables** (heads, temperatures, concentrations, etc.)
3. **Couplings are explicitly established** by the user
4. **Solution process control** remains in the user's hands

## Key Features

### Modular Architecture
- **Network Class**: Represents physical topology (nodes, elements, connectivity, geometry)
- **System Class**: Represents a physical phenomenon (fluid dynamics, thermal, chemical, etc.)
- **Flexible Coupling**: Multiple systems can interact on the same network
- **Add/Remove Systems**: Dynamically attach physical systems to the network

### Multi-Physics Coupling Capability
- **Cross-System Dependencies**: Systems can read variables from other systems
- **Custom Coupling Functions**: User-defined relationships between different physical phenomena
- **Sequential Coupling**: Iterative solution of interdependent systems

### Extensible System Design
- **User-Defined Properties**: Define physical properties as functions of network state or system variables
- **Element Conductance Functions**: General k-element calculation using full network state
- **Nonlinear Relationships**: Complex dependencies between variables
- **Thermophysical Properties**: Constant or state-dependent (temperature, pressure, concentration)

### Advanced Solvers
- **Linear Solvers**: Iterative methods (CG, GMRES, BiCGStab) and direct methods
- **Nonlinear Solvers**: Newton-Raphson, Newton-Krylov, least squares optimization
- **Automatic Partitioning**: Efficient separation of known/unknown nodes
- **Coupled System Solvers**: Sequential resolution of interdependent systems

### Physical System Types
- **Diffusive Systems**: Pipe flow, diffusion processes, conduction
- **Advective Systems**: Heat transfer, convective transport, species transport

## Mathematical Flexibility

NetSystems is built around the generalized system representation:

```
A · x = b
```

Where:
- **A**: Coefficient matrix (conductance, advection, diffusion terms)
- **x**: Nodal unknowns (heads, temperatures, concentrations)
- **b**: Source/sink terms (flows, heat, mass)

Users can define:
- **Element Conductance Functions**: k = f(network_state) - where network_state includes all system variables
- **Property Functions**: ρ = f(T, P, composition) - using variables from multiple systems
- **Coupling Functions**: Relationships between different system variables

## User Workflow

1. **Define Network**: Create nodes, elements, and connectivity
2. **Create Systems**: Instantiate fluid, thermal, or custom systems
3. **Set Properties**: Define element conductance and material properties
4. **Establish Couplings**: Specify how systems interact
5. **Solve**: Use built-in solvers or implement custom solution strategies
6. **Analyze**: Access results (matplotlib)

## Illustrated Example: Highly Coupled Model

```python
# Create a network with coupled fluid-thermal-chemical systems
net = Network()

...

# Add multiple interacting systems
fluid = System("fluid")
thermal = System("thermal")
chemical = System("concentration")

net.add_system(fluid)
net.add_system(thermal)
net.add_system(chemical)

...

# Define element conductivity as function of all system states
def coupled_conductivity(network):
    H = network.get("nodal_head")             # Access fluid heads
    T = network.get("nodal_temperature")      # Access temperatures
    C = network.get("concentration")          # Access concentrations
    
    # Complex relationship using all coupled variables
    k = f(H, T, C, network.element_diameters, ...)
    return k

fluid.element_conductance_function = coupled_conductivity

# Define temperature-dependent viscosity
def temperature_dependent_viscosity(network):
    T = network.get("temperature")
    return PropsSI("VISCOSITY", "T", T, ...)


# Solve coupled systems
solve_coupled_sequential(net, ["fluid", "thermal", "concentration"])
```

## Integration with Scientific Python Ecosystem

NetSystems seamlessly integrates with:
- **NumPy/SciPy**: Core numerical operations and solvers
- **CoolProp**: Real fluid properties for accurate simulations
- **Matplotlib/Pandas**: Visualization and data analysis
- **Custom Libraries**: User-defined property databases




## Advective Nodal Systems

### General description

An **advective system** models the transport of a nodal scalar quantity induced by a **previously computed flow field**, typically obtained as the solution of another physical system. In hydraulic and thermo–hydraulic networks, this formulation is used to describe the transport of **energy**, **enthalpy**, **concentration**, or other passive scalars once the **volumetric flow rate** is known.

In the proposed framework, advective systems are solved in a **partitioned manner**:

1. a diffusive (or hydraulic) system is solved first to compute the flow field,
2. the resulting flow is then used as input data for the assembly and solution of the advective system.

This explicit decoupling provides both numerical robustness and conceptual clarity, while allowing different physical models to interact through shared variables.

---

### Transported variable and advective property

Let $\phi$ denote the nodal scalar variable being transported (e.g., temperature). The advective transport is governed by the product of:

* the **volumetric flow rate** $Q_e$ on each element $e$,
* an **advective property** $\gamma$, which multiplies the flow.

For thermal transport, for instance:

$$
\gamma = \rho c_p
$$

so that the advective flux represents an enthalpy flow.

In the implementation, this property is defined **at nodes**, either as:

* a constant nodal vector `nodal_gamma`, or
* a state-dependent function `nodal_gamma_function(network)`.

The elemental value is obtained by nodal averaging:

$$
\gamma_e = \frac{\gamma_i + \gamma_j}{2}
$$

This is directly reflected in the assembly routine:

```python
gamma_e = 0.5 * (gamma_nodal[i] + gamma_nodal[j])
adv_coeff = abs(Q) * gamma_e
```

---

### Governing equation of an advective system

For a node $i$, the discrete advective balance can be written as:

$$
\sum_{e \in \mathcal{E}_i}
|Q_e| \gamma_e  (\phi_i - \phi_{\text{up}}) = S_i
$$

where:

* $\mathcal{E}_i$ is the set of elements connected to node $i$,
* $\phi_{\text{up}}$ is the upstream nodal value according to the flow direction,
* $S_i$ represents an equivalent nodal source term.

Enthalpy flow

$$\text{Flux} = Q_e \gamma_e (\phi_i - \phi_{\text{up}})$$

- $Q_e$ $\rightarrow$ flow rate
- $\gamma_e$ $\rightarrow$ $\rho c_p$
- $\phi$ $\rightarrow$ $T$

---

### Upwind discretization in the code

The upwind direction is determined exclusively from the sign of the elemental flow rate $Q_e$:

```python
if Q >= 0.0:
    up, down = i, j
else:
    up, down = j, i
```

The advective contribution is assembled only on the **downstream node**, yielding a positive diagonal term and a negative coupling with the upstream node:

```python
A[down, down] += adv_coeff
A[down, up]   -= adv_coeff
```

This corresponds exactly to the discrete term

$$
|Q_e| \gamma_e (\phi_{\text{down}} - \phi_{\text{up}})
$$

---

### Element-based source terms

The framework allows **element-level source terms** (`element_flux`)  , which are projected onto the downstream node:

```python
if system.element_flux is not None:
    b[down] += system.element_flux[e]
```

These terms may represent:

* heat generation per unit length,
* thermal exchange with the surroundings,
* localized energy injection.

Alternatively, a function

```python
system.element_flux_function(network)
```

can be provided, allowing sources to depend on the current network state (temperature, flow rate, iteration counter, etc.).

---

### Solution of the advective system

The assembled advective system is **linear** with respect to the transported nodal variable and is solved using a **partitioned approach with prescribed nodal values**.

Nodes with known values are eliminated from the linear system, leading to:

$$
A_{uu} \phi_u = b_u - A_{uk} , \phi_k
$$

which is solved using a sparse linear solver. The full solution vector is then reconstructed by combining known and unknown values.

This procedure is implemented in `solve_advective_system`.

---

### Role in coupled simulations

Within a **sequentially coupled solver**, the advective system:

* **does not compute the flow**,
* **uses the flow** obtained from a diffusive or hydraulic system,
* updates a transported variable that may, in turn, affect material properties or source terms in subsequent iterations.



## Sequential Solution of Coupled Diffusive–Advective Systems

### Overview of the coupling strategy

The function `solve_coupled_sequential` implements a **sequential (staggered) coupling** between two physical systems defined on the same network:

1. a **diffusive system** (e.g. hydraulic head and flow),
2. an **advective system** (e.g. thermal transport).

The key assumption is that **advection depends on the flow**, while the flow may depend (directly or indirectly) on the advected variable through material properties. As a result, the coupled problem is solved iteratively by alternating between both systems until convergence is achieved.

The solution order is fixed:

* the **diffusive system is solved first**,
* its resulting elemental flux is passed to the **advective system**.

---

### Initialization of the coupled iteration

Before entering the coupled loop, both systems must be initialized:

```python
net.iteration = 0
```

For the **advective system**, an initial guess for the transported variable is required:

```python
if adv.initial_x is None:
    adv.x = np.ones(net.n_nodes) * 300.
```

This corresponds to an initial condition:

$$
\boldsymbol{\phi}^{(0)} = \phi_0
$$

(e.g. uniform initial temperature).

Similarly, the diffusive system is initialized if no previous solution is available:

```python
if diff.initial_x is None:
    diff.x = np.ones(net.n_nodes)
```

---

### Step 1: Diffusive system solution

The diffusive system is typically **nonlinear** and is solved using Newton–Raphson iterations.

A suitable initial guess is obtained:

```python
guess = get_initial_guess(net, diffusive_system_name)
```

The nonlinear problem is then solved:

```python
res = solve_newton_raphson(
    net,
    system_name=diffusive_system_name,
    initial_guess=guess,
    tolerance=hydraulic_tol,
    max_iterations=hydraulic_max_iter
)
```

After convergence, the solution is finalized:

```python
finalize_solution(net, diffusive_system_name, res)
```

This step updates:

* nodal diffusive unknowns (e.g. heads),
* elemental fluxes $\mathbf{Q}$, stored as `diff.element_flux`.

---

### Step 2: Advective system solution

Once the flux field is known, the advective system becomes **linear** and is solved in a single step:

```python
solve_advective_system(
    net,
    advective_system_name,
    diff.element_flux
)
```

The solution is stored in:

```python
adv.x
```

---

### Convergence criterion

After the advective solve, convergence is assessed using the norm of the change in the transported variable:

```python
err = np.linalg.norm(T - T_old)
```

which corresponds to:

$$
|\boldsymbol{\phi}^{k+1} - \boldsymbol{\phi}^{k}|_2 < \varepsilon
$$

If the tolerance is satisfied:

```python
if err < tol:
    break
```

the coupled iteration terminates.

The iteration counter is updated:

```python
net.iteration += 1
```

allowing material properties or source terms to depend on the iteration index.



## Core Abstractions: `Network` and `System`

The framework is built around two fundamental abstractions:

* the **`Network`**, which represents the discrete spatial domain (nodes, elements, topology),
* the **`System`**, which represents a physical model defined on that network.

This separation allows multiple physical systems (hydraulic, thermal, transport, etc.) to coexist, interact, and be coupled in a modular way.

---

## The `Network` class

### Conceptual role

The `Network` class represents the **computational graph** on which all physical systems are defined. It contains:

* the **topology** (nodes and element connectivity),
* the **geometry** (node coordinates, element lengths),
* a **registry of physical systems**,
* a global **iteration counter** used in coupled simulations.

Mathematically, the network defines:

* a set of nodes $\mathcal{N}$,
* a set of elements $\mathcal{E} \subset \mathcal{N} \times \mathcal{N}$.

All systems share this same discretization.

---

### Topology and geometry

```python
self.n_nodes
self.n_elements
self.connectivity          # (n_elements, 2)
self.node_coordinates      # (n_nodes, dim)
self.element_lengths
```

Each element $e$ connects two nodes $(i, j)$, stored in `connectivity`. The element length is computed as:

$$
L_e = |\mathbf{x}_j - \mathbf{x}_i|
$$

which is implemented in the method of the class Network:

```python
def calculate_element_lengths(self):
    for e in range(self.n_elements):
        i, j = self.connectivity[e].astype(int)
        self.element_lengths[e] = np.linalg.norm(
            self.node_coordinates[j] - self.node_coordinates[i]
        )
```

---

### System registry and access

The network stores all physical systems in a dictionary:

```python
self.systems = {}
```

Systems are added via:

```python
net.add_system(system)
```

which also creates a **back-reference**:

```python
system.network = self
```

This bidirectional link allows:

* systems to query network data (geometry, other variables),
* the network to retrieve variables from systems.

---

### Unified variable access

The method:

```python
net.get(var_name)
```

allows retrieving variables **by physical name**, not by system. If exactly one system declares ownership of a variable name (e.g. `"temperature"`), the network returns it transparently.

This corresponds conceptually to:

$$
\text{Network}[\text{variable name}] \longrightarrow \text{System variable}
$$

This design is essential for multi-physics coupling, where systems need access to each other’s results without tight coupling.

---

## The `System` class

### Conceptual role

A `System` represents a **discrete physical model** defined on the network. Each system is either:

* **diffusive**, solving a balance of fluxes driven by gradients,
* **advective**, solving a transport problem driven by a known flow.

Formally, each system defines:

* a nodal unknown $\mathbf{x}$,
* a nodal source term $\mathbf{b}$,
* a system matrix $\mathbf{A}$.

---

### Mathematical variables

```python
self.x    # nodal unknown
self.b    # nodal source
self.A    # system matrix
self.initial_x
```

These correspond to the linear (or linearized) system:

$$
\mathbf{A}\mathbf{x} = \mathbf{b}
$$

`initial_x` is used to initialize iterative or coupled solvers.

---

### Physical metadata

```python
self.x_name
self.b_name
self.element_variable_name
```

These attributes define **semantic names** for:

* the nodal unknown (e.g. `"nodal_head"`, `"temperature"`),
* the nodal source (e.g. `"nodal_flow"`, `"heat"`),
* the element-wise variable (e.g. `"flow_rate"`).

This enables system-agnostic access through the network.

---

### Boundary conditions

```python
self.known_x_nodes
self.known_x_values
self.known_b_nodes
self.known_b_values
```

These arrays define Dirichlet-type conditions:

$$
x_i = \bar{x}_i 
$$

They are enforced during system assembly and solution via matrix partitioning.

---

### Element properties and constitutive behavior

```python
self.element_conductance
self.element_conductance_function
```

These represent **element-level constitutive coefficients**, which may be:

* constant,
* or computed dynamically from the network state.

For advective systems, the transported quantity is weighted by a **nodal advective property**:

```python
self.nodal_gamma
self.nodal_gamma_function
```

This corresponds to quantities such as $\rho c_p$ in thermal transport.


---

## Example: hydraulic–thermal coupled setup

### Diffusive system: hydraulics

```python
fluid = System("fluid", system_type="diffusive")
```

This system solves for **nodal hydraulic head**:

$$
\mathbf{h}
$$

with:

```python
fluid.x_name = "nodal_head"
fluid.b_name = "nodal_flow"
fluid.element_variable_name = "flow_rate"
```

Boundary conditions:

```python
fluid.known_x_nodes = [15]
fluid.known_x_values = [0]
fluid.known_b_nodes = [0]
fluid.known_b_values = [Q_vol]
```

which enforce:

$$
h_{15} = 0, \qquad q_0 = Q_{\text{vol}}
$$

The element conductivity is defined through:

```python
fluid.element_conductance_function = calculate_element_conductivity
```

which computes a nonlinear, Reynolds-dependent hydraulic conductance.

---

### Advective system: heat transport

```python
heat = System("heat", system_type="advective")
```

This system solves for **nodal temperature**:

$$
\boldsymbol{T}
$$

with:

```python
heat.x_name = "temperature"
heat.b_name = "heat"
heat.element_variable_name = "enthalpy_flow"
```

Boundary condition:

```python
heat.known_x_nodes = [0]
heat.known_x_values = [300]
```

which enforces:

$$
T_0 = 300;\text{K}
$$

The advective property is provided via:

```python
heat.nodal_gamma_function = gamma_function_heat
```

corresponding to (enthalpy flow):

$$
\gamma = \rho c_p
$$

Element heat sources are defined as:

```python
heat.element_property = q_w * pi * D * L
```

which represents distributed heat input per element.

---

### Registering systems in the network

Finally, both systems are attached to the network:

```python
net.add_system(fluid)
net.add_system(heat)
```

At this point, both systems share the same topology and geometry.


---

## General design logic

The combined `Network`–`System` design achieves:

* **Separation of concerns** (geometry vs physics),
* **Coupling through shared variables**,
* **Extensibility** to new physics without modifying existing solvers,


