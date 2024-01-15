# 2D Ising Model Simulation with the Metropolis Algorithm
The 2D Ising Model is a pivotal model in statistical mechanics, renowned for its ability to describe phase transitions in ferromagnetic materials.
This project aims to provide an interactive and visual simulation of the Ising Model, using the metropolis algorithm and to calculate key physical observables such as magnetization, susceptibility and energy.
## Theory
### Ising Model
We consider a set $\Lambda$ of lattice sites, each with a set of adjacent sites forming a two-dimensional lattice. For each lattice site $i \in \Lambda$ there is a classical spin $\sigma_i$ that can be in two states $\sigma_i= \pm 1$.
The interaction of two spin couples nearest neighbors is denoted as $\langle i j\rangle$.
For any two adjacent sites $i, j \in \Lambda$ there is an interaction $J_{ij}$. In addition, a site $j \in \Lambda$ has an external field $B_j$ interacting with it.
Then the energy of a certain configuration is given by the Hamiltonian
$$H=-\sum_{\langle i j\rangle} J_{ij} \sigma_i \sigma_j-\sum_i B_j \sigma_i$$
If we chose $J>0$ the system forms a ferromagnet, whereas $J<0$ implies a antiferromagnetic system.
To simplify matters, we set $J=1$ and $B=0$.
### Metropolis Algorithm
The system is modeled on a quadratic lattice $\Lambda$ with $N$ nodes in each dimension. This project implements periodic boundary conditions, i.e a spin at the left/lower boundary $(k, j=0) /(k=0, j)$ has its nearest neighbor at the opposite right/upper boundary $(k, j=N-1) /(k=N-1, j)$.
At time $t$, the system is in a state $x_t=\sigma_i \mid i \in \Lambda$. To move forward to time $t+1$, a new state $x^{\prime}$ is randomly proposed.
For the transition from the current state $x_t$ to the new state $x^{\prime}$, which is referred to as an update, a transition rate $w_{x_t x^{\prime}}$ is used within the time interval $dt$:

$$w_{x_t x^{\prime}} d t= \begin{cases}e^{-\beta \cdot\left(H\left(x^{\prime}\right)-H\left(x_t\right)\right)} & \text { if } H\left(x^{\prime}\right)-H\left(x_t\right) \geq 0 \\\ 1 & \text { otherwise }\end{cases}$$

where $\beta$ is the inverse temperature $\beta=\frac{1}{k_b T}$.
To maintain ergodicity, which means ensuring that the algorithm can reach and propose all possible states, a proposed state $x^{\prime}$ is only different from the current state by a single spin $\sigma_a$. 
If the energy of the new state $x^{\prime}$ is less than the energy of the current state $x_t$ (i.e., $H\left(x^{\prime}\right)-H\left(x_t\right) \leq 0$), the system adopts state $x^{\prime}$. 
f not, the system assumes state $x^{\prime}$ with a probability of $e^{-\beta \cdot\left(H\left(x^{\prime}\right)-H\left(x_t\right)\right)}$; otherwise, $x^{\prime}$ is rejected and the system remains in state xt for the next step ($x_t+1 = x_t$).
A random number $r$ within the range $[0, 1)$ determines whether $x^{\prime}$ is accepted or rejected:




## Module `animation`
There are a few custom classes implemented in the code that need a little bit of explanation. The module [animation](https://github.com/dantona02/projects/blob/main/animation.py) has the following classes:
- ### class `HaloPoint`
  An object of the class `HaloPoint` can be created like the following:
  ```python
  halopoint = HaloPoint(ax, mass, color_decay, color1, color2='white')
  ```
  - `ax` corresponds to the current instance of `Axes` of the animation.
  - `mass` is a parameter to change the size of `halopoint` according to the mass, even if that doesn't reflect the physical reality.
  - `color_decay` changes the decay of the halo. The larger the value of `color_decay`, the smaller the halo.
  - `color1` sets the color of the point.
  - `color2` should generally not be changed.
    
  **It is very important to note that the `get_artists()` method returns a list of `Artists`, which must then be merged into one list to be returned by the update-function.
    This applies to all other animation classes contained in this module.**