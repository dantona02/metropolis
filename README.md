# 2D Ising Model Simulation with the Metropolis Algorithm
The 2D Ising Model is a pivotal model in statistical mechanics, renowned for its ability to describe phase transitions in ferromagnetic materials.
This project aims to provide an interactive and visual simulation of the Ising Model, using the metropolis algorithm and to calculate key physical observables such as magnetization, susceptibility and energy.
## Theory
### Ising Model
We consider a set $\Lambda$ of lattice sites, each with a set of adjacent sites forming a two-dimensional lattice. For each lattice site $i \in \Lambda$ there is a classical spin $\sigma_i$ that can be in two states $\sigma_i= \pm 1$.
The interaction of two spin couples nearest neighbors is denoted as $\langle i j\rangle$.
For any two adjacent sites $i, j \in \Lambda$ there is an interaction $J_{ij}$. In addition, a site $j \in \Lambda$ has an external field $B_j$ interacting with it.
Then the energy of a certain configuration is given by the Hamiltonian
\begin{equation}
H=-\sum_{J_{ij}\langle i j\rangle} \sigma_i \sigma_j-\sum_i B_j \sigma_i .
\end{equation}
If we chose $J>0$ the systems forms a ferromagnet, whereas $J<0$ implies a antiferromagnetic system.
To simplify matters, we set $J=1$ and $B=0$.
### Metropolis Algorithm


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