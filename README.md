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
 if $r$ is less than or equal to $e^{-\beta \cdot\left(H\left(x^{\prime}\right)-H\left(x_t\right)\right)}$, then $x^{\prime}$ is accepted. In state $x^{\prime}$, the selected spin $\sigma_a$ is flipped: 
$x^{\prime}\left(\sigma_a\right)=-x_t\left(\sigma_a\right)$. This entire process of choosing a spin, proposing a new state, and then accepting or rejecting this proposal is referred to as an update.
After performing $N\times N$ metropolis updates, which equals one sweep, each spin was statistically updated once. To properly thermalize the system, the number of sweeps should not be to low.
Once the system reached its equilibrium state, we can approximate the observables as time average of these states.
For the magnetization we get
$$\left\langle M\right\rangle=\frac{1}{N^2}\sum_{i} \sigma_i$$
and for the susceptibility
$$\chi=\beta\left(\left\langle M^2\right\rangle-\left\langle M\right\rangle^2\right)$$


## Module `ising`
The module [ising](https://github.com/dantona02/metropolis/blob/main/ising.py) contains the class `Isingmodel`.
An object of this class can be initialized like the following:
### Initializing class object
An object of the class `Isingmodel` can be created like the following:
```python
model = Isingmodel(N=50, init=True)
```
- `N` corresponds to the size of the lattice. This is only needed if one wants to use the function `.plot()`.
- `init` if set to `True`all spins of the lattice are initialized up (i.e. $\sigma_i=1$). If set to `Fasle`, about 50 percent of the spins will be positiv and the other half will be negativ.
### Plotting lattice
To plot the lattice at a certain temperature after $N$ sweeps, we use the `.plot()` method.
First let's initialize the model and set the required parameters:
```python
N = 50
beta = 0.5
sweeps = 2500
model = Isingmodel(N=N, init=True)
```
The next step is to get the lattice and the total energy of the system:
```python
lattice = model.grid(N)
energy = model.get_energy(lattice)
```
- `.grid()` returns the two-dimensional lattice as an array with shape $N\times N$
- `.get_energy()` calculates the total energy of the system with its initial spin configuration

Now we can run the `.metropolis()` method and pass in all the parameters. This method returns an array for the current total spin and the current total energy of the system, each with shape `(1, sweeps)`,
an array with shape `(N, N)`, which represents the spin configuration of the lattice after a certain number of sweeps and an array with the squared total spin.
```python
spins, energies, equilibrium, spin_sq = model.metropolis(lattice, sweeps, beta, energy, N)
```
Running the `.metropolis()` method might take some time depending on the parameters passed to the method, even though it is already optimized with `numba`.
Once it finished running, we can plot the result using `.plot()`:
```python
model.plot(equilibrium=equilibrium, cmap='binary', times=sweeps, beta=beta, save=True)
```
- `equilibrium` takes the third return argument of the `metropolis()` method as the input.
- `cmap` can be chosen arbitrarily, but it needs to be supported by `matplotlib`.
- `times` represents the number of sweeps.
- `save` if set to `True` the method will save the image to the current directory.

### Animating the ising model
With the method `.animate()`of the Isingmodel class, it is possible to visualize und animate the ising model.
![400 by 400 lattice](images/ising400.gif)