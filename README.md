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
If the energy of the new state $x^{\prime}$ is less than the energy of the current state $x_t$ (i.e. $H\left(x^{\prime}\right)-H\left(x_t\right) \leq 0$), the system adopts state $x^{\prime}$. 
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
In a system approaching the thermodynamic limit and existing below the critical temperature, it's equally probable for the system to magnetize positively or negatively.
However, without an external magnetic field, it's impossible for the system to switch between positive and negative spontaneous magnetization.
In contrast, for systems of finite size, there is a specific time period, dependent on the lattice size, during which the system can transition from one state of spontaneous magnetization to the other.
Therefore, simply averaging the magnetization would result in an inaccurate value. For this reason, the absolute value of the magnetization, represented as $\langle|M|\rangle$, is used as the order parameter instead of $\left\langle M\right\rangle$ [1].

## Module `ising`
The module [ising](https://github.com/dantona02/metropolis/blob/main/ising.py) contains the class `Isingmodel`.
An object of this class can be initialized like the following:
### Initializing class object
An object of the class `Isingmodel` can be created like the following:
```python
model = Isingmodel(N=50, init=True)
```
- `N` corresponds to the size of the lattice. This is only needed if one wants to use the function `plot()`.
- `init` if set to `True` all spins of the lattice are initialized up (i.e. $\sigma_i=1$). If set to `Fasle`, about 50 percent of the spins will be positiv and the other half will be negativ.
### Plotting lattice
To plot the lattice at a certain temperature after $N$ sweeps, we use the `plot()` method.
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
- `grid()` returns the two-dimensional lattice as an array with shape $N\times N$
- `get_energy()` calculates the total energy of the system with its initial spin configuration

Now we can run the `metropolis()` method and pass in all the parameters. This method returns an array for the current total spin and the current total energy of the system, each with shape `(1, sweeps)`,
an array with shape `(N, N)`, which represents the spin configuration of the lattice after a certain number of sweeps and an array with the squared total spin.
```python
spins, energies, equilibrium, spin_sq = model.metropolis(lattice, sweeps, beta, energy, N)
```
Running the `metropolis()` method might take some time depending on the parameters passed to the method, even though it is already optimized with `numba`.
Once it finished running, we can plot the result using `plot()`:
```python
model.plot(equilibrium=equilibrium, cmap='binary', times=sweeps, beta=beta, save=True)
```
- `equilibrium` takes the third return argument of the `metropolis()` method as the input.
- `cmap` can be chosen arbitrarily, but it needs to be supported by `matplotlib`.
- `times` represents the number of sweeps.
- `save` if set to `True` the method will save the image to the current directory.

### Animating the ising model
It is possible to visualize und animate the ising model for different configurations.
The initialization of N, beta, sweeps, model, lattice and energy works analogously to the previous paragraph.
The method `animate()` returns an array of the shape `(sweeps, N, N)`:
```python
animation = model.animate(lattice, sweeps, beta, energy, N)
```
All method arguments operate as in the function `metropolis()`.
The animation array can then be passed into the `save()` method.
```python
model.save(animation, path='C:\Users\<username>\myproject\animation.gif', font_size=20, display_sweeps=True)
```
- `animation` return value of the `animate()` method.
- `path` determines the location of the animation.
- `font_size` sets the size of the text displayed in the animation.
- `display_sweeps` if set to `True`, the past sweeps are displayed.

**Important**: Before using the method `save()`, the path to the user's font must be changed to the corresponding location.
To do that, navigate to `ImageFont.truetype("/Library/Fonts/Arial Unicode.ttf", font_size)` in the `save()`method of `Isingmodel` and change to path to the one of the desired font.

The animation might look like the following:

<p align="center">
  <img src="images/ising400.gif" alt="400 by 400 lattice" width="300">
</p>

### Calculating magnetization and susceptibility
The two methods `get_magnetization()` and `get_susceptibility()` of the class `Isingmodel` allow to calculate the magnetization
and the susceptibility of a system for several different $\beta$ and lattice sizes.
As this computation takes a lot of resources, it is recommended to reduce the size of the system's lattice.
The magnetization can be computed like the following:
```python
betas = np.linspace(0.2, 0.65, 30)
N = [8, 10]
sweeps = [20000 for i in N]
sampleSize = [12000 for i in N]
iterations = 10
```
```python
model = Isingmodel()
magnetization = model.get_magnetization(N, betas, sweeps, sampleSize, iterations)
```
Both methods return an array of shape `(len(N), len(betas))`.
The array can then be visualized with `matplotlib`.

- `N` corresponds to the lattice sizes. This needs to be an array, even if `N` contains only one element.
- `betas` provides an array of values of the inverse temperature. An associated magnetization value is calculated for each value from `betas`.
- `sweeps` represents the number of sweeps.
- `sampleSize` corresponds to the number of sweeps to be averaged from the return values of the `metropolis()` method in order to reduce statistical fluctuations.
- `iterations` determines how often the calculations are performed per value from `betas`. This further reduces the statistical fluctuations and thus improves the quality of the calculated values.
Depending on the number of available processor cores, this in turn also increases the calculation time.

`get_susceptibility()` works exactly the same way as `get_magnetization()`.

### Required packages
`numpy`, `scipiy`, `matplotlib`, `pillow` and `numba`. They can all be installed via `pip`or `conda`.

### References
- [1] Müller-Krumbhaar, H. (2016). [The Hobbyhorse of Magnetic Systems: The Ising Model](https://www.researchgate.net/publication/304163856_The_Hobbyhorse_of_Magnetic_Systems_The_Ising_Model). ResearchGate.