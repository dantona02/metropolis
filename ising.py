import numba
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import convolve, generate_binary_structure
import matplotlib.patches as mpatches
from PIL import Image, ImageDraw, ImageFont
import multiprocessing


def calc_magnetization(N_array, betas, steps, sweeps):
    all_magnetizations = []
    for index1, i in enumerate(N_array):
        magnetization = np.zeros(len(betas))
        helpmodel = Isingmodel(N=i)
        grid = helpmodel.grid(i)
        energy = helpmodel.get_energy(grid)
        for index2, s in enumerate(betas):
            spins, energies, equilibrium, mag_squared = helpmodel.metropolis(grid, steps[index1], s, energy, i)
            sweep_index = int(sweeps[index1])
            magnetization[index2] = spins[-sweep_index:].mean() / i ** 2
        all_magnetizations.append(magnetization)
    return all_magnetizations


def calc_susceptibility(N_array, betas, steps, sweeps):
    all_sus = []
    for index1, i in enumerate(N_array):
        variance = np.zeros(len(betas))
        helpmodel = Isingmodel(N=i)
        grid = helpmodel.grid(i)
        energy = helpmodel.get_energy(grid)
        for index2, s in enumerate(betas):
            spins, energies, equilibrium, mag_squared = helpmodel.metropolis(grid, steps[index1], s, energy, i)
            sweep_index = int(sweeps[index1])
            variance[index2] = s * (
                    mag_squared[-sweep_index:].mean() - spins[-sweep_index:].mean() ** 2) / i ** 2
        all_sus.append(variance)
    return all_sus


class Isingmodel:
    def __init__(self, N=50, J=1, init=True):
        self.J = J
        self.N = N
        self.init = init

    def grid(self, N):
        init_random = np.random.random((N, N))
        if self.init:
            grid_n = np.ones((N, N))
        else:
            grid_n = np.zeros((N, N))
            grid_n[init_random >= .5] = 1
            grid_n[init_random < .5] = -1
        return grid_n

    def get_energy(self, grid):
        kernel = generate_binary_structure(2, 1)
        kernel[1][1] = False
        result = - self.J * grid * convolve(grid, kernel, mode='wrap')
        return result.sum()

    @staticmethod
    @numba.njit("Tuple((f8[:], f8[:], f8[:,:], f8[:]))(f8[:,:], i8, f8, f8, i8)", nopython=True, nogil=True)
    def metropolis(arr_spin, sweeps, beta, energy, N):
        times = sweeps * N ** 2
        arr_spin = arr_spin.copy()
        total_spin = np.zeros(sweeps)
        total_energy = np.zeros(sweeps)
        for t in range(0, times - 1):
            x = np.random.randint(0, N)
            y = np.random.randint(0, N)
            spin_t = arr_spin[x, y]
            spin_prime = -spin_t
            E_t = 0
            E_prime = 0
            neighbours = [(x - 1) % N, (x + 1) % N, (y - 1) % N, (y + 1) % N]
            for nx in [neighbours[0], neighbours[1]]:
                E_t += - spin_t * arr_spin[nx, y]
                E_prime += - spin_prime * arr_spin[nx, y]
            for ny in [neighbours[2], neighbours[3]]:
                E_t += - spin_t * arr_spin[x, ny]
                E_prime += - spin_prime * arr_spin[x, ny]
            dE = E_prime - E_t
            if (dE > 0) & (np.random.random() <= np.exp(-beta * dE)):
                arr_spin[x, y] = spin_prime
                energy += dE
            elif dE <= 0:
                arr_spin[x, y] = spin_prime
                energy += dE
            if t % (N ** 2) == 0:
                total_spin[t // (N ** 2)] = arr_spin.sum()
                total_energy[t // (N ** 2)] = energy
        mag_squared = total_spin ** 2
        return total_spin, total_energy, arr_spin, mag_squared

    def plot(self, equilibrium, cmap, times, beta, save=False):
        fig, ax = plt.subplots(dpi=150)
        plt.rcParams["font.family"] = "times"
        plt.rcParams["text.usetex"] = True
        plt.title(
            fr'{self.N}$\times${self.N}-lattice after {times:.1e} sweeps ($N^2$ updates)' + '\n' + rf'$\beta={beta}$',
            fontsize=12, pad=10)
        ax.imshow(equilibrium, cmap=cmap, interpolation='none')
        white_patch = mpatches.Patch(color='white', label=r'$\sigma_i=-1$')
        black_patch = mpatches.Patch(color='black', label=r'$\sigma_i=+1$')

        legend = ax.legend(handles=[white_patch, black_patch], loc="upper right", bbox_to_anchor=(1.4, 1.02),
                           fancybox=False, edgecolor='black', fontsize=12, facecolor='whitesmoke')
        legend.set_zorder(10)
        legend.get_frame().set_linewidth(0.5)
        plt.show()
        if save:
            fig.tight_layout()
            fig.savefig(f"/Users/danielmiksch/Downloads/{self.N}by{self.N}_grid.png")

    @staticmethod
    @numba.njit("f8[:,:,:](f8[:,:], i8, f8, f8, i8)", nopython=True, nogil=True)
    def animation(arr_spin, sweeps, beta, energy, N):
        times = sweeps * N ** 2
        arr_spin = arr_spin.copy()
        animation_arr = np.zeros((sweeps, N, N))
        for t in range(0, times - 1):
            x = np.random.randint(0, N)
            y = np.random.randint(0, N)
            spin_t = arr_spin[x, y]
            spin_prime = -spin_t
            E_t = 0
            E_prime = 0
            neighbours = [(x - 1) % N, (x + 1) % N, (y - 1) % N, (y + 1) % N]
            for nx in [neighbours[0], neighbours[1]]:
                E_t += - spin_t * arr_spin[nx, y]
                E_prime += - spin_prime * arr_spin[nx, y]
            for ny in [neighbours[2], neighbours[3]]:
                E_t += - spin_t * arr_spin[x, ny]
                E_prime += - spin_prime * arr_spin[x, ny]
            dE = E_prime - E_t
            if (dE > 0) & (np.random.random() <= np.exp(-beta * dE)):
                arr_spin[x, y] = spin_prime
                energy += dE
            elif dE <= 0:
                arr_spin[x, y] = spin_prime
                energy += dE
            if t % (N ** 2) == 0:
                animation_arr[t // (N ** 2)] = arr_spin.copy()
        return animation_arr

    @staticmethod
    def save(frame_array, path, font_size, display_sweeps=True):
        pil_frames = []
        for i, frame in enumerate(frame_array):
            img = Image.fromarray(np.uint8(frame))
            b, _ = img.size
            if img.mode != 'RGB':
                img = img.convert('RGB')
            draw = ImageDraw.Draw(img)
            if display_sweeps:
                text = f"Sweeps: {i}"
                font = ImageFont.truetype("/Library/Fonts/Arial Unicode.ttf", font_size)
                draw.text((b // 2, 10), text, fill="red", font=font, anchor='mt')
            pil_frames.append(img)

        first_image = pil_frames[0]
        first_image.save(path, format='GIF', append_images=pil_frames[1:], save_all=True, duration=20, loop=0)

    def magnetization(self, spins):
        return spins / self.N ** 2

    def print_energy(self):
        print(self.get_energy(self.grid(self.N)))

    @staticmethod
    def get_magnetization(N_array, betas, steps, sweeps, iterations):
        args_list = [(N_array, betas, steps, sweeps) for t in range(0, iterations)]
        with multiprocessing.Pool() as pool:
            flat_results = pool.starmap(calc_magnetization, args_list)
        magnetization_list = np.array(flat_results).reshape((iterations * len(N_array), len(betas)))
        magnetization_avg = np.zeros((len(N_array), len(betas)))
        for i in range(0, len(N_array)):
            group = range(i, len(magnetization_list), len(N_array))
            magnetization_avg[i] = np.mean(magnetization_list[group, :], axis=0)
        return magnetization_avg

    @staticmethod
    def get_susceptibility(N_array, betas, steps, sweeps, iterations):
        args_list = [(N_array, betas, steps, sweeps) for t in range(0, iterations)]
        with multiprocessing.Pool() as pool:
            flat_results = pool.starmap(calc_susceptibility, args_list)
        variance_list = np.array(flat_results).reshape((iterations * len(N_array), len(betas)))
        variance_avg = np.zeros((len(N_array), len(betas)))
        for i in range(0, len(N_array)):
            group = range(i, len(variance_list), len(N_array))
            variance_avg[i] = np.mean(variance_list[group, :], axis=0)
        return variance_avg
