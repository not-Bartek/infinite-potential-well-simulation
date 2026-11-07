import pandas as pd
import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation


plt.style.use('dark_background')
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['axes.edgecolor'] = '#333333'
plt.rcParams['axes.linewidth'] = 1.2

def animate_quantum_simulation(csv_filename="simulation.csv", L=10.0, x_min=0.0):

    data = pd.read_csv(csv_filename, header=None)
    
    t_vals = data.iloc[:, 0].values
    psi_sq = data.iloc[:, 1:].values
    
    num_steps, N = psi_sq.shape
    x_vals = np.linspace(x_min, L, N)
    
    fig, ax = plt.subplots(figsize=(10, 6), dpi=100)
    fig.patch.set_facecolor('#0b0b12')
    ax.set_facecolor('#0b0b12')

    line, = ax.plot(x_vals, psi_sq[0], color='#ff007f', lw=2.5, label=r'$|\Psi(x,t)|^2$')
    fill = ax.fill_between(x_vals, psi_sq[0], color='#7b2cbf', alpha=0.35)

    max_val = np.max(psi_sq) * 1.15
    ax.set_xlim(x_min, L)
    ax.set_ylim(0, max_val)

    ax.set_xlabel('Position $x$', fontsize=13, labelpad=10, color='#e0e0e0')
    ax.set_ylabel('Probability density $|\Psi(x,t)|^2$', fontsize=13, labelpad=10, color='#e0e0e0')
    title_text = ax.set_title(f'Czas $t = {t_vals[0]:.3f}$', fontsize=14, pad=15, color='#ffffff', weight='bold')

    # ax.axvline(0, color='#00f0ff', linestyle='--', linewidth=2, alpha=0.7, label='Ściany studni')
    # ax.axvline(L, color='#00f0ff', linestyle='--', linewidth=2, alpha=0.7)

    ax.grid(True, linestyle=':', alpha=0.2, color='#ffffff')
    ax.legend(loc='upper right', frameon=True, facecolor='#161625', edgecolor='#33334d', labelcolor='#ffffff')

    def update(frame):
        nonlocal fill
        line.set_ydata(psi_sq[frame])
        
        fill.remove()
        fill = ax.fill_between(x_vals, psi_sq[frame], color='#7b2cbf', alpha=0.35)
        
        title_text.set_text(f'Evolution: $t = {t_vals[frame]:.3f}$')
        return line, fill, title_text

    anim = FuncAnimation(fig, update, frames=num_steps, interval=30, blit=False)
    
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Quantum simulation visualizer")
    
    parser.add_argument("csv", nargs="?", default="simulation.csv", help="Path to results csv")
    parser.add_argument("L", nargs="?", type=float, default=10.0, help="Lenght of the infinite potential well (L)")
    parser.add_argument("x_min", nargs="?", type=float, default=0.0, help="Left bound of the simulation")

    
    args = parser.parse_args()
    
    animate_quantum_simulation(args.csv, L=args.L, x_min=args.x_min)