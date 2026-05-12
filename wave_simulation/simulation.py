import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve
import matplotlib.pyplot as plt

# --- 1. Parâmetros Físicos ---
f = 100.0
c0 = 1500.0
k0 = 2 * np.pi * f / c0
zs = -400.0       # Fonte a -400m
theta0 = 5.28     # Inclinado para baixo (+z)

# Grade Computacional (com espaço extra para a Sponge Layer)
z_plot_min, z_plot_max = -800.0, 200.0
z_grid_max = 600.0
dz = 0.5
r_max = 10000.0   # 10 km
dr = 5.0

z = np.arange(z_plot_min, z_grid_max + dz, dz)
r = np.arange(0, r_max + dr, dr)
Nz = len(z)

# --- 2. Perfil de Meio (Velocidade, Densidade e Atenuação) ---
c = np.zeros(Nz)
rho = np.ones(Nz)
alpha_sponge = np.zeros(Nz)

for i, zi in enumerate(z):
    if zi <= 0:
        # COLUNA DE ÁGUA
        c[i] = 1500.0 + 0.016 * (zi + 400.0) # 1500 na fonte, 1506.4 no fundo
        rho[i] = 1000.0
    else:
        # SEDIMENTO
        c[i] = 1655.0
        rho[i] = 1500.0
        # Sponge Layer: absorve energia profunda para evitar reflexão numérica
        if zi > 200:
            alpha_sponge[i] = 0.5 * ((zi - 200) / 400.0)**2 

n_complex = (c0 / c) + 1j * alpha_sponge

# --- 3. Montagem da Matriz IFD com Salto de Densidade ---
# Em vez de fórmulas complexas de interface, usamos o operador gradiente discreto correto:
rho_mid = (rho[:-1] + rho[1:]) / 2.0
V_up = rho[:-1] / rho_mid / dz**2
V_down = rho[1:] / rho_mid / dz**2

diag_V = np.zeros(Nz)
diag_V[:-1] -= V_up
diag_V[1:] -= V_down

# O termo 'u' agora contém a física da refração e da densidade
u_diag = k0**2 * (n_complex**2 - 1) + diag_V
gamma = 1j * dr / (4 * k0)

A = diags([-gamma * V_down, 1 - gamma * u_diag, -gamma * V_up], [-1, 0, 1], format='csc')
B = diags([ gamma * V_down, 1 + gamma * u_diag,  gamma * V_up], [-1, 0, 1], format='csc')

# --- 4. Feixe Inicial (Starting Field) ---
W = 30.0 # Largura do feixe (controla a dispersão)
tilt = np.deg2rad(theta0)
psi = np.exp(-((z - zs)**2) / W**2) * np.exp(1j * k0 * (z - zs) * np.sin(tilt))

# --- 5. Marcha de Alcance (Range Marching) ---
idx_min = np.argmin(np.abs(z - z_plot_min))
idx_max = np.argmin(np.abs(z - z_plot_max))
tl_map = np.zeros((idx_max - idx_min, len(r)))

for ir in range(len(r)):
    # Normalização cilíndrica e conversão para dB
    pressao = np.abs(psi) / np.sqrt(r[ir] + 1e-6)
    tl_map[:, ir] = -20 * np.log10(pressao[idx_min:idx_max] + 1e-12)
    
    # Resolver o sistema tridiagonal
    psi = spsolve(A, B.dot(psi))
    
    # Condição de contorno (Superfície Livre)
    psi[0] = 0 

# --- 6. Plotagem Idêntica ao Livro ---
R_mesh, Z_mesh = np.meshgrid(r / 1000.0, z[idx_min:idx_max])

plt.figure(figsize=(10, 4.5))
# Usamos contourf com níveis de 34 a 46 para recriar as "linhas" do gabarito
cont = plt.contourf(R_mesh, Z_mesh, tl_map, levels=np.arange(34, 48, 2), 
                    cmap='magma_r', extend='both')

# Inverter o eixo Y para colocar -800 no topo
plt.gca().invert_yaxis()

# Adicionar a estrela da fonte e a linha do fundo
plt.plot(0, zs, 'k*', markersize=10)
plt.axhline(0, color='black', linewidth=0.8)

plt.colorbar(cont, label='dB')
plt.xlabel('Range (km)')
plt.ylabel('Depth (m)')
plt.title('Simulação PE-IFD: Beam Splitting (F=100Hz, SD=-400m)')
plt.tight_layout()
plt.show()