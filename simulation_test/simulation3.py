import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve
import matplotlib.pyplot as plt

# ============================================================
# 1. Parâmetros físicos e da grade
# ============================================================
freq = 100.0
c0 = 1500.0
k0 = 2 * np.pi * freq / c0
grad = 0.016                # m/s/m
theta0_deg = 5.28           # inclinação descendente
theta0 = np.deg2rad(theta0_deg)

# Domínio vertical
z_surf = -800.0             # superfície livre
z_intf = 0.0                # interface água-fundo
z_plot_max = 200.0          # até onde plotar (mostrar sedimento raso)
z_max = 400.0               # domínio total (camada absorvente)
dz = 0.5
dr = 2.0
r_max = 10000.0

z = np.arange(z_surf, z_max + dz, dz)
Nz = len(z)
r = np.arange(0, r_max + dr, dr)
Nr = len(r)

# Índices importantes
idx_intf = np.argmin(np.abs(z - z_intf))
idx_plot_end = np.argmin(np.abs(z - z_plot_max))

# ============================================================
# 2. Meio: velocidade, densidade, atenuação (sponge)
# ============================================================
c = np.zeros(Nz)
rho = np.zeros(Nz)
alpha = np.zeros(Nz)

# Água (z ≤ 0)
mask_agua = z <= 0
c[mask_agua] = 1500.0 + 0.016 * (z[mask_agua] + 400.0)  # 1500 na fonte, 1506.4 no fundo
rho[mask_agua] = 1000.0

# Sedimento (z > 0)
mask_sed = z > 0
c[mask_sed] = 1655.0
rho[mask_sed] = 1500.0
# Sponge layer a partir de 200 m
sponge_start = 200.0
sponge_strength = 0.8
mask_sponge = z > sponge_start
alpha[mask_sponge] = sponge_strength * ((z[mask_sponge] - sponge_start) / (z_max - sponge_start)) ** 2

# Índice de refração complexo
n_sq = (c0 / c) ** 2 + 1j * alpha

# ============================================================
# 3. Coeficientes IFD (Claerbout wide-angle)
# ============================================================
a0, a1 = 1.0, 0.75
b0, b1 = 1.0, 0.25

w2 = b1 + (1j * k0 * dr / 2) * (a1 - b1)
w2c = np.conj(w2)

inv_w2 = 1.0 / w2
inv_w2c = 1.0 / w2c
k0dz2_half = (k0 * dz) ** 2 / 2.0

u = np.zeros(Nz, dtype=complex)
u_hat = np.zeros(Nz, dtype=complex)
v = np.zeros(Nz, dtype=complex)

for l in range(Nz):
    # Determinar densidades e índices acima/abaixo do ponto l
    if l == 0:
        rho1, rho2 = rho[0], rho[1] if Nz > 1 else rho[0]
        n1_sq, n2_sq = n_sq[0], n_sq[1] if Nz > 1 else n_sq[0]
    elif l == Nz - 1:
        rho1, rho2 = rho[-2], rho[-1]
        n1_sq, n2_sq = n_sq[-2], n_sq[-1]
    else:
        # Meio acima: ponto l-1, meio abaixo: ponto l+1 (ou o próprio l para a interface)
        # Na interface, a densidade logo acima é 1000, abaixo é 1500.
        if z[l] <= 0 and z[l+1] <= 0:
            rho_a, rho_b = 1000.0, 1000.0
        elif z[l] <= 0 and z[l+1] > 0:
            rho_a, rho_b = 1000.0, 1500.0
        elif z[l] > 0 and z[l-1] <= 0:
            rho_a, rho_b = 1000.0, 1500.0  # ponto da interface pertence ao sedimento
        else:
            rho_a, rho_b = 1500.0, 1500.0
        rho1, rho2 = rho_a, rho_b
        n1_sq = (c0 / c[l-1])**2 if z[l-1] <= 0 else (c0/1655.0)**2
        n2_sq = (c0 / c[l+1])**2 if z[l+1] <= 0 else (c0/1655.0)**2

    termo1 = (rho1 + rho2) / rho2
    v[l] = rho1 / rho2

    u[l] = termo1 * (k0dz2_half * inv_w2c - 1) + k0dz2_half * ((n1_sq - 1) + v[l] * (n2_sq - 1))
    u_hat[l] = termo1 * (k0dz2_half * inv_w2 - 1) + k0dz2_half * ((n1_sq - 1) + v[l] * (n2_sq - 1))

# ============================================================
# 4. Montagem das matrizes tridiagonais (pontos internos 1..N-2)
# ============================================================
N_int = Nz - 2
diag_u = u[1:-1]          # índices 1..N-2
diag_hat = u_hat[1:-1]
sub = np.ones(N_int - 1)  # subdiagonal = 1

# CORREÇÃO: superdiagonal = v[l] para l = 1..N-3 (linhas 1 a N-3)
# superdiagonal de tamanho N_int-1
vi_sup = v[1:-2]          # v[1], v[2], ..., v[N-3]   (comprimento N_int-1)

L = diags([sub, diag_u, vi_sup], [-1, 0, 1], format='csc')
R = diags([sub, diag_hat, vi_sup], [-1, 0, 1], format='csc')

# ============================================================
# 5. Campo inicial – feixe gaussiano com largura angular 2°
# ============================================================
sigma_theta = np.deg2rad(2.0)   # halfwidth
angulos = np.linspace(-15, 15, 301) * np.pi / 180
pesos = np.exp(-0.5 * ((angulos - theta0) / sigma_theta) ** 2)
pesos /= pesos.sum()            # normalização

psi = np.zeros(Nz, dtype=complex)
for ang, wgt in zip(angulos, pesos):
    psi += wgt * np.exp(1j * k0 * (z - (-400.0)) * np.sin(ang))
psi /= np.max(np.abs(psi))      # amplitude máxima = 1

# Condições de contorno
psi[0] = 0.0
psi[-1] = 0.0

# ============================================================
# 6. Marcha em alcance
# ============================================================
campo = np.zeros((Nr, Nz), dtype=complex)

for ir in range(Nr):
    campo[ir, :] = psi
    rhs = R.dot(psi[1:-1])
    psi_novo_interno = spsolve(L, rhs)
    psi[1:-1] = psi_novo_interno
    psi[0] = 0.0
    psi[-1] = 0.0

# ============================================================
# 7. Perda de transmissão (dB)
# ============================================================
amp = np.abs(campo)
amp_max = np.max(amp[0, :])
with np.errstate(divide='ignore'):
    tl = -20.0 * np.log10(amp / amp_max + 1e-12)

tl_clip = np.clip(tl, 34, 46)

# ============================================================
# 8. Visualização (profundidade negativa no topo)
# ============================================================
z_plot = z[:idx_plot_end+1]
tl_plot = tl_clip[:, :idx_plot_end+1].T   # transpor para (z, r)

R_km = r / 1000.0
plt.figure(figsize=(10, 4.5))
levels = np.arange(32, 40, 0.5)
cont = plt.contourf(R_km, z_plot, tl_plot, levels=levels, cmap='magma_r', extend='both')
plt.colorbar(cont, label='dB', ticks=np.arange(32, 40, 0.5))

plt.gca().invert_yaxis()
plt.axhline(0, color='black', linewidth=0.8, linestyle='--', label='Interface água‑fundo')
plt.plot(0, -400, 'k*', markersize=12, label='Fonte (100 Hz)')
plt.xlabel('Distância (km)')
plt.ylabel('Profundidade (m)')
plt.title('Beam Splitting — PE IFD Claerbout\nF=100 Hz, θ=5.28°, Grad=0.016 m/s/m, Fundo 1655 m/s')
plt.legend(loc='upper right')
plt.tight_layout()
plt.show()