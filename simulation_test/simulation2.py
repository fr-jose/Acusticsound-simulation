import numpy as np
from scipy.sparse import diags, eye
from scipy.sparse.linalg import spsolve
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

# ============================================================
# 1. Parâmetros físicos conforme o enunciado
# ============================================================
freq = 100.0                     # Hz
c0   = 1500.0                    # velocidade na fonte (m/s)
grad = 0.016                     # gradiente m/s por metro de profundidade
zs   = 400.0                     # profundidade da fonte (m)
H    = 800.0                     # profundidade do fundo (interface)
theta0_deg = 5.28                # ângulo de inclinação descendente
theta0 = np.deg2rad(theta0_deg)

# Número de onda de referência
k0 = 2 * np.pi * freq / c0

# ============================================================
# 2. Domínio discretizado (profundidade z e alcance r)
# ============================================================
dz = 1.0                         # resolução em profundidade (m)
z  = np.arange(0, H + dz, dz)    # 0 = superfície, H = fundo
Nz = len(z)

dr = 2.0                         # passo em alcance (m)
r_max = 10000.0                  # 10 km
r = np.arange(0, r_max + dr, dr)
Nr = len(r)

# ============================================================
# 3. Perfil de velocidade e índice de refração (apenas água)
# ============================================================
c_agua = c0 + grad * (z - zs)    # c(400)=1500, c(800)=1506.4
n = c0 / c_agua                  # índice de refração

# ============================================================
# 4. Condição inicial – feixe gaussiano inclinado
# ============================================================
w = 20.0                         # largura do feixe (m)
envoltoria = np.exp(-((z - zs) / w) ** 2)
fase = np.exp(1j * k0 * (z - zs) * np.sin(theta0))
psi = envoltoria * fase
psi[0] = 0.0                     # condição de superfície livre (p=0)

# ============================================================
# 5. Montagem das matrizes Crank‑Nicolson
# ============================================================
inv_dz2 = 1.0 / dz**2
diag_T = -2.0 * inv_dz2 + k0**2 * (n**2 - 1.0)
off_diag = np.ones(Nz - 1) * inv_dz2

T = diags([off_diag, diag_T, off_diag], [-1, 0, 1], format='csc')

a = 1j * dr / (4.0 * k0)
I = eye(Nz, format='csc')

L = I - a * T
R = I + a * T

# ----------------------------------------------------------
# Condições de contorno:
#   z = 0 (superfície): Dirichlet → psi(0) = 0
#   z = H (fundo)     : Neumann   → ∂ψ/∂z = 0 (refletor rígido)
# ----------------------------------------------------------
# Superfície: anula primeira linha e impõe identidade
L[0, :] = 0
L[0, 0] = 1.0
R[0, :] = 0
R[0, 0] = 1.0

# Fundo: impõe derivada zero usando fórmula de segunda ordem
#         ψ[M] = ψ[M-1]  (onde M = Nz-1)
L[-1, :] = 0
L[-1, -2] = -1.0
L[-1, -1] = 1.0
R[-1, :] = 0
R[-1, -2] = -1.0
R[-1, -1] = 1.0

# ============================================================
# 6. Marcha em alcance
# ============================================================
campo = np.zeros((Nr, Nz), dtype=np.complex128)

for idx in range(Nr):
    campo[idx, :] = psi
    psi = spsolve(L, R.dot(psi))
    psi[0] = 0.0      # mantém superfície livre a cada passo

# ============================================================
# 7. Perda de transmissão (TL) em dB
# ============================================================
amp = np.abs(campo)
amp_max = np.max(amp[0, :])
with np.errstate(divide='ignore'):
    tl = -20.0 * np.log10(amp / amp_max + 1e-12)

tl_min, tl_max = 30, 70
tl = np.clip(tl, tl_min, tl_max)

# ============================================================
# 8. Visualização
# ============================================================
fig, ax = plt.subplots(figsize=(11, 5))

cmap = plt.cm.jet_r
norm = mcolors.Normalize(vmin=tl_min, vmax=tl_max)

img = ax.imshow(tl.T, extent=[0, r_max/1000, H, 0],
                aspect='auto', cmap=cmap, norm=norm, origin='upper')

# Interface (fundo) e superfície
ax.axhline(y=H, color='white', linestyle='--', linewidth=1.2, alpha=0.9)
ax.text(r_max/1000 * 0.98, H - 25, 'Fundo (interface)', color='white',
        ha='right', fontsize=9, fontstyle='italic')
ax.axhline(y=0, color='lightgray', linewidth=1)

# Fonte
ax.plot(0, zs, marker='*', color='white', markersize=14, markeredgecolor='black', zorder=5)

# Rótulos
ax.set_xlabel('Distância (km)', fontsize=12)
ax.set_ylabel('Profundidade (m)', fontsize=12)
ax.set_title(f'Beam Splitting – Fonte {freq} Hz, θ = {theta0_deg}°\n'
             f'Gradiente = {grad} m/s/m, Fundo rígido (reflexão total)', fontsize=12, fontweight='bold')

cbar = plt.colorbar(img, ax=ax, label='Perda de Transmissão (dB)')
cbar.ax.invert_yaxis()          # cores quentes no topo = menor TL = sinal mais intenso

ax.invert_yaxis()               # profundidade cresce para baixo
plt.tight_layout()
plt.show()