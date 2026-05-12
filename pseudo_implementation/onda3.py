import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import spdiags
from scipy.sparse.linalg import spsolve

# --- 1. Configurações do Domínio ---
Profundidade = 800
Distancia = 5000
dz = 1 
dr = 20 
f = 100
c0 = 1500
k0 = 2 * np.pi * f / c0

z = np.arange(0, Profundidade + dz, dz)
r = np.arange(0, Distancia + dr, dr)
Nz = len(z)
Nr = len(r)

# --- 2. Escolha do Estilo do Feixe ---
# Opções: 'gaussian', 'greene', 'generalized'
tipo_feixe = 'generalized' 

zs = 400      # Profundidade da fonte (m)
theta_1 = 10  # Meia-largura do feixe (graus) - para o Generalized
theta_2 = 15  # Inclinação do feixe (graus) - 0 é reto, positivo aponta para baixo

def gerar_fonte():
    if tipo_feixe == 'gaussian':
        # Gaussiana clássica (Tappert)
        W = 20
        return np.exp(-(z - zs)**2 / W**2).astype(complex)

    elif tipo_feixe == 'greene':
        # Fonte de Greene otimizada para ângulo largo
        term1 = 1.4467 - 0.4201 * k0**2 * (z - zs)**2
        term2 = np.exp(-(k0**2 * (z - zs)**2) / 3.0512)
        return (np.sqrt(k0) * term1 * term2).astype(complex)

    elif tipo_feixe == 'generalized':
        # Gaussiana Generalizada (Permite inclinação/Tilt)
        th1_rad = np.radians(theta_1)
        th2_rad = np.radians(theta_2)
        amp = np.sqrt(k0) * np.tan(th1_rad)
        exp_gauss = np.exp(-(k0**2 / 2) * (z - zs)**2 * np.tan(th1_rad)**2)
        exp_tilt = np.exp(1j * k0 * (z - zs) * np.sin(th2_rad))
        return (amp * exp_gauss * exp_tilt).astype(complex)

psi = gerar_fonte()

# --- 3. Coeficientes de Greene (Otimizados para PE) ---
a0, a1, b0, b1 = 0.99987, 0.79624, 1.0, 0.30102
w1 = b0 + (1j * k0 * dr / 2) * (a0 - b0)
w1_star = b0 - (1j * k0 * dr / 2) * (a0 - b0)
w2 = b1 + (1j * k0 * dr / 2) * (a1 - b1)
w2_star = b1 - (1j * k0 * dr / 2) * (a1 - b1)

# --- 4. Simulação de Marcha ---
U_final = np.zeros((Nz, Nr))
fundo_perfil = np.zeros(Nr)

print(f"Calculando propagação com fonte: {tipo_feixe}...")

for m in range(Nr - 1):
    # Definição da rampa no fundo
    inicio_rampa = 1200 / dr
    prof_fundo = 800 - (m - inicio_rampa) * 3.5 if m > inicio_rampa else 800
    prof_fundo = max(prof_fundo, 250)
    fundo_perfil[m] = prof_fundo

    u_diag = np.zeros(Nz, dtype=complex)
    u_hat_diag = np.zeros(Nz, dtype=complex)
    comum = (k0**2 * dz**2) / 2

    for i in range(Nz):
        # n_atual = 1 na água, 1.6 + absorção no fundo
        n_atual = 1.0 if z[i] < prof_fundo else 1.6 + 0.3j
        
        u_diag[i] = 2 * (comum * (w1_star/w2_star) - 1) + 2 * comum * (n_atual**2 - 1)
        u_hat_diag[i] = 2 * (comum * (w1/w2) - 1) + 2 * comum * (n_atual**2 - 1)

    L = spdiags([np.ones(Nz), u_diag, np.ones(Nz)], [-1, 0, 1], Nz, Nz).tocsc()
    R = spdiags([np.ones(Nz), u_hat_diag, np.ones(Nz)], [-1, 0, 1], Nz, Nz).tocsc()
    
    # Próximo passo de distância
    psi = spsolve(L, (w2 / w2_star) * (R @ psi))
    U_final[:, m+1] = 20 * np.log10(np.abs(psi) + 1e-10)

# --- 5. Visualização ---
U_final -= np.max(U_final) # Normaliza para 0 dB no pico

plt.figure(figsize=(12, 6))
# Usando 'jet' ou 'magma' para contraste
plt.imshow(U_final, extent=[0, Distancia, Profundidade, 0], cmap='jet', aspect='auto', vmin=-35, vmax=0)
plt.plot(r, fundo_perfil, 'w--', lw=1.5, label='Perfil do Fundo')
plt.colorbar(label='Intensidade Relativa (dB)')

titulo = f'Fonte: {tipo_feixe.capitalize()}'
if tipo_feixe == 'generalized':
    titulo += f' (Inclinação: {theta_2}°, Abertura: {theta_1}°)'
plt.title(titulo)

plt.xlabel('Distância r (m)')
plt.ylabel('Profundidade z (m)')
plt.legend()
plt.show()