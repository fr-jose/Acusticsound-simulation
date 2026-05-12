import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import spdiags
from scipy.sparse.linalg import spsolve

# --- 1. Configurações do Domínio ---
Profundidade = 800
Distancia = 5000
dz = 2
dr = 20 
f = 100
c0 = 1500
k0 = 2 * np.pi * f / c0

z = np.arange(0, Profundidade + dz, dz)
r = np.arange(0, Distancia + dr, dr)
Nz = len(z)
Nr = len(r)

# --- 2. Coeficientes de Greene (Ângulo Largo) ---
a0, a1, b0, b1 = 0.99987, 0.79624, 1.0, 0.30102
w1 = b0 + (1j * k0 * dr / 2) * (a0 - b0)
w1_star = b0 - (1j * k0 * dr / 2) * (a0 - b0)
w2 = b1 + (1j * k0 * dr / 2) * (a1 - b1)
w2_star = b1 - (1j * k0 * dr / 2) * (a1 - b1)

# --- 3. Condição Inicial (Fonte) ---
zs = 400 
W = 20 
psi = np.exp(-(z - zs)**2 / W**2).astype(complex)

# --- 4. Preparação de Dados para o Gráfico ---
U_final = np.zeros((Nz, Nr))
fundo_perfil = np.zeros(Nr)

print("Iniciando processamento da marcha...")

# --- 5. Loop de Marcha (Execução Direta) ---
for m in range(Nr - 1):
    # Lógica da RAMPA
    inicio_rampa = 1000 / dr
    if m > inicio_rampa:
        prof_fundo = 800 - (m - inicio_rampa) * 2.5
        prof_fundo = max(prof_fundo, 300)
    else:
        prof_fundo = 800
    
    fundo_perfil[m] = prof_fundo

    # Montagem do Sistema
    u_diag = np.zeros(Nz, dtype=complex)
    u_hat_diag = np.zeros(Nz, dtype=complex)
    comum = (k0**2 * dz**2) / 2

    for i in range(Nz):
        if z[i] < prof_fundo:
            n_atual = 1.0 
        else:
            # n_atual ajustado para valor físico realista (sedimento)
            n_atual = 1.5 + 0.1j 

        u_diag[i] = 2 * (comum * (w1_star/w2_star) - 1) + 2 * comum * (n_atual**2 - 1)
        u_hat_diag[i] = 2 * (comum * (w1/w2) - 1) + 2 * comum * (n_atual**2 - 1)

    L = spdiags([np.ones(Nz), u_diag, np.ones(Nz)], [-1, 0, 1], Nz, Nz).tocsc()
    R = spdiags([np.ones(Nz), u_hat_diag, np.ones(Nz)], [-1, 0, 1], Nz, Nz).tocsc()
    
    # Passo de propagação
    psi = spsolve(L, (w2 / w2_star) * (R @ psi))
    
    # Armazena Intensidade em dB
    U_final[:, m+1] = 20 * np.log10(np.abs(psi) + 1e-10)

# Ajuste do fundo para o último ponto
fundo_perfil[-1] = fundo_perfil[-2]

# --- 6. Visualização Final ---
# Normalização para o máximo ser 0 dB
U_final -= np.max(U_final)

plt.figure(figsize=(12, 6))
img = plt.imshow(U_final, extent=[0, Distancia, Profundidade, 0], 
                 cmap='magma', aspect='auto', vmin=-40, vmax=0)

# Desenha a linha branca do fundo (montanha)
plt.plot(r, fundo_perfil, 'w--', lw=1.5, label='Perfil do Fundo')

plt.colorbar(img, label='Perda de Transmissão (dB)')
plt.title(f'Propagação Acústica Greene (Padé 1,1) - Rampa no Fundo')
plt.xlabel('Distância r (m)')
plt.ylabel('Profundidade z (m)')
plt.legend()
plt.tight_layout()
plt.show()