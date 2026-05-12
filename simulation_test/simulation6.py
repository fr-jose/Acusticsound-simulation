import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import solve_banded

def run_pe_simulation():
    print("Iniciando a simulação da Equação Parabólica (PE)...")
    
    # ---------------------------------------------------------
    # 1. Parâmetros Físicos e do Ambiente (Física do Capítulo 6)
    # ---------------------------------------------------------
    f = 50.0               # Frequência em Hz
    c0 = 1500.0            # Velocidade de referência da água (m/s)
    cb = 1700.0            # Velocidade do som no fundo (areia/rocha)
    rho_w = 1.0            # Densidade da água (g/cm^3)
    rho_b = 1.5            # Densidade do fundo (g/cm^3)
    omega = 2 * np.pi * f
    k0 = omega / c0        # Número de onda de referência
    zs = 200.0             # Profundidade da fonte (m)
    
    # Malha Espacial
    D = 1000.0             # Profundidade total do domínio (m)
    dz = 2.0               # Passo em z (m)
    z = np.arange(0, D + dz, dz)
    nz = len(z)
    
    Rmax = 20000.0         # Alcance máximo (m) -> 20 km
    dr = 10.0              # Passo em r (m)
    r = np.arange(0, Rmax + dr, dr)
    nr = len(r)
    
    # ---------------------------------------------------------
    # 2. Construindo o Monte Submarino (Seamount Batimetria)
    # ---------------------------------------------------------
    # O fundo começa a 800m, sobe para 400m no meio do caminho, e desce.
    bottom_depth = np.ones(nr) * 800.0
    seamount_center = 10000.0
    seamount_width = 3000.0
    seamount_height = 400.0
    
    for ir in range(nr):
        dist_from_center = abs(r[ir] - seamount_center)
        if dist_from_center < seamount_width:
            # Perfil gaussiano/cossenoide para o monte
            bottom_depth[ir] -= seamount_height * (1 + np.cos(np.pi * dist_from_center / seamount_width)) / 2

    # ---------------------------------------------------------
    # 3. Inicialização da Fonte (Gaussian Starter - Eq 6.100)
    # ---------------------------------------------------------
    psi = np.zeros((nz, nr), dtype=complex)
    fac = 1.0  # Alargamento da Gaussiana
    psi[:, 0] = np.sqrt(k0 / fac) * np.exp(-0.5 * (k0 / fac)**2 * (z - zs)**2)
    
    # ---------------------------------------------------------
    # 4. Marcha no Alcance (O motor do factor_Mod / backsub_Mod)
    # ---------------------------------------------------------
    # O método de Crank-Nicolson resolve: (I - 0.5*i*dr*A) * psi_next = (I + 0.5*i*dr*A) * psi_curr
    
    diag_indices = [0, 1, 2] # Para matriz tridiagonal no solve_banded do Scipy
    
    for ir in range(nr - 1):
        if ir % 500 == 0:
            print(f"Calculando passo {ir}/{nr-1} (Alcance: {r[ir]/1000:.1f} km)...")
            
        # Perfil de velocidade e densidade para este passo (r)
        c = np.where(z <= bottom_depth[ir], c0, cb)
        
        # O operador de refração n^2 - 1
        n = c0 / c
        k2 = (k0**2) * (n**2 - 1.0)
        
        # Absorção artificial no fundo para evitar reflexões espúrias do fim do grid
        sponge_start = 900.0
        alpha = np.where(z > sponge_start, 0.05 * (z - sponge_start), 0.0)
        k2 = k2 + 1j * alpha
        
        # Termos da derivada segunda (matriz tridiagonal)
        coef_sub_sup = 1.0 / (dz**2)
        coef_diag = -2.0 / (dz**2) + k2
        
        # Matrizes do passo explícito (Lado Direito) e implícito (Lado Esquerdo)
        coef_r = 0.5 * 1j * dr / k0
        
        R_diag = 1.0 + coef_r * coef_diag
        R_sub_sup = coef_r * coef_sub_sup
        
        L_diag = 1.0 - coef_r * coef_diag
        L_sub_sup = -coef_r * coef_sub_sup
        
        # Construindo o vetor do lado direito (B = R * psi)
        rhs = np.zeros(nz, dtype=complex)
        rhs[0] = R_diag[0] * psi[0, ir] + R_sub_sup * psi[1, ir]
        rhs[1:-1] = R_sub_sup * psi[:-2, ir] + R_diag[1:-1] * psi[1:-1, ir] + R_sub_sup * psi[2:, ir]
        rhs[-1] = R_sub_sup * psi[-2, ir] + R_diag[-1] * psi[-1, ir]
        
        # Configurando a matriz de banda para o Scipy (L)
        L_banded = np.zeros((3, nz), dtype=complex)
        L_banded[0, 1:] = L_sub_sup    # Sup-diagonal
        L_banded[1, :] = L_diag        # Diagonal principal
        L_banded[2, :-1] = L_sub_sup   # Sub-diagonal
        
        # Condições de Contorno de Dirichlet (Pressão zero na superfície e no fundo profundo)
        L_banded[1, 0] = 1.0
        L_banded[0, 1] = 0.0
        rhs[0] = 0.0
        
        L_banded[1, -1] = 1.0
        L_banded[2, -2] = 0.0
        rhs[-1] = 0.0
        
        # Resolução do sistema tridiagonal (Equivalente ao Factor/Backsub em Fortran)
        psi[:, ir + 1] = solve_banded((1, 1), L_banded, rhs)

    # ---------------------------------------------------------
    # 5. Processamento e Gráfico (Perda de Transmissão - TL)
    # ---------------------------------------------------------
    print("Gerando o gráfico...")
    
    # Espalhamento Cilíndrico (Hankel assintótico)
    # Ignoramos r=0 para evitar divisão por zero
    r_mat = np.tile(r, (nz, 1))
    r_mat[r_mat == 0] = 1e-10 
    
    pressure = psi * np.sqrt(2.0 / (np.pi * k0 * r_mat))
    TL = -20 * np.log10(np.abs(pressure) + 1e-10)
    
    # Criando a figura com o visual de artigo acústico
    plt.figure(figsize=(12, 6))
    
    # Plot do campo acústico
    plt.pcolormesh(r / 1000, z, TL, cmap='jet', vmin=40, vmax=90, shading='auto')
    
    # Plot da linha do fundo (Batimetria)
    plt.plot(r / 1000, bottom_depth, 'k-', linewidth=3, label='Batimetria (Seamount)')
    plt.fill_between(r / 1000, bottom_depth, D, color='saddlebrown', alpha=0.8)
    
    plt.gca().invert_yaxis()
    plt.colorbar(label='Perda de Transmissão (dB)')
    plt.xlabel('Alcance (km)')
    plt.ylabel('Profundidade (m)')
    plt.title('Propagação Acústica sobre Monte Submarino (Equação Parabólica)')
    plt.legend(loc='lower left')
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    run_pe_simulation()