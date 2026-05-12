import numpy as np

# Dados
frequencia = 100.0        # Hz
c0 = 1500.0               # Velocidade na fonte (z=0)
profundidade_max = 800.0  # Metros
distancia_max = 10000.0   # Metros 
gradiente = 0.016         # m/s por metro
c_fundo = 1655.0          # Velocidade no sedimento (m/s)
densidade_fundo = 1500    # kg/m^3
atenuacao_fundo = 0.0     # Sem atenuação
theta_0_deg = 5.28  
theta_0_rad = np.deg2rad(theta_0_deg)

# Constantes para a simulação
k0 = 2 * np.pi * frequencia / c0  # Número de onda de referência

def ref(z):
    "Termo de refração"
    cz = c0 + gradiente * z
    n_z = c0 / cz
    return k0**2 * (n_z**2 - 1)

# Exemplo de saída para conferência
print(f"Número de onda k0: {k0:.4f}")
print(f"Ângulo de inclinação: {theta_0_deg}°")