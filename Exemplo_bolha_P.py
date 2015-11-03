# -*- coding: utf-8 -*-
"""
Exemplo do uso do cálculo do ponto de bolha com saídas gráficas.
O exemplo foi uma aplicação do estudo de caso proposto em : V, E. S. P. B., & Oracz, P.
(1987). For binary mixtures of methanol, ethanol, 1-propanol, 35, 253–278.
"""
from matplotlib import use
use('Agg')

from Graficos import Graficos
from Conexao import Componente_Caracterizar, UNIQUAC, VIRIAL,SRK, PengRobinson
from VLE import VLE
from numpy import array

# Caracterizando os componentes
Comp1 = Componente_Caracterizar('Etanol',ConfigPsat=('Prausnitz4th',1),T=340.0)
Comp2 = Componente_Caracterizar('Benzeno',ConfigPsat=('Prausnitz4th',1),T=340.0)

Componentes = [Comp1,Comp2]

# Caracterizando os modelos e equações            
#model_vap = VIRIAL(Componentes)
#model_vap = SRK(Componentes)
model_vap = PengRobinson(Componentes)
model_liq = UNIQUAC(Componentes,313.15,1)

exemplo = VLE('PontoBolha_P',Componentes,model_liq,model_vap,z=[0.3703,1-0.3703],Temp=313.15,Pressao=1.013,tolAlg=1e-20)
exemplo.run()

print exemplo.Bolha.Pressao
print exemplo.Bolha.comp_molar
# Caracterização do gráfico
    # Caracterização do gráfico teórico
exemplo.Predicao('temperatura',313.15)
    # Caracterização dos pontos experimentais
comp_vap = [0.00000000001,0.1725,0.2630,0.3072,0.2746,0.3072,0.3703,0.3736,0.3872,0.4023,0.4215,0.4452,0.4778,0.5284,0.5737,0.5934,0.6596,0.7169,0.8029,0.8795,0.9999999999999] # MOLAR
comp_liq = [0.00000000001,0.0293,0.0823,0.1602,0.1942,0.3055,0.3703,0.4357,0.4981,0.5587,0.6218,0.6826,0.7438,0.8103,0.8516,0.8662,0.9052,0.9303,0.9586,0.9776,0.9999999999999] # MOLAR
P = 0.01*array([24.367,28.831,31.544,32.764,33.015,33.381,33.365,33.276,33.121,32.849,32.388,31.704,30.672,28.991,27.464,26.914,24.938,23.446,21.506,19.938,17.897])
P = P.tolist()
# Incertezas da temperatura
yerr = array([0.015]*21).tolist()

Graph = Graficos(x_experimentais = comp_liq, y_experimentais=comp_vap ,P_exp = P, P_incertezas = yerr)
Graph.P_x_y(exemplo,313.15)
