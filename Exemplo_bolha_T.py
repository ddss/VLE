# -*- coding: utf-8 -*-
"""
Exemplo do uso do cálculo do ponto de bolha com saídas gráficas.
O exemplo foi uma aplicação do estudo de caso proposto em : Ku, H.-C., & Tu, C.-H. (2005). 
Isobaric vapor–liquid equilibria for mixtures of acetone, ethanol, and 2,2,4-trimethylpentane
at 101.3kPa. Fluid Phase Equilibria, 231, 99–108. 
http://doi.org/10.1016/j.fluid.2005.01.007
"""
from Graficos import Graficos
from Conexao import Componente_Caracterizar, UNIQUAC, VIRIAL
from VLE import VLE
from numpy import array

# Caracterizando os componentes
Comp1 = Componente_Caracterizar('Acetona',ConfigPsat=('Prausnitz4th',1),T=340.0)
Comp2 = Componente_Caracterizar('Etanol',ConfigPsat=('Prausnitz4th',1),T=340.0)

Componentes = [Comp1,Comp2]

# Caracterizando os modelos e equações            
model_vap = VIRIAL(Componentes)
model_liq = UNIQUAC(Componentes,340.0,1)

exemplo = VLE('PontoBolha_T',Componentes,model_liq,model_vap,z=[0.044,1-0.044],Temp=340.0,Pressao=1.013,estphi=[1.0,1.0],estBeta = 0.5,tolAlg=1e-10,toleq=1e-4,maxiter=500,z_coordenacao=10.0)
exemplo.run()

# Saídas dos cálculos
print exemplo.Bolha.Temp
print exemplo.Bolha.comp_molar

# Caracterização do gráfico
    # Caracterização do gráfico teórico
exemplo.Predicao('pressao',1.013) 
    # Caracterização dos pontos experimentais
comp_vap = [0.00000000001,0.044,0.085,0.132,0.175,0.227,0.276,0.328,0.380,0.433,0.481,0.538,0.586,0.639,0.688,0.745,0.796,0.840,0.901,0.951,0.9999999999999] # MOLAR
comp_liq = [0.00000000001,0.129,0.232,0.325,0.403,0.463,0.517,0.571,0.609,0.651,0.685,0.719,0.749,0.777,0.803,0.835,0.867,0.898,0.929,0.965,0.9999999999999] # MOLAR
T = [351.44,348.91,346.65,344.53,342.69,341.17,339.81,338.39,337.32,336.24,335.37,334.51,333.78,333.06,332.43,331.81,331.20,330.66,330.17,329.69,329.26]
    # Incertezas da temperatura
yerr = array([0.05]*21).tolist()
Graph = Graficos(x_experimentais = comp_liq, y_experimentais=comp_vap ,T_exp = T, T_incertezas = yerr, P_exp = 1.013 )
Graph.T_x_y(exemplo,1.013)
    