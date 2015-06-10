# -*- coding: utf-8 -*-
"""
Exemplo do uso do cálculo do ponto de orvalho com saídas gráficas.
O exemplo foi uma aplicação do estudo de caso proposto em : Gupta, B. S., & Lee, M. J. (2013). 
Isobaric vapor-liquid equilibrium for binary systems of toluene+o-xylene, benzene+o-xylene, 
nonane+benzene and nonane+heptane at 101.3kPa. Fluid Phase Equilibria, 352, 86–92.
http://doi.org/10.1016/j.fluid.2013.05.016
"""
from Graficos import Graficos
from Conexao import Componente_Caracterizar, UNIQUAC, VIRIAL
from VLE import VLE
from numpy import array

# Caracterizando os componentes
Comp1 = Componente_Caracterizar('Tolueno',ConfigPsat=('Prausnitz4th',1),T=340.0)
Comp2 = Componente_Caracterizar('o-Xileno',ConfigPsat=('Prausnitz4th',1),T=340.0)
                
Componentes = [Comp1,Comp2]

# Caracterizando os modelos e equações            
model_vap = VIRIAL(Componentes)
model_liq = UNIQUAC(Componentes,340.0,1)

exemplo = VLE('PontoOrvalho_T',Componentes,model_liq,model_vap,z=[0.500,0.500],Temp=340.0,Pressao=1.013,estphi=[1.0,1.0],estBeta = 0.5,tolAlg=1e-10,toleq=1e-4,maxiter=500,z_coordenacao=10.0)
exemplo.run()

# Saídas dos cálculos
print exemplo.Bolha.Temp
print exemplo.Bolha.vapor.comp

# Caracterização do gráfico
    # Caracterização do gráfico teórico
exemplo.Predicao(P = 1.013) 
    # Caracterização dos pontos experimentais
comp_vap = array([0.9999999999999,0.981,0.976,0.961,0.930,0.886,0.871,0.825,0.765,0.693,0.626,0.526,0.436,0.269,0.196,0.147,0.00000000001]) # MOLAR
comp_liq = array([0.9999999999999,0.937,0.920,0.880,0.796,0.693,0.662,0.579,0.500,0.426,0.367,0.291,0.230,0.132,0.094,0.069,0.00000000001]) # MOLAR
comp = [comp_liq,comp_vap]
T = array([383.4,385.0,385.4,386.3,388.1,390.5,391.4,393.9,396.3,399.0,401.3,404.4,407.0,411.4,413.2,414.3,417.5])
    # Incertezas da temperatura
yerr = array([0.1]*17)

Graficos('Txy',Componentes, exemplo.Composicao_x1,exemplo.Temperatura_Ponto_Bolha, exemplo.Temperatura_Ponto_Orvalho, x_experimentais = comp ,y_experimentais = T, y_incertezas = yerr, P = 1.013 )
