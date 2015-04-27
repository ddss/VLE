# -*- coding: utf-8 -*-
'''
Exemplo de uso do VLE. A mistura utilizada será Acetona-Metanol, os modelos para as fases líquida e de vapor são, respectivamente
UNIQUAC e VIRIAL com a regra de Hayden O'Connel, a temperatura será 330.0 K e a pressão 1.013 bar.
'''

from Conexao import Componente_Caracterizar, UNIQUAC, VIRIAL
from VLE import VLE
T = 330.0
# Caracterizando os componentes
Comp1 = Componente_Caracterizar('Acetona',ConfigPsat=('Prausnitz4th',1),T=330.0)
Comp2 = Componente_Caracterizar('Metanol',ConfigPsat=('Prausnitz4th',1),T=330.0)
                
Componentes = [Comp1,Comp2]

# Caracterizando os modelos e equações            
model_vap = VIRIAL(Componentes)
model_liq = UNIQUAC(Componentes,330.0,1)

# Realizando os cálculos            
CalculoBolha = VLE('PontoBolha_P',Componentes,model_liq,model_vap,z=[0.95,0.05],Temp=330.0,Pressao=1.013,estgama=None,estphi=[1.0,1.0],estBeta = 0.5,tolAlg=1e-10,toleq=1e-4,maxiter=500,z_coordenacao=10.0)
CalculoOrvalho = VLE('PontoOrvalho_P',Componentes,model_liq,model_vap,z=[0.95,0.05],Temp=330.0,Pressao=1.013,estgama=None,estphi=None,estBeta = 0.5,tolAlg=1e-5,toleq=1e-4,maxiter=500,z_coordenacao=10.0)
Calculo_Coef_Atividade = VLE('Coeficiente_Atividade',Componentes,model_liq,model_vap,z=[0.95,0.05],Temp=330.0,Pressao=1.013,estgama=None,estphi=None,estBeta = 0.5,tolAlg=1e-5,toleq=1e-4,maxiter=500,z_coordenacao=10.0)
Calculo_Coef_Fugacidade = VLE('Coeficiente_Fugacidade',Componentes,model_liq,model_vap,z=[0.95,0.05],Temp=330.0,Pressao=1.013,estgama=None,estphi=None,estBeta = 0.5,tolAlg=1e-5,toleq=1e-4,maxiter=500,z_coordenacao=10.0)

# Saídas dos cálculos
print CalculoBolha.Bolha.Pressao
print CalculoOrvalho.Orvalho.Pressao
print Calculo_Coef_Atividade.coefAct
print Calculo_Coef_Fugacidade.coefFug

