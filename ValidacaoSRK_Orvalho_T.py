# -*- coding: utf-8 -*-
"""
Validação das rotinas que adicionam a equação SRK
Exemplo de calculo orvalho T para comparar com dados experimentais. Dados experimentais presente no  artigo:
van Winkle M.: Vapor-liquid equilibria. Ind.Eng.Chem. 48 (1956) 142-146
"""



from Conexao import Componente_Caracterizar, UNIQUAC, SRK
from VLE import VLE
from numpy import array,zeros
from xlwt import Workbook

# Caracterizando os componentes
Comp1 = Componente_Caracterizar('Acetona',ConfigPsat=('Prausnitz4th',1),T=328.15)
Comp2 = Componente_Caracterizar('Metanol',ConfigPsat=('Prausnitz4th',1),T=328.15)

Componentes = [Comp1,Comp2]

# Caracterizando os modelos e equações

model_vap = SRK(Componentes)


vetorComposicao=[0.08200,0.16100,0.25100,0.33600,0.42300,0.50000,0.58000,0.63900,0.70500,0.74500,0.80600,0.84300]
vetorTemperatura=[336.65,335.35,333.85,332.55,331.25,330.05,329.35,329.05,328.95,328.95,328.95,328.95]

a=array(zeros(len(vetorComposicao)))
b=array(zeros(len(vetorComposicao)))

for i in  range(len(vetorComposicao)):
    model_liq = UNIQUAC(Componentes, vetorTemperatura[i], 1)
    exemplo = VLE('PontoOrvalho_T', Componentes, model_liq, model_vap, [vetorComposicao[i], 1 - vetorComposicao[i]], Temp=vetorTemperatura[i], Pressao=1.013,tolAlg=1e-20)
    exemplo.run()
    a[i]=exemplo.Orvalho.Temp
    b[i]=exemplo.Orvalho.comp_molar[0]

#Impressão dos dados

planilha=Workbook(encoding='utf-8')

aba1=planilha.add_sheet('aba1')

for i in range(len(vetorComposicao)):
    aba1.write(i,0,a[i])
    aba1.write(i,1,b[i])

planilha.save('impressão.xls')