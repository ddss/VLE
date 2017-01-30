# -*- coding: utf-8 -*-
"""
Validação das rotinas que adicionam a equação SRK
Exemplo de calculo bolha P para comparar com dados experimentais. Dados experimentais presente no artigo:
Ho J.C.K.: Part III. System Ethanol - Cyclohexane at Atmospheric Pressure. J.Chem.Eng.Data 8 (1963) 549-559
"""



from Conexao import Componente_Caracterizar, UNIQUAC, SRK
from VLE import VLE
from numpy import array,zeros
from xlwt import Workbook

# Caracterizando os componentes
Comp1 = Componente_Caracterizar('Etanol',ConfigPsat=('Prausnitz4th',1),T=328.15)
Comp2 = Componente_Caracterizar('Benzeno',ConfigPsat=('Prausnitz4th',1),T=328.15)

Componentes = [Comp1,Comp2]

# Caracterizando os modelos e equações

model_vap = SRK(Componentes)
model_liq = UNIQUAC(Componentes,328.15,1)

vetorComposicao=[0.057,0.159,0.266,0.367,0.526,0.632,0.743,0.83,0.916]

a=array(zeros(len(vetorComposicao)))
b=array(zeros(len(vetorComposicao)))

for i in  range(len(vetorComposicao)):
    exemplo = VLE('PontoBolha_P', Componentes, model_liq, model_vap, [vetorComposicao[i], 1 - vetorComposicao[i]], Temp=328.15, Pressao=1.013,tolAlg=1e-20)
    exemplo.run()
    a[i]=exemplo.Bolha.Pressao
    b[i]=exemplo.Bolha.comp_molar[0]

#Impressão dos dados

planilha=Workbook(encoding='utf-8')

aba1=planilha.add_sheet('aba1')

for i in range(len(vetorComposicao)):
    aba1.write(i,0,a[i])
    aba1.write(i,1,b[i])

planilha.save('impressão.xls')
