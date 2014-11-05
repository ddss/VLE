# -*- coding: utf-8 -*-
"""
@author: Daniel
@Contribuição: Anderson e Caique
@ Data: 27/10/2013
@Atualizado em : 05/09/2014

Rotina para Cálculo do Equilíbrio de fases em Misturas.

Métodos:
    - Second_Virial_Coef: Cálculo do segundo coeficiente do Virial
    - Coeficiente_Atividade: Cálculo do Coeficiente de Atividade
    - Coeficiente_Fugacidade: Cálculo do Coeficiente de Fugacidade
    - Flash: Cáculo de um flash
    - PhiSat: Cálculo do coeficiente de fugacidade nas condições de saturação
    - PontoBolha: Cálculo do ponto de bolha (T conhecido) (Dado: x, T/K -> Cálcula: y, P/bar)
    - PontoOrvalho: Cálculo do ponto de orvalho (T conhecido) (Dado: y, T/K -> Cálcula: x, P/bar)

Referências:
- PRAUSNITZ, J. M. et al. Computer Calculations for multicomponent vapor-liquid and liquid-liquid equilibria. [s.l.] Prendice-Hall, 1980. p. 353
- HAYDEN, J. G.; O’CONNELL, J. P. A Generalized Method for Predicting Second Virial Coefficients. Industrial & Engineering Chemistry Process Design
  and Development, v. 14, n. 3, p. 209–216, jul. 1975.
- SMITH, J. M.; NESS, H. C. VAN; ABBOTT, M. M. Introduction to Chemical Engineering Thermodinamics. 7th. ed. [s.l.] Mc-Graw Hills, [s.d.]. 
"""

from numpy import log, exp, size, mean, abs, zeros

class Fase:
    '''
    ==============================================================================
            Definição:
    ==============================================================================    
    Classe para definição das características das fases líquida e vapor.
    
    ==============================================================================
            Entradas:
    ==============================================================================        
    *Composição (composicao):
        Inserido como *list*.
    *Coeficiente de fugacidade (coeffug):
        Inserido como *list*.
    *Coeficiciente de atividade (coefAct):
        Inserido como *list*.
    '''
    
    def __init__(self,composicao=None,coeffug=None,coefAct = None):
        self.comp      = composicao        
        self.coeffug   = coeffug
        self.coefAct   = coefAct

class Condicao:
    '''
    ==============================================================================
            Definição:
    ==============================================================================    
    Classe para identificar as condições termodinâmicas.
    
    ==============================================================================
            Entradas:
    ==============================================================================        
    *Pressão (P): 
        Dada em bar e inserido como *float*.
    *Temperatura (T): 
        Dada em Kelvin e inserido como *float*.
    *líquido:
        Características da fase líquida, acessada pela classe ``Fase``.
    *vapor:
        Características da fase de vapor, acessada pela classe ``Fase``.
    *Beta:
    '''
    def __init__(self,P,T,liquido,vapor,Beta):
        self.liquido = liquido
        self.vapor   = vapor
        self.Pressao = P
        self.Temp    = T
        self.Beta    = Beta

class VLE:        

    def __init__(self,Algoritmo,Componentes,z=None,Temp=None,Pressao=None,estgama=None,estphi=None,estBeta = 0.5, tolAlg=1e-10, toleq=1e-4, maxiter=100, model_liq=['UNIQUAC',[None]], model_vap=['VIRIAL',['Hayden_OConnel',None]]):    
        '''
        ==============================================================================
        Definição:
        ==============================================================================    
        Classe para realizar os cálculos do equilíbrio líquido-vapor, desde o cálculo de
        coeficiente de atividade aos cálculos de ponto de bolha e de orvalho.
        
        ==============================================================================
        Entradas:
        ==============================================================================        
        *Algoritmo: O método desejado, vide ``Métodos``. Inserido em forma de *string*;
        *Componentes:Nomes dos componentes a serem usados. Inserido em forma de uma lista que contém o
        método ``Componente_Caracterizar`` do ``Thermo_Data_Bank`` para cada um dos componentes.
        Ex.: Componentes = [TDB.Componente_Caracterizar('Acetone',T),TDB.Componente_Caracterizar('Metanol',T)]
        Onde TDB é Thermo_Data_Bank
        *Número de coordenação(z): Número de coordenação, inserido como número. Comumente é um valor inteiro, mas dado em *float*, em geral o valor usado é 10.0 e é um argumento que pode não ser inserido.
        *Temperatura do sistema (Temp): Dada em kelvin, é inserido como *float* e é um argumento que pode não ser inserido
        *Estimativa para gamma (estgama): Inserido em forma de *float* e é um argumento que pode não ser inserido.
        *Estimativa para phi (estphi): Inserido em forma de *float* e e é um argumento que pode não ser inserido.
        *Estimativa para Beta (estBeta): Inserido em forma de *float* e e é um argumento que pode não ser inserido.
        *Tolerância do algoritmo (tolAlg):Tolerância desejada para a operação dos métodos, inserido como número.
        *Tolerância do equilíbrio (toleq): Tolerância desejada para o equilíbrio, inserido como número.
        *Núemro máximo de iterações (maxiter): Número máximo de iterações desejadas para a operação dos métodos, inserido como número.
        *Modelo da fase líquida (model_liq): Modelo utilizado para a fase líquida, inserido em forma de lista, onde o primeiro elemento é o nome do modelo
        e o segundo é uma lista de parâmetros do modelo que são acessados no ``Thermo_Data_Bank``.
        Ex. model_liq = ['UNIQUAC', [a,z_UNIQUAC]]
        *Modelo da fase vapor (model_vap): Modelo utilizado para a fase líquida, inserido em forma de lista, onde o primeiro elemento é o nome do modelo
        e o segundo é uma lista contendo a regra utilizada e os parâmetros do modelo acessados no ``Thermo_Data_Bank``.
        Ex.: model_vap=['VIRIAL',['Hayden_OConnel',Eta]]
        
        ==============================================================================
        Métodos:
        ==============================================================================        
        *Second_Virial_Coef:  Calcula o segundo coeficiente Virial.
        *Coeficiente_Atividade:  Cálculo do coeficiente de atividade.
        *Coeficiente_Fugacidade: Cálculo do coeficietne de fugacidade.
        *Flash: Cálculo de flash.
        *Phisat: Cálculo de phisat (coeficiente de fugacidade nas condições de saturação).
        *PontoBolha: Cálculo do ponto de bolha.
        *PontoOrvalho:  Cálculo do ponto de orvalho.
        '''
        # Definindo o __init__ de VLE:
        self.Algoritmo = Algoritmo  # Algoritmo a ser utilizado: Flash, PontoBolha, PontoOrvalho, Coeficiente_Atividade, Coeficiente_Fugacidade
        self.NC = size(Componentes) # Número de componentes da mistura

        self.Componente = Componentes

        # Definição dos modelos termodinâmicos:
        self.model_liq       = model_liq[0] # Modelo para a fase líquida: UNIQUAC ou NRTL
        self.model_liq_param = model_liq[1] # Parâmetros para o modelo da fase líquida: UNIQUAC(a,coordnumber), NRTL(t,G)
        self.model_vap       = model_vap[0] # Modelo para a fase vapor: VIRIAl
        self.model_vap_param = model_vap[1] # Parâmetros para o modelo da fase vapor: VIRIAL(Hayden_OConnel, Eta)

        # Condições do VLE
        self.z       = z       # Composição global
        self.Pressao = Pressao # Pressão em bar
        self.Temp    = Temp    # Temperatura em K
        
        # Estimativas iniciais
        if estgama == None:        
            self.estgama = [1.0]*len(self.Componente) # estimativa para o coeficiente de atividade
        else:
            self.estgama = estgama
            
        if estphi == None:
            self.estphi  = [1.0]*len(self.Componente)  # estimativa para o coeficiente de fugacidade
        else:
            self.estphi = estphi
            
        self.estBeta = estBeta # estimativa para a fração entre líquido e vapor
        
        # Tolerância
        self.tolAlg = tolAlg   # Tolerância do algortimo
        self.toleq   = toleq   # Tolerância do equilíbrio
        self.maxiter = maxiter # Número máximo de iterações

        # Calcular Métodos gerais:
        # Cálculo do segundo coeficiente da equação do VIRIAL
        if self.model_vap == 'VIRIAL':
            self.Second_Virial_Coef()
        # Cálculo do coeficiente de fugacidade nas condições de saturação
        self.PhiSat()


        if self.Algoritmo == 'Flash':
            self.Flash()
        elif self.Algoritmo == 'PontoBolha':
            self.PontoBolha()

        elif self.Algoritmo == 'PontoOrvalho':
            self.PontoOrvalho()   

        elif self.Algoritmo == 'Coeficiente_Atividade':    
            self.coefAct = self.Coeficiente_Atividade(self.z,self.Temp)

        elif self.Algoritmo == 'Coeficiente_Fugacidade':
            self.coefFug = self.Coeficiente_Fugacidade(self.z,self.Pressao,self.Temp)

    def Second_Virial_Coef(self):
        '''
        ==============================================================================
                Definição:
        ==============================================================================    
        Módulo para calcular o segundo coeficiente da equação Viral.
    
        ==============================================================================
                Entradas:
        ==============================================================================        
        Algoritmo sem entradas.
        
        ==============================================================================
                Saídas:
        ==============================================================================
        Coeficiente B Virial em forma de *list*.
        '''
        parametros  = self.model_vap_param
        T           = self.Temp
        Componente  = self.Componente
        NC          = self.NC

        if parametros[0] == 'Hayden_OConnel':
            # T     = Temperatura / K
            # ek    = energia característica da interação i-j, K
            # sigma = tamanho molecular , A
            # mi    = momento dipolo
            # Eta   = Parâmetro de associação (i=j), parâmetro de solvatação (i != j)
            # w     = Fator acêntrico não polar
            
            Eta = parametros[1]

            # -----------------------------------------------------------------------------------------------------------------------------------------------------------
            # Parâmetros independentes da temperatura:
         
            ek     = [[0.0 for j in xrange(NC)] for i in xrange(NC)]
            sigma  = [[0.0 for j in xrange(NC)] for i in xrange(NC)]
            w      = [[0.0 for j in xrange(NC)] for i in xrange(NC)]
            ekl    = [[0.0 for j in xrange(NC)] for i in xrange(NC)]
            sigmal = [[0.0 for j in xrange(NC)] for i in xrange(NC)]
          
            # Parâmetros Puros
            for i in xrange(NC):
                for j in xrange(NC):
                    if i == j:
                        w[i][j]      = 0.006026*Componente[i].radius_giration + 0.02096*(Componente[i].radius_giration**2.0) - 0.001366*(Componente[i].radius_giration**3.0)
                        sigmal[i][j] = (2.4507 - w[i][j])*(1.0133*Componente[i].Tc/Componente[i].Pc)**(1.0/3.0)
                        ekl[i][j]    = Componente[i].Tc*(0.748 + 0.91*w[i][j] - 0.4*Eta[i][j]/(2.0+20.0*w[i][j]))
                       
                        if Componente[i].dipole_moment < 1.45:
                            Xi = 0.0
                        else:
                            den1 = 2.882 - 1.882*w[i][j]/(0.03 + w[i][j])
                            den2 = Componente[i].Tc*(sigmal[i][j]**6.0)*ekl[i][j]
                            Xi   = (1.7941*10**7.0)*(Componente[i].dipole_moment**4.0)/(den1*den2)
                          
                        c1 = (16.0+400.0*w[i][j])/(10.0+400.0*w[i][j])
                        c2 = 3.0/(10.0 + 400.0*w[i][j])
                          
                        ek[i][j]     = ekl[i][j]*(1.0-Xi*c1*(1.0-Xi*(1.0+c1)/2.0))                            
                        sigma[i][j]  = sigmal[i][j]*(1+Xi*c2)**(1.0/3.0) 
            
            # Parâmetros cruzados
            for i in xrange(NC):
                for j in xrange(NC):
                    if i != j:
                        w[i][j]      = (0.5)*(w[i][i]+w[j][j])
                        sigmal[i][j] = (sigma[i][i]*sigma[j][j])**0.5
                        ekl[i][j]    = 0.7*((ek[i][i]*ek[j][j])**0.5) + 0.6/(1.0/ek[i][i] + 1.0/ek[j][j])
                            
                        if Componente[i].dipole_moment >= 2.0 and Componente[j].dipole_moment == 0.0:
                            Xil = (Componente[i].dipole_moment**2.0) * ek[j][j]**(2.0/3.0) * (sigma[j][j]**4.0) / (ekl[i][j]*sigmal[i][j]**6.0)
                        elif Componente[j].dipole_moment >= 2.0 and Componente[i].dipole_moment == 0.0:
                            Xil = (Componente[j].dipole_moment**2.0) * ek[i][i]**(2.0/3.0) * (sigma[i][i]**4.0) / (ekl[i][j]*sigmal[i][j]**6.0)
                        else:
                            Xil = 0.0
                         
                        c1l = (16.0+400.0*w[i][j])/(10.0+400.0*w[i][j])
                        c2l = 3.0/(10.0 + 400.0*w[i][j])
                      
                        ek[i][j]     = ekl[i][j]*(1.0+Xil*c1l)                            
                        sigma[i][j]  = sigmal[i][j]*(1.0-Xil*c2l)**(1.0/3.0)
                      
            mi_ast = [[7243.8*self.Componente[i].dipole_moment*Componente[j].dipole_moment/(ek[i][j]*(sigma[i][j]**3.0)) for j in xrange(NC)] for i in xrange(NC)]
          
            b0 = [[1.26184*(sigma[i][j]**3.0) for j in xrange(NC)] for i in xrange(NC)]
         
            mi_astl = [[0.0 for j in xrange(NC)] for i in xrange(NC)] 
            for i in xrange(NC):
                for j in xrange(NC):
                    if mi_ast[i][j] < 0.04:
                        mi_astl[i][j] = mi_ast[i][j]
                    elif mi_ast[i][j] < 0.25 and mi_ast[i][j] >= 0.04:
                        mi_astl[i][j] = 0.0
                    elif mi_ast[i][j] >= 0.25:
                        mi_astl[i][j] = mi_ast[i][j] - 0.25
          
            A      = [[-0.3 - 0.05*mi_ast[i][j]        for j in xrange(NC)] for i in xrange(NC)]         
            deltah = [[1.99 + 0.20*(mi_ast[i][j]**2.0) for j in xrange(NC)] for i in xrange(NC)]
                 
            E = [[0.0 for j in xrange(NC)] for i in xrange(NC)] 
            for i in xrange(NC):
                for j in xrange(NC):
                    if Eta[i][j] < 4.5:
                        E[i][j] = exp(Eta[i][j]*(650.0/(ek[i][j]+300.0) - 4.27))
                    elif Eta[i][j] >= 4.5:
                        E[i][j] = exp(Eta[i][j]*(42800.0/(ek[i][j]+22400.0) - 4.27))
                  
            # -----------------------------------------------------------------------------------------------------------------------------------------------------------            
            # Parâmetros dependentes da temperatura:
            T_ast   = [[T/ek[i][j]                  for j in xrange(NC)] for i in xrange(NC)]
            T_astll = [[1/T_ast[i][j] - 1.6*w[i][j] for j in xrange(NC)] for i in xrange(NC)]
            
            #Cálculos dos BF's:
            BFnonpolar = [[ b0[i][j]*(0.94 - 1.47*T_astll[i][j] - 0.85*(T_astll[i][j]**2.0) + 1.015*(T_astll[i][j]**3.0))           for j in xrange(NC)] for i in xrange(NC)]
            BFpolar    = [[-b0[i][j]*mi_astl[i][j]*(0.74 - 3.0*T_astll[i][j] + 2.1*(T_astll[i][j]**2.0) + 2.1*(T_astll[i][j]**3.0)) for j in xrange(NC)] for i in xrange(NC)]
         
            # Cálculos para BD:
            Bmetastable_Bbound = [[b0[i][j]*A[i][j]*exp(deltah[i][j]/T_ast[i][j]) for j in xrange(NC)] for i in xrange(NC)]
            Bchemical          = [[b0[i][j]*E[i][j]*(1 - exp(1500.0*Eta[i][j]/T)) for j in xrange(NC)] for i in xrange(NC)]
            # -----------------------------------------------------------------------------------------------------------------------------------------------------------
         
            BD            = [[Bmetastable_Bbound[i][j] + Bchemical[i][j] for j in xrange(NC)] for i in xrange(NC)] # D bound or dimerizes molecules (Chemical forces)
            BF            = [[BFnonpolar[i][j]         + BFpolar[i][j]   for j in xrange(NC)] for i in xrange(NC)] # Free molecules
            self.Bvirial  = [[BF[i][j]                 + BD[i][j]        for j in xrange(NC)] for i in xrange(NC)]

    def Coeficiente_Atividade(self,x,T):
        '''
        ==============================================================================
                Definição:
        ==============================================================================    
        Módulo para calcular o coeficiente de atividade (gamma).
    
        ==============================================================================
                Entradas:
        ==============================================================================        
        *Composições dos componentes da fase líquida (x):
            Inserido como *list*.
        *Temperatura (T):
            Dada em Kelvin e inserido como *float*
        ==============================================================================
                Saídas:
        ==============================================================================
        Coeficiente de atividade em forma de *list*, onde o primeiro elemento é o coeficiente
        de atividade para o primeiro componente da lista ``Componentes``, inserida em VLE. E o
        segundo elemento é o coeficiente de atividade para o segundo componente da mesma lista.
        '''                        
        # Modelo UNIQUAC
        if self.model_liq == 'UNIQUAC':
            self.a           = self.model_liq_param[0]
            self.coordnumber = self.model_liq_param[1]
            #print 'COORD',self.coordnumber

            tau   = [[exp(-self.a[i][j]/T) for j in xrange(2)] for i in xrange(2)]
            #print 'a',self.a,'\n'
            #print 'tau',tau,'\n'
            phi   = [self.Componente[i].r  * x[i] / sum([self.Componente[j].r  * x[j] for j in xrange(self.NC)])        for i in xrange(self.NC)]
            teta  = [self.Componente[i].q  * x[i] / sum([self.Componente[j].q  * x[j] for j in xrange(self.NC)])        for i in xrange(self.NC)]
            tetal = [self.Componente[i].ql * x[i] / sum([self.Componente[j].ql * x[j] for j in xrange(self.NC)])        for i in xrange(self.NC)]
            l     = [(self.coordnumber/2.0)*(self.Componente[i].r - self.Componente[i].q) - (self.Componente[i].r - 1)  for i in xrange(self.NC)]
             
            #print 'phi',phi,'\n'
            #print 'teta',teta, '\n'
            #print 'tetal',tetal,'\n'
            #print  'l',l,'\n'

            A = [self.Componente[i].ql*sum([tetal[j]*tau[i][j]/sum([tetal[k]*tau[k][j] for k in xrange(self.NC)]) for j in xrange(self.NC)]) for i in xrange(self.NC)]
        
            #print 'A', A, '\n'
            Combinatorial = [log(phi[i]/x[i]) + (self.coordnumber/2.0)*self.Componente[i].q*log(teta[i]/phi[i]) + l[i] - (phi[i]/x[i]) * sum([x[j]*l[j] for j in xrange(self.NC)]) for i in xrange(self.NC)]
            Residual      = [-self.Componente[i].ql*log(sum([tetal[j]*tau[j][i] for j in xrange(self.NC)])) + self.Componente[i].ql - A[i]      for i in xrange(self.NC)] 
            #print 'Comb',Combinatorial,'\n'
            #print 'Res',Residual,'\n'

            Coeficiente_Atividade = [exp(Combinatorial[i]+Residual[i]) for i in xrange(self.NC)]
        
        # Modelo NRTL
        if self.model_liq == 'NRTL':
            R = 83.144621            
            self.g           = self.model_liq_param[0]            
            self.alpha       = self.model_liq_param[1]                      
            
            tau    =  [[self.g[i][j]/(R*T) for j in xrange(self.NC)] for i in xrange(self.NC)]
            G  = [[exp(-self.alpha[i][j]*tau[i][j]) for j in xrange(self.NC)] for i in xrange(self.NC)]

            Coeficiente_Atividade = []
            # Cálculo do coeficiente de atividade

            for i in xrange(self.NC):
    
                soma1 = 0
                for j in xrange(self.NC):
                    soma1 = tau[j][i]*G[j][i]*x[j] + soma1
    
                soma2 = 0
                for k in xrange(self.NC):
                    soma2 = G[k][i]*x[k] + soma2
        
                soma5 = 0          
                for j in xrange(self.NC):
                    soma3 = 0
                    soma4 = 0
                    for k in xrange(self.NC):
                        soma3 = G[k][j]*x[k] + soma3
                        soma4 = x[k]*tau[k][j]*G[k][j]+soma4
                    soma5 = ( (x[j]*G[i][j])/soma3 )*( tau[i][j]-(soma4/soma3) ) + soma5
                    
                Coeficiente_Atividade_aux = exp( soma1/soma2 + soma5 ) # |Eq.3| Equação do modelo NRTL para multicomponentes
                Coeficiente_Atividade.append(Coeficiente_Atividade_aux)
                
                
       # Modelo Van Laar
        if self.model_liq == 'Van_Laar':
            R = 83.144621                        
            self.A_VL       =   self.model_liq_param[0]
            self.B_VL       =   self.model_liq_param[1]
            
            self.p_VL = [self.A_VL,self.B_VL]          
            
            Coeficiente_Atividade_aux   =    [[exp( ( self.p_VL[i]/(R*T) )*( 1+(self.p_VL[i]/self.p_VL[j])*\
            (x[i]/x[j]) )**-2 ) for j in xrange(self.NC) if i!=j] for i in xrange(self.NC)]
            
            Coeficiente_Atividade       =    [Coeficiente_Atividade_aux[0][0],Coeficiente_Atividade_aux[1][0]]
        
        # Modelo de Wilson                      
        if self.model_liq == 'Wilson':
            self.A          =   self.model_liq_param[0]
            
            Coeficiente_Atividade = []

            for k in xrange(self.NC):
    
                soma1 = 0
                for j in xrange(self.NC):
                    soma1 = self.A[k][j]*x[j] + soma1
        
                soma2 = 0          
                for i in xrange(self.NC):
                    soma3 = 0
                    for j in xrange(self.NC):
                        soma3 = self.A[i][j]*x[j] + soma3
                    soma2 = (x[i]*self.A[i][k])/soma3 + soma2
        
                Coeficiente_Atividade_aux = exp( -log(soma1)+1-soma2 ) # |Eq.2| Equação do modelo NRTL para multicomponentes
                Coeficiente_Atividade.append(Coeficiente_Atividade_aux)
                
        return Coeficiente_Atividade

    def Coeficiente_Fugacidade(self,y,P,T):
        '''
        ==============================================================================
                Definição:
        ==============================================================================    
        Módulo para calcular o coeficiente de fugacidade (phi).
    
        ==============================================================================
                Entradas:
        ==============================================================================        
        *Composições dos componentes da fase de vapor (y):
            Inserido como *list*.
        *Temperatura (T):
            Dada em Kelvin e inserido como *float*
        *Pressão do sistema (P):
            Dada em bar e inserido como *float*.
        
        ==============================================================================
                Saídas:
        ==============================================================================
        Coeficiente de fugacidade em forma de *list*, onde o primeiro elemento é o coeficiente
        de fugacidade para o primeiro componente da lista ``Componentes``, inserida em VLE. E o
        segundo elemento é o coeficiente de fugacidade para o segundo componente da mesma lista.
        '''   
        if self.model_vap == 'VIRIAL':
            NC = size(y)
            R = 83.144621 # cm3.bar/(K.mol)
            B = self.Bvirial
            Bmixture = sum([sum([y[i]*y[j]*B[i][j] for j in xrange(NC)])      for i in xrange(NC)])
            A        = [ 2*sum([y[j]*B[i][j] for j in xrange(NC)]) - Bmixture for i in xrange(NC)]
            phi      = [exp(A[i]*P/(R*T)) for i in xrange(NC)]
            
            # Validação grosseira das condições de pressão:
            P_lim = (T/2.0)*sum([y[i]*self.Componente[i].Pc for i in xrange(NC)])/sum([y[i]*self.Componente[i].Tc for i in xrange(NC)])
            if P <= P_lim:
                self.ValVirial = 'A pressão do sistema é inferior à da validação'
            else:
                self.ValVirial = 'Warning: A pressão do sistema é SUPERIOR à da validação - Verificar uso da equacao do VIRIAL'

        return phi

    def Flash(self):
        '''
        ==============================================================================
                Definição:
        ==============================================================================    
        Módulo para realizar o calculo de flash
    
        ==============================================================================
                Entradas:
        ==============================================================================        
        Algoritmo sem entradas.
        
        ==============================================================================
                Saídas:
        ==============================================================================
        ?
        '''   

        self.PontoBolha()
        self.PontoOrvalho()
        
        if (self.Pressao <= self.Bolha.Pressao) and (self.Pressao >= self.Orvalho.Pressao):
            interp = (self.Pressao - self.Orvalho.Pressao)/(self.Bolha.Pressao - self.Orvalho.Pressao)
            coefAct = [(self.Bolha.liquido.coefAct[i] - self.Bolha.liquido.coefAct[i])*interp + self.Orvalho.liquido.coefAct[i]  for i in xrange(self.NC) ]
            coefFug = [(self.Bolha.vapor.coefFug[i]   - self.Bolha.vapor.coefFug[i])*interp   + self.Orvalho.vapor.coefFug[i]    for i in xrange(self.NC) ]
            Beta    = [(self.Bolha.Pressao - self.Pressao)/(self.Bolha.Pressao - self.Orvalho.Pressao)]

            cont = 0; deltabeta = 1e10; Valeq = 1e10; tBeta = False
            while ((deltabeta>self.tolAlg or mean(abs(deltaeq))>self.toleq) and (cont<(self.maxiter))) or (tBeta == True):
                print cont
		# Início das listas K, x e y
		K = range(self.NC)
		x = range(self.NC)
		y = range(self.NC)
		
		# Cálculo de K
		K = [(coefAct[i]*self.Componente[i].Psat) / (coefFug[i]*self.Pressao) for i in xrange(self.NC)]
		
		# Cálculo de F e sua derivada em relação à fração de vapor
		F     = sum([self.z[i]*(K[i]-1)    / (1+Beta[cont]*(K[i]-1))    for i in xrange(self.NC)])
		diffF = sum([self.z[i]*(K[i]-1)**2 / (1+Beta[cont]*(K[i]-1))**2 for i in xrange(self.NC)])
		
		# Aplicação do método de Newton
		Beta.append(Beta[cont]+0.05*F/diffF)
		
		if Beta[cont+1] > 1.0 or Beta[cont+1] < 0.0:
		   Beta[cont+1] = 0.5
		   tBeta        = True

		# Cálculo das composições e subsequente normalização:
		x    = [self.z[i] / (1 + Beta[cont+1]*(K[i]-1)) for i in xrange(self.NC)]
		norm = sum(x)
		x    = [x[i]/norm for i in xrange(self.NC)]
		
		y    = [K[i]*x[i] for i in xrange(self.NC)]
		norm = sum(y)
		y    = [y[i]/norm for i in xrange(self.NC)]
	
		# Cálculo dos coeficientes de atividade:
		coefAct    = self.Coeficiente_Atividade(x,self.Temp)
		coefFug    = self.Coeficiente_Fugacidade(y,self.Pressao,self.Temp)
		
		coefFugsat = [self.Coeficiente_Fugacidade([y[i]],self.Componente[i].Psat,self.Temp) for i in xrange(2)]
	
		# Validação do equilíbrio líquido-vapor
		# Líquido
		Eqliq = [x[i]*coefAct[i]*self.Componente[i].Psat for i in xrange(self.NC)]
		# Vapor
		Eqvap = [y[i]*coefFug[i]*self.Pressao for i in xrange(self.NC)]    
		# Diferença entre Eqliq e Eqvap (Eles deveriam ser iguais) para cada componente       
		deltaeq   = [Eqliq[i] - Eqvap[i] for i in xrange(self.NC)]

		# Diferença entre valores de beta
		deltabeta = Beta[cont] - Beta[cont-1]
	
		cont+=1
	
	    if deltabeta>self.tolAlg:
		 self.wBeta = "Convergência em Beta NÃO alcançada"
	    else:
		self.wBeta = "Convergência em Beta alcançada"
		
	    if mean(abs(deltaeq))>self.toleq:
		self.wVal = "Convergência no equilíbrio NÃO alcançada"
	    else:
		self.wVal = "Convergência no equilíbrio alcançada"

	    if cont>=(self.maxiter):
		self.witer = "Número máximo de iterações atingido"
	    else:
		self.witer = "Número máximo de iterações NÃO atingido"
	    
	    liquido = Fase(composicao=x,coeffug=None   ,coefAct = coefAct)
	    vapor   = Fase(composicao=y,coeffug=coefFug,coefAct = None)
	    Beta    = Beta
            self.Flash = Condicao(self.Pressao,self.Temp,liquido,vapor,Beta)
	else:
	    self.wFlash = 'Não é possível executar o Flash'
            self.Flash = None

    def PhiSat(self):
        '''
        ==============================================================================
                Definição:
        ==============================================================================    
        Módulo para calcular o coeficiente de fugacidade nas condições de saturação (phisat).
    
        ==============================================================================
                Entradas:
        ==============================================================================        
        Algoritmo sem entradas.
        
        ==============================================================================
                Saídas:
        ==============================================================================
        Coeficiente de fugacidade em forma de *list*, onde o primeiro elemento é o coeficiente
        de fugacidade para o primeiro componente da lista ``Componentes``, inserida em VLE. E o
        segundo elemento é o coeficiente de fugacidade para o segundo componente da mesma lista.
        ''' 
        # Cálculo de phisat (coeficiente de fugacidade nas condições de saturação)
        # Utilizado o mesmo modelo de phi, quando x->1.
        phisat = []
        for i in xrange(self.NC):
            comp    = zeros((1,self.NC)) + 0.00000001
            comp    = comp.tolist()[0]
            comp[i] = 0.99999            
            fator  = (1-comp[i])/((self.NC-1)*0.00000001)
            for j in xrange(self.NC):
                if i != j:
                    comp[j] = fator*0.00000001
            phisat.append(self.Coeficiente_Fugacidade(comp,self.Componente[i].Psat,self.Temp)[i])
        self.phisat = phisat
    def PontoBolha(self):
        '''
        ==============================================================================
                Definição:
        ==============================================================================    
        Módulo para calcular o ponto de bolha, quando a temperatura é conhecida.
    
        ==============================================================================
                Entradas:
        ==============================================================================        
        Algoritmo sem entradas.
        
        ==============================================================================
                Saídas:
        ==============================================================================
        Pressão e composição da fase de vapor, onde a pressão retorna na forma de *float*
        e a composição na forma de *list*.
        ''' 
        
        

        P        = []
        coefAct  = self.Coeficiente_Atividade(self.z,self.Temp)
        liquido  = Fase(composicao=self.z,coeffug=None,coefAct = coefAct)
        coeffug  = self.estphi
        
        cont   = 0; deltaP = 10000
        while (deltaP > self.tolAlg) and (cont<self.maxiter+1):
            P.append(sum([liquido.comp[i]*liquido.coefAct[i]*self.Componente[i].Psat*self.phisat[i]/coeffug[i]        for i in xrange(self.NC)]))
#            print liquido.comp[0], liquido.coefAct[0], self.Componente[1].Psat
            y       = [liquido.comp[i]*liquido.coefAct[i]*self.Componente[i].Psat*self.phisat[i]/(coeffug[i]*P[cont]) for i in xrange(self.NC)]
            coeffug = self.Coeficiente_Fugacidade(y,P[cont],self.Temp)
            deltaP = abs(P[cont] - P[cont-1])
            cont+=1
        vapor =  Fase(composicao=y,coeffug=coeffug,coefAct = None)
        self.Bolha = Condicao(P[cont-1],self.Temp,liquido,vapor,0.0)

    def PontoOrvalho(self):
        '''
        ==============================================================================
                Definição:
        ==============================================================================    
        Módulo para calcular o ponto de orvalho, quando a temperatura é conhecida.
    
        ==============================================================================
                Entradas:
        ==============================================================================        
        Algoritmo sem entradas.
        
        ==============================================================================
                Saídas:
        ==============================================================================
        Pressão e composição da fase líquida, onde a pressão retorna na forma de *float*
        e a composição na forma de *list*.
        ''' 
        P        = [self.Pressao]
        coeffug  = [1.0,1.0]
        coefAct  = [1.0,1.0]        
        cont   = 1; deltaP = 10000
        while (deltaP > self.tolAlg) and (cont<self.maxiter+1):
            P.append(1/sum([self.z[i]*coeffug[i]/(coefAct[i]*self.Componente[i].Psat*self.phisat[i])    for i in xrange(self.NC)]))
            x       = [self.z[i]*coeffug[i]*P[cont]/(coefAct[i]*self.Componente[i].Psat*self.phisat[i]) for i in xrange(self.NC)]
            coeffug = self.Coeficiente_Fugacidade(self.z,P[cont],self.Temp)
            coefAct = self.Coeficiente_Atividade(x,self.Temp)
            
            deltaP = abs(P[cont] - P[cont-1])
            cont+=1

        liquido = Fase(composicao=x,coeffug=None,coefAct = coefAct)
        vapor   = Fase(composicao=self.z,coeffug=coeffug,coefAct = None)
        self.Orvalho = Condicao(P[cont-1],self.Temp,liquido,vapor,1.0)
