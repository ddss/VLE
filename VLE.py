# -*- coding: utf-8 -*-
"""
@author: Daniel
@Contribuição:Caique
@ Data: 27/10/2013

Rotinas para Cálculo do Equilíbrio de fases em Misturas através da abordagem gamma-phi

Modelos termodinâmicos:
    - UNIQUAC
    - VIRIAL
    - NRTL
    - WILSON
    - Van Laar

Métodos:
    - Second_Virial_Coef: Cálculo do segundo coeficiente do Virial
    - Coeficiente_Atividade: Cálculo do Coeficiente de Atividade
    - Coeficiente_Fugacidade: Cálculo do Coeficiente de Fugacidade
    - Flash: Cáculo de um flash (Incompleto)
    - PhiSat: Cálculo do coeficiente de fugacidade nas condições de saturação
    - PontoBolha_P: Cálculo do ponto de bolha (T conhecido) (Dado: x, T/K -> Cálcula: y, P/bar)
    - PontoOrvalho_P: Cálculo do ponto de orvalho (T conhecido) (Dado: y, T/K -> Cálcula: x, P/bar)
    - PontoBolha_T: Cálculo do ponto de bolha (P conhecido) (Dado: x, P/bar -> Cálcula: y, T/K)
    - PontoOrvalho_T: Cálculo do ponto de orvalho (P conhecido) (Dado: y, P/bar -> Cálcula: x, T/K)

Referências:
[1] PRAUSNITZ, J. M. et al. Computer Calculations for multicomponent vapor-liquid and liquid-liquid equilibria. [s.l.] Prendice-Hall, 1980. p. 353
[2] HAYDEN, J. G.; O’CONNELL, J. P. A Generalized Method for Predicting Second Virial Coefficients. Industrial & Engineering Chemistry Process Design
    and Development, v. 14, n. 3, p. 209–216, jul. 1975.
[3] SMITH, J. M.; NESS, H. C. VAN; ABBOTT, M. M. Introduction to Chemical Engineering Thermodinamics. 7th. ed. [s.l.] Mc-Graw Hills, [s.d.]. 
"""
# Sistema
import sys
reload(sys)
sys.setdefaultencoding("utf-8") # Forçar o sistema utilizar o coding utf-8

from threading import Thread
from warnings import warn
from numpy import log, exp, size, abs, zeros, linspace

class Condicao:
    
    def __init__(self,P,T,composicao,coeffug,coefAct,**kwargs):
        '''
        Rotina que caracteriza as condições do equilíbrio líquido-vapor.
        
        =====================
        Entradas obrigatórias
        =====================
        
        * P (float): Pressao em bar;
        * T (float): Temperatura em Kelvin;
        * composicao (list): Composição da mistura na fase;
        * coeffug (list): Os valores dos coeficientes de fugacidade;
        * coefAct (list): Os valores dos coeficientes de atividade.
        
        ==============================
        Keywords (Entradas opcionais):
        ==============================
        
        * ``Numero_componentes``(int): Número de componentes utilizado
        * ``Massa_molar``       (list): lista com as massas molares dos componentes;
        * ``beta``              (float): relação entre vapor e líquido.

        =========
        Atributos
        =========
        
        Os atributos desta classe possuem as mesmas propriedades das entradas da classe.
        
        * ``Pressao`` (float): Pressao em bar;
        * ``Temp`` (float): Temperatura em Kelvin;
        * ``comp`` (list): Composição da mistura na fase;
        * ``coeffug`` (list): Os valores dos coeficientes de fugacidade;
        * ``coefAct`` (list): Os valores dos coeficientes de atividade;
        '''
        # ----------------------------------------------------
        # VALIDAÇÃO
        # ----------------------------------------------------
        # Keywords disponíveis
        __keywordsEntrada = ['massa_molar','numero_componentes','beta']

        # Validação se houve keywords digitadas incorretamente:
        keyincorreta  = [key for key in kwargs.keys() if not key in __keywordsEntrada]

        if len(keyincorreta) != 0:
            raise NameError(u'keyword(s) incorretas: '+', '.join(keyincorreta)+'.'+u' Keywords disponíveis: '+', '.join(self.__keywordsEntrada)+'.')

        # ----------------------------------------------------
        # EXECUÇÃO
        # ----------------------------------------------------

        self.Pressao = P
        self.Temp    = T
        self.coeffug = coeffug
        self.coefAct = coefAct        
        self.comp_molar = composicao        
        self.beta     = kwargs.get('beta')

        # Caracterização dos kwargs
        mm_comp = kwargs.get(__keywordsEntrada[0])
        NC = kwargs.get(__keywordsEntrada[1])

        if NC and mm_comp is not None:
            # Cálculo da composição mássica
            mm_medio = sum([self.comp_molar[0]*mm_comp[0],self.comp_molar[1]*mm_comp[1]]) # em g/mol
            self.comp_massica = [mm_comp[i]*self.comp_molar[i]/mm_medio for i in xrange(NC)]




class VLE(Thread):        

    def __init__(self,Algoritmo,Componentes,model_liq, model_vap,z=None,Temp=None,Pressao=None,estgama=None,estphi=None, estBeta = 0.5, tolAlg=1e-10, toleq=1e-4, maxiter=100, z_coordenacao = 10.0 ):    
        '''
        ************************
        Vapor-Liquid Equilibrium
        ************************        
        
        Rotina para realizar os cálculos do equilíbrio líquido-vapor, desde o cálculo de
        coeficiente de atividade aos cálculos de ponto de bolha e de orvalho.
        
        ========
        Entradas
        ========  
        
        *  Algoritmo (str): O nome do cálculo que se deseja fazer. Os algoritmos de cálculo disponíveis são: 'Coeficiente_Atividade', 'Coeficiente_Fugacidade', 'PontoBolha', 'PontoOrvalho' e 'Flash'.
        
            * Exemplos da entrada Algoritmo: ::
                
                Algoritmo = 'PontoBolha_P'
                Algoritmo = 'Coeficiente_Atividade'
                Algoritmo = 'Flash'   
                
        * Componentes (list): É uma lista de objetos ``Componente_Caracterizar``, vide documentação da dessa classe;
        
            * A entrada pode ser dada da seguinte forma: ::
            
                >>> Comp1 = Componente_Caracterizar('Metano',ConfigPsat=('Prausnitz4th',1),T=100.0)
                >>> Comp2 = Componente_Caracterizar('Etano',ConfigPsat=('Prausnitz4th',1),T=289.9)
                
                >>> Componentes = [Comp1,Comp2]
                
        * model_liq (list): É um objeto do grupo de classes de modelos da rotina ``Conexao``. Vide documentação dos modelos presentes nesta rotina. 
          Os modelos para a fase líquida disponíveis na rotina são: UNIQUAC[1] , NRTL[2], WILSON[3] e Van_Laar[4];

            * Exemplo da entrada model_liq: ::
            
                >>> from Conexao import Van_Laar
                
                >>> Comp1 = Componente_Caracterizar('Metano',ConfigPsat=('Prausnitz4th',1),T=100.0)
                >>> Comp2 = Componente_Caracterizar('Etano',ConfigPsat=('Prausnitz4th',1),T=289.9)
                
                >>> Componentes = [Comp1,Comp2]
                
                >>> model_liq = Van_Laar(Componentes)
                
        * model_vap (list): É um objeto do grupo de classes de modelos da rotina ``Conexao``. Vide documentação dos modelos presentes nesta rotina.
          Os modelos para a fase de vapor disponíveis na rotina são: VIRIAL[5];
            
            * Exemplo da entrada model_vap: ::
            
                >>> from Conexao import VIRIAL
                
                >>> Comp1 = Componente_Caracterizar('Metano',ConfigPsat=('Prausnitz4th',1),T=100.0)
                >>> Comp2 = Componente_Caracterizar('Etano',ConfigPsat=('Prausnitz4th',1),T=289.9)
                
                >>> Componentes = [Comp1,Comp2]
                
                >>> model_vap = VIRIAL(Componentes,'Tsonopoulos')        
        
        * z (list): Composição global da mistura. Onde o primeiro elemento é a composição do componente(1) e o segundo é a composição do componente(2);
        * Temp (float): Temperatura do sistema dada em kelvin;
        * estgama (float): Estimativa para o parametro gamma;
        * estphi (float): Estimativa para o parametro phi;
        * estBeta (float): Estimativa para o parametro Beta;    
        * tolAlg (float): Tolerância do algoritmo, a tolerância desejada para a operação dos métodos;
        * toleq (float): Tolerância do equilíbrio, a tolerância desejada para o equilíbrio;
        * maxiter (int): Número máximo de iterações desejadas para a operação dos métodos;
        * z_coordenacao (float): Número de coordenação do componente.
        
        
        ===============
        Valores padrões
        ===============

        As seguintes entradas possuem valores padrões.
        
            * estphi: Pode ser None;
            * estgama: Pode ser None;
            * estBeta = 0.5;
            * tolAlg = 1e-10;
            * toleq = 1e-4;
            * maxiter = 100;
            * z_coordenacao = 10.0.
        
        =========
        Atributos
        =========
        
        Os atributos dependem do algoritmo de cálculo escolhido. Além dos Métodos que podem ser acessados como atributos, vide documentação dos métodos,
        os atributos para cada algoritmo de cálculo são:
        
        * ``Coeficiente_Atividade``:
            * ``coefAct`` (list): O atributo gerado representa os coeficientes de atividade dos componentes baseado nas condições de entrada da classe VLE.
        * ``Coeficiente_Fugacidade``:        
            * ``coefFug`` (list): O atributo gerado representa os coeficientes de atividade dos componentes baseado nas condições de entrada da classe VLE;
            * ``Bvirial`` (list): Uma lista de listas contendo os valores para os componentes puros e cruzados do segundo coeficiente Virial em unidade de volume.
        * ``PontoBolha_P`` & ``PontoBolha_T`` :
            * ``Bolha`` : Objeto da classe ``Condicao``, vide documentação da classe.
            * ``liquido`` : Objeto da classe ``Condicao``, vide documentação da classe.
            * ``Bvirial`` (list): Uma lista de listas contendo os valores para os componentes puros e cruzados do segundo coeficiente Virial em unidade de volume.
            * ``phisat`` (list): Os coeficientes de fugacidade nas condições de saturação.
        * ``PontoOrvalho_P`` & ``PontoOrvalho_T`` :
            * ``Orvalho`` : Objeto da classe ``Condicao``, vide documentação da classe.
            * ``vapor`` : Objeto da classe ``Condicao``, vide documentação da classe.
            * ``Bvirial`` (list): Uma lista de listas contendo os valores para os componentes puros e cruzados do segundo coeficiente Virial em unidade de volume.
            * ``phisat`` (list): Os coeficientes de fugacidade nas condições de saturação.
        * ``Flash``:
            * ``Bolha`` : Objeto da classe ``Condicao``, vide documentação da classe.
            * ``Orvalho`` : Objeto da classe ``Condicao``, vide documentação da classe.
            * ``Beta`` (float): Valor da coeficiente Beta calculado no Flash
            
        =======
        Métodos
        =======
        
        * ``Second_Virial_Coef``:
            * Método que realiza o cálculo do segundo coeficiente Virial, de acordo com as regras de Hayden O'Connel[6] e Tsonopoulos[7], vide documentação do método.
        * ``Coeficiente_Atividade``:
            * Método para cálcular do coeficiente de atividade, vide documentação do método.
        * ``Coeficiente_Fugacidade``:
            * Método para cálcular do coeficiente de fugacidade, vide documentação do método.
        * ``Flash``:
            * Método para realizar o cálculo de flash, vide documentação do método.
        * ``Phisat``:
            Método que realiza o cálculo de phisat (coeficiente de fugacidade nas condições de saturação).
        * ``PontoBolha_P``:
            Método para realizar o cálculo do ponto de bolha dado x e T.
        * ``PontoBolha_T``:
            Método para realizar o cálculo do ponto de bolha dado x e P.
        * ``PontoOrvalho_P``:
            Método para realizar o cálculo do ponto de orvalho dado y e T.
        * ``PontoOrvalho_T``:
            Método para realizar o cálculo do ponto de orvalho dado y e P.
            
        =========
        Exemplo 1
        =========
        
        A priori, é necessário utilizar a rotina ``Conexao`` para caracterizar os componentes e os modelos.
        Neste exemplo, será calculado o ponto de bolha utilizando o método correspondente. A mistura utilizada será Acetona-Metanol, os modelos para as fases líquida e de vapor são, respectivamente
        UNIQUAC e VIRIAL com a regra de Hayden O'Connel, a temperatura será 330.0 K e a pressão 1.013 bar. ::
        
            >>> from Conexao import Componente_Caracterizar, UNIQUAC, VIRIAL
            
            >>> Comp1 = Componente_Caracterizar('Acetona',ConfigPsat=('Prausnitz4th',1),T=330.0)
            >>> Comp2 = Componente_Caracterizar('Metanol',ConfigPsat=('Prausnitz4th',1),T=330.0)
            
            >>> Componentes = [Comp1,Comp2]
            
            >>> model_vap = VIRIAL(Componentes)
            >>> model_liq = UNIQUAC(Componentes,330.0,1)
            
            >>> CalculoBolha = VLE('PontoBolha_P',Componentes,model_liq,model_vap,z=[0.95,0.05],Temp=330.0,Pressao=1.013,estgama=None,estphi=None,estBeta = 0.5,tolAlg=1e-5,toleq=1e-4,maxiter=500,z_coordenacao=10.0)
            
        O acesso do valor da pressão do ponto de bolha é da seguinte forma:
        
            >>> CalculoBolha.Bolha.Pressao
        
        =========
        Exemplo 2
        =========
        
        Os passos deste exemplo são idênticos aos do Exemplo 1. Neste exemplo, será calculado o ponto de orvalho utilizando o método correspondente.
        A mistura, as condições da mesma e os modelos são os mesmos do Exemplo 1. ::
        
            >>> from Conexao import Componente_Caracterizar, UNIQUAC, VIRIAL
            
            >>> Comp1 = Componente_Caracterizar('Acetona',ConfigPsat=('Prausnitz4th',1),T=330.0)
            >>> Comp2 = Componente_Caracterizar('Metanol',ConfigPsat=('Prausnitz4th',1),T=330.0)
            
            >>> Componentes = [Comp1,Comp2]
            
            >>> model_vap = VIRIAL(Componentes)
            >>> model_liq = UNIQUAC(Componentes,330.0,1)
            
            >>> CalculoOrvalho = VLE('PontoOrvalho_P',Componentes,model_liq,model_vap,z=[0.95,0.05],Temp=330.0,Pressao=1.013,estgama=None,estphi=None,estBeta = 0.5,tolAlg=1e-5,toleq=1e-4,maxiter=500,z_coordenacao=10.0)
            
        O acesso do valor da pressão do ponto de orvalho é da seguinte forma:
        
            CalculoOrvalho.Orvalho.Pressao
            
        =========
        Exemplo 3
        =========
        
        Os passos deste exemplo são idênticos aos do Exemplo 1. Neste exemplo, serão calculados os coeficientes de atividade utilizando o método correspondente.
        A mistura, as condições da mesma e os modelos são os mesmos do Exemplo 1. ::
        
            >>> from Conexao import Componente_Caracterizar, UNIQUAC, VIRIAL
            
            >>> Comp1 = Componente_Caracterizar('Acetona',ConfigPsat=('Prausnitz4th',1),T=330.0)
            >>> Comp2 = Componente_Caracterizar('Metanol',ConfigPsat=('Prausnitz4th',1),T=330.0)
            
            >>> Componentes = [Comp1,Comp2]
            
            >>> model_vap = VIRIAL(Componentes)
            >>> model_liq = UNIQUAC(Componentes,330.0,1)
            
            >>> Calculo_Coef_Atividade = VLE('Coeficiente_Atividade',Componentes,model_liq,model_vap,z=[0.95,0.05],Temp=330.0,Pressao=1.013,estgama=None,estphi=None,estBeta = 0.5,tolAlg=1e-5,toleq=1e-4,maxiter=500,z_coordenacao=10.0)
            
        O acesso dos valores dos coeficientes de atividade é da seguinte forma:
        
            >>> Calculo_Coef_Atividade.coefAct
            
        =========
        Exemplo 4
        =========
        
        Os passos deste exemplo são idênticos aos do Exemplo 1. Neste exemplo, serão calculados os coeficientes de fugacidade utilizando o método correspondente.
        A mistura, as condições da mesma e os modelos são os mesmos do Exemplo 1. ::
        
            >>> from Conexao import Componente_Caracterizar, UNIQUAC, VIRIAL
            
            >>> Comp1 = Componente_Caracterizar('Acetona',ConfigPsat=('Prausnitz4th',1),T=330.0)
            >>> Comp2 = Componente_Caracterizar('Metanol',ConfigPsat=('Prausnitz4th',1),T=330.0)
                 
            >>> Componentes = [Comp1,Comp2]
           
            >>> model_vap = VIRIAL(Componentes)
            >>> model_liq = UNIQUAC(Componentes,330.0,1)
            
            >>> Calculo_Coef_Fugacidade = VLE('Coeficiente_Fugacidade',Componentes,model_liq,model_vap,z=[0.95,0.05],Temp=330.0,Pressao=1.013,estgama=None,estphi=None,estBeta = 0.5,tolAlg=1e-5,toleq=1e-4,maxiter=500,z_coordenacao=10.0)
            
        O acesso dos valores dos coeficientes de fugacidade é da seguinte forma:
        
            >>> Calculo_Coef_Fugacidade.coefFug
        
        ===========
        Referências
        ===========
        
        [1] ABRAMS, D. S.; PRAUSNITZ, J. M. Statistical thermodynamics of liquid
        mixtures: A new expression for the excess Gibbs energy of partly or completely
        miscible systems. AIChE Journal, v. 21, n. 1, p. 116–128, jan. 1975. ISSN
        0001-1541. Disponível em: <http://doi.wiley.com/10.1002/aic.690210115>
        
        [2] RENON, H.; PRAUSNITZ, J. M. Local compositions in thermodynamic excess
        functions for liquid mixtures. AIChE Journal, v. 14, n. 1, p. 135–144, jan. 1968.
        ISSN 0001-1541. Disponível em: <http://doi.wiley.com/10.1002/aic.690140124>.
        
        [3] WILSON, G. M. Vapor-Liquid Equilibrium. XI. A New Expression for the Excess
        Free Energy of Mixing. Journal of the American Chemical Society, v. 86, n. 2, p.
        127–130, jan. 1964. ISSN 0002-7863. Disponível em: <http://pubs.acs.org/doi/abs-
        /10.1021/ja01056a002>
        
        [4] VAN LAAR, J. J. The Vapor pressure of binary mixtures. Z. Phys.
        Chem. 1910, 72, 723−751.        
        
        [5] MASON, E. A.; SPURLING, T. H. The Virial Equation of State; The
        International Encyclopedia of Physical Chemistry and Chemical
        Physics, Topic 10: The Fluid State, Vol. 2; Pergamon Press: New
        York, 1969; p 297
        
        [6] HAYDEN, J. G.; O’CONNELL, J. P. A Generalized Method for Predicting
        Second Virial Coefficients. Industrial & Engineering Chemistry Process Design and
        Development, v. 14, n. 3, p. 209–216, jul. 1975. ISSN 0196-4305.
        
        [7] TSONOPOULOS, C.; HEIDMAN, J.L. From the Virial to the cubic equation of state. 
        Fluid Phase Equilib. 57 (1990) 261–276.

        '''
        Thread.__init__(self)
        # Definindo o __init_
        
        self.Algoritmo = Algoritmo  # Algoritmo a ser utilizado: Flash, PontoBolha, PontoOrvalho, Coeficiente_Atividade, Coeficiente_Fugacidade
        
        self.NC = size(Componentes) # Número de componentes da mistura

        self.Componente = Componentes

        # Definição dos modelos termodinâmicos:
        self.model_liq       = model_liq  # Modelo para a fase líquida
        self.model_vap       = model_vap  # Modelo para a fase vapor
        
        # Condições do VLE
        self.z           = z       # Composição global
        self.Pressao     = Pressao # Pressão em bar
        self.Temp        = Temp    # Temperatura em K
        
        # Estimativas iniciais

        if estgama is None:
            self.estgama = [1.0]*len(self.Componente) # estimativa para o coeficiente de atividade
        else:
            self.estgama = estgama
            
        if estphi is None:
            self.estphi  = [1.0]*len(self.Componente)  # estimativa para o coeficiente de fugacidade
        else:
            self.estphi = estphi

        if self.model_liq.nome_modelo == 'UNIQUAC':                
            self.coordnumber     = z_coordenacao # Número de coordenação do componente               
            
        self.estBeta = estBeta # estimativa para a fração entre líquido e vapor
        self.toleq   = toleq   # Tolerância do equilíbrio
        self.tolAlg = tolAlg   # Tolerância do algortimo            
        self.maxiter = maxiter # Número máximo de iterações
            
    def Second_Virial_Coef(self):
        '''
        Módulo para calcular o segundo coeficiente da equação Viral de acordo com as regras disponíveis.
        Estas são: Hayden O'Connel[1] e Tsonopoulos[2].
    
        ======
        Saídas
        ======
        
        * ``B.virial``: Atributo em forma de lista de listas contendo os valores do segundo coeficiente virial dos componentes puros e cruzados.
        
        ===========
        Referências
        ===========
        
        [1] HAYDEN, J. G.; O’CONNELL, J. P. A Generalized Method for Predicting
        Second Virial Coefficients. Industrial & Engineering Chemistry Process Design and
        Development, v. 14, n. 3, p. 209–216, jul. 1975. ISSN 0196-4305.
        
        [2] TSONOPOULOS, C.; HEIDMAN, J.L. From the Virial to the cubic equation of state. 
        Fluid Phase Equilib. 57 (1990) 261–276.
        
        '''
            
        if self.model_vap.regra_mistura == 'Hayden_o_Connel':
            
            R = 83.144621 # em cm3.bar/ K.mol
                    
            # T     = Temperatura / K
            # ek    = energia característica da interação i-j, K
            # sigma = tamanho molecular , A
            # mi    = momento dipolo
            # Eta   = Parâmetro de associação (i=j), parâmetro de solvatação (i != j)
            # w     = Fator acêntrico não polar
            
            Eta = self.model_vap.coef_solv

            # -----------------------------------------------------------------------------------------------------------------------------------------------------------
            # Parâmetros independentes da temperatura:
         
            ek     = [[0.0 for j in xrange(self.NC)] for i in xrange(self.NC)]
            sigma  = [[0.0 for j in xrange(self.NC)] for i in xrange(self.NC)]
            w      = [[0.0 for j in xrange(self.NC)] for i in xrange(self.NC)]
            ekl    = [[0.0 for j in xrange(self.NC)] for i in xrange(self.NC)]
            sigmal = [[0.0 for j in xrange(self.NC)] for i in xrange(self.NC)]
          
            # Parâmetros Puros
            for i in xrange(self.NC):
                for j in xrange(self.NC):
                    if i == j:
                        w[i][j]      = 0.006026*self.Componente[i].radius_giration + 0.02096*(self.Componente[i].radius_giration**2.0) - 0.001366*(self.Componente[i].radius_giration**3.0)
                        sigmal[i][j] = (2.4507 - w[i][j])*(1.0133*self.Componente[i].Tc/self.Componente[i].Pc)**(1.0/3.0)
                        ekl[i][j]    = self.Componente[i].Tc*(0.748 + 0.91*w[i][j] - 0.4*Eta[i][j]/(2.0+20.0*w[i][j]))
                       
                        if self.Componente[i].dipole_moment < 1.45:
                            Xi = 0.0
                        else:
                            den1 = 2.882 - 1.882*w[i][j]/(0.03 + w[i][j])
                            den2 = self.Componente[i].Tc*(sigmal[i][j]**6.0)*ekl[i][j]
                            Xi   = (1.7941*10**7.0)*(self.Componente[i].dipole_moment**4.0)/(den1*den2)
                          
                        c1 = (16.0+400.0*w[i][j])/(10.0+400.0*w[i][j])
                        c2 = 3.0/(10.0 + 400.0*w[i][j])
                          
                        ek[i][j]     = ekl[i][j]*(1.0-Xi*c1*(1.0-Xi*(1.0+c1)/2.0))                            
                        sigma[i][j]  = sigmal[i][j]*(1+Xi*c2)**(1.0/3.0) 
            
            # Parâmetros cruzados
            for i in xrange(self.NC):
                for j in xrange(self.NC):
                    if i != j:
                        w[i][j]      = (0.5)*(w[i][i]+w[j][j])
                        sigmal[i][j] = (sigma[i][i]*sigma[j][j])**0.5
                        ekl[i][j]    = 0.7*((ek[i][i]*ek[j][j])**0.5) + 0.6/(1.0/ek[i][i] + 1.0/ek[j][j])
                            
                        if self.Componente[i].dipole_moment >= 2.0 and self.Componente[j].dipole_moment == 0.0:
                            Xil = (self.Componente[i].dipole_moment**2.0) * ek[j][j]**(2.0/3.0) * (sigma[j][j]**4.0) / (ekl[i][j]*sigmal[i][j]**6.0)
                        elif self.Componente[j].dipole_moment >= 2.0 and self.Componente[i].dipole_moment == 0.0:
                            Xil = (self.Componente[j].dipole_moment**2.0) * ek[i][i]**(2.0/3.0) * (sigma[i][i]**4.0) / (ekl[i][j]*sigmal[i][j]**6.0)
                        else:
                            Xil = 0.0
                         
                        c1l = (16.0+400.0*w[i][j])/(10.0+400.0*w[i][j])
                        c2l = 3.0/(10.0 + 400.0*w[i][j])
                      
                        ek[i][j]     = ekl[i][j]*(1.0+Xil*c1l)                            
                        sigma[i][j]  = sigmal[i][j]*(1.0-Xil*c2l)**(1.0/3.0)
                      
            mi_ast = [[7243.8*self.Componente[i].dipole_moment*self.Componente[j].dipole_moment/(ek[i][j]*(sigma[i][j]**3.0)) for j in xrange(self.NC)] for i in xrange(self.NC)]
          
            b0 = [[1.26184*(sigma[i][j]**3.0) for j in xrange(self.NC)] for i in xrange(self.NC)]
         
            mi_astl = [[0.0 for j in xrange(self.NC)] for i in xrange(self.NC)] 
            for i in xrange(self.NC):
                for j in xrange(self.NC):
                    if mi_ast[i][j] < 0.04:
                        mi_astl[i][j] = mi_ast[i][j]
                    elif mi_ast[i][j] < 0.25 and mi_ast[i][j] >= 0.04:
                        mi_astl[i][j] = 0.0
                    elif mi_ast[i][j] >= 0.25:
                        mi_astl[i][j] = mi_ast[i][j] - 0.25
          
            A      = [[-0.3 - 0.05*mi_ast[i][j]        for j in xrange(self.NC)] for i in xrange(self.NC)]         
            deltah = [[1.99 + 0.20*(mi_ast[i][j]**2.0) for j in xrange(self.NC)] for i in xrange(self.NC)]
                 
            E = [[0.0 for j in xrange(self.NC)] for i in xrange(self.NC)] 
            for i in xrange(self.NC):
                for j in xrange(self.NC):
                    if Eta[i][j] < 4.5:
                        E[i][j] = exp(Eta[i][j]*(650.0/(ek[i][j]+300.0) - 4.27))
                    elif Eta[i][j] >= 4.5:
                        E[i][j] = exp(Eta[i][j]*(42800.0/(ek[i][j]+22400.0) - 4.27))
                  
            # Parâmetros dependentes da temperatura:
            T_ast   = [[self.Temp/ek[i][j]                  for j in xrange(self.NC)] for i in xrange(self.NC)]
            T_astll = [[1/T_ast[i][j] - 1.6*w[i][j] for j in xrange(self.NC)] for i in xrange(self.NC)]
            
            #Cálculos dos BF's:
            BFnonpolar = [[ b0[i][j]*(0.94 - 1.47*T_astll[i][j] - 0.85*(T_astll[i][j]**2.0) + 1.015*(T_astll[i][j]**3.0))           for j in xrange(self.NC)] for i in xrange(self.NC)]
            BFpolar    = [[-b0[i][j]*mi_astl[i][j]*(0.74 - 3.0*T_astll[i][j] + 2.1*(T_astll[i][j]**2.0) + 2.1*(T_astll[i][j]**3.0)) for j in xrange(self.NC)] for i in xrange(self.NC)]
         
            # Cálculos para BD:
            Bmetastable_Bbound = [[b0[i][j]*A[i][j]*exp(deltah[i][j]/T_ast[i][j]) for j in xrange(self.NC)] for i in xrange(self.NC)]
            Bchemical          = [[b0[i][j]*E[i][j]*(1 - exp(1500.0*Eta[i][j]/self.Temp)) for j in xrange(self.NC)] for i in xrange(self.NC)]
            # -----------------------------------------------------------------------------------------------------------------------------------------------------------
         
            BD            = [[Bmetastable_Bbound[i][j] + Bchemical[i][j] for j in xrange(self.NC)] for i in xrange(self.NC)] # D bound or dimerizes molecules (Chemical forces)
            BF            = [[BFnonpolar[i][j]         + BFpolar[i][j]   for j in xrange(self.NC)] for i in xrange(self.NC)] # Free molecules
            self.Bvirial  = [[BF[i][j]                 + BD[i][j]        for j in xrange(self.NC)] for i in xrange(self.NC)]
            
        elif self.model_vap.regra_mistura == 'Tsonopoulos' :
            
            R = 82.05746 # cm3.atm.K−1.mol−1            
            
            # Caracteristica da mistura quanto à polaridade. Ex.: mistura Polar-Polar.
            Caracteristica_mistura = self.Componente[0].polaridade + '-' +self.Componente[1].polaridade
                 
            w               = [[0.0 for j in xrange(self.NC)] for i in xrange(self.NC)]
                
            Tc              = [[0.0 for j in xrange(self.NC)] for i in xrange(self.NC)]
            Tr              = [[0.0 for j in xrange(self.NC)] for i in xrange(self.NC)]
            Pc              = [[0.0 for j in xrange(self.NC)] for i in xrange(self.NC)] 
            Vc              = [0.0 for j in xrange(self.NC)]
                
            mi              = [0.0 for j in xrange(self.NC)]
            mi_reduzido     = [0.0 for j in xrange(self.NC)]
                
            parametro_a     = [[0.0 for j in xrange(self.NC)] for i in xrange(self.NC)]
            parametro_b     = [[0.0 for j in xrange(self.NC)] for i in xrange(self.NC)] 
                
            # Funções da corelação de Tsonopoulos. Todas em função de Tr.
                
            parametro_f0    = [[0.0 for j in xrange(self.NC)] for i in xrange(self.NC)]
            parametro_f1    = [[0.0 for j in xrange(self.NC)] for i in xrange(self.NC)]
            parametro_f2    = [[0.0 for j in xrange(self.NC)] for i in xrange(self.NC)]
                 
            # Parâmetro de interação binária
                 
            k_int_binaria   = self.model_vap.k_int_binaria
            
            # PARÂMETROS PUROS
        
            for i in xrange(self.NC):
                for j in xrange(self.NC):
                    
                    if i == j:  
              
                        Tc[i][j]       = self.Componente[j].Tc           
                        Pc[i][j]       = self.Componente[j].Pc / 1.01325  # Deve ser em atm
                        w[i][j]        = self.Componente[j].w
                        Vc[j]          = self.Componente[j].Vc                
                        mi[j]          = self.Componente[j].dipole_moment    
                        
                        mi_reduzido[j] = (10**5)*(mi[j]**2)*Pc[i][j]/(Tc[i][j])**2
                        Tr[i][j]       = self.Temp/Tc[i][j]
                        
            # PARÂMETROS DE ASSOCIAÇÃO
            
            for i in xrange(self.NC):
                for j in xrange(self.NC):
                    
                    if self.Componente[j].grupo_funcional in ['Cetona','Aldeido','Alquil_nitrila','Eter','Acido_carboxilico', 'Ester']:
                        if i == j:
                            
                            parametro_a[i][j]        = -2.14e-4*mi_reduzido[j]-4.308e-21*(mi_reduzido[j])**8
                            parametro_b[i][j]        = 0     
                            
                    elif self.Componente[j].grupo_funcional in ['Haleto_organico', 'Mercaptan','Sulfeto', 'Dissulfeto']:
                        if i == j:
                            
                            parametro_a[i][j]        = -2.188e-11*(mi_reduzido[j])**4 - 7.831e-21*(mi_reduzido[j])**8
                            parametro_b[i][j]        = 0     
                
                    elif self.Componente[j].grupo_funcional == 'Alcool':
                        if self.Componente[j].nome != 'Metanol':
                            if i == j:
                                
                                parametro_a[i][j]    =  0.0878
                                parametro_b[i][j]    =  0.00908 + 0.0006957*mi_reduzido[j]
                        
                    elif self.Componente[j].nome == 'Metanol':
                        if i == j:
                            
                            parametro_a[i][j]        = 0.0878
                            parametro_b[i][j]        = 0.0525
                    
                    elif self.Componente[j].nome == 'Agua':
                        if i == j:
                            
                            parametro_a[i][j]        = -0.0109
                            parametro_b[i][j]        = 0.0

                    else:  
                        if i == j:
                            
                            parametro_a[i][j]        = 0.0
                            parametro_b[i][j]        = 0.0

            # PARÂMETROS CRUZADOS
            
            for i in xrange(self.NC):
                for j in xrange(self.NC):
                    
                    if Caracteristica_mistura == 'Polar-Polar':
                        if i!=j:
                            
                            parametro_a[i][j]  = 0.5*(parametro_a[i][i]+parametro_a[j][j])
                            parametro_b[i][j]  = 0.5*(parametro_b[i][i]+parametro_b[j][j])
                
                            w[i][j]            = 0.5*(w[i][i]+w[j][j])
                            
                            Tc[i][j]           = (Tc[i][i]*Tc[j][j])**0.5*(1-k_int_binaria[i][j])
                            Tr[i][j]           = self.Temp/Tc[i][j]
                            
                            Pc[i][j]           = 4*( Tc[i][j]*(Pc[i][i]*Vc[i]/Tc[i][i] + Pc[j][j]*Vc[j]/Tc[j][j])/((Vc[i]**(1.0/3.0) + Vc[j]**(1.0/3.0))**3) )
                    
                    elif Caracteristica_mistura == 'Apolar-Apolar': 
                        if i!=j:
                            
                            parametro_a[i][j]  = 0.0
                            parametro_b[i][j]  = 0.0
                            
# Estimativa para kij       k_int_binaria[i][j]= 1 - ( (2*(Vc[i]*Vc[j])**1.0/6.0)/(Vc[i]**(1.0/3.0) + Vc[j]**(1.0/3.0))  )**3
                            
                            w[i][j]            = 0.5*(w[i][i]+w[j][j])
                            
                            Tc[i][j]           = (Tc[i][i]*Tc[j][j])**0.5*(1-k_int_binaria[i][j])
                            Tr[i][j]           = self.Temp/Tc[i][j]
                            Pc[i][j]           = 4*( Tc[i][j]*(Pc[i][i]*Vc[i]/Tc[i][i] + Pc[j][j]*Vc[j]/Tc[j][j])/(Vc[i]**(1.0/3.0) + Vc[j]**(1.0/3.0) )**3 )
        
                    elif Caracteristica_mistura in ['Apolar-Polar','Polar-Apolar']: 
                        if i!=j:
                            
                            parametro_a[i][j]  = 0.0
                            parametro_b[i][j]  = 0.0
                            
                            w[i][j]            = 0.5*(w[i][i]+w[j][j])
                            
                            Tc[i][j]           = (Tc[i][i]*Tc[j][j])**0.5*(1-k_int_binaria[i][j])
                            Tr[i][j]           = self.Temp/Tc[i][j]
                            Pc[i][j]           = 4*( Tc[i][j]*(Pc[i][i]*Vc[i]/Tc[i][i] + Pc[j][j]*Vc[j]/Tc[j][j])/( Vc[i]**(1.0/3.0) + Vc[j]**(1.0/3.0) ) ** 3.0 )
 
       
            # PARAMETROS CÁLCULO Bij
                   
            for i in xrange(self.NC):
                for j in xrange(self.NC):
                    
                    # Funções da corelação de Tsonopoulos. Todas em função de Tr.
                    parametro_f0[i][j]    = 0.1445 - 0.330/Tr[i][j] - 0.1385/(Tr[i][j])**2 - 0.0121/(Tr[i][j])**3 - 0.000607/(Tr[i][j])**8
                    parametro_f1[i][j]    = 0.0637 + 0.331/(Tr[i][j])**2 - 0.423/(Tr[i][j])**3 - 0.008/(Tr[i][j])**8
                    parametro_f2[i][j]    = parametro_a[i][j]/(Tr[i][j])**6 - parametro_b[i][j]/(Tr[i][j])**8
            
            # CÁLCULO DO Bij
            
            self.Bvirial = [[ (Tc[i][j]*R/Pc[i][j])*(parametro_f0[i][j] + w[i][j]*parametro_f1[i][j] + parametro_f2[i][j]) for j in xrange(self.NC)] for i in xrange(self.NC)]



    def Coeficiente_Atividade(self,x,T):
        '''
        Módulo para calcular o coeficiente de atividade de acordo com os modelos disponíveis.
        Estes são: UNIQUAC[1] , NRTL[2], Wilson[3] e Van Laar[4].
        
        ========
        Entradas
        ========
        
        * x (list): Composições dos componentes da fase líquida;
        * T (float): Temperatura em Kelvin;
        
        ======
        Saídas
        ======
        
        * O método retorna os valores dos coeficientes de atividade dos componentes em forma de lista.
        
        ===========
        Referências
        ===========
        
        [1] ABRAMS, D. S.; PRAUSNITZ, J. M. Statistical thermodynamics of liquid
        mixtures: A new expression for the excess Gibbs energy of partly or completely
        miscible systems. AIChE Journal, v. 21, n. 1, p. 116–128, jan. 1975. ISSN
        0001-1541. Disponível em: <http://doi.wiley.com/10.1002/aic.690210115>
        
        [2] RENON, H.; PRAUSNITZ, J. M. Local compositions in thermodynamic excess
        functions for liquid mixtures. AIChE Journal, v. 14, n. 1, p. 135–144, jan. 1968.
        ISSN 0001-1541. Disponível em: <http://doi.wiley.com/10.1002/aic.690140124>.
        
        [3] WILSON, G. M. Vapor-Liquid Equilibrium. XI. A New Expression for the Excess
        Free Energy of Mixing. Journal of the American Chemical Society, v. 86, n. 2, p.
        127–130, jan. 1964. ISSN 0002-7863. Disponível em: <http://pubs.acs.org/doi/abs-
        /10.1021/ja01056a002>
        
        [4] VAN LAAR, J. J. The Vapor pressure of binary mixtures. Z. Phys.
        Chem. 1910, 72, 723−751.  
                
        '''                        
        # Modelo UNIQUAC
        if self.model_liq.nome_modelo == 'UNIQUAC':
            
            if self.model_liq.formaEq == 1:
#                R = 8.314462 # em J · K−1 · mol−1
                # Caso tenha em mãos a diferenca do parametro a. A formaEq 1 é a mais recorrente.
                a     = self.model_liq.parametro_int
    
                tau   = [[exp(-a[i][j]/(T)) for j in xrange(self.NC)] for i in xrange(self.NC)] 
                
                phi   = [self.Componente[i].r  * x[i] / sum([self.Componente[j].r  * x[j] for j in xrange(self.NC)])        for i in xrange(self.NC)]
                
                teta  = [self.Componente[i].q  * x[i] / sum([self.Componente[j].q  * x[j] for j in xrange(self.NC)])        for i in xrange(self.NC)]
                tetal = [self.Componente[i].ql * x[i] / sum([self.Componente[j].ql * x[j] for j in xrange(self.NC)])        for i in xrange(self.NC)]
                l     = [(self.coordnumber/2.0)*(self.Componente[i].r - self.Componente[i].q) - (self.Componente[i].r - 1)  for i in xrange(self.NC)]
    
                A     = [self.Componente[i].ql*sum([tetal[j]*tau[i][j]/sum([tetal[k]*tau[k][j] for k in xrange(self.NC)]) for j in xrange(self.NC)]) for i in xrange(self.NC)]
                
                Combinatorial = [log(phi[i]/x[i]) + (self.coordnumber/2.0)*self.Componente[i].q*log(teta[i]/phi[i]) + l[i] - (phi[i]/x[i]) * sum([x[j]*l[j] for j in xrange(self.NC)]) for i in xrange(self.NC)]
                Residual      = [-self.Componente[i].ql*log(sum([tetal[j]*tau[j][i] for j in xrange(self.NC)])) + self.Componente[i].ql - A[i]      for i in xrange(self.NC)] 
    
                Coeficiente_Atividade = [exp(Combinatorial[i]+Residual[i]) for i in xrange(self.NC)]
                
            elif self.model_liq.formaEq == 2:
                # Caso tenha em mãos o parametro tau 
                tau     = self.model_liq.parametro_int
    
                phi   = [self.Componente[i].r  * x[i] / sum([self.Componente[j].r  * x[j] for j in xrange(self.NC)])        for i in xrange(self.NC)]
                teta  = [self.Componente[i].q  * x[i] / sum([self.Componente[j].q  * x[j] for j in xrange(self.NC)])        for i in xrange(self.NC)]
                tetal = [self.Componente[i].ql * x[i] / sum([self.Componente[j].ql * x[j] for j in xrange(self.NC)])        for i in xrange(self.NC)]
                l     = [(self.coordnumber/2.0)*(self.Componente[i].r - self.Componente[i].q) - (self.Componente[i].r - 1)  for i in xrange(self.NC)]
    
                A     = [self.Componente[i].ql*sum([tetal[j]*tau[i][j]/sum([tetal[k]*tau[k][j] for k in xrange(self.NC)]) for j in xrange(self.NC)]) for i in xrange(self.NC)]
                
                Combinatorial = [log(phi[i]/x[i]) + (self.coordnumber/2.0)*self.Componente[i].q*log(teta[i]/phi[i]) + l[i] - (phi[i]/x[i]) * sum([x[j]*l[j] for j in xrange(self.NC)]) for i in xrange(self.NC)]
                Residual      = [-self.Componente[i].ql*log(sum([tetal[j]*tau[j][i] for j in xrange(self.NC)])) + self.Componente[i].ql - A[i]      for i in xrange(self.NC)] 
    
                Coeficiente_Atividade = [exp(Combinatorial[i]+Residual[i]) for i in xrange(self.NC)]

            elif self.model_liq.formaEq == 3:
                # Caso tenha em mãos o parametro a do componente puro
                a_ij  = self.model_liq.parametro_int
    
                tau   = [[exp(-(a_ij[i][j]-a_ij[j][j])/T) for j in xrange(self.NC)] for i in xrange(self.NC)]
                phi   = [self.Componente[i].r  * x[i] / sum([self.Componente[j].r  * x[j] for j in xrange(self.NC)])        for i in xrange(self.NC)]
                teta  = [self.Componente[i].q  * x[i] / sum([self.Componente[j].q  * x[j] for j in xrange(self.NC)])        for i in xrange(self.NC)]
                tetal = [self.Componente[i].ql * x[i] / sum([self.Componente[j].ql * x[j] for j in xrange(self.NC)])        for i in xrange(self.NC)]
                l     = [(self.coordnumber/2.0)*(self.Componente[i].r - self.Componente[i].q) - (self.Componente[i].r - 1)  for i in xrange(self.NC)]
    
                A     = [self.Componente[i].ql*sum([tetal[j]*tau[i][j]/sum([tetal[k]*tau[k][j] for k in xrange(self.NC)]) for j in xrange(self.NC)]) for i in xrange(self.NC)]
                
                Combinatorial = [log(phi[i]/x[i]) + (self.coordnumber/2.0)*self.Componente[i].q*log(teta[i]/phi[i]) + l[i] - (phi[i]/x[i]) * sum([x[j]*l[j] for j in xrange(self.NC)]) for i in xrange(self.NC)]
                Residual      = [-self.Componente[i].ql*log(sum([tetal[j]*tau[j][i] for j in xrange(self.NC)])) + self.Componente[i].ql - A[i]      for i in xrange(self.NC)] 
    
                Coeficiente_Atividade = [exp(Combinatorial[i]+Residual[i]) for i in xrange(self.NC)]                
        #==============================================================================
        # Modelo NRTL
        #==============================================================================
        
        if self.model_liq.nome_modelo == 'NRTL':
            
            R = 83.144621# em cm3.bar/ K.mol
            alpha       = self.model_liq.alpha            
            
            if self.model_liq.formaEq == 1:
                # Caso tenha em mãos a diferenca do parametro g. A formaEq 1 é a mais recorrente.
                g      = self.model_liq.parametro_int                                
                
                tau    =  [[g[i][j]/(R*T) for j in xrange(self.NC)] for i in xrange(self.NC)]
                G      = [[exp(-alpha[i][j]*tau[i][j]) for j in xrange(self.NC)] for i in xrange(self.NC)]
    
                parte1 = [(sum([tau[j][i]*G[j][i]*x[j] for j in xrange(self.NC)]))/(sum([G[k][i]*x[k] for k in xrange(self.NC)])) for i in xrange(self.NC)]
                parte2 = [sum( [( x[j]*G[i][j]/sum([G[k][j]*x[k] for k in xrange(self.NC)]) )*(tau[i][j] -sum([x[k]*tau[k][j]*G[k][j] for k in xrange(self.NC)])/sum([G[k][j]*x[k] for k in xrange(self.NC)]) )   for j in xrange(self.NC)] ) for i in xrange(self.NC)]
    
                Coeficiente_Atividade = [exp( parte1[i] + parte2[i] ) for i in xrange(self.NC)]             
                
            elif self.model_liq.formaEq == 2:
                # Caso tenha em mãos o parametro tau 
                tau    = self.model_liq.parametro_int                                
            
                G      = [[exp(-alpha[i][j]*tau[i][j]) for j in xrange(self.NC)] for i in xrange(self.NC)]
                
                parte1 = [(sum([tau[j][i]*G[j][i]*x[j] for j in xrange(self.NC)]))/(sum([G[k][i]*x[k] for k in xrange(self.NC)])) for i in xrange(self.NC)]
                parte2 = [sum( [( x[j]*G[i][j]/sum([G[k][j]*x[k] for k in xrange(self.NC)]) )*(tau[i][j] -sum([x[k]*tau[k][j]*G[k][j] for k in xrange(self.NC)])/sum([G[k][j]*x[k] for k in xrange(self.NC)]) )   for j in xrange(self.NC)] ) for i in xrange(self.NC)]
    
                Coeficiente_Atividade = [exp( parte1[i] + parte2[i] ) for i in xrange(self.NC)]  
                
            elif self.model_liq.formaEq == 3:
                # Caso tenha em mãos o parametro g do componente puro 
                g_ij      = self.model_liq.parametro_int                                
                
                tau    =  [[(g_ij[i][j]-g_ij[j][j])/(R*T) for j in xrange(self.NC)] for i in xrange(self.NC)]
                G      = [[exp(-alpha[i][j]*tau[i][j]) for j in xrange(self.NC)] for i in xrange(self.NC)]
                    
                parte1 = [(sum([tau[j][i]*G[j][i]*x[j] for j in xrange(self.NC)]))/(sum([G[k][i]*x[k] for k in xrange(self.NC)])) for i in xrange(self.NC)]
                parte2 = [sum( [( x[j]*G[i][j]/sum([G[k][j]*x[k] for k in xrange(self.NC)]) )*(tau[i][j] -sum([x[k]*tau[k][j]*G[k][j] for k in xrange(self.NC)])/sum([G[k][j]*x[k] for k in xrange(self.NC)]) )   for j in xrange(self.NC)] ) for i in xrange(self.NC)]
    
                Coeficiente_Atividade = [exp( parte1[i] + parte2[i] ) for i in xrange(self.NC)]
                
        #==============================================================================
        # Modelo de Wilson                               
        #==============================================================================
        
        if self.model_liq.nome_modelo == 'Wilson':
            
            if self.model_liq.formaEq == 1:
            # Caso tenha em mãos o parametro LAMBDA. 
                A          =   self.model_liq.parametro_int
            
                Coeficiente_Atividade = [exp( -log( sum([x[j]*A[i][j] for j in xrange(self.NC)]) ) + 1.0 - sum([ x[k]*A[k][i]/(sum([x[j]*A[k][j] for j in xrange(self.NC)])) for k in xrange(self.NC) ]) ) for i in xrange(self.NC) ]
               
       # Modelo Van Laar
        if self.model_liq.nome_modelo == 'Van Laar':
            
            R  = 83.144621 # em cm3.bar/ K.mol
            
            p_VL = self.model_liq.Parametro
            # Onde A é o parâmetro p_VL[0][1] e B é o parâmetro p_VL[1][0]
            
            Coeficiente_Atividade_aux   =    [[exp( ( p_VL[i][j]/(R*T) )*( 1+(p_VL[i][j]/p_VL[j][i])*\
            (x[i]/x[j]) )**-2 ) for j in xrange(self.NC) if i!=j] for i in xrange(self.NC)]
                        
            Coeficiente_Atividade       =    [Coeficiente_Atividade_aux[0][0],Coeficiente_Atividade_aux[1][0]]
                
        return Coeficiente_Atividade

    def Coeficiente_Fugacidade(self,y,P,T):
        '''
        Módulo para calcular o coeficiente de fugacidade de acordo com as equações de estado disponíveis.
        Estas são: Virial[1].
    
        ========
        Entradas
        ========
        
        * y (list): Composições dos componentes da fase de vapor.
        * T (float): Temperatura em Kelvin.
        * P (float): Pressão do sistema em bar.
        
        ======
        Saídas
        ======
        
        * O método retorna os valores dos coeficientes de fugacidade dos componentes em forma de lista.
        
        ===========
        Referências
        ===========
        
        [1] MASON, E. A.; SPURLING, T. H. The Virial Equation of State; The
        International Encyclopedia of Physical Chemistry and Chemical
        Physics, Topic 10: The Fluid State, Vol. 2; Pergamon Press: New
        York, 1969; p 297
        
        '''   
        R = 83.144621 # em cm3.bar/ K.mol
        if self.model_vap.nome_modelo == 'Virial':
            NC = size(y)
            self.Second_Virial_Coef()
            B = self.Bvirial
            Bmixture = sum([sum([y[i]*y[j]*B[i][j] for j in xrange(NC)])      for i in xrange(NC)])
            A        = [ 2*sum([y[j]*B[i][j] for j in xrange(NC)]) - Bmixture for i in xrange(NC)]
            phi      = [exp(A[i]*P/(R*T)) for i in xrange(NC)]
            
            # Validação grosseira das condições de pressão:
            P_lim = (T/2.0)*sum([y[i]*self.Componente[i].Pc for i in xrange(NC)])/sum([y[i]*self.Componente[i].Tc for i in xrange(NC)])
            if P <= P_lim:
                warn(u'A pressão do sistema é inferior à da validação da equação VIRIAL, vide documentação da mesma.')
            else:
                warn(u'A pressão do sistema é superior à da validação da equação VIRIAL, vide documentação da mesma.')
                
        return phi
        
    def PhiSat(self,T):
        '''
        Módulo para calcular o coeficiente de fugacidade nas condições de saturação segundo [1] e [2].
        
        ======
        Saídas
        ======
        
        * O método retorna os valores dos coeficientes de fugacidade dos componentes em forma de lista.
        
        ===========
        Referências
        ===========
        
        [1] PRAUSNITZ, J. M. et al. Computer Calculations for multicomponent vapor-liquid and liquid-liquid equilibria. [s.l.] Prendice-Hall, 1980. p. 353
        
        [2] SMITH, J. M.; NESS, H. C. VAN; ABBOTT, M. M. Introduction to Chemical Engineering Thermodinamics. 7th. ed. [s.l.] Mc-Graw Hills, [s.d.]. 
        
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
            phisat.append(self.Coeficiente_Fugacidade(comp,self.Componente[i].Pvap_Prausnitz_4th(self.Componente[i].VPA,self.Componente[i].VPB,self.Componente[i].VPC,self.Componente[i].VPD,T,self.Componente[i].Tc,self.Componente[i].Pc),T)[i])
        self.phisat = phisat
        
    
    def PontoBolha_P(self,x,T):
        '''
        Módulo para calcular o ponto de bolha segundo [1] e [2], quando a temperatura e composição são conhecidas. 
        
        ========
        Entradas
        ========
        
        * T (float): Temperatura em Kelvin;
        * x (list): Composição da fase líquida.
            
        ======
        Saídas
        ======
        
        As seguintes saídas são em forma de atributos.
        
        * ``liquido``: Um objeto da classe ``Condicao``, vide documentação da classe;
        * ``Bolha``: Um objeto da classe ``Condicao``, vide documentação da classe.
        
        ===========
        Referências
        ===========
        
        [1] PRAUSNITZ, J. M. et al. Computer Calculations for multicomponent vapor-liquid and liquid-liquid equilibria. [s.l.] Prendice-Hall, 1980. p. 353
        
        [2] SMITH, J. M.; NESS, H. C. VAN; ABBOTT, M. M. Introduction to Chemical Engineering Thermodinamics. 7th. ed. [s.l.] Mc-Graw Hills, [s.d.].                 
        ''' 
        #==============================================================================
        #------------------ Estimativas iniciais --------------------------------------
        #==============================================================================                
        x = [x[i]/(sum([x[i] for i in xrange(self.NC)])) for i in xrange(self.NC)] # Normalização das composições
        self.PhiSat(T)
        P        = []
        coefAct  = self.Coeficiente_Atividade(x,T)
        # Caracterização da fase
        coeffug  = self.estphi
        
        cont   = 0; deltaP = 10000
        while (deltaP > self.tolAlg) and (cont<self.maxiter+1):
            # Atualização do valor de P por VLE
            P.append(sum([x[i]*coefAct[i]*self.Componente[i].Pvap_Prausnitz_4th(self.Componente[i].VPA,self.Componente[i].VPB,self.Componente[i].VPC,self.Componente[i].VPD,T,self.Componente[i].Tc,self.Componente[i].Pc)*self.phisat[i]/coeffug[i]        for i in xrange(self.NC)]))
            # Cálculo de y por VLE
            y       = [x[i]*coefAct[i]*self.Componente[i].Pvap_Prausnitz_4th(self.Componente[i].VPA,self.Componente[i].VPB,self.Componente[i].VPC,self.Componente[i].VPD,T,self.Componente[i].Tc,self.Componente[i].Pc)*self.phisat[i]/(coeffug[i]*P[cont]) for i in xrange(self.NC)]
            # Normalização do valor de y
            y        = [y[i]/(sum([y[i] for i in xrange(self.NC)])) for i in xrange(self.NC)]
            # Atualização de phi por EoS
            coeffug = self.Coeficiente_Fugacidade(y,P[cont],T)
            
            deltaP = abs(P[cont] - P[cont-1])
            cont+=1
            
        # Caracterização das fases
        self.Bolha   = Condicao(P[cont-1],T,y,coeffug,None)
        self.liquido = Condicao(P[cont-1],T,x,None,coefAct)

    def PontoBolha_T(self,x,P):
        ''' 
        Módulo para calcular o ponto de bolha segundo [1] e [2], quando a pressão e composição são conhecidas.

        ========
        Entradas
        ========
        
        
        * P (float): Pressão em bar;
        * x (list): Composição da fase líquida.
        
        ======
        Saídas
        ======
        
        As seguintes saídas são em forma de atributos.
        
        * ``liquido``: Um objeto da classe ``Condicao``, vide documentação da classe;
        * ``Bolha``: Um objeto da classe ``Condicao``, vide documentação da classe.
        
        ===========
        Referências
        ===========
        
        [1] PRAUSNITZ, J. M. et al. Computer Calculations for multicomponent vapor-liquid 
        and liquid-liquid equilibria. [s.l.] Prendice-Hall, 1980. p. 353
        
        [2] SMITH, J. M.; NESS, H. C. VAN; ABBOTT, M. M. Introduction to Chemical 
        Engineering Thermodinamics. 7th. ed. [s.l.] Mc-Graw Hills, [s.d.]. 
        '''
        #==============================================================================
        #------------------ Estimativas iniciais --------------------------------------
        #==============================================================================
        
        # Normalização do valor de x
        x = [x[i]/(sum([x[i] for i in xrange(self.NC)])) for i in xrange(self.NC)]
        T = [sum([self.Componente[i].Tsat_Prausnitz_4th(self.Componente[i].VPA,self.Componente[i].VPB,self.Componente[i].VPC,self.Componente[i].VPD,P,self.Componente[i].Tc,self.Componente[i].Pc)*x[i] for i in xrange(self.NC)])]
        
        coeffug  = self.estphi
        cont   = 0; deltaT = 10        
        while (deltaT > self.tolAlg) and (cont<self.maxiter+1):
            # cálculo da pressão de saturação P_i^(sat) por Prausnitz
            psat_ini     = [self.Componente[i].Pvap_Prausnitz_4th(self.Componente[i].VPA,self.Componente[i].VPB,self.Componente[i].VPC,self.Componente[i].VPD,T[cont],self.Componente[i].Tc,self.Componente[i].Pc) for i in xrange(self.NC)]            
            # Cálculo de phisat
            self.PhiSat(T[cont])
            # Cálculo de gamma por modelos termodinâmicos
            coefAct  = self.Coeficiente_Atividade(x,T[cont])            
            # Predição de y por VLE
            y        = [x[i]*coefAct[i]*psat_ini[i]*self.phisat[i]/(coeffug[i]*P) for i in xrange(self.NC)]            
            # Normalização do valor de y
            y        = [y[i]/(sum([y[i] for i in xrange(self.NC)])) for i in xrange(self.NC)]
            # Atualização de phi por EoS
            coeffug  = self.Coeficiente_Fugacidade(y,P,T[cont])
            # Cálculo da pressão de saturação P_j^(sat) por VLE
            psat     = [P/(sum([(x[i]*coefAct[i]*self.phisat[i]*psat_ini[i])/(psat_ini[j]*coeffug[i]) for i in xrange(self.NC)])) for j in xrange(self.NC)]            
            # Atualização do valor de T por Prausnitz
            T.append(self.Componente[1].Tsat_Prausnitz_4th(self.Componente[1].VPA,self.Componente[1].VPB,self.Componente[1].VPC,self.Componente[1].VPD,psat[1],self.Componente[1].Tc,self.Componente[1].Pc))
            # Atualização do valor de deltaT
            deltaT  = abs((T[cont+1] - T[cont])/T[cont])
            cont+=1
            
        # Caracterização das fases
        self.liquido = Condicao(P,T[cont],x,None,coefAct)
        self.Bolha   = Condicao(P,T[cont],y,coeffug,None)
        
    def PontoOrvalho_P(self,y,T):
        ''' 
        Módulo para calcular o ponto de orvalho segundo [1] e [2], quando a temperatura e composição são conhecidas.

        ========
        Entradas
        ========
        
        
        * T (float): Temperatura em Kelvin;
        * y (list): Composição da fase vapor.
        
        ======
        Saídas
        ======
        
        As seguintes saídas são em forma de atributos.
        
        * ``vapor``: Um objeto da classe ``Condicao``, vide documentação da classe;
        * ``Orvalho``: Um objeto da classe ``Condicao``, vide documentação da classe.
        
        ===========
        Referências
        ===========
        
        [1] PRAUSNITZ, J. M. et al. Computer Calculations for multicomponent vapor-liquid 
        and liquid-liquid equilibria. [s.l.] Prendice-Hall, 1980. p. 353
        
        [2] SMITH, J. M.; NESS, H. C. VAN; ABBOTT, M. M. Introduction to Chemical 
        Engineering Thermodinamics. 7th. ed. [s.l.] Mc-Graw Hills, [s.d.]. 
        ''' 
        #==============================================================================
        #------------------ Estimativas iniciais --------------------------------------
        #==============================================================================        
        # Normalização do valor de y
        y        = [y[i]/(sum([y[i] for i in xrange(self.NC)])) for i in xrange(self.NC)]
        self.PhiSat(T)
        P        = [self.Pressao]
        coeffug  = self.estphi
        coefAct  = self.estgama        
        cont   = 0; deltaP = 10000
        
        while (deltaP > self.tolAlg) and (cont<self.maxiter+1):
            # Atualização do valor de P por VLE
            P.append(1/sum([y[i]*coeffug[i]/(coefAct[i]*self.Componente[i].Pvap_Prausnitz_4th(self.Componente[i].VPA,self.Componente[i].VPB,self.Componente[i].VPC,self.Componente[i].VPD,T,self.Componente[i].Tc,self.Componente[i].Pc)*self.phisat[i])    for i in xrange(self.NC)]))
            # Cálculo de x por VLE
            x       = [y[i]*coeffug[i]*P[cont]/(coefAct[i]*self.Componente[i].Pvap_Prausnitz_4th(self.Componente[i].VPA,self.Componente[i].VPB,self.Componente[i].VPC,self.Componente[i].VPD,T,self.Componente[i].Tc,self.Componente[i].Pc)*self.phisat[i]) for i in xrange(self.NC)]
            # Normalização de x
            x       = [x[i]/(sum([x[i] for i in xrange(self.NC)])) for i in xrange(self.NC)]
            # Cálculo de phi por EoS
            coeffug = self.Coeficiente_Fugacidade(y,P[cont],T)
            # Cálculo de gamma por modelos termodinamicos
            coefAct = self.Coeficiente_Atividade(x,T)
            
            deltaP = abs(P[cont] - P[cont-1])
            cont+=1

        self.vapor   = Condicao(P[cont-1],T,y,coeffug,None)
        self.Orvalho = Condicao(P[cont-1],T,x,None,coefAct)

        
    def PontoOrvalho_T(self,y,P):
        ''' 
        Módulo para calcular o ponto de orvalho segundo [1] e [2], quando a pressão e composição são conhecidas.

        ========
        Entradas
        ========
        
        * P (float): Pressão em bar;
        * y (list): Composição da fase vapor.
        
        ======
        Saídas
        ======
        
        As seguintes saídas são em forma de atributos.
        
        * ``vapor``: Um objeto da classe ``Condicao``, vide documentação da classe;
        * ``Orvalho``: Um objeto da classe ``Condicao``, vide documentação da classe.
        
        ===========
        Referências
        ===========
        
        [1] PRAUSNITZ, J. M. et al. Computer Calculations for multicomponent vapor-liquid 
        and liquid-liquid equilibria. [s.l.] Prendice-Hall, 1980. p. 353
        
        [2] SMITH, J. M.; NESS, H. C. VAN; ABBOTT, M. M. Introduction to Chemical 
        Engineering Thermodinamics. 7th. ed. [s.l.] Mc-Graw Hills, [s.d.]. 
        ''' 
        #==============================================================================
        #------------------ Estimativas iniciais --------------------------------------
        #==============================================================================
        # Normalização do valor de y
        y        = [y[i]/(sum([y[i] for i in xrange(self.NC)])) for i in xrange(self.NC)]
        T = [sum([self.Componente[i].Tsat_Prausnitz_4th(self.Componente[i].VPA,self.Componente[i].VPB,self.Componente[i].VPC,self.Componente[i].VPD,P,self.Componente[i].Tc,self.Componente[i].Pc)*y[i] for i in xrange(self.NC)])]
        
        coeffug  = [1.0,1.0]
        coefAct  = [[1.0,1.0]]
        cont   = 0; deltaT = 10
        
        while (deltaT > self.tolAlg) and (cont<self.maxiter+1):
            # cálculo da pressão de saturação P_i^(sat) por Prausnitz
            psat_ini     = [self.Componente[i].Pvap_Prausnitz_4th(self.Componente[i].VPA,self.Componente[i].VPB,self.Componente[i].VPC,self.Componente[i].VPD,T[cont],self.Componente[i].Tc,self.Componente[i].Pc) for i in xrange(self.NC)]            
            # Cálculo de phisat
            self.PhiSat(T[cont])            
            # Atualização de phi por EoS
            coeffug  = self.Coeficiente_Fugacidade(y,P,T[cont])
            # Loop interno para o cálculo de gamma
            contador = 0 ; delta_gamma = 10
            while (delta_gamma > self.tolAlg) and (contador<self.maxiter+1):
                # Predição de x por equilíbrio
                x = [y[i]*coeffug[i]*P/(coefAct[contador][i]*psat_ini[i]*self.phisat[i]) for i in xrange(self.NC)]
                # Normalização do valor de x
                x = [x[i]/(sum([x[i] for i in xrange(self.NC)])) for i in xrange(self.NC)]
                # Atualização do valor de gamma por modelos termodinâmicos
                coefAct.append(self.Coeficiente_Atividade(x,T[cont]))
                # Atualização do valor de delta_gamma
                delta_gamma  = abs((coefAct[contador+1][0] - coefAct[contador][0])/coefAct[contador][0])
                contador+=1 
                
            # Cálculo da pressão de saturação P_j^(sat) por VLE
            psat     = [P*(sum([(y[i]*coeffug[i]*psat_ini[j])/(coefAct[contador][i]*self.phisat[i]*psat_ini[i]) for i in xrange(self.NC)])) for j in xrange(self.NC)]            
            # Atualização do valor de T por Prausnitz
            T.append(self.Componente[1].Tsat_Prausnitz_4th(self.Componente[1].VPA,self.Componente[1].VPB,self.Componente[1].VPC,self.Componente[1].VPD,psat[1],self.Componente[1].Tc,self.Componente[1].Pc))
            # Atualização do valor de deltaT
            deltaT  = abs((T[cont+1] - T[cont])/T[cont])
            cont+=1
            
        # Caracterização da fase vapor
        self.vapor   = Condicao(P,T[cont],y,coeffug,None)
        self.Orvalho = Condicao(P,T[cont],x,None,coefAct)
        
    def Flash(self,z,T,P):
        '''        
        Módulo para realizar o calculo de flash segundo [1] e [2], dada pressão, composições globais e temperatura.    

        ========
        Entradas
        ========
        
        * z (list): Composição global da mistura;
        * T (float): Temperatura em Kelvin;
        * P (float): Pressão em bar.
        
        ======
        Saídas
        ======
        
        As seguintes saídas são em forma de atributos.
        
        * ``Bolha``: Um objeto da classe ``Condicao``, vide documentação da classe;
        * ``Orvalho``: Um objeto da classe ``Condicao``, vide documentação da classe;
        * ``Beta``: Um objeto da classe ``Condicao``, vide documentação da classe.
        
        ===========
        Referências
        ===========
        
        [1] PRAUSNITZ, J. M. et al. Computer Calculations for multicomponent vapor-liquid 
        and liquid-liquid equilibria. [s.l.] Prendice-Hall, 1980. p. 353
        
        [2] SMITH, J. M.; NESS, H. C. VAN; ABBOTT, M. M. Introduction to Chemical 
        Engineering Thermodinamics. 7th. ed. [s.l.] Mc-Graw Hills, [s.d.]. 
        '''
        # TDODO: documentaao
        z    = [z[i]/(sum([z[i] for i in xrange(self.NC)])) for i in xrange(self.NC)]
        Psat = [self.Componente[i].Pvap_Prausnitz_4th(self.Componente[i].VPA,self.Componente[i].VPB,self.Componente[i].VPC,self.Componente[i].VPD,T,self.Componente[i].Tc,self.Componente[i].Pc) for i in xrange(self.NC)]            
        self.PontoBolha_P(z,T)
        self.PontoOrvalho_P(z,T)

        if (P < self.Bolha.Pressao) and (P > self.Orvalho.Pressao):
            
            interp = (P - self.Orvalho.Pressao)/(self.Bolha.Pressao - self.Orvalho.Pressao)
            coefAct = [(self.liquido.coefAct[i] - self.Orvalho.coefAct[i])*interp + self.Orvalho.liquido.coefAct[i]  for i in xrange(self.NC) ]
            coeffug = [(self.Bolha.coeffug[i]   - self.vapor.coeffug[i])*interp   + self.Orvalho.vapor.coefFug[i]    for i in xrange(self.NC) ]
            
            self.PhiSat(T)
            
            V    = [(self.Bolha.Pressao - self.Pressao)/(self.Bolha.Pressao - self.Orvalho.Pressao)]
            x    = [z]
            y    = [z]
            cont = 0; deltaV = 1e10; delta_x = 1e4; delta_y = 1e4            
            while (deltaV>self.tolAlg) and (delta_x>self.toleq) and (delta_y>self.toleq) and (cont<(self.maxiter)) :
                
                # Cálculo de K
                
                K = [(coefAct[i]*Psat[i]*self.phisat[i]) / (coeffug[i]*P) for i in xrange(self.NC)]
		
                # Cálculo de F e sua derivada em relação à fração de vapor
                F     = sum([z[i]*(K[i]-1)    / (1+V[cont]*(K[i]-1))    for i in xrange(self.NC)])
                dFdV = (-1)*sum([z[i]*(K[i]-1)**2 / (1+V[cont]*(K[i]-1))**2 for i in xrange(self.NC)])
		
                # Aplicação do método de Newton
                V.append(V[cont] - F/dFdV)

                # Cálculo das composições e subsequente normalização:
                x_aux = [z[i] / (1 + V[cont+1]*(K[i]-1)) for i in xrange(self.NC)]
                x.append([x_aux[i]/(sum([x_aux[i] for i in xrange(self.NC)])) for i in xrange(self.NC)])
            
                y_aux = [K[i]*x[i] for i in xrange(self.NC)]
                y.append([y_aux[i]/(sum([y_aux[i] for i in xrange(self.NC)])) for i in xrange(self.NC)])

                # Cálculo das fugacidades
                coefAct = self.Coeficiente_Atividade(x[cont+1],T)
                coeffug = self.Coeficiente_Fugacidade(y[cont+1],P,T)
                
                # Diferença entre valores de x
                delta_x = abs((x[cont+1][0] - x[cont][0])/x[cont][0])
                
                # Diferença entre valores de y
                delta_y = abs((y[cont+1][0] - y[cont][0])/y[cont][0])
                
                # Diferença entre valores de beta
                deltaV = abs((V[cont+1] - V[cont])/V[cont])
                
                cont+=1
                
          # Configuração das fases
	    self.Orvalho = Condicao(P,T,x[cont],None,coefAct) # configuração da fase líquida
	    self.Bolha   = Condicao(P,T,y[cont],coeffug,None) # configuração da fase vapor
	    self.Beta     = V[cont]
	else:
	    
	    raise ValueError(u'Não é possível realizar o cálculo de Flash, dado que a condição de equilíbrio não é satisfeita.') 
        
    def Predicao(self,Constante,Valor_cte):
        '''
        Metodo para caracterização dos eixos Ox e Oy para a realização dos gráficos.
        
        ========
        Entradas
        ========
        
        * Constante (str): Nome da variável que será mantida constante: Temperatura ou Pressao;
        * Valor_cte (float): Valor da constante de acordo com a variável inserida em *Constante*.
        
        ======
        Saídas
        ======
        
        As seguintes saídas são em forma de atributos.
        
        * ``Bolha``: Um objeto da classe ``Condicao``, vide documentação da classe;
        * ``Orvalho``: Um objeto da classe ``Condicao``, vide documentação da classe.
        
        '''
        # Validação da string Constante
        if Constante not in ['Pressao','Temperatura']:
            
            raise NameError(u'A constante escolhida deve ser a temperatura ou pressão.')  # Emite um erro com a mensagem inserida no método
        
        
        if Constante == 'Temperatura':        
            
            T = Valor_cte
            # Criação do eixo X para fazer os gráficos
            z_1 = linspace(1e-13,0.1,1000).tolist() # Devido à união das pontas, o passo nas extremidades é menor
            z_2 = linspace(0.1,0.9,500).tolist()
            z_3 = linspace(0.9,0.9999999999999,1000).tolist() # Devido à união das pontas, o passo nas extremidades é menor
            
            z = z_1+z_2+z_3 # Forma~çao do eixo X, eixo das composições, completo

            # Criação das listas vazias
            y_1 = []
            y_2 = []
            x_1 = []
            x_2 = []
            Pressao_Ponto_Bolha   = []
            Pressao_Ponto_Orvalho = []
            
            # Realiza o cálculo do ponto de bolha e de orvalho de cada par de concetrações
            for i in xrange(len(z)):
            
                # Cálculos
                self.PontoBolha_P([z[i],1-z[i]],T)
                self.PontoOrvalho_P([z[i],1-z[i]],T)
                          
                # Preenchimento das listas          
                y_1.append(self.Bolha.comp_molar[0])
                y_2.append(self.Bolha.comp_molar[1])
                x_1.append(self.Orvalho.comp_molar[0])
                x_2.append(self.Orvalho.comp_molar[1])
                Pressao_Ponto_Bolha.append(self.Bolha.Pressao)
                Pressao_Ponto_Orvalho.append(self.Orvalho.Pressao)

            # caracterização das fases
            self.Bolha   = Condicao(Pressao_Ponto_Bolha,T,[y_1,y_2],None,None)
            self.Orvalho = Condicao(Pressao_Ponto_Orvalho,T,[x_1,x_2],None,None)
        
        elif Constante == 'Pressao':
            
            P = Valor_cte
            # Criação do eixo X para fazer os gráficos
            z_1 = linspace(1e-13,0.1,1000).tolist() # Devido à união das pontas, o passo nas extremidades é menor
            z_2 = linspace(0.1,0.9,500).tolist()
            z_3 = linspace(0.9,0.9999999999999,1000).tolist() # Devido à união das pontas, o passo nas extremidades é menor
            
            z = z_1+z_2+z_3 # Forma~çao do eixo X, eixo das composições, completo

            # Criação das listas vazias
            y_1 = []
            y_2 = []
            x_1 = []
            x_2 = []
            Temperatura_Ponto_Bolha   = []
            Temperatura_Ponto_Orvalho = []
            
            # Realiza o cálculo do ponto de bolha e de orvalho de cada par de concetrações
            for i  in xrange(len(z)):
                
                # Cálculos
                self.PontoBolha_T([z[i],1-z[i]],P)
                self.PontoOrvalho_T([z[i],1-z[i]],P)
                          
                # Preenchimento das listas    
                y_1.append(self.Bolha.comp_molar[0])
                y_2.append(self.Bolha.comp_molar[1])
                x_1.append(self.Orvalho.comp_molar[0])
                x_2.append(self.Orvalho.comp_molar[1])
                Temperatura_Ponto_Bolha.append(self.Bolha.Temp)
                Temperatura_Ponto_Orvalho.append(self.Orvalho.Temp)
            
            # Caracterização das fases
            self.Bolha   = Condicao(P,Temperatura_Ponto_Bolha,[y_1,y_2],None,None)
            self.Orvalho = Condicao(P,Temperatura_Ponto_Orvalho,[x_1,x_2],None,None)

    def run(self):
        
        # Validação do algoritmo
        algoritmos_disponiveis = ['Coeficiente_Fugacidade','Coeficiente_Atividade','PontoBolha_P','PontoBolha_T','PontoOrvalho_P','PontoOrvalho_T','Flash']
        if self.Algoritmo not in algoritmos_disponiveis:
            
            raise NameError(u'O algoritmo escolhido não consta na lista de algoritmo disponíveis: '+', '.join(algoritmos_disponiveis)+'.')
            
        if self.Algoritmo == 'Coeficiente_Fugacidade':
            
            self.coefFug = self.Coeficiente_Fugacidade(self.z,self.Pressao,self.Temp)                        
            
        elif self.Algoritmo == 'Coeficiente_Atividade':
                        
            self.coefAct = self.Coeficiente_Atividade(self.z,self.Temp)
        
        elif self.Algoritmo == 'PontoBolha_P':
            
            self.PontoBolha_P(self.z,self.Temp)

        elif self.Algoritmo == 'PontoBolha_T':
            
            self.PontoBolha_T(self.z,self.Pressao)
            
        elif self.Algoritmo == 'PontoOrvalho_P':
                        
            self.PontoOrvalho_P(self.z,self.Temp)
        
        elif self.Algoritmo == 'PontoOrvalho_T':
            
            self.PontoOrvalho_T(self.z,self.Pressao)

        elif self.Algoritmo == 'Flash':
            
            self.Flash(self.z,self.Temp,self.Pressao)