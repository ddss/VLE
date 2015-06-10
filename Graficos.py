# -*- coding: utf-8 -*-
"""
Created on Wed Feb 04 11:01:10 2015

Rotina para a plotagem de gráficos. Os gráficos possíveis são:

- Diagrama de (P) x (x,y);
- Diagrama de (T) x (x,y);
- Diagrama de (x) x (y).

Vide [1] para maiores informações sobre os diagramas.

Referências:

[1] SMITH, J. M.; NESS, H. C. V.; ABBOTT, M. M. Introduction to 
    Chemical Engineering Thermodinamics. 7. ed. [S.l.]: Mc-Graw Hills, 2004.        

@author: CaiqueFerreira
"""

from matplotlib.pyplot import plot, xlabel, ylabel, title, xlim, legend, grid, show, figure, errorbar

class Graficos:

    def __init__(self,Diagrama,VLE,x_experimentais = None,y_experimentais = None, y_incertezas = None, x_incertezas = None, T = None, P = None, P_incertezas = None, T_incertezas = None):
        '''
        Algoritmo para criação das saídas gráficas dos cálculos envolvendo o equilíbrio líquido vapor, vide [1].
        
        ========
        Entradas
        ========

        * Diagrama (str): O tipo de diagrama ou gráficor que se pretende realizar;
        * VLE (objeto): objeto VLE após execuação do método Predicao executado
        * x_experimentais (list):
        * y_experimentais (array):
        * y_incertezas (array):
        * T (float): A temperatura da mistura em Kelvin;
        * P (float): A pressão em bar.
        
        =======
        Métodos
        =======
        
        * ``P_x_y``:
            * Método para criação do diagrama Pxy, pressão em função das composições, vide [1].
        * ``T_x_y``:
            * Método para criação do diagrama Txy, temperatura em função das composições, vide [1].
        * ``x_y``:
            * Método para criação do diagrama xy, composição da fase de vapor em função da composicção da fase líquida,vide [1].
        
        =========       
        Exemplo 1
        =========
        
        A priori, é necessário utilizar a rotina ``Conexao`` para caracterizar os componentes e os modelos.
        Neste exemplo, será calculado o ponto de bolha utilizando o método correspondente. A mistura utilizada será Acetona-Metanol, os modelos para as fases líquida e de vapor são, respectivamente
        UNIQUAC e VIRIAL com a regra de Hayden O'Connel, a temperatura será 330.0 K e a pressão 1.013 bar. ::
        
            from Conexao import Componente_Caracterizar, UNIQUAC, VIRIAL
            from VLE import VLE
            
            Comp1 = Componente_Caracterizar('Acetona',ConfigPsat=('Prausnitz4th',1),T=330.0)
            Comp2 = Componente_Caracterizar('Metanol',ConfigPsat=('Prausnitz4th',1),T=330.0)
                
            Componentes = [Comp1,Comp2]
            
            model_vap = VIRIAL(Componentes)
            model_liq = UNIQUAC(Componentes,330.0,1)
            
            CalculoBolha = VLE('PontoBolha',Componentes,model_liq,model_vap,z=[0.95,0.05],Temp=330.0,Pressao=1.013,estgama=None,estphi=None,estBeta = 0.5,tolAlg=1e-5,toleq=1e-4,maxiter=500,z_coordenacao=10.0)
        
        Após este passo, é acessado o método de predição, para caracterizar os gráficos. ::
        
            CalculoBolha.Predicao(330.0)
            
        Em seguida, pode plotado o diagrama Pxy, conforme consta em [1].  ::
        
            Graficos('Pxy',Componentes, CalculoBolha.Composicao_x1,CalculoBolha.Pressao_Ponto_Bolha, CalculoBolha.Pressao_Ponto_Orvalho,T = 330.0)
        
        =========
        Exemplo 2        
        =========
        
        Os passos deste exemplo são idênticos aos do Exemplo 1. ::
        
            from Conexao import Componente_Caracterizar, UNIQUAC, VIRIAL
            from VLE import VLE
            
            Comp1 = Componente_Caracterizar('Acetona',ConfigPsat=('Prausnitz4th',1),T=330.0)
            Comp2 = Componente_Caracterizar('Metanol',ConfigPsat=('Prausnitz4th',1),T=330.0)
                
            Componentes = [Comp1,Comp2]
            
            model_vap = VIRIAL(Componentes)
            model_liq = UNIQUAC(Componentes,330.0,1)
            
            CalculoBolha = VLE('PontoBolha',Componentes,model_liq,model_vap,z=[0.95,0.05],Temp=330.0,Pressao=1.013,estgama=None,estphi=None,estBeta = 0.5,tolAlg=1e-5,toleq=1e-4,maxiter=500,z_coordenacao=10.0)
        
        Após este passo, é acessado o método de predição, para caracterizar os gráficos. ::
        
            CalculoBolha.Predicao(330.0)
        
        Em seguida, pode plotado o diagrama xy para os componentes, conforme consta em [1].  ::
        
            Graficos('xy',Componentes, CalculoBolha.Composicao_x1,CalculoBolha.Composicao_y1,CalculoBolha.Composicao_y2,CalculoBolha.Composicao_x2,T = 330.0)  
        
        Os gráficos são gerados automaticamente.
        
        ===========
        Referências
        ===========
        
        [1] SMITH, J. M.; NESS, H. C. V.; ABBOTT, M. M. Introduction to 
        Chemical Engineering Thermodinamics. 7. ed. [S.l.]: Mc-Graw Hills,p. 254, 2004.        
        '''
        # TODO: converter para keywargs
        #==============================================================================
        #         Definição das variáveis
        #==============================================================================
        self.T  = T
        self.P_exp  = P
        self.P_incertezas = P_incertezas
        self.x_exp = x_experimentais
        self.y_exp = y_experimentais
        self.y_incertezas = y_incertezas
        
        if Diagrama == 'Pxy':
            
            self.P_x_y(VLE)
        
        elif Diagrama == 'Txy':
             
            self.T_x_y()
        
        elif Diagrama == 'xy':
            
            self.x_y()
            
    def P_x_y(self,VLE):
        '''
        Método para criação do diagrama Pxy, pressão em função das composições, vide [1].
        
        ======
        Saídas
        ======
        
        * Gera o diagrama Pxy.
        
        ===========
        Referências
        ===========
        
        [1] SMITH, J. M.; NESS, H. C. V.; ABBOTT, M. M. Introduction to 
        Chemical Engineering Thermodinamics. 7. ed. [S.l.]: Mc-Graw Hills,p. 254, 2004.

        '''
        #==============================================================================
        #         Plotagem dos pontos calculados
        #==============================================================================
        plot(VLE.Bolha.comp[0],VLE.Bolha.Pressao,'-',color ='green')
        plot(VLE.Orvalho.comp[0],VLE.Orvalho.Pressao,'-',color ='blue')
        
        if self.y_exp is not None:
            # TODO: usar quando definido
            print 'aqui'
            # ==============================================================================
            #         Plotagem dos pontos experimentais
            # ==============================================================================
            #plot(self.x_exp,self.P_exp, ':',markeredgecolor ='yellow', marker="o", markerfacecolor="w")
            #plot(self.y_exp,self.P_exp,':',markeredgecolor ='cyan', marker="o", markerfacecolor="k")
            # ==============================================================================
            #         Plotagem das incertezas
            # ==============================================================================
            errorbar(self.x_exp,self.P_exp,yerr=self.P_incertezas,fmt='o', ecolor='g',markeredgecolor ='magenta', marker="o", markerfacecolor="w")
            errorbar(self.y_exp,self.P_exp,yerr=self.P_incertezas,fmt='o', ecolor='b',markeredgecolor ='magenta', marker="o", markerfacecolor="w")
        
        xlabel(u'x,y')
        ylabel(u'Pressão /bar')
#        title('Diagrama Pxy para {:s}-{:s} a {:.1f}K'.format(self.Componentes[0].nome,self.Componentes[1].nome,self.T))
        xlim(0,1)
        legend(['Dew Pressure','Bubble Pressure'],loc='best')
        grid() # Adiciona a grade ao gráfico
        show()

    def T_x_y(self):
        '''
        Método para criação do diagrama Txy, temperatura em função das composições, vide [1].
        
        ======
        Saídas
        ======
        
        * Gera o diagrama Txy.
        
        ===========
        Referências
        ===========
        
        [1] SMITH, J. M.; NESS, H. C. V.; ABBOTT, M. M. Introduction to 
        Chemical Engineering Thermodinamics. 7. ed. [S.l.]: Mc-Graw Hills,p. 254, 2004.

        '''
        #==============================================================================
        #         Plotagem dos pontos calculados
        #==============================================================================
        plot(self.x1,self.y1,'-',color ='green')
        plot(self.x1,self.y2,'-',color ='blue')
        
        if self.y_exp != None:
            #==============================================================================
            #         Plotagem dos pontos experimentais
            #==============================================================================
            plot(self.x_exp[0],self.y_exp, ':',markeredgecolor ='yellow', marker="o", markerfacecolor="w")
            plot(self.x_exp[1],self.y_exp,':',markeredgecolor ='cyan', marker="o", markerfacecolor="k")
            #==============================================================================
            #         Plotagem das incertezas
            #==============================================================================
            errorbar(self.x_exp[0],self.y_exp,self.y_incertezas,fmt='o', ecolor='g')
            errorbar(self.x_exp[1],self.y_exp,self.y_incertezas,fmt='o', ecolor='b')
        
        xlabel(u'x,y')
        ylabel(u'Temperatura /K')
        title('Diagrama Txy para {:s}-{:s} a {:.3f} bar'.format(self.Componentes[0].nome,self.Componentes[1].nome,self.P))
        xlim(0,1)
        legend(['Bubble Temperature','Dew Temperature'],loc='best')
        grid() # Adiciona a grade ao gráfico
        show()
        
    def x_y(self):
        '''
        Método para criação do diagrama xy, composição da fase de vapor em função 
        da composicção da fase líquida,vide [1].
        
        ======
        Saídas
        ======
        
        * Gera o diagrama xy.
        
        ===========
        Referências
        ===========
        
        [1] SMITH, J. M.; NESS, H. C. V.; ABBOTT, M. M. Introduction to 
        Chemical Engineering Thermodinamics. 7. ed. [S.l.]: Mc-Graw Hills,p. 254, 2004.

        '''
        fig = figure() # Adiciona a figura
        fig.add_subplot(111)
        plot(self.x1,self.y1,'-') 
        xlabel(u'Composição de {:s} na fase líquida'.format(self.Componentes[0].nome))
        ylabel(u'Composição de {:s} na fase de vapor'.format(self.Componentes[0].nome))
        title('Diagrama xy para {:s} a {:.1f}K'.format(self.Componentes[0].nome,self.T))
        xlim(0,1) # Dimensiona a imagem do gráfico
        
        grid() # Adiciona a grade ao gráfico

        fig = figure() # Adiciona a figura
        fig.add_subplot(111)
        plot(self.x2,self.y2,'-')
        xlabel(u'Composição de {:s} na fase líquida'.format(self.Componentes[1].nome))
        ylabel(u'Composição de {:s} na fase de vapor'.format(self.Componentes[1].nome))
        title('Diagrama xy para {:s} a {:.2f}K'.format(self.Componentes[1].nome,self.T))
        xlim(0,1) # Dimensiona a imagem do gráfico
        
        grid() # Adiciona a grade ao gráfico
        show() # Mostra os gráficos
