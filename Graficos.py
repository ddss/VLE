# -*- coding: utf-8 -*-
"""
Created on Wed Feb 04 11:01:10 2015

@author: labpic-14
"""

from matplotlib.pyplot import plot, xlabel, ylabel, title, xlim, legend, grid, show, figure

class Graficos:

    def __init__(self,Diagrama,T,P,Componentes,Eixo_X1,Eixo_Y1,Eixo_Y2 = None, Eixo_X2 = None):
        '''
        Algoritmo para criação das saídas gráficas dos cálculos envolvendo o equilíbrio líquido vapor, vide [1].
        
        ========
        Entradas
        ========

        * Diagrama (str): O tipo de diagrama ou gráficor que se pretende realizar;
        * T (float): A temperatura da mistura em Kelvin;
        * P (float): A pressão em bar;
        * Componentes (list): É uma lista de objetos ``Componente_Caracterizar``, vide documentação da dessa classe;
        * Eixo_X1 (array): Contém os valores do eixo x, entrada obrigratória;
        * Eixo_Y1 (array): Contém os valores do eixo y, entrada obrigratória;
        * Eixo_Y2 (array): Contém valores que podem ser implementados no eixo y, entrada opcional;
        * Eixo_X2 (array): Contém valores que podem ser implementados no eixo x, entrada opcional.
        
        =======
        Métodos
        =======
        
        * ``P_x_y``:
            * Método para criação do diagrama Pxy, pressão em função das composições, vide [1].
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
        
            Graficos('Pxy',330.0,Componentes, CalculoBolha.Composicao_x1,CalculoBolha.Pressao_Ponto_Bolha, CalculoBolha.Pressao_Ponto_Orvalho)
        
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
        
            Graficos('xy',330.0,Componentes, CalculoBolha.Composicao_x1,CalculoBolha.Composicao_y1,CalculoBolha.Composicao_y2,CalculoBolha.Composicao_x2)  
        
        Os gráficos são gerados automaticamente.
        
        ===========
        Referências
        ===========
        
        [1] SMITH, J. M.; NESS, H. C. V.; ABBOTT, M. M. Introduction to 
        Chemical Engineering Thermodinamics. 7. ed. [S.l.]: Mc-Graw Hills,p. 254, 2004.        
        '''
        self.T  = T
        self.P  = P
        self.x1 = Eixo_X1
        self.x2 = Eixo_X2
        self.y1 = Eixo_Y1
        self.y2 = Eixo_Y2        
        self.Componentes   = Componentes
        
        if Diagrama == 'Pxy':
            
            self.P_x_y()
        
        elif Diagrama == 'Txy':
            
            self.T_x_y()
        
        elif Diagrama == 'xy':
            
            self.x_y()
            
    def P_x_y(self):
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
        plot(self.x1,self.y1,'-')
        plot(self.x1,self.y2,'-')
        
        xlabel(u'x,y')
        ylabel(u'Pressão(bar)')
        title('Diagrama Pxy para {:s}-{:s} a {:.1f}K'.format(self.Componentes[0].nome,self.Componentes[1].nome,self.T))
        xlim(0,1)
        legend(['Bubble Temperature','Dew Temperature'],loc='best')
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
        plot(self.x1,self.y1,'-')
        plot(self.x1,self.y2,'-')
        
        xlabel(u'x,y')
        ylabel(u'Pressão(bar)')
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
