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
# Construção de gráficos
from matplotlib.pyplot import figure, close, plot, xlabel, ylabel, title, xlim, legend, grid, show, figure, errorbar

# Pacotes do sistema operacional
from os import getcwd, sep, mkdir, path

class Graficos:

    def __init__(self,Diagrama,VLE,**kwargs):

        '''
        Algoritmo para criação das saídas gráficas dos cálculos envolvendo o equilíbrio líquido vapor, vide [1].
        
        =====================
        Entradas obrigatórias
        =====================

        * Diagrama (str): O tipo de diagrama ou gráficor que se pretende realizar;
        * VLE (objeto): objeto VLE após execuação do método Predicao executado
        
        ==============================
        Keywords (Entradas opcionais):
        ==============================
        
        * x_experimentais (list): Valores experimentais das composições da fase líquida;
        * y_experimentais (list): Valores experimentais das composições da fase de vapor;
        * T_exp (float or list): A(s) temperatura(s) da mistura em Kelvin, será float caso for constante e list caso for variável;
        * P_exp (float or list): As pressões em bar, será float caso for constante e list caso for variável;
        * x_incertezas (list): Valores das incertezas das composições da fase líquida;
        * y_incertezas (list): Valores das incertezas das composições da fase de vapor;
        * T_incertezas (float): Valores das incertezas das temperaturas da mistura em Kelvin;
        * P_incertezas (float): Valores das incertezas das pressões em bar.
        
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
        
            >>> from Conexao import Componente_Caracterizar, UNIQUAC, VIRIAL
            >>> from VLE import VLE
            
            >>> Comp1 = Componente_Caracterizar('Acetona',ConfigPsat=('Prausnitz4th',1),T=330.0)
            >>> Comp2 = Componente_Caracterizar('Metanol',ConfigPsat=('Prausnitz4th',1),T=330.0)
                
            >>> Componentes = [Comp1,Comp2]
            
            >>> model_vap = VIRIAL(Componentes)
            >>> model_liq = UNIQUAC(Componentes,330.0,1)
            
            >>> CalculoBolha = VLE('PontoBolha_P',Componentes,model_liq,model_vap,z=[0.95,0.05],Temp=330.0,Pressao=1.013,estgama=None,estphi=None,estBeta = 0.5,tolAlg=1e-5,toleq=1e-4,maxiter=500,z_coordenacao=10.0)
        
        Após este passo, é acessado o método de predição, para caracterizar os gráficos. ::
        
            >>> CalculoBolha.Predicao(330.0)
            
        Em seguida, pode plotado o diagrama Pxy, conforme consta em [1].  ::
        
            >>> Graficos('Pxy',Componentes, CalculoBolha.Composicao_x1,CalculoBolha.Pressao_Ponto_Bolha, CalculoBolha.Pressao_Ponto_Orvalho,T = 330.0)
        
        =========
        Exemplo 2        
        =========
        
        Os passos deste exemplo são idênticos aos do Exemplo 1. ::
        
            >>> from Conexao import Componente_Caracterizar, UNIQUAC, VIRIAL
            >>> from VLE import VLE
            
            >>> Comp1 = Componente_Caracterizar('Acetona',ConfigPsat=('Prausnitz4th',1),T=330.0)
            >>> Comp2 = Componente_Caracterizar('Metanol',ConfigPsat=('Prausnitz4th',1),T=330.0)
                
            >>> Componentes = [Comp1,Comp2]
            
            >>> model_vap = VIRIAL(Componentes)
            >>> model_liq = UNIQUAC(Componentes,330.0,1)
            
            >>> CalculoBolha = VLE('PontoBolha',Componentes,model_liq,model_vap,z=[0.95,0.05],Temp=330.0,Pressao=1.013,estgama=None,estphi=None,estBeta = 0.5,tolAlg=1e-5,toleq=1e-4,maxiter=500,z_coordenacao=10.0)
        
        Após este passo, é acessado o método de predição, para caracterizar os gráficos. ::
        
            >>> CalculoBolha.Predicao(330.0)
        
        Em seguida, pode plotado o diagrama xy para os componentes, conforme consta em [1].  ::
        
            >>> Graficos('xy',Componentes, CalculoBolha.Composicao_x1,CalculoBolha.Composicao_y1,CalculoBolha.Composicao_y2,CalculoBolha.Composicao_x2,T = 330.0)  
        
        Os gráficos são gerados automaticamente.
        
        ===========
        Referências
        ===========
        
        [1] SMITH, J. M.; NESS, H. C. V.; ABBOTT, M. M. Introduction to 
        Chemical Engineering Thermodinamics. 7. ed. [S.l.]: Mc-Graw Hills,p. 254, 2004.        
        '''
        
        self.__validacaoArgumentosEntrada(kwargs)
        #==============================================================================
        #         Definição das variáveis
        #==============================================================================
        self.x_exp = kwargs.get(self.__keywordsEntrada[0])
        self.y_exp = kwargs.get(self.__keywordsEntrada[1])
        self.T_exp = kwargs.get(self.__keywordsEntrada[2])
        self.P_exp = kwargs.get(self.__keywordsEntrada[3])
        self.x_incertezas = kwargs.get(self.__keywordsEntrada[4])
        self.y_incertezas = kwargs.get(self.__keywordsEntrada[5])
        self.T_incertezas = kwargs.get(self.__keywordsEntrada[6])
        self.P_incertezas = kwargs.get(self.__keywordsEntrada[7])
        
        # Cria caminho para salvar os gráficos
        self.base_path = getcwd() + sep +'Graficos termodinamicos'+sep
        # Cria o dirétorio do caminho  
        if path.exists(self.base_path) is False:
            mkdir(self.base_path)
        
        if Diagrama == 'Pxy':
            
            self.P_x_y(VLE)
        
        elif Diagrama == 'Txy':
             
            self.T_x_y(VLE)
        
        elif Diagrama == 'xy':
            
            self.x_y()
            
    def __validacaoArgumentosEntrada(self,keywargs):
        
        self.__keywordsEntrada = ['x_experimentais','y_experimentais','T_exp','P_exp','x_incertezas','y_incertezas','T_incertezas','P_incertezas']
        
        #==============================================================================
        #         # Validação se houve keywords digitadas incorretamente:
        #==============================================================================
        keyincorreta  = [key for key in keywargs.keys() if not key in self.__keywordsEntrada]
        
        if len(keyincorreta) != 0:
            raise NameError(u'keyword(s) incorretas: '+', '.join(keyincorreta)+'.'+u' Keywords disponíveis: '+', '.join(self.__keywordsEntrada)+'.')

        #==============================================================================
        #  Validação se houve a presença das incertezas sem os pontos experimentais
        #==============================================================================
        # Validação da composição da fase líquida
        if keywargs.get(self.__keywordsEntrada[0]) is None:
            if keywargs.get(self.__keywordsEntrada[4]) is not None:
                raise ValueError(u'Faz-se necessário os pontos experimentais como argumentos de entrada.')
        # Validação da composição da fase de vapor
        if keywargs.get(self.__keywordsEntrada[1]) is None:
            if keywargs.get(self.__keywordsEntrada[5]) is not None:
                raise ValueError(u'Faz-se necessário os pontos experimentais como argumentos de entrada.')
        # Validação da temperatura
        if keywargs.get(self.__keywordsEntrada[2]) is None:
            if keywargs.get(self.__keywordsEntrada[6]) is not None:
                raise ValueError(u'Faz-se necessário os pontos experimentais como argumentos de entrada.')
        # Validação da pressão
        if keywargs.get(self.__keywordsEntrada[3]) is None:
            if keywargs.get(self.__keywordsEntrada[7]) is not None:
                raise ValueError(u'Faz-se necessário os pontos experimentais como argumentos de entrada.')

        #==============================================================================
        #  Validação se houve a presença dos pares das composiões e de suas incertezas
        #==============================================================================
        # Validação da composição da fase líquida
        if keywargs.get(self.__keywordsEntrada[0]) is None: # Caso não for inserido a composição da fase líquida
            if keywargs.get(self.__keywordsEntrada[1]) is not None: # E for inserido a composição da fase de vapor
                raise ValueError(u'Faz-se necessário que todas as composições sejam argumentos de entrada.')
        else: # Caso for inserido a composição da fase líquida
            if keywargs.get(self.__keywordsEntrada[1]) is None: # E não for inserido a composição da fase de vapor
                raise ValueError(u'Faz-se necessário que todas as composições sejam argumentos de entrada.')
            else: # Caso também for inserido a composição da fase líquida
                # Validação das incertezas, caso for inserido as duas composições
                if keywargs.get(self.__keywordsEntrada[4]) is None: # Caso não for inserida a incerteza da composição da fase líquida
                    if keywargs.get(self.__keywordsEntrada[5]) is not None: # E for inserida a incerteza da composição da fase de vapor
                        raise ValueError(u'Faz-se necessário que todas incertezas das composições sejam argumentos de entrada.')
                else: # Caso for inserida a incerteza da composição da fase líquida
                    if keywargs.get(self.__keywordsEntrada[5]) is None: # E não for inserida a incerteza da composição da fase de vapor
                        raise ValueError(u'Faz-se necessário que todas incertezas das composições sejam argumentos de entrada.')
                                
                
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
        fig = figure()
        fig.add_subplot(1,1,1) 
        #==============================================================================
        #         Plotagem dos pontos calculados
        #==============================================================================
        plot(VLE.Bolha.comp_molar[0],VLE.Bolha.Pressao,'-',color ='green')
        plot(VLE.Orvalho.comp_molar[0],VLE.Orvalho.Pressao,'-',color ='blue')

        #==============================================================================
        #         Plotagem dos pontos experimentais
        #==============================================================================
        
        if self.x_exp is not None and self.P_exp is not None and self.y_exp is not None:

            errorbar(self.x_exp,self.P_exp,yerr=self.P_incertezas,xerr=self.x_incertezas,fmt='o', ecolor='g',markeredgecolor ='magenta', marker="o", markerfacecolor="w")
            errorbar(self.y_exp,self.P_exp,yerr=self.P_incertezas,xerr=self.y_incertezas,fmt='o', ecolor='b',markeredgecolor ='magenta', marker="o", markerfacecolor="w")
                                                                    
        xlabel(u'x,y')
        ylabel(u'Pressão /bar')
        title('Diagrama Pxy para {:s}-{:s} a {:.1f}K'.format(VLE.Componente[0].nome,VLE.Componente[1].nome,self.T_exp))
        xlim(0,1)
        legend(['Dew Pressure','Bubble Pressure'],loc='best')
        grid() # Adiciona a grade ao gráfico
        show()
        fig.savefig(self.base_path+'Diagrama_P_x_y.png')
        close()

    def T_x_y(self,VLE):
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
        fig = figure()
        fig.add_subplot(1,1,1) 
        #==============================================================================
        #         Plotagem dos pontos calculados
        #==============================================================================
        plot(VLE.Bolha.comp_molar[0],VLE.Bolha.Temp,'-',color ='green')
        plot(VLE.Orvalho.comp_molar[0],VLE.Orvalho.Temp,'-',color ='blue')

        #==============================================================================
        #         Plotagem dos pontos experimentais
        #==============================================================================
        
        if self.x_exp and self.T_exp is not None:
            plot(self.x_exp,self.T_exp,ls ='None',markeredgecolor ='magenta', marker="o", markerfacecolor="w")
            plot(self.y_exp,self.T_exp,ls ='None',markeredgecolor ='magenta', marker="o", markerfacecolor="w")
            
            if self.T_incertezas is not None:
                # ==============================================================================
                #         Plotagem das incertezas de P
                # ==============================================================================
                errorbar(self.x_exp,self.T_exp,yerr=self.T_incertezas,fmt='o', ecolor='g',markeredgecolor ='magenta', marker="o", markerfacecolor="w")
                errorbar(self.y_exp,self.T_exp,yerr=self.T_incertezas,fmt='o', ecolor='g',markeredgecolor ='magenta', marker="o", markerfacecolor="w")
                
                if self.x_incertezas is not None:
                    # ==============================================================================
                    #         Plotagem das incertezas das composições
                    # ==============================================================================
                    errorbar(self.x_exp,self.T_exp,yerr=self.T_incertezas,xerr=self.x_incertezas,fmt='o', ecolor='g',markeredgecolor ='magenta', marker="o", markerfacecolor="w")
                    errorbar(self.y_exp,self.T_exp,yerr=self.T_incertezas,xerr=self.y_incertezas,fmt='o', ecolor='b',markeredgecolor ='magenta', marker="o", markerfacecolor="w")
        
        xlabel(u'x,y')
        ylabel(u'Temperatura /K')
        title('Diagrama Txy para {:s}-{:s} a {:.3f} bar'.format(VLE.Componente[0].nome,VLE.Componente[1].nome,self.P_exp))
        xlim(0,1)
        legend(['Bubble Temperature','Dew Temperature'],loc='best')
        grid() # Adiciona a grade ao gráfico
        show()
        fig.savefig(self.base_path+'Diagrama_T_x_y.png')
        close()
        
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
        fig.add_subplot(1,1,1)
        plot(self.x1,self.y1,'-') 
        xlabel(u'Composição de {:s} na fase líquida'.format(self.Componentes[0].nome))
        ylabel(u'Composição de {:s} na fase de vapor'.format(self.Componentes[0].nome))
        title('Diagrama xy para {:s} a {:.1f}K'.format(self.Componentes[0].nome,self.T))
        xlim(0,1) # Dimensiona a imagem do gráfico
        
        grid() # Adiciona a grade ao gráfico

        fig = figure() # Adiciona a figura
        fig.add_subplot(1,1,1)
        plot(self.x2,self.y2,'-')
        xlabel(u'Composição de {:s} na fase líquida'.format(self.Componentes[1].nome))
        ylabel(u'Composição de {:s} na fase de vapor'.format(self.Componentes[1].nome))
        title('Diagrama xy para {:s} a {:.2f}K'.format(self.Componentes[1].nome,self.T))
        xlim(0,1) # Dimensiona a imagem do gráfico
        
        grid() # Adiciona a grade ao gráfico
        show() # Mostra os gráficos
        close()