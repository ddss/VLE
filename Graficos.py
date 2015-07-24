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

    def __init__(self,**kwargs):

        '''
        Algoritmo para criação das saídas gráficas dos cálculos envolvendo o equilíbrio líquido vapor, vide [1].
        
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
        
            >>> CalculoBolha.Predicao('temperatura',313.15)
            
        Em seguida, pode plotado o diagrama Pxy, conforme consta em [1].  ::
        
            >>> Graph = Graficos()
            >>> Graph.P_x_y(CalculoBolha,313.15)
        
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
        
            >>> CalculoBolha.Predicao('temperatura',313.15)
        
        Em seguida, pode plotado o diagrama xy para os componentes, conforme consta em [1].  ::
        
            >>> Graph = Graficos()  
            >>> Graph.x_y(CalculoBolha,313.15)
        
        Os gráficos são salvos em na pasta 'Gráficos termodinâmicos' automaticamente.
        
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

    def __validacaoArgumentosEntrada(self,keywargs):
        
        self.__keywordsEntrada = ['x_experimentais','y_experimentais','T_exp','P_exp','x_incertezas','y_incertezas','T_incertezas','P_incertezas']
        
        #==============================================================================
        # Validação se houve keywords digitadas incorretamente:
        #==============================================================================
        keyincorreta  = [key for key in keywargs.keys() if not key in self.__keywordsEntrada]
        
        if len(keyincorreta) != 0:
            raise NameError(u'keyword(s) incorretas: '+', '.join(keyincorreta)+'.'+u' Keywords disponíveis: '+', '.join(self.__keywordsEntrada)+'.')

        #==============================================================================
        # Validação se houve a presença dos pares das composiões e de suas incertezas
        #==============================================================================
        # Caso definido alguma composição experimental, é necessário definir a outra
        if self.__keywordsEntrada[0] in keywargs.keys() or self.__keywordsEntrada[1] in keywargs.keys():
            # Se as keywords 'x_experimentais','y_experimentais' não estiverem em keywargs.keys, erro.
            if not {self.__keywordsEntrada[0],self.__keywordsEntrada[1]}.issubset(keywargs.keys()):
                raise ValueError(u'Caso seja definido alguma composição experimental, faz-se necessário definir todas, tanto do líquido quanto do vapor como argumentos de entrada.')

        # Caso definido alguma incerteza experimental, é necessário definir a outra, bem como as composições
        if self.__keywordsEntrada[4] in keywargs.keys() or self.__keywordsEntrada[5] in keywargs.keys():
            # Se as keywords 'x_experimentais','y_experimentais', 'x_incertezas','y_incertezas' não estiverem em keywargs.keys, erro.
            if not {self.__keywordsEntrada[0],self.__keywordsEntrada[1],self.__keywordsEntrada[4],self.__keywordsEntrada[5]}.issubset(keywargs.keys()):
                raise ValueError(u'Caso seja definido alguma inerteza, faz-se neessário definir as composição experimentais para as fases líquida e vapor, bem como suas respectivas incertezas.')

        #==============================================================================
        # Validação dos tipos de variáveis
        #==============================================================================


        if keywargs.get(self.__keywordsEntrada[0]) is not None:
            if not isinstance(keywargs.get(self.__keywordsEntrada[0]),list):
                 raise TypeError(u'As composições experimentais da fase líquida devem ser iseridas numa LISTA.')

        if keywargs.get(self.__keywordsEntrada[1]) is not None:
            if not isinstance(keywargs.get(self.__keywordsEntrada[1]),list):
                 raise TypeError(u'As composições experimentais da fase vapor devem ser iseridas numa LISTA.')

        if keywargs.get(self.__keywordsEntrada[2]) is not None:
            if not isinstance(keywargs.get(self.__keywordsEntrada[2]),list):
                 raise TypeError(u'As temperaturas devem ser iseridas numa LISTA.')

        if keywargs.get(self.__keywordsEntrada[3]) is not None:
            if not isinstance(keywargs.get(self.__keywordsEntrada[3]),list):
                 raise TypeError(u'As pressões devem ser inseridas numa LISTA.')
                 
        if keywargs.get(self.__keywordsEntrada[4]) is not None:
            if not isinstance(keywargs.get(self.__keywordsEntrada[4]),list):
                 raise TypeError(u'As incertezas das composições da fase líquida devem ser iseridas numa LISTA.')

        if keywargs.get(self.__keywordsEntrada[5]) is not None:
            if not isinstance(keywargs.get(self.__keywordsEntrada[5]),list):
                 raise TypeError(u'As incertezas das composições da fase vapor devem ser iseridas numa LISTA.')

        if keywargs.get(self.__keywordsEntrada[6]) is not None:
            if not isinstance(keywargs.get(self.__keywordsEntrada[6]),list):
                 raise TypeError(u'As incertezas das temperaturas devem ser iseridas numa LISTA.')

        if keywargs.get(self.__keywordsEntrada[7]) is not None:
            if not isinstance(keywargs.get(self.__keywordsEntrada[7]),list):
                 raise TypeError(u'As incertezas das pressões devem ser inseridas numa LISTA.')

        #==============================================================================
        # Validação se houve a presença das incertezas sem os pontos experimentais
        #==============================================================================

        # Validação da temperatura
        if keywargs.get(self.__keywordsEntrada[2]) is None:
            if keywargs.get(self.__keywordsEntrada[6]) is not None:
                raise ValueError(u'Foi inserida informação para a incerteza da temperatura, entretanto não foram inseridos dados de temperatura.')
        # Validação da pressão
        if keywargs.get(self.__keywordsEntrada[3]) is None:
            if keywargs.get(self.__keywordsEntrada[7]) is not None:
                raise ValueError(u'Foi inserida informação para a incerteza da pressão, entretanto não foram inseridos dados de pressão.')

    def P_x_y(self,VLE,T,unidT='K'):
        '''
        Método para criação do diagrama Pxy, pressão em função das composições, vide [1].

        ========
        Entradas
        ========

        * VLE (objeto)  : Objeto VLE que executou o método Predicao
        * T (float)     : temperatura para a qual foi executado
        * unidT (string): unidade da temperatura

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
        title('Diagrama Pxy para {:s}-{:s} a {:.1f} {:s}'.format(VLE.Componente[0].nome,VLE.Componente[1].nome,T,unidT))
        xlim(0,1)
        legend([u'Pressão de orvalho',u'Pressão de bolha'],loc='best')
        grid() # Adiciona a grade ao gráfico
        fig.savefig(self.base_path+'Diagrama_P_x_y.png')
        close()

    def T_x_y(self,VLE,P,unidP='bar'):
        '''
        Método para criação do diagrama Txy, temperatura em função das composições, vide [1].
        
        ========
        Entradas
        ========

        * VLE (objeto)  : Objeto VLE que executou o método Predicao
        * P (float)     : temperatura para a qual foi executado
        * unidP (string): unidade da temperatura

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
        
        if self.x_exp is not None and self.T_exp is not None and self.y_exp is not None:
            
            errorbar(self.x_exp,self.T_exp,yerr=self.T_incertezas,xerr=self.x_incertezas,fmt='o', ecolor='g',markeredgecolor ='magenta', marker="o", markerfacecolor="w")
            errorbar(self.y_exp,self.T_exp,yerr=self.T_incertezas,xerr=self.y_incertezas,fmt='o', ecolor='b',markeredgecolor ='magenta', marker="o", markerfacecolor="w")
        
        xlabel(u'x,y')
        ylabel(u'Temperatura /K')
        title('Diagrama Txy para {:s}-{:s} a {:.3f} {:s}'.format(VLE.Componente[0].nome,VLE.Componente[1].nome,P,unidP))
        xlim(0,1)
        legend([u'Temperatura de orvalho',u'Temperatura de bolha'],loc='best')
        grid() # Adiciona a grade ao gráfico
        fig.savefig(self.base_path+'Diagrama_T_x_y.png')
        close()
        
    def x_y(self,VLE,T,unidT='K'):
        '''
        Método para criação do diagrama xy, composição da fase de vapor em função 
        da composicção da fase líquida, a uma temperatura constante, vide [1].
        
        ========
        Entradas
        ========

        * VLE (objeto)  : Objeto VLE que executou o método Predicao;
        * T (float)     : temperatura para a qual foi executado;
        * unidT (string): unidade da temperatura.
        
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
        plot(VLE.Orvalho.comp_molar[0],VLE.Bolha.comp_molar[0],'-') 
        xlabel(u'Composição de {:s} na fase líquida'.format(VLE.Componentes[0].nome))
        ylabel(u'Composição de {:s} na fase de vapor'.format(VLE.Componentes[0].nome))
        title('Diagrama xy para {:s} a {:.1f} {:s}'.format(VLE.Componentes[0].nome,T,unidT))
        xlim(0,1) # Dimensiona a imagem do gráfico        
        grid() # Adiciona a grade ao gráfico

        fig = figure() # Adiciona a figura
        fig.add_subplot(1,1,1)
        plot(VLE.Orvalho.comp_molar[1],VLE.Bolha.comp_molar[1],'-')
        xlabel(u'Composição de {:s} na fase líquida'.format(VLE.Componentes[1].nome))
        ylabel(u'Composição de {:s} na fase de vapor'.format(VLE.Componentes[1].nome))
        title('Diagrama xy para {:s} a {:.2f} {:s}'.format(VLE.Componentes[1].nome,T,unidT))
        xlim(0,1) # Dimensiona a imagem do gráfico
        grid() # Adiciona a grade ao gráfico
        
        fig.savefig(self.base_path+'Diagrama_x_y.png')
        close()