# -*- coding: utf-8 -*-

from sqlite3 import connect
from warnings import warn
from scipy import exp, log
from scipy.optimize import root
from numpy import zeros

class Componente_Caracterizar:
    
    def __init__(self,Componente,T,ConfigPsat=('Prausnitz4th',None)):
        u'''
        Objeto para definir as propriedades de componentes puros.
        
        ========
        Entradas
        ========
        
        * Componente ('string')     : Nome do componente conforme consta no Banco de dados
        
            * A lista de componentes disponíveis no Banco de dados pode ser acessada da seguinte forma: ::
                
                >>> Componente_Caracterizar(None) 
        
        * T: Temperatura em Kelvin;
        * P: A pressão dada em bar;
        * ConfigPsat ('tuple' or 'list'): São as opções para o cálculo da pressão de vapor. Primeiro argumento são as referências disponíveis ('Prausnitz4th') e o segundo o número da equação conforme consta no Banco de dados.
        
        =========
        Atributos
        =========
        
        As propriedades estão disponíveis na forma de atributos:
                  
            * ``Tc``: Temperatura crítica em Kelvin;
            * ``Pc``: Pressão crítica em bar;
            * ``w``: Fator acêntrico (Adimensional);
            * ``MM``: Massa molar em g.mol-1;
            * ``radius_giration``: Raio médio de giração (Adimensional) ;
            * ``Zc``: Fator de compressibilidade crítico (Adimensional) ;
            * ``Vc``: Volume molar crítico (Adimensional) ;
            * ``dipole_moment``: Momento dipolo (Adimensional)
            * ``r``, ``q`` & ``ql``: Parâmetros do modelo UNIQUAC (Adimensionais);
            * ``VPA``, ``VPB`` , ``VPC`` & ``VPD``: Parâmetros para o cálculo de Psat, vide [1];
            * ``TminPsat``: Temperatura mínima para a aplicaçao da fórmula do cálculo de Psat em Kelvin;
            * ``TmaxPsat``: Temperatura máxima para a aplicaçao da fórmula do cálculo de Psat em Kelvin;
            * ``ro``: Densidade do líquido em g/cm3;
            * ``Td``: Temperatura da densidade do líquido em Kelvin;
            * ``grupo_funcional``: Grupo funcional do componente.        
            
        =======
        Métodos     
        =======
        
        Os métodos dispníveis desta classe são:
            * ``Pvap_Prausnitz_4th``:
                * Método para o cálculo da pressão de vapor. Vide documentação do método.            
            * ``Propriedade``:        
                * Método para a busca no Banco de dados das propriedades puras dos componentes. Vide documentação do método.
        
        =======
        Exemplo
        =======
        
        * A classe é acessada de forma usual, como pode ser visto no caso a seguir: ::
        
            Componente = Componente_Caracterizar('Metano',ConfigPsat=('Prausnitz4th',1),T=100.0)
        
        ===========
        Referências
        ===========
        
        [1] REID, R.C.; PRAUSNITZ, J.M.; POLING, B.E. The properties of Gases and Liquids, 4th edition, McGraw-Hill, 1987.
        '''
        
        #==============================================================================
        #         CONEXÃO COM O BANCO DE DADOS
        #==============================================================================
        conector     = connect('THERMO_DATA_BANK_EXEMPLO.db') # Conecta a rotina ao banco de dados
        self.__cursor  = conector.cursor()                      # Permite a navegação no banco de dados a partir da rotina

        #==============================================================================
        #         LISTAGEM DE MÉTODOS DISPONÍVEIS
        #==============================================================================
        self.__lista_EqPsat = ['Prausnitz4th'] # Lista dos métodos utilizados para o cálculo da pressão de vapor
    
        #==============================================================================
        #         LISTA DE COMPONENTES - CASO SOLICITADO
        #==============================================================================
        # Mostrar lista de componentes para o usuário, caso solicitado
        if Componente == None:
            
            print u'Os seguintes componentes estão disponíveis no banco de dados: '+', '.join(self.lista_componentes())+'.'
        
        else:
            #==============================================================================
            #             VALIDAÇÃO         
            #==============================================================================
            # Nome dos componentes
            self.Validacao_Nome(Componente)  # Verificar se o nome do componente está no banco
            
            #==============================================================================
            #             CONSULTA AO BANCO E DEFINIÇÃO DE VARIÁVEIS
            #==============================================================================
            self.nome = Componente # Nome do componente
            
            # ID da variável
            self.Busca_ID() # Busca do ID do componente e cria o atributo ID.
    
            # GRUPO FUNCIONAL
            self.Busca_grupo() # Busca o grupo funcional e cria atributo
            
            #==============================================================================
            #             DEFINIÇÃO DE OUTRAS VARIÁVEIS
            #==============================================================================
            self.T    = T # Temperatura
            
            # PRESSÃO DE VAPOR
            self.eqPsat   = ConfigPsat[0] # Nome da equação para cálculo da pressão de vapor
            self.nEqPsat  = ConfigPsat[1] # Número da equação de Psat (Vide referência(?))
            
            self.Validacao_e_Default_de_EqPsat() # Validação dos atributos de pressão de vapor
    
            # PROPRIEDADES DO COMPONENTE PURO
            self.Propriedade() # Criação dos atributos contendo as propriedades do componentes puro
        
        #==============================================================================
        #         ENCERRAR CONEXÃO COM O BANCO
        #==============================================================================
        conector.close()
        
    def lista_componentes(self):
        u'''
        Método para gerar uma lista contendo os nomes dos componentes disponíveis no Banco de dados.
        
        ======
        Saídas
        ======
        
        * Retorna uma lista contendo os nomes dos componentes.
        '''        
        
        self.__cursor.execute('SELECT Nome FROM Componentes')  # Seleciona a coluna Nome da tabela Componentes do banco de dados
        row = self.__cursor.fetchall()                         # Retorna a coluna selecionada em forma de lista de tuplas
        return [i[0] for i in row]                             # Transforma a lista de tuplas em uma lista com o contéudo das tuplas (os nomes dos componentes)
        
    def Validacao_Nome(self,Nome):
        u'''
        Método para validar a entrada ``Componente``. Comparando esta
        entrada com os nomes obtidos pelo método ``lista_componentes``.
        Se a entrada de ``Componente`` não passam pela validação, há a emissão
        de um erro como saída deste método. O erro emitido interrompe a execuçao do programa.
        '''
        # Validação do nome do componente
        if Nome not in self.lista_componentes():
            raise NameError(u'O nome do componente não consta no Banco de dados. Os seguintes componentes estão disponíveis no banco de dados: '+', '.join(self.lista_componentes())+'.')  # Emite um erro com a mensagem inserida no método

    def Busca_ID(self):
        u'''
        Método para busca do ID, chave primária, da tabela Componentes no Banco de dados.
        
        ======
        Saídas
        ======
        
        * Gera o atributo ID em forma de número inteiro
        '''
        
        self.__cursor.execute('SELECT ID FROM Componentes WHERE  Nome=?',(self.nome,)) # Busca do ID do componente na tabela Componentes no banco de dados
        row = self.__cursor.fetchall()                                                 # Retorna a linha contendo o ID em forma de lista de tupla
        self.ID = row[0][0]                                                            # Cria o atributo ID

    def Busca_grupo(self):
        u'''
        Método para busca do grupo funcional do componente.
        
        ======
        Saídas
        ======
        
        * O método gera o atributo ``grupo_funcional``. Este é uma ``string`` representando o grupo funcional do componente requerido.
        '''

        #==============================================================================
        #         BUSCA ID_GRUPO NA TABELA COMPONENTES
        #==============================================================================
        self.__cursor.execute('SELECT ID_grupo FROM Componentes WHERE ID=?',(self.ID,))        
        row = self.__cursor.fetchall()        
        ID_grupo = row[0][0]
        
        #==============================================================================
        #         BUSCA NOME DO GRUPO PELA ID_GRUPO NA TABELA GRUPO
        #==============================================================================
        self.__cursor.execute('SELECT Nome FROM Grupo WHERE ID=?',(ID_grupo,))        
        row = self.__cursor.fetchall()
        self.grupo_funcional = row[0][0]
        

    def Busca_FormaEqPsat(self):        
        u'''
        Método para buscar a forma da equação para o cálculo de  Psat no Banco de dados.
                
        ======
        Saídas
        ======
        
        * Este método retorna uma lista contendo as formas de equações, que são números inteiros, para o componente requerido.
        '''
        
        if self.eqPsat == self.__lista_EqPsat[0]: 
            
            self.__cursor.execute('SELECT ID_forma FROM Parametros_Psat_Prausnitz_4th_edition WHERE ID_componente=?',(self.ID,)) # Busca as possíveis formas de equação do cálculo da pressão de vapor
            row = self.__cursor.fetchall()         # Cria uma lista de tuplas contendo as formas de equações
            
        return [i[0] for i in row]                 # Transforma a lista de tuplas em uma lista com o contéudo das tuplas (Os marcadores das formas de equações)

    def Validacao_e_Default_de_EqPsat(self):
        u'''
        Este método valida as entradas de ``ConfigPsat``. Assim são 
        validadas as entradas do método de cálculo de Psat e as formas de equações dos métodos. 
        Se as entradas de ``ConfigPsat`` não passam pela validação, há a emissão de um erro como 
        saída deste método. O erro emitido interrompe a execuçao do programa.
        '''
        #==============================================================================
        #         VALIDAÇÃO DO NOME DO MÉTODO DE CÁLCULO DE PRESSÃO DE VAPOR
        #==============================================================================
        
        # Caso o método inserido não constar na lista de métodos disponíveis
        if self.eqPsat not in self.__lista_EqPsat:                      
            raise NameError(u'O método de cálculo de Psat não consta no Banco de dados. Métodos disponíveis: '+', '.join(self.__lista_EqPsat)+'.')
        
        # Caso o método inserido conste na lista de métodos disponíveis:
        if self.eqPsat == self.__lista_EqPsat[0]:
            dadosbanco = self.Busca_FormaEqPsat()   # Busca as formas de equações disponíveis do banco de dados

            if self.nEqPsat != None:                # Caso a forma da equação seja inserida pelo usuário
                if self.nEqPsat not in dadosbanco:  # Caso a forma da equação inserida não conste no banco
                    raise ValueError(u'Foi escolhida a equação %s para calcular a pressão de saturação. Contudo foi inserido uma forma de equação que não consta no Banco de dados (Vide documentação do Banco de dados). Formas de equações disponíveis: '%(self.eqPsat,)+', '.join(str(forma) for forma in dadosbanco)+'.') 
            
            if self.nEqPsat == None:                # Caso a forma da equação não seja inserida
                if len(dadosbanco) == 1:            # Caso haja apenas uma forma de equação, dados banco tem tamanho 1
                    self.nEqPsat = dadosbanco[0]    # Como não seja inserida a forma de equação, o programa usará a única forma disponível
                else: # Caso haja mais de uma forma de equação disponível no banco de dados
                    raise ValueError(u'Foi escolhida a equação %s para calcular a pressão de saturação. No entanto, é necessário informar expressamente a forma da equaçao (Vide documentação do Banco de dados), visto que há diferentes opções disponíveis: '%(self.eqPsat,)+', '.join(str(forma) for forma in dadosbanco)+'.') 
        
    def warnings(self):
        u'''        
        Método para verificar se a temperatura inserida está dentro ou não da 
        faixa de aplicação das fórmulas do calculo da pressão de vapor. 
        Caso não esteja dentro da faixa de aplicabilidade, o programa gerará um aviso. 
        Contudo a execução não será interrompida.
        '''
        
        #==============================================================================
        #         Verificação se a temperatura inserida está dentro ou não da faixa de aplicabilidade das fórmulas do cálculo de Psat
        #==============================================================================
        if self.T < self.__TminPsat or self.T > self.__TmaxPsat:
            warn(u'A temperatura especificada está fora da faixa de aplicabilidade da equaçao de Psat. A temperatura de pertencer ao intervalo: (%f, '%self.__TminPsat+'%f).'%self.__TmaxPsat) # Emite um aviso. Contudo o programa continua a rodar
 
        #==============================================================================
        #        Validar se a temperatura especificada é menor do que Tc, se self.nEqPsat = 1.
        #==============================================================================
        if self.T > self.__TmaxPsat:
            # Sendo a equação de Psat Prausnitz4th
            if self.eqPsat == self.__lista_EqPsat[0]:  
                # sendo a primeira forma disponível
                if self.nEqPsat == 1:
                    if self.T > self.Tc:
                        raise ValueError(u'Para a equação de pressão de vapor escolhida para o método de cálculo: %s, é necessário que a temperatura esteja abaixo da temperatura crítica,'%(self.eqPsat,)+' Tc = %f.'%self.Tc)


    def Pvap_Prausnitz_4th(self,VPA,VPB,VPC,VPD,T,nEq,Tc=None,Pc=None,Pvp_ini=101325,tol=1e-10):
        u'''
        Método para cálculo da pressão de vapor de componentes puros, conforme [1].
        
        ========
        Entradas
        ========
        
        * VPA (float): Parâmetro VPA;
        * VPB (float): Parâmetro VPB;
        * VPC (float): Parâmetro VPC;
        * VPD (float): Parâmetro VPD;
        * T (float): Temperatura em Kelvin;
            
            * O valor de T inserido deve ser menor do que o valor de Tc.
            
        * nEq (int): Número de identificação da forma (ou tipo) de equação;
        * Tc (float): Temperatura crítica em Kelvin; 
        * Pc (float): Pressão crítica em bar
        * Pvp_ini: Estimativa inicial para a pressão de vapor, quando nEq = 2 (Equação implícta).
        * tol: teolerância para o cálculo da raiz da equação nEq = 2 (Equação implícta)
        
        ================
        Valores default 
        ================        
        
        Valores utilizados apenas quando nEq = 2.
        
        * Pvp_ini = 101325 bar
        * tol     = 1e-10
        
        ======
        Saídas
        ======
        
        * Retorna a pressão de vapor em bar
        
        ===========
        Referências
        ===========
        
        [1] REID, R.C.; PRAUSNITZ, J.M.; POLING, B.E. The properties of Gases and Liquids, 4th edition, McGraw-Hill, 1987.
        '''
    
        # Equação implícita para cálculo de Psat para nEq = 2
        def Eq2(Pvp,VPA,VPB,VPC,VPD,T): 
            Res = VPA - VPB/T + VPC*log(T) + VPD*Pvp/(T**2.0)-log(Pvp) # Vide [1]
            return Res
    
        if self.nEqPsat == 1: # Cálculo de Psat quando nEq = 1
            Pc  = Pc # Bar
            x   = 1 - T/Tc
            Pvp = exp((VPA*x+VPB*(x**1.5)+VPC*(x**3.0)+VPD*(x**6.0))/(1.0-x))*Pc # Vide [1]
            # P.s.: Caso o T inserido seja maior que o valor de Tc, haverá um erro.
            
        elif self.nEqPsat == 2:  # Cálculo de Psat quando nEq = 2
            Pvp_ini = Pvp_ini # Bar, Estimativa inicial
            Resul   = root(Eq2,Pvp_ini,args=(VPA,VPB,VPC,VPD,T),tol=tol) # Determinação das raízes da equação implícia
            return Resul.x # Retorno
    
        elif self.nEqPsat == 3: # Cálculo de Psat quando nEq = 3
           Pvp = exp(VPA-VPB/(T+VPC)) # Vide [1]
    
        # Todos os Pvp são dados em bar
        return Pvp
        
    def Tsat_Prausnitz_4th(self,VPA,VPB,VPC,VPD,P,nEq,Tc=None,Pc=None,Tsat_ini=500,tol=1e-10):
        u'''
        Método para cálculo da temperatura de componentes puros, conforme [1].
        
        ========
        Entradas
        ========
        
        * VPA (float): Parâmetro VPA;
        * VPB (float): Parâmetro VPB;
        * VPC (float): Parâmetro VPC;
        * VPD (float): Parâmetro VPD;
        * P (float): Pressão em bar;
        * nEq (int): Número de identificação da forma (ou tipo) de equação;
        * Tc (float): Temperatura crítica em Kelvin; 
        * Pc (float): Pressão crítica em bar
        * Pvp_ini: Estimativa inicial para a pressão de vapor, quando nEq = 2 (Equação implícta).
        * tol: teolerância para o cálculo da raiz da equação nEq = 2 (Equação implícta)
        
        ================
        Valores default 
        ================        
        
        * Tsat_ini = 500 K
        * tol     = 1e-10
        
        ======
        Saídas
        ======
        
        * Retorna a temperatura de saturação em K
        
        ===========
        Referências
        ===========
        
        [1] REID, R.C.; PRAUSNITZ, J.M.; POLING, B.E. The properties of Gases and Liquids, 4th edition, McGraw-Hill, 1987.
        '''
    
        # Equações implícitas
        def Eq1(T,Pc,Tc,VPA,VPB,VPC,VPD,P):
            Pc  = Pc # Bar
            x   = 1 - T/Tc
            Res = exp((VPA*x+VPB*(x**1.5)+VPC*(x**3.0)+VPD*(x**6.0))/(1.0-x))*Pc - P
            return Res
        
        def Eq2(T,VPA,VPB,VPC,VPD,P): 
            Res = VPA - VPB/T + VPC*log(T) + VPD*P/(T**2.0)-log(P) # Vide [1]
            return Res
        
        def Eq3(T,VPA,VPB,VPC,P):
            Res = exp(VPA-VPB/(T+VPC)) - P # Vide [1]
            return Res
            
        if self.nEqPsat == 1: # Cálculo de Psat quando nEq = 1
            Tsat_ini = Tsat_ini # K, Estimativa inicial
            Resul   = root(Eq1,Tsat_ini,args=(Pc,Tc,VPA,VPB,VPC,VPD,P),tol=tol) # Determinação das raízes da equação implícia
            return Resul.x[0] # Retorno
            
        elif self.nEqPsat == 2:  # Cálculo de Psat quando nEq = 2
            Tsat_ini = Tsat_ini # K, Estimativa inicial
            Resul   = root(Eq2,Tsat_ini,args=(VPA,VPB,VPC,VPD,P),tol=tol) # Determinação das raízes da equação implícia
            return Resul.x[0] # Retorno
    
        elif self.nEqPsat == 3: # Cálculo de Psat quando nEq = 3
            Tsat_ini = Tsat_ini # K, Estimativa inicial
            Resul   = root(Eq3,Tsat_ini,args=(VPA,VPB,VPC,P),tol=tol) # Determinação das raízes da equação implícia
            return Resul.x[0] # Retorno
    
        # Todos os Tsat são em Kelvin


    def Propriedade(self):
        u'''
        Algoritmo para busca das propriedades dos componentes puros.
        
        ======        
        Saídas
        ======
        
        * Todas as propriedades que sao buscadas retornam em forma de atributos do tipo float. As propriedades que são selecionadas por este algoritmo são:
            
            * ``Tc``: Temperatura crítica em Kelvin;
            * ``Pc``: Pressão crítica em bar;
            * ``w``: Fator acêntrico (Adimensional);
            * ``MM``: Massa molar em g.mol-1;
            * ``radius_giration``: Raio médio de giração (Adimensional) ;
            * ``Zc``: Fator de compressibilidade crítico (Adimensional) ;
            * ``Vc``: Volume molar crítico (Adimensional) ;
            * ``dipole_moment``: Momento dipolo (Adimensional)
            * ``r``, ``q`` & ``ql``: Parâmetros do modelo UNIQUAC (Adimensionais);
            * ``VPA``, ``VPB`` , ``VPC`` & ``VPD``: Parâmetros para o cálculo de Psat, vide [1];
            * ``__TminPsat``: Temperatura mínima para a aplicaçao da fórmula do cálculo de Psat em Kelvin;
            * ``__TmaxPsat``: Temperatura máxima para a aplicaçao da fórmula do cálculo de Psat em Kelvin;
            * ``ro``: Densidade do líquido em g/cm3;
            * ``Td``: Temperatura da densidade do líquido em Kelvin. 
            
        ===========
        Referências
        ===========
        
        [1] REID, R.C.; PRAUSNITZ, J.M.; POLING, B.E. The properties of Gases and Liquids, 4th edition, McGraw-Hill, 1987.
        '''
        
        self.__cursor.execute('SELECT * FROM Propriedades_puras WHERE  ID_componente=?',(self.ID,))
        row = self.__cursor.fetchall() # linha do banco de dados para o ID

        #==============================================================================
        # PROPRIEDADES DA SUBSTÂNCIA (Tc,Pc,Fator acêtrico(w),Massa molar (MM),radius_giration,Fator de Compressibilidade crítico(Zc)
        # ,Volume molar crítico(Vc),dipole_moment,Parâmetros do UNIQUAC(r,q & ql),d & Td)
        #==============================================================================
        
        self.Tc              = row[0][2]  # Temperatura crítica / K
        self.Pc              = row[0][3]  # Pressão crítica / bar
        self.w               = row[0][4]  # Fator acêtrico / admnesional
        self.MM              = row[0][5]  # Massa Molar
        self.radius_giration = row[0][6]  # mean radius of gyration
        self.Zc              = row[0][7]  # Fator de Compressibilidade crítico
        self.Vc              = row[0][8]  # Volume molar crítico
        self.dipole_moment   = row[0][9]  # Momento dipolo
        self.r               = row[0][10] # Parâmetro r do UNIQUAC
        self.q               = row[0][11] # Parâmetro q do UNIQUAC
        self.ql              = row[0][12] # Parâmetro ql do UNIQUAC (Correção para álcool)
        self.ro              = row[0][13] # Densidade líquido
        self.Td              = row[0][14] # Temp_Densidade líquido
        self.polaridade      = row[0][15] # A polaridade do componente

        #==============================================================================
        #  PARÂMETROS DO CÁLCULO DE PRESSÃO DE VAPOR (The Properties of Gases & Liquids, 4th edition)
        #==============================================================================
        
        if self.eqPsat == self.__lista_EqPsat[0]:
            
            self.__cursor.execute('SELECT * FROM Parametros_Psat_Prausnitz_4th_edition WHERE  ID_componente=? AND ID_forma=?',(self.ID,self.nEqPsat))
            row = self.__cursor.fetchall() 
            
            #==============================================================================
            # PARÂMETROS PARA O CÁLCULO DE PSAT FORNECIDO PELO Prausnitz_4th_edition 
            #==============================================================================
            
            self.VPA             = row[0][3]
            self.VPB             = row[0][4]
            self.VPC             = row[0][5]
            self.VPD             = row[0][6]
            
            #==============================================================================
            # FAIXA DE TEMPERATURA NA QUAL PODE-SE APLICAR OS CÁLCULOS DE PSAT
            #==============================================================================
            
            self.__TminPsat = row[0][7]  #  Temperatura mínima de Psat  
            self.__TmaxPsat = row[0][8]  #  Temperatura máxima de Psat
            
            #==============================================================================
            # VALIDAÇAO DA TEMPERATURA, SE A MESMA PERTENCE AO INTERVALO DE APLICAÇAO
            #==============================================================================
            
            self.warnings()    
            self.Psat = self.Pvap_Prausnitz_4th(self.VPA,self.VPB,self.VPC,self.VPD,self.T,self.nEqPsat,self.Tc,self.Pc)
        
        
class Modelo:

    def __init__(self,Componentes):
        u'''
        Classe auxiliar para facilitar as buscas e validações dos parâmetros dos modelos.
        
        ========
        Entradas
        ========
        
        * Componentes (list): É uma lista de objetos ``Componente_Caracterizar``, vide documentação da dessa classe.
        
        =========
        Atributos        
        =========
        
        Os atributos presentes nesta classe são:

            * ``__ID_Componentes`` (list): Uma lista contendo os ID's dos componentes. Os ID's são buscados como atributos da classe ``Componente_Caracterizar``, vide documentação dessa classe;
            * ``lista_forma_eq`` (list): Uma lista contendo as formas de equações disponíveis no Banco de dados para determinada mistura e modelo;
            * ``FormaEq`` (int): ID da forma da equação utilizada, conforme consta no Banco de dados.
                
        =======
        Métodos
        =======
        
        Os métodos disponíveis desta classe são:
        
        * ``Parametros``: 
            * Método utilizado para busca dos parâmetros dos modelos. Vide documentação do método.
        * ``Busca_e_Validacao_da_faixa_Temp``:
            * Método para buscar a faixa de temperatura de aplicação do modelo no Banco de dados e validar se a temperatura inserida está dentro dessa faixa. Vide documentação do método.
        * ``FormaEquacao``:
            * Método usado para realizar a busca das possíveis formas de equação para determinada mistura e determinado modelo. Vide documentação do método.
        * ``ValidacaoFormaEq``:
            * Método utilizado para validar a forma de equação inserida. Vide documentação do método.
        '''
        
        #==============================================================================
        #         CONEXAO COM O BANCO DE DADOS
        #==============================================================================
        self.__conector = connect('THERMO_DATA_BANK_EXEMPLO.db')  # Conecta a rotina ao banco de dados
        self.__cursor = self.__conector.cursor()

        #==============================================================================
        #         VALIDA SE A ENTRADA Componentes É UM OBJETO DA CLASSE Componente_Caracterizar        
        #==============================================================================
        teste                 = [isinstance(elemento,Componente_Caracterizar) for elemento in Componentes]
        if False in teste:
            raise NameError(u'A entrada Componentes não é um objeto da classe Componente_Caracterizar (Vide documentação da classe).')
            
        #==============================================================================
        #         BUSCA DOS ID'S    
        #==============================================================================
        self.__ID_Componentes = [Componente.ID for Componente in Componentes] # Criação da lista com as ID's dos componentes
        
    def Busca_Parametros(self,tabela,coluna,IDFORMA=False):
        u'''
        Método utilizado para busca dos parâmetros dos modelos.
        
        ========
        Entradas
        ========
        
        * tabela (str): Nome da tabela, conforme consta no Banco de dados, na qual a coluna se encontra;
        * coluna (str): Nome da coluna, conforme consta no Banco de dados, da qual deseja-se buscar o parâmetro;
        * IDFORMA (int): O ID da forma da equação desejada, vide documentação do Banco de dados. Caso o modelo inserido não possua diferentes formas de equações, esta entrada não é necessária.
                
        ======
        Saídas
        ======
        
        * Este método retorna uma lista de listas, contendo os parâmetros, que são números inteiros. 
        
        =========
        Exemplo 1
        =========
        
        * Para os modelos: NRTL, UNIQUAC e Wilson, a IDFORMA deve ser inserida. Como pode ser visto no exemplo a seguir: ::
        
            parametro_int = Parametros('ParametroInteracao','UNIQUAC_parametros_interacao_binaria',1)    
            
        * O retorno deste método é da seguinte forma: ::
        
            [[0.0, 124.5], [365.78, 0.0]]
        
        =========
        Exemplo 2
        =========
        
        * Para equações como a equação Virial, não se faz necessário a entrada de uma forma de equação. Assim, o método é acessado da seguinte forma: ::
        
            coef_solv = Parametros('CoeficienteSolvatacao','Propriedade_mistura')   
            
        * O retorno deste método é da seguinte forma: ::
        
            [[5.67, 9.8], [8.75, 4.68]]
    
        '''
        
        #==============================================================================
        #         Criação da matriz de zeros em forma de lista de listas
        #==============================================================================
        retorno = zeros((len(self.__ID_Componentes),len(self.__ID_Componentes))).tolist() 
        
        #==============================================================================
        #         Substituição dos elementos da matriz criada
        #==============================================================================
        for i,ID_i in enumerate(self.__ID_Componentes):
            for j,ID_j in enumerate(self.__ID_Componentes):
                if IDFORMA == False:
                    selecao = 'SELECT '+coluna+' FROM '+tabela+' WHERE ID_componente_i=? AND ID_componente_j=? '
                else:
                    selecao = 'SELECT '+coluna+' FROM '+tabela+' WHERE ID_forma=%d AND ID_componente_i=? AND ID_componente_j=? '%IDFORMA
                self.__cursor.execute(selecao,(ID_i,ID_j))
                row = self.__cursor.fetchall() 
                retorno[i][j] = row[0][0]
                
        return retorno
    
    def Busca_e_Validacao_da_faixa_Temp(self,tabela,T,FormaEq):
        u'''
        Método para buscar a faixa de temperatura de aplicação do modelo no Banco de dados e validar se a temperatura inserida está dentro dessa faixa.
        
        ========
        Entradas
        ========
        
        * tabela (str): Nome da tabela, conforme consta no Banco de dados, na qual deseja-se buscar a faixa de temperatura de aplicação do modelo;
        * T (float): Temperatura em Kelvin na qual deseja-se trabalhar;
        * FormaEq (int): ID da forma da equação utilizada, conforme consta no Banco de dados.

        '''
        
        #==============================================================================
        #         Busca da faixa de temperatura no banco de dados        
        #==============================================================================        
        selecao    =  'SELECT TempMin, TempMax FROM '+tabela+' WHERE ID_componente_i=? AND ID_componente_j=? AND ID_forma=?'      
        self.__cursor.execute(selecao,(self.__ID_Componentes[0],self.__ID_Componentes[1],FormaEq)) 
        faixa      =  self.__cursor.fetchall()

        #==============================================================================
        #         Validação da temperatura
        #==============================================================================
        if T < faixa[0][0] or T > faixa[0][1]:
            warn(u'A temperatura especificada está fora da faixa de aplicabilidade da mistura utilizada para o modelo desejado. A temperatura deve pertencer ao intervalo: (%f, '%faixa[0][0]+'%f).'%faixa[0][1]) # Emite um aviso quando a temperatura está fora da faixa de aplicação. Contudo o programa continua a rodar

        
    def FormaEquacao(self,tabela):
        u'''
        Método usado para realizar a busca das possíveis formas de equação para determinada mistura e determinado modelo.
        
        ========
        Entradas
        ========
        
        * tabela (str): Nome da tabela, conforme consta no Banco de dados, na qual deseja-se buscar as formas de equações
        
        ======
        Saídas
        ======
        
        * Gera o atributo ``lista_forma_eq``. Este é uma lista contendo as formas de equações disponíveis no Banco de dados para determinada mistura e modelo;
        
        '''
        self.tabela = tabela # Criação do atributo tabela
        
        self.__cursor.execute('SELECT ID_forma FROM '+self.tabela+' WHERE ID_componente_i=? AND ID_componente_j=?',(self.__ID_Componentes[0],self.__ID_Componentes[1]))
        row                 = self.__cursor.fetchall() # linha contendo as formas de equações disponíveis em forma de lista de tupla..
        self.lista_forma_eq = [i[0] for i in row] # Criação do atriubto lista_forma_eq em forma de lista de inteiros.
        
         
         
    def ValidacaoFormaEq(self):
        u'''
        Método utilizado para validar a forma de equação inserida. 
        A validação ocorre com dois casos maiores: caso a forma de equação for inserida e 
        caso não for inserida. No primeiro caso, se a forma de equação inserida não constar
        em ``lista_forma_eq`` haverá um erro. 
        No segundo, se a forma da equação não for inserida
        e houver apenas uma forma de equação disponível para o modelo e mistura desejados, esta
        forma de equação será usada. Se houver mais de uma forma de equação disponível para o modelo
        e mistura desejados, há a necessidade de informar expressamente a forma da equação a ser usada.
        '''
        
        if self.formaEq == None:                         # Caso a forma da equação não seja inserida
            if len(self.lista_forma_eq) == 1:            # Caso haja apenas uma forma de equação, dados banco tem tamanho 1
                self.formaEq = self.lista_forma_eq[0]    # Como não seja inserida a forma de equação, o programa usará a única forma disponível
            else: # Caso haja mais de uma forma de equação disponível no banco de dados
                raise ValueError(u'Foi escolhido o modelo %s. No entanto, é necessário informar expressamente a forma da equaçao (Vide documentação do Banco de dados), visto que há diferentes opções disponíveis para a mistura desejada: '%(self.tabela,)+', '.join(str(model) for model in self.lista_forma_eq)+'.')         
        
        if self.formaEq not in self.lista_forma_eq: # Caso a forma da equação inserida não conste no banco de dados
            raise ValueError(u'A forma de equação inserida não consta no Banco de dados para a mistura e o modelo desejados (Vide documentação do Banco de dados). Para o caso requerido, as formas de equações disponíveis são: '+', '.join(str(model) for model in self.lista_forma_eq)+'.')
                
class VIRIAL(Modelo):
   
    def __init__(self,Componentes,regra_mistura='Hayden_o_Connel',parametro_int=None):
        u'''
        Rotina para busca dos parâmetros da equação Virial, vide [1].
        
        ========
        Entradas
        ========
        
        * Componentes (list): É uma lista de objetos ``Componente_Caracterizar``, vide documentação da dessa classe;
        * regra_mistura (str): Nome da regra utilizada para as correlações da equação Virial. Esta entrada possui a regra 'Hayden_o_Connel', vide [2], como valor padrão.
        
            * A lista de regras disponíveis pode ser acessada da seguinte forma: ::
            
                VIRIAL(None)
        
        =========
        Atributos
        =========
        
        * ``coef_solv``: Uma lista de listas contendo os coeficientes de solvatação e associação da mistura desejada.
        
        
        =======
        Exemplo 
        =======
        
        Como já foi citado, a entrada ``Componentes`` é uma lista de objetos ``Componente_Caracterizar``, portanto
        o primeiro passo para usar esta classe é acessar a classe ``Componente_Caracterizar``, vide documentação da
        classe, da seguinte forma: ::
        
            Comp1 = Componente_Caracterizar('Metano',ConfigPsat=('Prausnitz4th',1),T=100.0)
            Comp2 = Componente_Caracterizar('Etano',ConfigPsat=('Prausnitz4th',1),T=289.9)

        Em seguida, a classe ``VIRIAL`` pode ser acessada do seguinte modo: ::

            VIRIAL([Comp1,Comp2])
            
        * Como não foi inserido regra de mistura, a rotina utilizará o valor padrão ('Hayden_o_Connel').
        
        ===========
        Referências
        ===========
        
        [1] MASON, E. A.; SPURLING, T. H. The Virial Equation of State; The
        International Encyclopedia of Physical Chemistry and Chemical
        Physics, Topic 10: The Fluid State, Vol. 2; Pergamon Press: New
        York, 1969; p 297
        
        [2] HAYDEN, J. G.; O’CONNELL, J. P. A Generalized Method for Predicting
        Second Virial Coefficients. Industrial & Engineering Chemistry Process Design and
        Development, v. 14, n. 3, p. 209–216, jul. 1975. ISSN 0196-4305.
        
        [3] TSONOPOULOS, C.; HEIDMAN, J.L. From the Virial to the cubic equation of state. Fluid Phase Equilib. 57 (1990) 261–276.
        '''
        
        #==============================================================================
        #         NOME DO MODELO
        #==============================================================================
        self.nome_modelo = 'Virial' # Atributo útil para a rotina VLE
        
        #==============================================================================
        #         REGRA DE MISTURA
        #==============================================================================
        self.regra_mistura = regra_mistura
        self.ValidacaoREGRA()
        
        #==============================================================================
        #         MOSTRAR REGRAS DISPONÍVEIS
        #==============================================================================
        if Componentes == None:
            print u'As seguintes regras de mistura estão disponíveis: '+', '.join(self.__regras_mistura_disponiveis)+'.'
            raise NameError(u'Insira uma das regras mostradas na lista.')       
    
        #==============================================================================
        #         BUSCA ID NA CLASSE MÃE (CLASSE MODELO)  
        #==============================================================================
        Modelo.__init__(self,Componentes) 
        
        #==============================================================================
        #         REGRA PARA O CÁLCULO DO SEGUNDO COEFICIENTE     
        #==============================================================================
        if self.regra_mistura == 'Tsonopoulos':
            # Parâmetro utilizado em Hayden O'Connel
            if parametro_int == None:
                self.k_int_binaria = self.Busca_Parametros('Tsonopoulos','kij')
            else:
                self.k_int_binaria  = parametro_int
        elif self.regra_mistura == 'Hayden_o_Connel':
            # Parâmetro utilizado em Hayden O'Connel
            if parametro_int == None:
                self.coef_solv  = self.Busca_Parametros('Propriedade_mistura','CoeficienteSolvatacao') # Trasforma a def Parametros da classe Modelo em atributo da classe VIRIAL
            else:
                self.coef_solv  = parametro_int
        
        self._Modelo__conector.close()
            
    def ValidacaoREGRA(self):
        u'''
        Método para validação da regra de mistura da equação Virial. Caso a regra inserida não constar
        entre as regras disponíveis, ocorrerá um erro.
        '''
        self.__regras_mistura_disponiveis = ['Hayden_o_Connel','Tsonopoulos'] 
        if self.regra_mistura not in self.__regras_mistura_disponiveis:
            raise NameError(u'A regra de mistura inserida não está disponível. Regras disponíveis: '+'%, '.join(self.__regras_mistura_disponiveis)+'.')
    

class UNIQUAC(Modelo):
   
    def __init__(self,Componentes,T,FormaEqUNIQUAC = None,parametro_int=None):
        u'''
        Rotina para busca dos parâmetros do modelo UNIQUAC, vide [1].
        
        ========
        Entradas
        ========
        
        * Componentes (list): É uma lista de objetos ``Componente_Caracterizar``, vide documentação da dessa classe;
        * T (float): Temperatura, em Kelvin, na qual deseja-se trabalhar;
        * FormaEqUNIQUAC (int): ID da forma da equação utilizada para o modelo UNIQUAC, conforme consta no Banco de dados.
        
        =========
        Atributos
        =========
        
        * ``parametro_int``: Uma lista de listas contendo os parâmetros de interação da mistura desejada.
                
        =======
        Exemplo 
        =======
        
        Como já foi citado, a entrada ``Componentes`` é uma lista de objetos ``Componente_Caracterizar``, portanto
        o primeiro passo para usar esta classe é acessar a classe ``Componente_Caracterizar``, vide documentação da
        classe, da seguinte forma: ::
        
            Comp1 = Componente_Caracterizar('Metano',ConfigPsat=('Prausnitz4th',1),T=100.0)
            Comp2 = Componente_Caracterizar('Etano',ConfigPsat=('Prausnitz4th',1),T=289.9)

        Em seguida, a classe ``UNIQUAC`` pode ser acessada do seguinte modo: ::

            modelo = UNIQUAC([Comp1,Comp2],278.45,1)
            
        A entrada ``FormaEqUNIQUAC`` pode ser ``None``. ::
        
            modelo = UNIQUAC([Comp1,Comp2],278.45,None)
        
        * Dessa forma, caso conste apenas uma forma de equação para a mistura desejada, a escolha da forma de equação fica por conta do Banco de dados. Caso conste mais de uma forma de equação no Banco de dados, é necessário informar uma forma de equação.
        
        ===========
        Referências
        ===========
        
        [1] ABRAMS, D. S.; PRAUSNITZ, J. M. Statistical thermodynamics of liquid
        mixtures: A new expression for the excess Gibbs energy of partly or completely
        miscible systems. AIChE Journal, v. 21, n. 1, p. 116–128, jan. 1975. ISSN
        0001-1541. Disponível em: <http://doi.wiley.com/10.1002/aic.690210115>
        '''
        
        #==============================================================================
        #         NOME DO MODELO
        #==============================================================================
        self.nome_modelo = 'UNIQUAC' # Atributo útil para a rotina VLE
        
        #==============================================================================
        #         BUSCA ID NA CLASSE MÃE (CLASSE MODELO)        
        #==============================================================================
        Modelo.__init__(self,Componentes)
        
        #==============================================================================
        #         FORMA DE EQUACAO 
        #==============================================================================
        self.FormaEquacao('UNIQUAC_parametros_interacao_binaria')
        self.formaEq       = FormaEqUNIQUAC
        self.ValidacaoFormaEq()
        
        #==============================================================================
        #         FAIXA DE TEMPERATURA
        #==============================================================================
        self.Busca_e_Validacao_da_faixa_Temp('UNIQUAC_parametros_interacao_binaria',T,self.formaEq)
    
        #==============================================================================
        #         PARAMETROS DO MODELO
        #==============================================================================
        if parametro_int == None:
            self.parametro_int = self.Busca_Parametros('UNIQUAC_parametros_interacao_binaria','ParametroInteracao',self.formaEq) # Trasforma a def Parametros da classe Modelo em atributo da classe UNIQUAC
        else:
            self.parametro_int = parametro_int

        #==============================================================================
        #         ENCERRAR CONEXÃO COM O BANCO
        #==============================================================================         
        self._Modelo__conector.close()
         
class NRTL(Modelo):
   
    def __init__(self,Componentes,T,FormaEqNRTL = None,parametro_int=None,alpha = None):
        u'''
        Rotina para busca dos parâmetros do modelo NRTL, vide [1].
        
        ========
        Entradas
        ========
        
        * Componentes (list): É uma lista de objetos ``Componente_Caracterizar``, vide documentação da dessa classe;
        * T (float): Temperatura, em Kelvin, na qual deseja-se trabalhar;
        * FormaEqNRTL (int): ID da forma da equação utilizada para o modelo NRTL, conforme consta no Banco de dados.
        
        =========
        Atributos
        =========
        
        * ``parametro_int``: Uma lista de listas contendo os parâmetros de interação da mistura desejada;
        * ``alpha``: Um lista de listas contendo os parâmetros da mistura desejada.
        
        =======
        Exemplo 
        =======
        
        Como já foi citado, a entrada ``Componentes`` é uma lista de objetos ``Componente_Caracterizar``, portanto
        o primeiro passo para usar esta classe é acessar a classe ``Componente_Caracterizar``, vide documentação da
        classe, da seguinte forma: ::
        
            Comp1 = Componente_Caracterizar('Metano',ConfigPsat=('Prausnitz4th',1),T=100.0)
            Comp2 = Componente_Caracterizar('Etano',ConfigPsat=('Prausnitz4th',1),T=289.9)

        Em seguida, a classe ``NRTL`` pode ser acessada do seguinte modo: ::

            modelo = NRTL([Comp1,Comp2],278.45,1)
            
        A entrada ``FormaEqNRTL`` pode ser ``None``. ::
        
            modelo = NRTL([Comp1,Comp2],278.45,None)
        
        * Dessa forma, caso conste apenas uma forma de equação para a mistura desejada, a escolha da forma de equação fica por conta do Banco de dados. Caso conste mais de uma forma de equação no Banco de dados, é necessário informar uma forma de equação.
               
        ===========
        Referências
        ===========
        
        [1] RENON, H.; PRAUSNITZ, J. M. Local compositions in thermodynamic excess
        functions for liquid mixtures. AIChE Journal, v. 14, n. 1, p. 135–144, jan. 1968.
        ISSN 0001-1541. Disponível em: <http://doi.wiley.com/10.1002/aic.690140124>.
        '''
        
        #==============================================================================
        #         NOME DO MODELO
        #==============================================================================
        self.nome_modelo = 'NRTL' # Atributo útil para a rotina VLE
        
        #==============================================================================
        #         BUSCA ID NA CLASSE MÃE (CLASSE MODELO)
        #==============================================================================
        Modelo.__init__(self,Componentes)
 
        #==============================================================================
        #         FORMA DE EQUACAO 
        #==============================================================================
        self.FormaEquacao('NRTL')                
        self.formaEq       = FormaEqNRTL        
        self.ValidacaoFormaEq()
        
        #==============================================================================
        #         FAIXA DE TEMPERATURA
        #==============================================================================
        self.Busca_e_Validacao_da_faixa_Temp('NRTL',T,self.formaEq)
        
        #==============================================================================
        #         PARAMETROS DO MODELO
        #==============================================================================
        if parametro_int == None:
            self.parametro_int = self.Busca_Parametros('NRTL','ParametroInteracao',self.formaEq) # Trasforma a def Parametros da classe Modelo em atributo da classe NRTL
        else:
            self.parametro_int = parametro_int
        
        if alpha == None:
            self.alpha   = self.Busca_Parametros('NRTL','alphaij',self.formaEq)              
        else:
            self.alpha   = alpha

        #==============================================================================
        #         ENCERRAR CONEXÃO COM O BANCO
        #==============================================================================        
        self._Modelo__conector.close()

class WILSON(Modelo):
    
    def __init__(self,Componentes,T,FormaEqWILSON = None,parametro_int=None):
        u'''
        Rotina para busca dos parâmetros do modelo de Wilson, vide [1].
        
        ========
        Entradas
        ========
        
        * Componentes (list): É uma lista de objetos ``Componente_Caracterizar``, vide documentação da dessa classe;
        * T (float): Temperatura, em Kelvin, na qual deseja-se trabalhar;
        * FormaEqWILSON (int): ID da forma da equação utilizada para o modelo de Wilson, conforme consta no Banco de dados.
        
        =========
        Atributos
        =========
        
        * ``parametro_int``: Uma lista de listas contendo os parâmetros de interação da mistura desejada.
        
        =======
        Exemplo 
        =======
        
        Como já foi citado, a entrada ``Componentes`` é uma lista de objetos ``Componente_Caracterizar``, portanto
        o primeiro passo para usar esta classe é acessar a classe ``Componente_Caracterizar``, vide documentação da
        classe, da seguinte forma: ::
        
            Comp1 = Componente_Caracterizar('Metano',ConfigPsat=('Prausnitz4th',1),T=100.0)
            Comp2 = Componente_Caracterizar('Etano',ConfigPsat=('Prausnitz4th',1),T=289.9)

        Em seguida, a classe ``WILSON`` pode ser acessada do seguinte modo: ::

            modelo = WILSON([Comp1,Comp2],278.45,1)
            
        A entrada ``FormaEqWILSON`` pode ser ``None``. ::
        
            modelo = WILSON([Comp1,Comp2],278.45,None)
        
        * Dessa forma, caso conste apenas uma forma de equação para a mistura desejada, a escolha da forma de equação fica por conta do Banco de dados. Caso conste mais de uma forma de equação no Banco de dados, é necessário informar uma forma de equação.        
        
        ===========
        Referências
        ===========
        
        [1] WILSON, G. M. Vapor-Liquid Equilibrium. XI. A New Expression for the Excess
        Free Energy of Mixing. Journal of the American Chemical Society, v. 86, n. 2, p.
        127–130, jan. 1964. ISSN 0002-7863. Disponível em: <http://pubs.acs.org/doi/abs-
        /10.1021/ja01056a002>
        '''
        
        #==============================================================================
        #         NOME DO MODELO
        #==============================================================================
        self.nome_modelo = 'Wilson' # Atributo útil para a rotina VLE
        
        #==============================================================================
        #         BUSCA ID NA CLASSE MÃE (CLASSE MODELO)        
        #==============================================================================
        Modelo.__init__(self,Componentes) 
       
        #==============================================================================
        #         FORMA DE EQUACAO        
        #==============================================================================
        self.FormaEquacao('Wilson')                
        self.formaEq       = FormaEqWILSON                                
        self.ValidacaoFormaEq()
        
        #==============================================================================
        #         FAIXA DE TEMPERATURA
        #==============================================================================
        self.Busca_e_Validacao_da_faixa_Temp('Wilson',T,self.formaEq)

        #==============================================================================
        #         PARAMETROS DO MODELO
        #==============================================================================
        if parametro_int == None:
                    self.parametro_int = self.Busca_Parametros('Wilson','ParametroInteracao',self.formaEq) # Trasforma a def Parametros da classe Modelo em atributo da classe WILSON
        else:
            self.parametro_int = parametro_int

        #==============================================================================
        #         ENCERRAR CONEXÃO COM O BANCO
        #==============================================================================
        self._Modelo__conector.close()
        
class Van_Laar(Modelo):
    
    def __init__(self,Componentes,parametro=None):
        u'''
        Rotina para busca dos parâmetros do modelo de Wilson, vide [1].
        
        ========
        Entradas
        ========
        
        * Componentes (list): É uma lista de objetos ``Componente_Caracterizar``, vide documentação da dessa classe;
        
        =========
        Atributos
        =========
        
        * ``parametro``: Uma lista de listas contendo os parâmetros do modelo de Van Laar para a mistura desejada.
        
        =======
        Exemplo 
        =======
        
        Como já foi citado, a entrada ``Componentes`` é uma lista de objetos ``Componente_Caracterizar``, portanto
        o primeiro passo para usar esta classe é acessar a classe ``Componente_Caracterizar``, vide documentação da
        classe, da seguinte forma: ::
        
            Comp1 = Componente_Caracterizar('Metano',ConfigPsat=('Prausnitz4th',1),T=100.0)
            Comp2 = Componente_Caracterizar('Etano',ConfigPsat=('Prausnitz4th',1),T=289.9)

        Em seguida, a classe ``Van_Laar`` pode ser acessada do seguinte modo: ::

            modelo = Van_Larr([Comp1,Comp2])
        
        ===========
        Referências
        ===========
        
        [1] VAN LAAR, J. J. The Vapor pressure of binary mixtures. Z. Phys.
        Chem. 1910, 72, 723−751.
        '''
        #==============================================================================
        #         NOME DO MODELO
        #==============================================================================
        self.nome_modelo = 'Van Laar' # Atributo útil para a rotina VLE
        
        #==============================================================================
        #         BUSCA ID NA CLASSE MÃE (CLASSE MODELO)        
        #==============================================================================
        Modelo.__init__(self,Componentes) 

        #==============================================================================
        #         PARAMETROS DO MODELO
        #==============================================================================
        if parametro == None:
            self.parametro = self.Busca_Parametros('Van_Laar','Parametro') # Trasforma a def Parametros da classe Modelo em atributo da classe Van_Laar
        else:
            self.parametro = parametro
            
        #==============================================================================
        #         ENCERRAR CONEXÃO COM O BANCO
        #==============================================================================
        self._Modelo__conector.close()