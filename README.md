# SIR-mono
Códigos utilizados na monografia.
Esses códigos foram utilizados para gerar os dados da monografia. 

O código modelosirRodovia.c foi utilizado para gerar dados para a simulação em uma dimensão com condição de contorno fechada e para uma única população. Foi construído na linguagem de programação C e pode ser executado com a linha de comando 'gcc modelosirRodovia.c -L/user/local/lib64 -lm -o t2.out' utilizando o compilador gcc.
As constantes: 
#define QUANT_AMOSTRAS 100
#define POPULACOES 30
#define N 1000000
contém a quantidade de amostras que serão geradas, a quantidade de populações consideradas na simulação e o número de indivíduos em cada população, respectivamente.
Na linha 26 deste código é possível alterar a constante 'a', a qual altera a taxa de difusão entre as cidades.
Nas linhas 53, 54 e 56, são definidos parâmetros epidêmicos e o parâmetro 'f', que reduz a taxa de contatos entre os indivíduos. Taxa esta que foi definida nas linhas 80 e 83.
Ao final da execução do código, são retornados os arquivos com: a média das densidades de infectados e removidos no tempo; a densidade de infectados e removidos por população; o pico de prevalência e o tempo para alcançar o pico de prevalência; a largura das curvas de infectados por população.

O código modelosirPlano.c foi utilizado para gerar os dados para a simulação em duas dimensões de uma rede quadrada com condição de contorno fechada. Foi construído na linguagem de programação C e pode ser executado com a linha de comando 'gcc modelosirPlano.c -L/user/local/lib64 -lm -o t2.out' utilizando o compilador gcc.
Esse código é similar ao código modelosirRodovia.c, a sua principal diferença está na linha 6 com a constante 'L', que determina quantas populações terão em cada linha e em cada coluna. Ele também retorna os mesmos resultados.

O código runge-kutta.c foi utilizado para a obtenção dos dados determinísticos do modelo SIR a partir das taxas de difuão entre as populações. Pode ser utilizado também para uma única população. Foi construído na linguagem de programação C e pode ser executado com a linha de comando 'gcc runge-kutta.c -L/user/local/lib64 -lm -o t2.out' utilizando o compilador gcc. Esse código possui os mesmos parâmetros do código para o 'modelosirRodovia.c', a diferença está na maneira como os dados são obtidos. Ele também retorna os mesmos resultados.

O código difusaoMG.c foi utilizado para obtenção dos dados para o estado de Minas Gerais e necessida dos dados 'dadosMunicipioMG.dat' e 'difusaoMunicipiosMG.dat'. As taxas de difusão e a população dos municípios são importadas desses arquivos. Os dados de cada município estão organizados utilizando o código IBGE desses municípios. Os parâmetros epidêmicos e a taxa de contatos estão organizados de maneira similar ao código 'modelosirRodovia.c'. Ele também retorna os mesmos resultados.
