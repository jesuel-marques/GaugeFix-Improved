// Paramêtros da fixação de gauge

//  Códigos para os algoritmos de fixação

#define LA 00   
#define SR 10
#define SRE 20
#define MD 30

#define tolerance 1e-12
#define sweeps_to_measurement_e2 10   //  Quantidade de varreduras para fazer medição de e2

#define maxhits 2   //  Hits de iteração para maximização de Tr[w(n)]

#define omega_OR 1.83    //  Parâmetro omega da sobrerelaxação
#define pSRE 0.78   //  Parâmetro p da sobrerelaxação estocástica
#define alfaMD 0.222    //  Parâmetro alfa do método de Cornell ou máximo declive