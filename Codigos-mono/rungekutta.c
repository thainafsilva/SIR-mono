#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#define N 1000000
#define POPULACOES 30
#define QUANT_DADOS 100000
#define TEMPO_MAXIMO 1000

double h;
double dados[2*POPULACOES+1][QUANT_DADOS], media[QUANT_DADOS];
double beta, mu;
double k[POPULACOES];
double tempo=0.0;
double picoEpidemia[2][POPULACOES];
double Ni[POPULACOES][POPULACOES], Ns[POPULACOES][POPULACOES], Nr[POPULACOES][POPULACOES];
double d[POPULACOES][POPULACOES];
double NiTotal[POPULACOES], NsTotal[POPULACOES], NrTotal[POPULACOES], NiTotalSoma, NsTotalSoma;
double f;
double d_tot;
double d0;
double a=8.0;
double soma[2][POPULACOES], larguraCurva[POPULACOES];

void inicializacao(){
  int i, j;
  i=0;
  tempo=0.0;
  beta=0.06;
  mu=0.25;

  f=0.0;
  h=0.0001;

  memset(Ni, 0, POPULACOES*POPULACOES*sizeof(Ni[0][0]));
  memset(Ns, 0, POPULACOES*POPULACOES*sizeof(Ns[0][0]));
  memset(Nr, 0, POPULACOES*POPULACOES*sizeof(Nr[0][0]));
  memset(NiTotal, 0, POPULACOES*sizeof(NiTotal[0]));
  memset(NsTotal, 0, POPULACOES*sizeof(NsTotal[0]));
  memset(NrTotal, 0, POPULACOES*sizeof(NrTotal[0]));
  memset(larguraCurva, 0, POPULACOES*sizeof(larguraCurva[0]));
  memset(soma, 0, 2*POPULACOES*sizeof(soma[0][0]));

  for (i=0; i<POPULACOES; i++){
      if (i==0){
        Ni[i][i]=N*0.0001;
      }
      NiTotal[i]+=Ni[i][i];
      NiTotalSoma+=Ni[i][i];
      NsTotalSoma+=Ni[i][i];
      Ns[i][i]=N-Ni[i][i];
      NsTotal[i]+=Ns[i][i];
      NsTotalSoma+=Ns[i][i];
      k[i]=13;


      k[i]=(1-f)*k[i];
  }

  d_tot=0.0;
  d0=0.0;
  for(i=1; i<POPULACOES; i++){
      d_tot+=0.03/(pow(i,1));
      d0+=1.0/pow(i,a);
  }
  d0=d_tot/d0;

  d_tot=0.0;
  for (i=0; i<POPULACOES; i++){
    for (j=0; j<POPULACOES; j++){
      if (i!=j){
        d_tot+=d0/(pow(abs(i-j), a));
        d[i][j]=d0/(pow(abs(i-j),a));
        printf("%d %d: %lf\n", i, j, d[i][j]);
      }
    }
  }

}

void F(double x, double ki[][POPULACOES][POPULACOES], double ki1[][POPULACOES][POPULACOES]){
  int i, j, l;
  double populacaoTotal[POPULACOES];

  memset(populacaoTotal, 0, POPULACOES*sizeof(populacaoTotal[0]));
  if (POPULACOES>1){
    for (i=0; i<POPULACOES; i++){
      for (j=0; j<POPULACOES;j++){

        populacaoTotal[i]+=Ni[i][j]+x*ki1[0][i][j] + Ns[i][j] + x*ki1[1][i][j] + Nr[i][j] + x*ki1[2][i][j];


        //printf("Aqui: %lf %lf %lf %lf %lf %lf\n", Ni[i][j], ki1[0][i][j], Ns[i][j], x*ki1[1][i][j], Nr[i][j], x*ki1[2][i][j]);
      }
    }
  }
  else{
    populacaoTotal[0]=N;
  }

  for (i=0; i<POPULACOES; i++){
    ki[0][i][i]=-mu*(Ni[i][i]+x*ki1[0][i][i])+beta*k[i]*(Ni[i][i]+x*ki1[0][i][i])*(Ns[i][i]+x*ki1[1][i][i])/populacaoTotal[i];
    ki[1][i][i]=-beta*k[i]*(Ni[i][i]+x*ki1[0][i][i])*(Ns[i][i]+x*ki1[1][i][i])/populacaoTotal[i];
    ki[2][i][i]=mu*(Ni[i][i]+x*ki1[0][i][i]);
    for (j=0; j<POPULACOES; j++){
      if (j!=i){
        ki[0][i][i]+=beta*k[i]*(Ni[i][j]+x*ki1[0][i][j])*(Ns[i][i]+x*ki1[1][i][i])/populacaoTotal[i];
        ki[1][i][i]-=beta*k[i]*(Ni[i][j]+x*ki1[0][i][j])*(Ns[i][i]+x*ki1[1][i][i])/populacaoTotal[i];

        ki[0][i][i]+=-(Ni[i][i]+x*ki1[0][i][i])*d[i][j]+0.5*Ni[j][i]+x*ki1[0][j][i];
        ki[1][i][i]+=-(Ns[i][i]+x*ki1[1][i][i])*d[i][j]+0.5*Ns[j][i]+x*ki1[1][j][i];
        ki[2][i][i]+=Nr[j][i]+x*ki1[2][j][i];

        ki[0][i][j]=-mu*(Ni[i][j]+x*ki1[0][i][j])-0.5*(Ni[i][j]+x*ki1[0][i][j])+d[j][i]*(Ni[j][j]+x*ki1[0][j][j]);
        ki[1][i][j]=-0.5*(Ns[i][j]+x*ki1[1][i][j])+d[j][i]*(Ns[j][j]+x*ki1[1][j][j]);
        ki[2][i][j]=mu*(Ni[i][j]+x*ki1[0][i][j])-(Nr[i][j]+x*ki1[2][i][j]);



        for (l=0; l<POPULACOES; l++){
          ki[0][i][j]+=beta*k[i]*(Ni[i][l]+x*ki1[0][i][l])*(Ns[i][j]+x*ki1[1][i][j])/populacaoTotal[i];
          ki[1][i][j]-=beta*k[i]*(Ni[i][l]+x*ki1[0][i][l])*(Ns[i][j]+x*ki1[1][i][j])/populacaoTotal[i];

        }
      }
    }
  }
}

void rungeKutta(){
  double k1[3][POPULACOES][POPULACOES], k2[3][POPULACOES][POPULACOES], k3[3][POPULACOES][POPULACOES], k4[6][POPULACOES][POPULACOES];
  int i, j;


  memset(k1, 0, 3*POPULACOES*POPULACOES*sizeof(k1[3][0][0]));
  memset(k2, 0, 3*POPULACOES*POPULACOES*sizeof(k2[3][0][0]));
  memset(k3, 0, 3*POPULACOES*POPULACOES*sizeof(k3[3][0][0]));
  memset(k4, 0, 3*POPULACOES*POPULACOES*sizeof(k4[3][0][0]));

  F(0, k1, k1);
  F(0.5*h, k2, k1);
  F(0.5*h, k3, k2);
  F(h, k4, k3);

  NiTotalSoma=0;
  NsTotalSoma=0;
  for (i=0; i<POPULACOES; i++){
    NiTotal[i]=0;
    NsTotal[i]=0;
    NrTotal[i]=0;
    for (j=0; j<POPULACOES; j++){
      Ni[i][j]+=h*(k1[0][i][j]+2*k2[0][i][j]+2*k3[0][i][j]+k4[0][i][j])/6.0;
      Ns[i][j]+=h*(k1[1][i][j]+2*k2[1][i][j]+2*k3[1][i][j]+k4[1][i][j])/6.0;
      Nr[i][j]+=h*(k1[2][i][j]+2*k2[2][i][j]+2*k3[2][i][j]+k4[2][i][j])/6.0;
      NiTotal[i]+=Ni[i][j];
      NsTotal[i]+=Ns[i][j];
      NrTotal[i]+=Nr[i][j];
    }
    NsTotalSoma+=NsTotal[i];

    NiTotalSoma+=NiTotal[i];

  //  printf("Infectados: %lf %lf\n", NiTotal[i], NiTotalSoma);
  }
}

void largura_da_curva(){
  int i, t;
  double normalizacaoAux, aux;
  for (t=0; t<TEMPO_MAXIMO; t++){
    for (i=0; i<POPULACOES; i++){
      soma[0][i]+=(1.0*dados[2*i+1][t]);
    }
  }
  for (i=0; i<POPULACOES; i++){
    aux=soma[0][i];
    soma[0][i]=0.0;
    for (t=0; t<TEMPO_MAXIMO; t++){
      normalizacaoAux=(1.0*dados[2*i+1][t])/aux;
      soma[0][i]+=pow(t, 2)*normalizacaoAux;
      soma[1][i]+=t*normalizacaoAux;
    }
    larguraCurva[i]=sqrt(soma[0][i]-pow(soma[1][i],2));
  }

}

int main(){
  FILE *arq, *arq2, *arq3, *arq4, *arq5, *arq6;
  char nome[50]="";
  int n, i, j, tempoViagem;
  double t, dt;
  double tempoIntervalo;
  int auxTempoIntervalo;
  dt=0.1;
  auxTempoIntervalo=0;
  tempoIntervalo=0.1;

  inicializacao();

  memset(dados,0,(3)*(QUANT_DADOS)*sizeof(dados[0][0]));
  sprintf (nome, "mediaCidades-%d-tempo%d-%lf.dat", POPULACOES, (int)TEMPO_MAXIMO, a);
  arq = fopen (nome, "w");

  sprintf (nome, "pico-epidemia-%d-%d-%lf.dat", POPULACOES, (int)TEMPO_MAXIMO, a);
  arq2 = fopen (nome, "w");

  memset(dados,0,(2*POPULACOES+1)*(QUANT_DADOS)*sizeof(dados[0][0]));
  memset(picoEpidemia,0,(2)*(POPULACOES)*sizeof(dados[0][0]));
  memset(media, 0, QUANT_DADOS*sizeof(media[0]));

  sprintf (nome, "NumeroRecuperador-%d-%d%lf.dat", POPULACOES, (int)TEMPO_MAXIMO, a);
  arq6 = fopen (nome, "w");

  n=0;
  do{

      if (tempo>=t){

        media[n]+=1.0*NiTotalSoma/(1.0*N*POPULACOES);
        for (i=1, j=0; i<(POPULACOES*2+1); i+=2, j++){
          dados[i][n]+=(1.0*NiTotal[j])/(NiTotal[j]+NsTotal[j]+NrTotal[j]);
          dados[i+1][n]+=1.0*NrTotal[j]/(1.0*(NiTotal[j]+NsTotal[j]+NrTotal[j]));
          //if (picoEpidemia[1][j]<(1.0*NiTotal[j])/(NiTotal[j]+NsTotal[j]+NrTotal[j])){
            //picoEpidemia[0][j]=tempo;
          //  picoEpidemia[1][j]=(1.0*NiTotal[j])/(NiTotal[j]+NsTotal[j]+NrTotal[j]);

          //}
          if (picoEpidemia[1][j]<1.0*NiTotal[j]/(NiTotal[j]+NsTotal[j]+NrTotal[j])){
            picoEpidemia[0][j]=tempo;
            picoEpidemia[1][j]=(1.0*NiTotal[j])/(NiTotal[j]+NsTotal[j]+NrTotal[j]);

          }
        }




        t+=dt;
        n++;


      }
      rungeKutta();
      tempo+=h;
    }while (tempo<TEMPO_MAXIMO && NiTotalSoma>0.0000000001 && n<QUANT_DADOS);

    double mediaRecuperados[POPULACOES];
    memset(mediaRecuperados, 0, POPULACOES*sizeof(mediaRecuperados[0]));
    for (i=0; i<POPULACOES; i++){
      mediaRecuperados[i]+=NrTotal[i];
    }

    t=0.0;
    double auxiliar[2][QUANT_DADOS];


    sprintf (nome, "infectados-%d-cidades-tempo%d-%lf.dat", POPULACOES, (int)TEMPO_MAXIMO, a);
    arq3 = fopen (nome, "w");

    sprintf (nome, "larguraCurva-%d-cidades-tempo%d-%lf.dat", POPULACOES, (int)TEMPO_MAXIMO, a);
    arq4 = fopen (nome, "w");
    n=0;
    for (i=0;t<TEMPO_MAXIMO;i++){
      dados[0][i]= t;
      if(auxiliar[0][n]==0 && n!=0){
        break;
      }
      auxiliar[0][i]=0;
      auxiliar[1][i]=0;
      fprintf(arq3, "%lf ", dados[0][i]);
      for (j=1; j<(POPULACOES*2+1); j+=2){
        fprintf(arq3, "%lf %lf ", dados[j][i], dados[j+1][i]);
        auxiliar[0][i]+=dados[j][i];
        auxiliar[1][i]+=dados[j+1][i];
      }
      fprintf (arq3, "\n");
      t+=dt;
      n=i;

    }


    largura_da_curva();

    for (i=0; i<POPULACOES; i++){
    //  printf("%lf\n", picoEpidemia[1][i]);

    //  printf("%lf\n", picoEpidemia[1][i]);
      //getchar();
      fprintf(arq6, "%d %lf \n ", i, 1.0*mediaRecuperados[i]);
      fprintf(arq2, "%d %lf %lf\n", i, picoEpidemia[0][i], picoEpidemia[1][i]);
      fprintf(arq4, "%d %lf\n", i, larguraCurva[i]);
    }

    for (i=0; i<=(n); i++){
      fprintf(arq, "%lf %lf %lf %lf", dados[0][i], auxiliar[0][i]/(1.0*POPULACOES), media[i], auxiliar[1][i]/(1.0*POPULACOES));
      fprintf(arq, "\n");
    }


    return 0;
}
