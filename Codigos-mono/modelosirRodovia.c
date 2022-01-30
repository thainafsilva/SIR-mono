#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#define QUANT_AMOSTRAS 100
#define POPULACOES 30
#define N 1000000
#define QUANT_DADOS 10000
#define TEMPO_MAXIMO 1000

int u;
double dados[2*POPULACOES+1][QUANT_DADOS], media[QUANT_DADOS];
double lambda, mi;
double k[POPULACOES];
double tempo=0.0;
double taxaInfeccao[POPULACOES], taxaCura[POPULACOES];
double phi[TEMPO_MAXIMO]; //inserir quantidade de dias
double taxaInfeccaoTotal, taxaCuraTotal, sumTaxas;
double picoEpidemia[2][POPULACOES];
int Ni[POPULACOES][POPULACOES], Ns[POPULACOES][POPULACOES], Nr[POPULACOES];
double d[POPULACOES][POPULACOES];
int NiTotal[POPULACOES], NsTotal[POPULACOES], NiTotalSoma, NsTotalSoma;
double f;
double d_tot;
double d0;
double a=20.0;
double soma[2][POPULACOES], larguraCurva[POPULACOES];
int lista[2][POPULACOES];

int im=1073741824;
int a1 = 65539;
int z=2431413;
int w=521288629;
int jsr=362436069;
int jcong=916191069;

double rand01_kiss(){
  int kiss;
  z = 69069 * z + 1327217885;
  jsr = jsr ^ (jsr << 13);
  jsr = jsr ^ (jsr >> 17);
  jsr = jsr ^ (jsr << 5);
  w = 18000 * (w & 65535) + (w >> 16);
  jcong = 30903 * (jcong & 65535) + (jcong >> 16);
  kiss = z + jsr + (w << 16) + jcong;
  return kiss * (1.0/(4.0*im+4.0)) + 0.5;
}

void inicializacao(){
  int i, j;
  i=0;
  tempo=0.0;
  lambda=0.06;
  mi=0.25;
  f=0.0; //fator que reduz o número de contatos diários


  memset(lista, 0, 2*POPULACOES*sizeof(lista[0][0]));
  memset(lista, -1, 2*POPULACOES*sizeof(lista[0][0]));

  memset(Ni, 0, POPULACOES*POPULACOES*sizeof(Ni[0][0]));
  memset(Ns, 0, POPULACOES*POPULACOES*sizeof(Ns[0][0]));
  memset(Nr, 0, POPULACOES*sizeof(Nr[0]));
  memset(NiTotal, 0, POPULACOES*sizeof(NiTotal[0]));
  memset(NsTotal, 0, POPULACOES*sizeof(NsTotal[0]));
  memset(soma, 0, 2*POPULACOES*sizeof(soma[0][0]));
  memset(larguraCurva, 0, POPULACOES*sizeof(larguraCurva[0]));
  NiTotalSoma=0;
  NsTotalSoma=0;
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

  NiTotalSoma=0;
  NsTotalSoma=0;

  taxaCuraTotal=0;
  taxaInfeccaoTotal=0;



  for (i=0; i<POPULACOES; i++){
    NiTotal[i]=0;
    NsTotal[i]=0;
    for (j=0; j<POPULACOES; j++){
      NiTotal[i]+=Ni[i][j];
      NsTotal[i]+=Ns[i][j];
    }
    NiTotalSoma+=NiTotal[i];
    NsTotalSoma+=NsTotal[i];
    taxaInfeccao[i]=lambda*k[i]*NiTotal[i]*NsTotal[i]/(NiTotal[i]+NsTotal[i]+Nr[i]);
    taxaCura[i]=mi*NiTotal[i];
    taxaCuraTotal+=taxaCura[i];
    taxaInfeccaoTotal+=taxaInfeccao[i];
  }

  sumTaxas=taxaCuraTotal+taxaInfeccaoTotal;
}

void curar(){
  double zr, aux;
  int cidade1, cidade2;
  int i, j;

  zr=rand01_kiss()*taxaCuraTotal;
  aux=0;

  for (i=0; i<POPULACOES; i++){
    aux+=taxaCura[i];
    if (zr<aux){
      aux-=taxaCura[i];
      aux+= Ni[i][i]*mi;
      if (zr<aux){
        Ni[i][i]--;
        Nr[i]++;
        cidade1=i;
        cidade2=i;
        i=POPULACOES;
      }
      else{
        for(j=0; j<POPULACOES;j++){
          if (i!=j){
            aux+=Ni[i][j]*mi;
          }
          if (zr<aux){
            Ni[i][j]--;
            Nr[j]++;
            cidade1=i;
            cidade2=j;
            i=POPULACOES;
            j=POPULACOES;
          }
        }
      }
    }
  }
  //calcular taxas que mudaram cura ocorre na população j que está em i

  NiTotal[cidade1]--;
  NiTotalSoma--;


  taxaInfeccaoTotal-=taxaInfeccao[cidade1];
  taxaCuraTotal-=taxaCura[cidade1];



  if (cidade1!=cidade2){
    taxaInfeccaoTotal-=taxaInfeccao[cidade2];
    taxaInfeccao[cidade2]=lambda*k[cidade2]*NiTotal[cidade2]*NsTotal[cidade2]/(NiTotal[cidade2]+NsTotal[cidade2]+Nr[cidade2]);
    taxaInfeccaoTotal+=taxaInfeccao[cidade2];
  }


  taxaInfeccao[cidade1]=lambda*k[cidade1]*NiTotal[cidade1]*NsTotal[cidade1]/(NiTotal[cidade1]+NsTotal[cidade1]+Nr[cidade1]);
  taxaCura[cidade1]=mi*NiTotal[cidade1];


  taxaInfeccaoTotal+=taxaInfeccao[cidade1];
  taxaCuraTotal+=taxaCura[cidade1];



  sumTaxas=taxaCuraTotal+taxaInfeccaoTotal;
}

void infectar(){
  double zr, aux, aux2;
  int cidade1, cidade2;
  int i, j;
  zr=rand01_kiss()*taxaInfeccaoTotal;
  aux=0;

  for (i=0; i<POPULACOES; i++){
    aux+=taxaInfeccao[i];
    aux2=aux;
    if (zr<aux){
      zr=(rand01_kiss())*NsTotal[i];
      aux=Ns[i][i];
      if (zr<aux){
        if (lista[0][i]==-1){
          aux2-=taxaInfeccao[i];
          for (j=0; j<POPULACOES; j++){
            aux2+=lambda*k[i]*Ni[i][j]*NsTotal[i]/(NiTotal[i]+NsTotal[i]+Nr[i]);
            if(zr<aux2 && i!=j){
              lista[0][i]=j;
              lista[1][i]=2; //indice de pegou na cidade
              j=POPULACOES;
            }
          }
        }
        Ni[i][i]++;
        Ns[i][i]--;
        cidade1=i;
        cidade2=i;
        i=POPULACOES;
      }
      else{
        zr=rand01_kiss()*(NsTotal[i]-Ns[i][i]);
        aux=0;
        for(j=0; j<POPULACOES;j++){
          if (i!=j){
            aux+=Ns[i][j];
          }
          if (zr<aux){
            Ni[i][j]++;
            Ns[i][j]--;
            if (lista[0][j]==-1){
              lista[0][j]=i;
              lista[1][j]=1; //indice de pegou em viagem
            }
            cidade1=i;
            cidade2=j;
            i=POPULACOES;
            j=POPULACOES;
          }
        }
      }
    }
  }
  //cidade1 onde ocorreu a infecção da população da cidade2
  NiTotal[cidade1]++;
  NsTotal[cidade1]--;
  NiTotalSoma++;
  NsTotalSoma--;

  taxaInfeccaoTotal-=taxaInfeccao[cidade1];
  taxaCuraTotal-=taxaCura[cidade1];

  taxaInfeccao[cidade1]=lambda*k[cidade1]*NiTotal[cidade1]*NsTotal[cidade1]/(NiTotal[cidade1]+NsTotal[cidade1]+Nr[cidade1]);
  taxaCura[cidade1]=mi*NiTotal[cidade1];

  taxaInfeccaoTotal+=taxaInfeccao[cidade1];
  taxaCuraTotal+=taxaCura[cidade1];

  sumTaxas=taxaCuraTotal+taxaInfeccaoTotal;
}

void viajar(){
  int i, j, aux1, aux2;
  double zr;
  for (i=0; i<POPULACOES; i++){
    aux1=Ni[i][i];
    aux2=Ns[i][i];
    for(j=0; j<POPULACOES;j++){
      if (i!=j){
        Ni[j][i]=aux1*d[i][j];
        Ns[j][i]=aux2*d[i][j];
        Ni[i][i]-=Ni[j][i];
        Ns[i][i]-=Ns[j][i];
      }
    }
    for (j=0; j<POPULACOES; j++){
      if (i!=j){
        zr=rand01_kiss();
        if (zr<=(aux1*d[i][j]-floor(aux1*d[i][j])) && (Ni[i][i]-1)>=0){
          Ni[j][i]++;
          Ni[i][i]--;
        }

        zr=rand01_kiss();
        if (zr<=(aux2*d[i][j]-floor(aux2*d[i][j])) && (Ns[i][i]-1)>=0){
          Ns[j][i]++;
          Ns[i][i]--;
        }
        if (Ni[i][i]<0 || Ns[i][i]<0){
          printf("i j: %d %d: %d %d %d %d\n", i, j, Ni[i][i], NiTotal[i], Ns[i][i], Ns[j][i]);
          getchar();
        }
      }
    }

  }
  taxaCuraTotal=0;
  taxaInfeccaoTotal=0;
  NiTotalSoma=0;
  NsTotalSoma=0;

  for (i=0; i<POPULACOES; i++){
    NiTotal[i]=0;
    NsTotal[i]=0;
    for (j=0; j<POPULACOES; j++){
      NiTotal[i]+=Ni[i][j];
      NsTotal[i]+=Ns[i][j];

      //printf("POPULACOES %d %d Nij %d Sij %d\n", i, j, Ni[i][j], Ns[i][j]);
    }
    NiTotalSoma+=NiTotal[i];
    NsTotalSoma+=NsTotal[i];

    taxaInfeccao[i]=lambda*k[i]*NiTotal[i]*NsTotal[i]/(NiTotal[i]+NsTotal[i]+Nr[i]);
    taxaCura[i]=mi*NiTotal[i];
    taxaCuraTotal+=taxaCura[i];
    taxaInfeccaoTotal+=taxaInfeccao[i];
    //printf("%d: V %lf I %lf C %lf \n\n", i, taxaViagem[i], taxaInfeccao[i], taxaCura[i]);
  }
  //printf("%lf\n", taxaViagemTotal);
  sumTaxas=taxaCuraTotal+taxaInfeccaoTotal;


}

void retornar(){
  int i, j;

  for (i=0; i<POPULACOES; i++){
    for(j=0; j<POPULACOES;j++){
      if (i!=j){
        Ni[j][j]+=Ni[i][j];
        Ni[i][j]=0;
        Ns[j][j]+=Ns[i][j];
        Ns[i][j]=0;
      }
    }
  }
  NiTotalSoma=0;
  NsTotalSoma=0;
  taxaCuraTotal=0;
  taxaInfeccaoTotal=0;

  for (i=0; i<POPULACOES; i++){

    NiTotal[i]=Ni[i][i];
    NsTotal[i]=Ns[i][i];

    NiTotalSoma+=NiTotal[i];
    NsTotalSoma+=NsTotal[i];

    taxaInfeccao[i]=lambda*k[i]*NiTotal[i]*NsTotal[i]/(NiTotal[i]+NsTotal[i]+Nr[i]);
    taxaCura[i]=mi*NiTotal[i];
    taxaCuraTotal+=taxaCura[i];
    taxaInfeccaoTotal+=taxaInfeccao[i];
    //printf("%d: V %lf I %lf C %lf \n\n", i, taxaViagem[i], taxaInfeccao[i], taxaCura[i]);
  }
  //getchar();
  //printf("%lf\n", taxaViagemTotal);
  sumTaxas=taxaCuraTotal+taxaInfeccaoTotal;
}

void passo_tempo_sir(){
  double zr;
  int aux, aux1, i, j, NiAux, NsAux;
  double propInfectados, propCura, propViagem;


  if (sumTaxas>0){
    propCura=taxaCuraTotal/(taxaCuraTotal+taxaInfeccaoTotal);
    propInfectados=taxaInfeccaoTotal/(taxaCuraTotal+taxaInfeccaoTotal);



    zr=rand01_kiss();
    if (zr<propCura){
      curar();
    }
    else{
      infectar();
    }
  }
  else{
    printf("%d %lf %lf %lf %d %d\n", u, tempo, taxaCuraTotal, taxaInfeccaoTotal, NiTotalSoma, NsTotalSoma);
    getchar();
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

int main() {
  FILE *arq, *arq2, *arq3, *arq4, *arq5, *arq6;
  char nome[50]="";
  int n, i, j, tempoViagem;
  double t, dt, auxPicoEpidemia[2][POPULACOES], auxDViagem;
  int auxViagem=0;
  double tempoIntervalo;
  int auxTempoIntervalo, auxPOPULACAO;
  auxTempoIntervalo=0;
  tempoIntervalo=0.1;
  int mediaRecuperados[POPULACOES];

  dt=1.0; //passo de tempo para salvar os resultados

  //encontrar d0 em função de a
  d_tot=0.0;
  d0=0.0;
  for(i=1; i<POPULACOES; i++){
      d_tot+=0.03/(pow(i,1));
      d0+=1.0/pow(i,a);
  }
  //d0=d_tot/(d0*10.0);
  d0=d_tot/d0;

  //fim
  d_tot=0.0;
  for (i=0; i<POPULACOES; i++){
    //auxDViagem=0;
    for (j=0; j<POPULACOES; j++){
      if (i!=j){

        d_tot+=d0/(pow(abs(i-j), a));
        d[i][j]=d0/(pow(abs(i-j),a));
        printf("%d %d: %lf\n", i, j, d[i][j]);
        //auxDViagem+=d[i][j];
      }
    }
    //printf("%lf\n", auxDViagem);
  }
  //d[0][POPULACOES-1]=d[0][1];
  //d[POPULACOES-1][0]=d[0][1];
  //getchar();
  inicializacao();

  memset(dados,0,(3)*(QUANT_DADOS)*sizeof(dados[0][0]));
  sprintf (nome, "mediaCidades-%d-tempo%d-TESTE.dat", POPULACOES, (int)TEMPO_MAXIMO);
  arq = fopen (nome, "w");

  sprintf (nome, "pico-epidemia-%d-tempo%d.dat", POPULACOES, (int)TEMPO_MAXIMO);
  arq2 = fopen (nome, "w");
  printf("Nome: %s\n", nome);
  u=0;

  memset(dados,0,(2*POPULACOES+1)*(QUANT_DADOS)*sizeof(dados[0][0]));
  memset(picoEpidemia,0,(2)*(POPULACOES)*sizeof(dados[0][0]));
  memset(media, 0, QUANT_DADOS*sizeof(media[0]));
  memset(mediaRecuperados, 0, POPULACOES*sizeof(mediaRecuperados[0]));
  sprintf (nome, "NumeroRecuperador-%d-%d.dat", POPULACOES, (int)TEMPO_MAXIMO);
  arq6 = fopen (nome, "w");
  while (u<QUANT_AMOSTRAS){
    memset(auxPicoEpidemia,0,2*(POPULACOES)*sizeof(auxPicoEpidemia[0][0]));
    inicializacao();
    t=0.0;
    n=0;
    tempoViagem=0;
    auxViagem=0;
    printf("%d\n", u);
    do{
      if (tempo>=t){


        if (u==0){
          media[n]=1.0*NiTotalSoma/(1.0*N*POPULACOES);
          for (i=1, j=0; i<(POPULACOES*2+1); i+=2, j++){
            dados[i][n]=(1.0*NiTotal[j])/(NiTotal[j]+NsTotal[j]+Nr[j]);
            dados[i+1][n]=1.0*Nr[j]/(1.0*(NiTotal[j]+NsTotal[j]+Nr[j]));
            if (auxPicoEpidemia[1][j]<(1.0*NiTotal[j])/(NiTotal[j]+NsTotal[j]+Nr[j])){
              auxPicoEpidemia[0][j]=tempo;
              auxPicoEpidemia[1][j]=(1.0*NiTotal[j])/(NiTotal[j]+NsTotal[j]+Nr[j]);


            }
          }
        }
        else{
          media[n]+=1.0*NiTotalSoma/(1.0*N*POPULACOES);
          for (i=1, j=0; i<(POPULACOES*2+1); i+=2, j++){
            dados[i][n]+=(1.0*NiTotal[j])/(NiTotal[j]+NsTotal[j]+Nr[j]);
            dados[i+1][n]+=1.0*Nr[j]/(1.0*(NiTotal[j]+NsTotal[j]+Nr[j]));
            if (auxPicoEpidemia[1][j]<(1.0*NiTotal[j])/(NiTotal[j]+NsTotal[j]+Nr[j])){
              auxPicoEpidemia[0][j]=tempo;
              auxPicoEpidemia[1][j]=(1.0*NiTotal[j])/(NiTotal[j]+NsTotal[j]+Nr[j]);

            }
          }
        }
        t+=dt;
        n++;

      }
      if (tempoViagem<=(tempo-1) && auxViagem==1){
        tempoViagem++;
        viajar();
        //printf("Viajar: %lf\n", tempo);
        auxViagem=0;
      }
      if (tempoViagem<=(tempo-0.5) && auxViagem==0){
        retornar();
        //printf("Retornar: %lf\n", tempo);
        auxViagem=1;

      }

      passo_tempo_sir();
      if ((sumTaxas)>0){
        tempo+=1/sumTaxas;
      }
      else{
        //getchar();
      }
      //printf("%d %d\n", NiTotalSoma, NsTotalSoma);
    } while (tempo<TEMPO_MAXIMO && NiTotalSoma>0);

    for (i=0; i<POPULACOES; i++){
      picoEpidemia[0][i]+=auxPicoEpidemia[0][i];
      picoEpidemia[1][i]+=auxPicoEpidemia[1][i];
      mediaRecuperados[i]+=Nr[i];
    }
  //  getchar();
    u++;
  }


  t=0.0;
  double auxiliar[2][QUANT_DADOS];


  sprintf (nome, "infectados-%d-cidades-tempo%d.dat", POPULACOES, (int)TEMPO_MAXIMO);
  arq3 = fopen (nome, "w");

  sprintf (nome, "larguraCurva-%d-cidades-tempo%d.dat", POPULACOES, (int)TEMPO_MAXIMO);
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
      dados[j][i]/=QUANT_AMOSTRAS;
      dados[j+1][i]/=QUANT_AMOSTRAS;
      fprintf(arq3, "%lf %lf ", dados[j][i], dados[j+1][i]);
      auxiliar[0][i]+=dados[j][i];
      auxiliar[1][i]+=dados[j+1][i];
    }
    fprintf (arq3, "\n");
    t+=dt;
    n=i;
    printf("%d\n", i);
  }

  largura_da_curva();

  sprintf (nome, "lista-%d-tempo%d.dat", POPULACOES, (int)TEMPO_MAXIMO);
  arq5 = fopen (nome, "w");


  for (i=0; i<POPULACOES; i++){
    picoEpidemia[0][i]/=QUANT_AMOSTRAS;
  //  printf("%lf\n", picoEpidemia[1][i]);
    picoEpidemia[1][i]/=QUANT_AMOSTRAS;
  //  printf("%lf\n", picoEpidemia[1][i]);
    //getchar();
    fprintf(arq6, "%d %lf \n ", i, 1.0*mediaRecuperados[i]/1.0*QUANT_AMOSTRAS);
    fprintf(arq2, "%d %lf %lf\n", i, picoEpidemia[0][i], picoEpidemia[1][i]);
    fprintf(arq4, "%d %lf\n", i, larguraCurva[i]);
    fprintf(arq5, "%d %d %d\n", i, lista[0][i], lista[1][i]);
  }

  for (i=0; i<=(n); i++){
    fprintf(arq, "%lf %lf %lf %lf", dados[0][i], auxiliar[0][i]/(1.0*POPULACOES), media[i]/(1.0*QUANT_AMOSTRAS), auxiliar[1][i]/(1.0*POPULACOES));
    fprintf(arq, "\n");
  }

  fclose(arq);
  fclose(arq2);

  return 0;
}
