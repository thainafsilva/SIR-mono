#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#define QUANT_AMOSTRAS 10
#define L 21
#define N 1000000
#define QUANT_DADOS 10000
#define TEMPO_MAXIMO 1000

int u, distFocoIni[L*L];
double dados[2*(L*L)+1][QUANT_DADOS], media[QUANT_DADOS];
double lambda, mi;
double k[L*L];
double tempo=0.0;
double taxaInfeccao[L*L], taxaCura[L*L];
double phi[TEMPO_MAXIMO]; //inserir quantidade de dias
double taxaInfeccaoTotal, taxaCuraTotal, sumTaxas;
double picoEpidemia[2][L*L];
int Ni[L*L][L*L], Ns[L*L][L*L], Nr[L*L];
int NiTotal[L*L], NsTotal[L*L], NiTotalSoma, NsTotalSoma;
double f;
double d0;
double a=5.0;
double soma[2][L*L], larguraCurva[L*L];
int matriz[L][L];


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
  int i, j, l, focoIni;
  focoIni=L/2;
  i=0;
  tempo=0.0;
  lambda=0.06;
  mi=0.25;
  l=0;
  for (i=0; i<L; i++){
    for(j=0; j<L; j++){
      matriz[i][j]=l;
      distFocoIni[l]=abs(i-focoIni)+abs(j-focoIni);
      l++;

    }
  }


  memset(Ni, 0, L*L*L*L*sizeof(Ni[0][0]));
  memset(Ns, 0, L*L*L*L*sizeof(Ns[0][0]));
  memset(Nr, 0, L*L*sizeof(Nr[0]));
  memset(NiTotal, 0, L*L*sizeof(NiTotal[0]));
  memset(NsTotal, 0, L*L*sizeof(NsTotal[0]));
  memset(soma, 0, 2*L*L*sizeof(soma[0][0]));
  memset(larguraCurva, 0, L*L*sizeof(larguraCurva[0]));
  NiTotalSoma=0;
  NsTotalSoma=0;
  for (i=0; i<L*L; i++){
      if (i==focoIni){
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



  for (i=0; i<L*L; i++){
    NiTotal[i]=0;
    NsTotal[i]=0;
    for (j=0; j<L*L; j++){
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
//  printf("Soma das taxas Inicio: %lf %lf %lf\n", sumTaxas, taxaCuraTotal, taxaInfeccaoTotal);
}

void curar(){
  double zr, aux;
  int cidade1, cidade2;
  int i, j;

  zr=rand01_kiss()*taxaCuraTotal;
  aux=0;

  for (i=0; i<L*L; i++){
    aux+=taxaCura[i];
    if (zr<aux){
      aux-=taxaCura[i];
      aux+= Ni[i][i]*mi;
      if (zr<aux){
        Ni[i][i]--;
        Nr[i]++;
        cidade1=i;
        cidade2=i;
        i=L*L;
      }
      else{
        for(j=0; j<L*L;j++){
          if (i!=j){
            aux+=Ni[i][j]*mi;
          }
          if (zr<aux){
            Ni[i][j]--;
            Nr[j]++;
            cidade1=i;
            cidade2=j;
            i=L*L;
            j=L*L;
          }
        }
      }
    }
  }
  //calcular taxas que mudaram
  taxaCuraTotal-=taxaCura[cidade1];
  NiTotal[cidade1]--;
  NiTotalSoma--;
  taxaCura[cidade1]=mi*NiTotal[cidade1];
  taxaCuraTotal+=taxaCura[cidade1];
  sumTaxas=taxaCuraTotal+taxaInfeccaoTotal;

//  printf("Soma das taxas Cura: %lf %lf %lf\n", sumTaxas, taxaCuraTotal, taxaInfeccaoTotal);
}

void infectar(){
  double zr, aux;
  int cidade1, cidade2, auxC=0;
  int i, j;
  zr=rand01_kiss()*taxaInfeccaoTotal;
  aux=0;

  for (i=0; i<L*L; i++){
    aux+=taxaInfeccao[i];

    if (zr<aux){
      zr=(rand01_kiss())*NsTotal[i];
      aux=Ns[i][i];
    //  if (aux<0)
      //  printf("%lf %d %d \n", aux, Ns[i][i], i);
      if (zr<aux && Ns[i][i]>0){
        Ni[i][i]++;
        Ns[i][i]--;
        cidade1=i;
        cidade2=i;
        i=L*L;
        auxC=1;
      }
      else{
        zr=rand01_kiss()*(NsTotal[i]-Ns[i][i]);
        aux=0;
        for(j=0; j<L*L;j++){
          if (i!=j){
            aux+=Ns[i][j];
          //  if (aux<0)
        //      printf("%lf %d %d %d\n", aux, Ns[i][j], i, j);
          }
          if (zr<aux && Ns[i][j]>0){
            Ni[i][j]++;
            Ns[i][j]--;
            cidade1=i;
            cidade2=j;
            i=L*L;
            j=L*L;
            auxC=1;
          }
        }
      }
    }
  }
  //cidade1 onde ocorreu a infecção da população da cidade2
  if (aux<0 || zr<0){
    printf("%d %lf %lf\n", cidade1, aux, zr);
//    getchar();
  }
  taxaInfeccaoTotal-=taxaInfeccao[cidade1];
  NiTotal[cidade1]++;
  NsTotal[cidade1]--;
  NiTotalSoma++;
  NsTotalSoma--;
  taxaInfeccao[cidade1]=lambda*k[cidade1]*NiTotal[cidade1]*NsTotal[cidade1]/(NiTotal[cidade1]+NsTotal[cidade1]+Nr[cidade1]);
  taxaInfeccaoTotal+=taxaInfeccao[cidade1];
  sumTaxas=taxaCuraTotal+taxaInfeccaoTotal;
  if (auxC==0){
    printf("aqui\n");
  //  getchar();
  }
//  printf("Soma das taxas infectar: %lf %lf %lf\n", sumTaxas, taxaCuraTotal, taxaInfeccaoTotal);
}


void viajar(){
  int i, j, l, m, aux1, aux2, teste, teste2;
  double zr, Dijlm;

  for (i=0; i<L; i++){
    for (j=0; j<L; j++){
      aux1=Ni[matriz[i][j]][matriz[i][j]];
      aux2=Ns[matriz[i][j]][matriz[i][j]];
      for (l=0; l<L; l++){
        for (m=0; m<L; m++){
          if(i!=l || j!=m){

            Dijlm=d0/pow((abs(i-l)+abs(j-m)),a);

            Ni[matriz[l][m]][matriz[i][j]]=aux1*Dijlm;
            Ns[matriz[l][m]][matriz[i][j]]=aux2*Dijlm;
            Ni[matriz[i][j]][matriz[i][j]]-=Ni[matriz[l][m]][matriz[i][j]];
            Ns[matriz[i][j]][matriz[i][j]]-=Ns[matriz[l][m]][matriz[i][j]];
            if (Ni[matriz[i][j]][matriz[i][j]]<0 || Ns[matriz[i][j]][matriz[i][j]]<0){
              printf("AQUI: %d %d %lf %lf %d %d %d %d\n", matriz[i][j], matriz[l][m], d0, Dijlm, aux1, aux2, Ni[matriz[i][j]][matriz[i][j]], Ns[matriz[i][j]][matriz[i][j]]);
            //  getchar();
            }
          }
        }
      }
      for (l=0; l<L; l++){
        for (m=0; m<L; m++){
          if(i!=l && j!=m){
            zr=rand01_kiss();
            Dijlm=d0/pow((abs(i-l)+abs(j-m)),a);
            if (zr<=(aux1*Dijlm-floor(aux1*Dijlm)) && (Ni[matriz[i][j]][matriz[i][j]]-1)>=0){
              Ni[matriz[l][m]][matriz[i][j]]++;
              Ni[matriz[i][j]][matriz[i][j]]--;
            }
            zr=rand01_kiss();
            if (zr<=(aux2*Dijlm-floor(aux2*Dijlm)) && (Ns[matriz[i][j]][matriz[i][j]]-1)>=0){
              Ns[matriz[l][m]][matriz[i][j]]++;
              Ns[matriz[i][j]][matriz[i][j]]--;
            }
          }
        }
      }
    }
  }

  taxaCuraTotal=0;
  taxaInfeccaoTotal=0;
  NiTotalSoma=0;
  NsTotalSoma=0;
//  printf("aqui\n" );


  for (i=0; i<L*L; i++){
    NiTotal[i]=0;
    NsTotal[i]=0;
    for (j=0; j<L*L; j++){
        NiTotal[i]+=Ni[i][j];
        NsTotal[i]+=Ns[i][j];

    }
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

//  printf("Soma das taxas Viajar: %lf %lf %lf\n", sumTaxas, taxaCuraTotal, taxaInfeccaoTotal);

}

void retornar(){
  int i, j;

  for (i=0; i<L*L; i++){
    for(j=0; j<L*L;j++){
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

  for (i=0; i<L*L; i++){

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

//  printf("Soma das taxas Retornar: %lf %lf %lf\n", sumTaxas, taxaCuraTotal, taxaInfeccaoTotal);
}

void passo_tempo_sir(){
  double zr;

  double propInfectados, propCura;


  if (sumTaxas>0){
    propCura=taxaCuraTotal/(taxaCuraTotal+taxaInfeccaoTotal);
    propInfectados=taxaInfeccaoTotal/(taxaCuraTotal+taxaInfeccaoTotal);



    zr=rand01_kiss();
    if (zr<propCura){
      //printf("Cura\n");
      curar();
    }
    else{
      //printf("infecta\n" );
      infectar();
    }
  }
  else{
    printf("%d %lf %lf %lf %d %d\n", u, tempo, taxaCuraTotal, taxaInfeccaoTotal, NiTotalSoma, NsTotalSoma);
  //  getchar();
  }
}

void largura_da_curva(){
  int i, t;
  double normalizacaoAux, aux;
  for (t=0; t<TEMPO_MAXIMO; t++){
    for (i=0; i<L*L; i++){
      soma[0][i]+=(1.0*dados[2*i+1][t]);
    }
  }
  for (i=0; i<L*L; i++){
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
  FILE *arq, *arq2, *arq3, *arq4;
  char nome[50]="";
  int n, i, j, tempoViagem;
  double t, p, auxPicoEpidemia[2][L*L], d_tot;
  int auxViagem=0;

  p=1;
  f=0.4;

  //encontrar d0 em função de a
  d_tot=0.0;
  d0=0.0;
  for(i=0; i<L; i++){
    for (j=0;j<L;j++){
      if (i!=j){
        d_tot+=0.01/pow(abs(i-L/2)+abs(j-L/2),1);
        d0+=1.0/pow(abs(i-L/2)+abs(j-L/2),a);
      }
    }
  }
  d0=d_tot/d0;

  //fim


  memset(dados,0,(3)*(QUANT_DADOS)*sizeof(dados[0][0]));
  sprintf (nome, "mediaCidades-%lf-tempo%d-TESTE.dat", a, (int)TEMPO_MAXIMO);
  arq = fopen (nome, "w");

  sprintf (nome, "pico-epidemia-%lf-tempo%d.dat", a, (int)TEMPO_MAXIMO);
  arq2 = fopen (nome, "w");
  printf("Nome: %s\n", nome);
  u=0;

  memset(dados,0,(2*L*L+1)*(QUANT_DADOS)*sizeof(dados[0][0]));
  memset(picoEpidemia,0,(2)*(L*L)*sizeof(picoEpidemia[0][0]));
  memset(media, 0, QUANT_DADOS*sizeof(media[0]));
  while (u<QUANT_AMOSTRAS){
    memset(auxPicoEpidemia,0,2*(L*L)*sizeof(auxPicoEpidemia[0][0]));
    inicializacao();
    t=0.0;
    n=0;
    tempoViagem=0;
    auxViagem=0;
    printf("%d\n", u);
    do{
      if (tempo>=t){

        if (u==0){
          media[n]=1.0*NiTotalSoma/(1.0*N*L*L);
          i=1;
          for (j=0; j<L*L; j++, i+=2){
            dados[i][n]=(1.0*NiTotal[j])/(NiTotal[j]+NsTotal[j]+Nr[j]);
            dados[i+1][n]=1.0*Nr[j]/(1.0*(NiTotal[j]+NsTotal[j]+Nr[j]));
            if (auxPicoEpidemia[1][j]<(1.0*NiTotal[j])/(NiTotal[j]+NsTotal[j]+Nr[j])){
              auxPicoEpidemia[0][j]=tempo;
              auxPicoEpidemia[1][j]=(1.0*NiTotal[j])/(NiTotal[j]+NsTotal[j]+Nr[j]);
            }
          }
        }
        else{
          media[n]+=1.0*NiTotalSoma/(1.0*N*L*L);
          i=1;
          for (j=0; j<L*L; j++, i+=2){
            dados[i][n]+=(1.0*NiTotal[j])/(NiTotal[j]+NsTotal[j]+Nr[j]);
            dados[i+1][n]+=1.0*Nr[j]/(1.0*(NiTotal[j]+NsTotal[j]+Nr[j]));
            if (auxPicoEpidemia[1][j]<(1.0*NiTotal[j])/(NiTotal[j]+NsTotal[j]+Nr[j])){
              auxPicoEpidemia[0][j]=tempo;
              auxPicoEpidemia[1][j]=(1.0*NiTotal[j])/(NiTotal[j]+NsTotal[j]+Nr[j]);
            }
          }
        }
        t+=p;
        n++;

      }
      if (tempoViagem<=(tempo-1) && auxViagem==1){
        tempoViagem++;
        //printf("viagem\n");
        viajar();
        //printf("Viajar: %lf\n", tempo);
        auxViagem=0;
      }
      if (tempoViagem<=(tempo-0.5) && auxViagem==0){
        //printf("retorno\n");
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

    for (i=0; i<L*L; i++){
      picoEpidemia[0][i]+=auxPicoEpidemia[0][i];
      picoEpidemia[1][i]+=auxPicoEpidemia[1][i];
    }
  //  getchar();
    u++;
  }


  t=0.0;
  double auxiliar[2][QUANT_DADOS];


  sprintf (nome, "infectados-%lf-cidades-tempo%d.dat", a, (int)TEMPO_MAXIMO);
  arq3 = fopen (nome, "w");

  sprintf (nome, "larguraCurva-%lf-cidades-tempo%d.dat", a, (int)TEMPO_MAXIMO);
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
    for (j=1; j<(L*L*2+1); j+=2){
      dados[j][i]/=QUANT_AMOSTRAS;
      dados[j+1][i]/=QUANT_AMOSTRAS;
      fprintf(arq3, "%lf %lf ", dados[j][i], dados[j+1][i]);
      auxiliar[0][i]+=dados[j][i];
      auxiliar[1][i]+=dados[j+1][i];
    }
    fprintf (arq3, "\n");
    t+=p;
    n=i;
  //  printf("%d\n", i);
  }

  largura_da_curva();



  for (i=0; i<L*L; i++){
    picoEpidemia[0][i]/=QUANT_AMOSTRAS;
    //  printf("%lf\n", picoEpidemia[1][i]);
    picoEpidemia[1][i]/=QUANT_AMOSTRAS;
    //  printf("%lf\n", picoEpidemia[1][i]);
    //getchar();
    fprintf(arq2, "%d %lf %lf\n", distFocoIni[i], picoEpidemia[0][i], picoEpidemia[1][i]);
    fprintf(arq4, "%d %lf\n", distFocoIni[i], larguraCurva[i]);
  }

  for (i=0; i<=(n); i++){
    fprintf(arq, "%lf %lf %lf %lf", dados[0][i], auxiliar[0][i]/(1.0*L), media[i]/(1.0*QUANT_AMOSTRAS), auxiliar[1][i]/(1.0*L));
    fprintf(arq, "\n");
  }

  fclose(arq);
  fclose(arq2);

  return 0;
}
