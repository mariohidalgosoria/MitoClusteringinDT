#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#define Nds 21 //number of rows in each file  id.dat
char  nom[50];//name of each file sprintf        
FILE *tm0;

int main(void){
   
  int id,i,Ndat;
  double Ser_X[Nds],xser;
  double Ser_t[Nds];
  
  int ini,Nexp;
   
  ini  =  1;
  Nexp = 10;//number of time series
  for(i=0;i<Nds;i++){
    Ser_X[i]=0;
    Ser_t[i]=0;
  }

   
  for(id=ini;id<=Nexp;id++){          
   
    i=0;
    //Captura de datos:
    sprintf(nom,"XDEUG15T%d.dat",id);//coloca el id, xid.dat en una cadena 'nom'
    printf("%s\n",nom);//imprime en la pantalla xid.dat
    tm0 = fopen(nom,"r");//abre xid.dat para lectura, creando el puntero tm0
    while((!feof(tm0))){//mientras no sea el final de xid.dat
      fscanf(tm0,"%lf",&Ser_t[i]);//lee el primer elemento de xid.dat un float y la guarda en Ser_t[i]
      //printf("t=%lf\n",Ser_t[i]);
      fscanf(tm0,"%lf",&xser);//lee el segunda elemento de xid.dat un float y la guarda en xser
      //printf("MDP=%lf\n",xser);
      //fscanf(tm0,"%lf",&dser);//lee el tercer elemento de xid.dat un int y la guarda en mser
      //printf("MDD=%lf\n",dser);
      //fscanf(tm0,"%lf",&dtser);//lee el cuarto elemento de xid.dat un int y la guarda en mser
      //printf("MDDT=%lf\n",dtser);
      Ser_X[i]+=xser;
      i++;
    }// final del while
    fclose(tm0);
             
  }//final del for
    
  Ndat = i; //da el numero de filas en cada archivo leido
  printf("Numero de datos: %d\n", Ndat);
      
  sprintf(nom,"EXDEUG15.dat");//la cadena nom='x_mean.dat'
  tm0 = fopen(nom,"w"); //se crea el archivo x_mean.dat
      
  for(i=0;i<Ndat-1;i++){
    fprintf(tm0,"%lf   %lf\n",Ser_t[i],Ser_X[i]/Nexp);
  }
  printf("Numero de datos: %d\n ", i);
  fclose(tm0);     
  return 0;
} 

