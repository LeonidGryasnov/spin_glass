/*Programka, modeliruet 3D Isingovskoe spinovoe steklo,
algoritmom Metropolisa
*/
//Versiya 1.03_no_Object_MPI
#include <stdio.h>
#include <stdlib.h>
//#include <conio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <mpi.h>
//Parametri modeli
const int L=6;//Lineyniy razmer sistemi
const double D0=1;//Koncentraciya spinov orientirovannih "vverh"
const int N=L*L*L;//Kolichestvo spinov v reshetke
const int NK=3000;//Kolichestvo shagov Monte-Karlo
const int NN=3000;//Kolichestvo nachalnopropuskaemih shagov
const int NM=2;//Kolichestvo usredneniy po odnoy temperature
//const double T0=1.15;//Nachalnaya temperatura
//const double dT=0.01;//Shag po temperature
//const double Tmax=1.16;//Konechnaya temperatura
const double T1 = 5;
const double T2 = 1.15;
const double P=0.5;//Koncentraciya vzaimodeystviya
//1 - ferromagnetik, 0 - antiferromagnetik, 0.5 - spin glass
//-------------------------------------
double D[2][3000];//Distanciya Hamminga
double D1[2][3000];//Distanciya Hamminga
//double D[NT][NK];//Distanciya Hamminga
//-------------------------------------
//const int NT=(int)((Tmax-T0)/dT);//Kolichestvo shagov po temperature
const int NT = 1;
//-------------------------------------
//Peremennie dlya generatora sluchaynih chisel
int mult1, mult2;
int modm1, modm2;
double rn[256];
int ibm1,ibm2;
int xx1,xx2;
//-------------------------------------
int Gener_Init()
//Inicializatiya generatora sluchaynih chisel
{
     int i;
     //Initializiruem zerna generatora
     xx1=1235*time(NULL);
     xx2=6356;
     mult1=16807;
     mult2=65539;
     modm1=2147483647;
     modm2=2147483647;

     ibm1=2*xx1+1;

     for (i=0;i<256;i++)
     {
          ibm1*=mult1;
          if (ibm1<0) ibm1+=(modm1+1);
          rn[i]=ibm1/(modm2*1.0);
     }
     return 0;
}
//-------------------------------------

//-------------------------------------
///*
int Init_Dinamic()
//Inicializatiya dinamicheskih peremennih
{
     int i,j;

     for (i=0;i<NT;i++)
     {
          for (j=0;j<NK;j++)
          {
               D[i][j]=0;
          }
     }

     return 0;
}
//*/
//-------------------------------------
int J[L][L][L][L][L][L];//Matrica vzaimodeystviya
//-------------------------------------
int Gran(int i)
//Periodicheskie granichnie ucloviya
{
     if (i==(-1)) i=L-1;
     if (i==L) i=0;
     return i;
};
//-------------------------------------
int Spin[2][L][L][L];//Reshetka spinov
//-------------------------------------
double Gener()
//Generatiya sluchaynogo chisla
{
     double r;
     int k;
     ibm2=2*xx2+1;
     ibm2*=mult2;
     if (ibm2<0) ibm2+=(modm1+1);
     k=ibm2/8388608+1;
     r=rn[k];
     ibm1*=mult1;
     if (ibm1<0) ibm1+=(modm1+1);
     rn[k]=ibm1/(modm2*1.0);
     return r;
}
//-------------------------------------
int Init_J()
//Initializatiya matrici vzaimodeystviya
{
     int i,j,k,i1,j1,k1;

     //Obnulyaem matricu obmennogo vzaimodeystviya
     for (i=0;i<L;i++)
     {
          for (j=0;j<L;j++)
          {
               for (k=0;k<L;k++)
               {
                    for (i1=0;i1<L;i1++)
                    {
                         for (j1=0;j1<L;j1++)
                         {
                              for (k1=0;k1<L;k1++)
                              {
                                   J[i][j][k][i1][j1][k1]=0;
                                   //J[i][j][k][i1][j1][k1]=-1;//Proverka - sie 100% ferromagnetik
                              }
                         }
                    }
               }
          }
     }

     //Probegaem vse spini
     for (i=0;i<L;i++)
     {
          for (j=0;j<L;j++)
          {
               for (k=0;k<L;k++)
               {
                    //Ischem sosedey u spina
                    //Sosed "sleva"
                    if (J[i][j][k][Gran(i-1)][j][k]==0)
                    {
                         if (Gener()>P)
                         {
                              J[i][j][k][Gran(i-1)][j][k]=+1;
                         }
                         else
                         {
                              J[i][j][k][Gran(i-1)][j][k]=-1;
                         };
                         J[Gran(i-1)][j][k][i][j][k]=J[i][j][k][Gran(i-1)][j][k];
                    };
                    //Sosed "sprava"
                    if (J[i][j][k][Gran(i+1)][j][k]==0)
                    {
                         if (Gener()>P)
                         {
                              J[i][j][k][Gran(i+1)][j][k]=+1;
                         }
                         else
                         {
                              J[i][j][k][Gran(i+1)][j][k]=-1;
                         }
                         J[Gran(i+1)][j][k][i][j][k]=J[i][j][k][Gran(i+1)][j][k];
                    };
                    //Sosed "speredi"
                    if (J[i][j][k][i][Gran(j-1)][k]==0)
                    {
                         if (Gener()>P)
                         {
                              J[i][j][k][i][Gran(j-1)][k]=+1;
                         }
                         else
                         {
                              J[i][j][k][i][Gran(j-1)][k]=-1;
                         }
                         J[i][Gran(j-1)][k][i][j][k]=J[i][j][k][i][Gran(j-1)][k];

                    };
                    //Sosed "szadi"
                    if (J[i][j][k][i][Gran(j+1)][k]==0)
                    {
                         if (Gener()>P)
                         {
                              J[i][j][k][i][Gran(j+1)][k]=+1;
                         }
                         else
                         {
                              J[i][j][k][i][Gran(j+1)][k]=-1;
                         };
                         J[i][Gran(j+1)][k][i][j][k]=J[i][j][k][i][Gran(j+1)][k];
                    };
                    //Sosed "snizu"
                    if (J[i][j][Gran(k-1)][i][j][k]==0)
                    {
                         if (Gener()>P)
                         {
                              J[i][j][Gran(k-1)][i][j][k]=+1;
                         }
                         else
                         {
                              J[i][j][Gran(k-1)][i][j][k]=-1;
                         };
                         J[i][j][k][i][j][Gran(k-1)]=J[i][j][Gran(k-1)][i][j][k];
                    };
                    //Sosed "Sverhu"
                    if (J[i][j][Gran(k+1)][i][j][k]==0)
                    {
                         if (Gener()>P)
                         {
                              J[i][j][Gran(k+1)][i][j][k]=+1;
                         }
                         else
                         {
                              J[i][j][Gran(k+1)][i][j][k]=-1;
                         };
                         J[i][j][k][i][j][Gran(k+1)]=J[i][j][Gran(k+1)][i][j][k];
                    };
               }
          }
     }
     return 0;
}
//-------------------------------------
int Init_Spin()
{
     //Initializaton reshetki spinov
     {
          int i,j,k,z,q;
          z=0;

          //Po vsem reshetkam
          for(i=0;i<L;i++)
          {
               for(j=0;j<L;j++)
               {
                    for(k=0;k<L;k++)
                    {
                         if (D0==1)
                         //Odna reshetka zapolnena spinami "vverh", dr. spinami "vniz"
                              {
                                   Spin[0][i][j][k]=+1;
                                   Spin[1][i][j][k]=-1;
                                   if ((i==0)&&(j==0)&&(k==0)) z=1;
                              }
                         if (D0==0.5)
                         //Obe reshetki zapolneni spinami "vverh"-"vniz" ravnoveroyatnim obrazom
                              {
                                   for (q=0;q<2;q++)
                                   {
                                        if (Gener()<0.5) Spin[q][i][j][k]=+1;
                                        else Spin[q][i][j][k]=-1;
                                   }
                                   if ((i==0)&&(j==0)&&(k==0)) z=1;
                              }
                         if (D0==0)
                         //Obe reshetki zapolneni identichno, za isklucheniem centralnogo spina
                         {
                              Spin[0][i][j][k]=+1;
                              Spin[1][i][j][k]=+1;
                         }
                    }
               }
          }
          if (D0==0)
          {
               i=L/2;
               j=L/2;
               k=L/2;
               Spin[0][i][j][k]=-Spin[1][i][j][k];
               z=1;
          }
          if (z!=1)
          {
               printf("\n-Devushka, a, devushka. A mogno s vami poznokomitsya. Menya Baba-Vanya zovut\n");
               printf("-Nu i dura!\n");
               return 1;
          }
          else return 0;
     }
}
//-------------------------------------
//Staticheskie peremennit
double M,M2;//Namagnichennost & ee kvadrat
double Q,Q2;//Parametr poryadka Parizi & ee kvadrat
double E,E2;//Energiya & ee kvadrat
//-------------------------------------
int Init_Static(int im, int it, int myid)
//Inicializatiya staticheskih peremennih
{
     FILE *fp;

     //Obnulyaem peremennie
     M=0;
     M2=0;
     Q=0;
     Q2=0;
     E=0;
     E2=0;

     //Fail s nachalnimi usloviyami
     if ((im==0) && (im==0) && (myid==0))
     //Proverka fayla na suchestvovanie
     {
          //Stali bit fayla net
          fp=fopen("Experiment.xls","w");//Sozdayem noviy fail
          fprintf(fp,"L\t%i\n",L);
          fprintf(fp,"D0\t%lf\n",D0);
          fprintf(fp,"NK\t%i\n",NK);
          fprintf(fp,"NN\t%i\n",NN);
          fprintf(fp,"NM\t%i\n",NM);
          fprintf(fp,"T0\t%.3lf\n",T0);
          fprintf(fp,"dT\t%.3lf\n",dT);
          fprintf(fp,"Tmax\t%.3lf\n",Tmax);
          fprintf(fp,"P\t%.1lf\n",P);
          fclose(fp);
     }
     return 0;
}
//-------------------------------------
int Sum(int z, int i, int j, int k)
//Summa sosedey s uchetom obmennogo integrala
{
     int s=0;
     //Sosed "sleva"
     s+=Spin[z][Gran(i-1)][j][k]*J[Gran(i-1)][j][k][i][j][k];
     //Sosed "sprava"
     s+=Spin[z][Gran(i+1)][j][k]*J[Gran(i+1)][j][k][i][j][k];
     //Sosed "speredi"
     s+=Spin[z][i][Gran(j-1)][k]*J[i][Gran(j-1)][k][i][j][k];
     //Sosed "szadi"
     s+=Spin[z][i][Gran(j+1)][k]*J[i][Gran(j+1)][k][i][j][k];
     //Sosed "snizu"
     s+=Spin[z][i][j][Gran(k-1)]*J[i][j][Gran(k-1)][i][j][k];
     //Sosed "sverhu"
     s+=Spin[z][i][j][Gran(k+1)]*J[i][j][Gran(k+1)][i][j][k];

     return s;
};
//-------------------------------------
int FlipSpin(double T)
//Perevorot spina
{
     int i,j,k,z,is;
     int dE;//Energiya ot perevorota spina
     double w;//Veroyatnost perevorota

     for (is=0;is<N;is++)
     //Daem vozmognost perevernutsya
     {
          //Vibiraem sluchayniy spin
          i=(int) ((L)*Gener());
          j=(int) ((L)*Gener());
          k=(int) ((L)*Gener());

          if (i == L) i--;
          if (j == L) j--;
          if (k == L) k--;

          //Probuem ego oprokinut
          for (z=0;z<2;z++)
          //Po oboim reshetkam
          {
               dE=-2*Spin[z][i][j][k]*Sum(z,i,j,k);//Energiya ot perevorota odnogo spina
               w=Gener();
               if (dE<0)
               //Esli perevorot energeticheski vigoden, to perevorachivaem spin
               {
                    Spin[z][i][j][k]=-Spin[z][i][j][k];
               }
               else
               {
                    if (w<exp(-dE/(T*1.0)) )
                    {
                         Spin[z][i][j][k]=-Spin[z][i][j][k];
                    }
               }
          }
     }

     return 0;
}
//-------------------------------------
///*
int Dinamic_Value(int it, int in)
//Podschet dinamicheskih velichin
{
     double dValue=0;
     int i,j,k;

     for (i=0;i<L;i++)
     {
          for (j=0;j<L;j++)
          {
               for (k=0;k<L;k++)
               {
                    dValue+=(Spin[1][i][j][k]*Spin[0][i][j][k])/(2.0*N);
               }
          }
     }
     D[it][in]=N*((dValue*dValue)/NK);
     return 0;
};
//*/
//-------------------------------------
int Static_Value()
//Poschet staticheskih peremennih
{
     int i,j,k;
     double dValue=0;

     dValue=0;
     //Podschet namagnichennosti
     for (i=0;i<L;i++)
     {
          for (j=0;j<L;j++)
          {
               for (k=0;k<L;k++)
               {
                    dValue+=Spin[0][i][j][k];
               }
          }
     }

     M+=fabs(dValue);
     M2+=dValue*dValue;

     dValue=0;

     //Podschet parametra poriyadka Parizi
     for (i=0;i<L;i++)
     {
          for (j=0;j<L;j++)
          {
               for (k=0;k<L;k++)
               {
                    dValue+=Spin[1][i][j][k]*Spin[0][i][j][k];
               }
          }
     }
     Q+=fabs(dValue);
     Q2+=dValue*dValue;

     dValue=0;

     //Podschet energii sistemi
     for (i=0;i<L;i++)
     {
          for (j=0;j<L;j++)
          {
               for (k=0;k<L;k++)
               {
                    dValue+=Spin[0][i][j][k]*Sum(0,i,j,k);
               }
          }
     }
     E+=dValue;
     E2+=dValue*dValue;

     dValue=0;
     return 0;
}
//-------------------------------------
int Otvet_Static(int im, int it,int numproc,int myid,int n_proc,double T)
//Vivod otveta v fail
{
     FILE *fp;
     int i;
     char name[25];
     char sqrname[30];
	 double Xi;//Vospriimchivost, ili teploemkost dlya dannoy velichini

	 //----------------------------------------------------------
	 //M

     strcpy(name,"M(T,L).xls");
     strcpy(sqrname,"sqr");
     strcat(sqrname,name);

     //Usrednenie velichin
     M=M/(1.0*(NK-NN)*N);
     M2=M2/(1.0*(NK-NN)*N);

     if ((im==n_proc*myid)&&(it==0))
     {
          fp=fopen(name,"w");
          fprintf(fp,"im\\T");
          //Zabivaem shapku
          for(i=0;i<NT;i++)
          {
               fprintf(fp,"\t%lf",T0+dT*i);
          }
     }
     else fp=fopen(name,"a");

     if (it==0)
     {
          fprintf(fp,"\n%i",im+1);
     }
     //Zapolnyaem tablicu
     fprintf(fp,"\t%lf",M);

     fclose(fp);
	 //---------------------------
	 //Xi_M

	 strcpy(name,"Xi_M(T,L).xls");

	 Xi=0;
	 Xi=(M2 - M*M*N)/(1.0*T);

	 if ((im==n_proc*myid)&&(it==0))
     {
          fp=fopen(name,"w");
          fprintf(fp,"im\\T");
          //Zabivaem shapku
          for(i=0;i<NT;i++)
          {
               fprintf(fp,"\t%lf",T0+dT*i);
          }
     }
     else fp=fopen(name,"a");

     if (it==0)
     {
          fprintf(fp,"\n%i",im+1);
     }
     //Zapolnyaem tablicu
     fprintf(fp,"\t%lf",Xi);

     fclose(fp);

     //---------------------------
	 //sqrM

     if ((im==n_proc*myid)&&(it==0))
     {
          fp=fopen(sqrname,"w");
          fprintf(fp,"im\\T");
          //Zabivaem shapku
          for(i=0;i<NT;i++)
          {
               fprintf(fp,"\t%lf",T0+dT*i);
          }
     }
     else fp=fopen(sqrname,"a");

     if (it==0)
     {
          fprintf(fp,"\n%i",im+1);
     }
     //Zapolnyaem tablicu
     fprintf(fp,"\t%lf",M2);

     fclose(fp);

     M=0;
     M2=0;

     //-----------------------------------------------------------
	 //Q

     strcpy(name,"Q(T,L).xls");
     strcpy(sqrname,"sqr");
     strcat(sqrname,name);

     //Usrednenie velichin
     Q=Q/(1.0*(NK-NN)*N);
     Q2=Q2/(1.0*(NK-NN)*N);

     if ((im==n_proc*myid)&&(it==0))
     {
          fp=fopen(name,"w");
          fprintf(fp,"im\\T");
          //Zabivaem shapku
          for(i=0;i<NT;i++)
          {
               fprintf(fp,"\t%lf",T0+dT*i);
          }
     }
     else fp=fopen(name,"a");

     if (it==0)
     {
          fprintf(fp,"\n%i",im+1);
     }
     //Zapolnyaem tablicu
     fprintf(fp,"\t%lf",Q);

     fclose(fp);

	 //---------------------------
	 //Xi_Q

	 strcpy(name,"Xi_Q(T,L).xls");

	 Xi=0;
	 Xi=(Q2 - Q*Q*N)/(1.0*T);

	 if ((im==n_proc*myid)&&(it==0))
     {
          fp=fopen(name,"w");
          fprintf(fp,"im\\T");
          //Zabivaem shapku
          for(i=0;i<NT;i++)
          {
               fprintf(fp,"\t%lf",T0+dT*i);
          }
     }
     else fp=fopen(name,"a");

     if (it==0)
     {
          fprintf(fp,"\n%i",im+1);
     }
     //Zapolnyaem tablicu
     fprintf(fp,"\t%lf",Xi);

     fclose(fp);

     //---------------------------
	 //sqrQ

     if ((im==n_proc*myid)&&(it==0))
     {
          fp=fopen(sqrname,"w");
          fprintf(fp,"im\\T");
          //Zabivaem shapku
          for(i=0;i<NT;i++)
          {
               fprintf(fp,"\t%lf",T0+dT*i);
          }
     }
     else fp=fopen(sqrname,"a");

     if (it==0)
     {
          fprintf(fp,"\n%i",im+1);
     }
     //Zapolnyaem tablicu
     fprintf(fp,"\t%lf",Q2);

     fclose(fp);

     Q=0;
     Q2=0;

     //-----------------------------------------------------------
	 //E

     strcpy(name,"E(T,L).xls");
     strcpy(sqrname,"sqr");
     strcat(sqrname,name);

     //Usrednenie velichin
     E=E/(1.0*(NK-NN)*N);
     E2=E2/(1.0*(NK-NN)*N);

     if ((im==n_proc*myid)&&(it==0))
     {
          fp=fopen(name,"w");
          fprintf(fp,"im\\T");
          //Zabivaem shapku
          for(i=0;i<NT;i++)
          {
               fprintf(fp,"\t%lf",T0+dT*i);
          }
     }
     else fp=fopen(name,"a");

     if (it==0)
     {
          fprintf(fp,"\n%i",im+1);
     }
     //Zapolnyaem tablicu
     fprintf(fp,"\t%lf",E);

     fclose(fp);

	 //---------------------------
	 //C

	 strcpy(name,"C(T,L).xls");

	 Xi=0;
	 Xi=(E2 - E*E*N)/(1.0*T*T);

	 if ((im==n_proc*myid)&&(it==0))
     {
          fp=fopen(name,"w");
          fprintf(fp,"im\\T");
          //Zabivaem shapku
          for(i=0;i<NT;i++)
          {
               fprintf(fp,"\t%lf",T0+dT*i);
          }
     }
     else fp=fopen(name,"a");

     if (it==0)
     {
          fprintf(fp,"\n%i",im+1);
     }
     //Zapolnyaem tablicu
     fprintf(fp,"\t%lf",Xi);

     fclose(fp);

     //---------------------------
	 //sqrE

     if ((im==n_proc*myid)&&(it==0))
     {
          fp=fopen(sqrname,"w");
          fprintf(fp,"im\\T");
          //Zabivaem shapku
          for(i=0;i<NT;i++)
          {
               fprintf(fp,"\t%lf",T0+dT*i);
          }
     }
     else fp=fopen(sqrname,"a");

     if (it==0)
     {
          fprintf(fp,"\n%i",im+1);
     }
     //Zapolnyaem tablicu
     fprintf(fp,"\t%lf",E2);

     fclose(fp);

     E=0;
     E2=0;

     return 0;
};
//-------------------------------------
///*
int Dinamic_Otvet(int im,int myid)
{
     FILE *fp;
     int i,j;

     //char name[25];
     char s[30];

     //strcpy(name,"D(T,L).xls");
     //itoa(im*(myid+1),s,10);
	 sprintf(s,"%d - D(T,L).dat",im*(myid+1));

     //strcat(s," - ");
     //strcat(s,name);

     //Usrednenie velichin

     fp=fopen(s,"w");
     //Zaplonyaem shapku
     fprintf(fp,"in\\T");
     //for (i=0;i<NT;i++)
     //{
          fprintf(fp,"\t%lf", T2 /*T0+dT*i*/);
     //}
     fprintf(fp,"\n");
     //Zapolnyaem tablicu
     for (j=0;j<NK;j++)
     {
          fprintf(fp,"%i",j+1);
          for (i=0;i<NT;i++)
          {
               fprintf(fp,"\t%lf",D1[i][j]);
          }
          fprintf(fp,"\n");
     }

     fclose(fp);

    /* for (i=0;i<NT;i++)
     {
          for (j=0;j<NK  ;j++)
          {
               D[i][j]=0;
          }
     } */

     return 0;
}
//*/
//-------------------------------------
int main(int argc, char **argv)
{
	//-----------------------------------------

	int numproc;//Kolichestvo procov
	int myid;//Identificator procov
	int n_proc;//Kolichestvo usredneniy na proc

	MPI_Status status();//Poneslas mocha po trubam
	MPI_Init(&argc,&argv);

	MPI_Comm_size (MPI_COMM_WORLD,&numproc);//Schitaem proci
	MPI_Comm_rank (MPI_COMM_WORLD,&myid);//Daem kagdomu klichku "Motya"

	n_proc=NM/numproc;//Vichislyaem skolko dolgen poschitat kagdiy proc

	if (n_proc*numproc<NM)
	//Na vsyakiy pogarniy
	{
		n_proc=n_proc+1;
	}
	//----------------------------------------

     int it,im,in;
     //double T;//Temperatura
      for (auto i=0;i<NT;i++)
        {
            for (auto j=0;j<NK  ;j++)
            {
                D1[i][j]=0;
            }
        }

     Gener_Init();//Inicializatiya generatora sluchaynih chisel

     for (im=n_proc*myid;im<n_proc*(myid+1);im++)
     //Po usredneniyam
     {
          Init_Dinamic();//Inicializatiya dinamicheskih peremennih

          //Inicializiruem velichini
         // for (it=0;it<NT;it++)
          //Po temperaturam
          //{
               Init_J();//Inicializiruem matricu obmennogo vzaimodeystviya
               Init_Spin();//Inicializiruem nachalnoe pologenie spinov

              // T=T0+dT*it;
              // if (myid==0) printf("\tT=%.3lf\n",T);

              // Init_Static(im,it,myid);

               for (in=0;in<(NK+NN);in++)
               //Shagi Monte-Karlo
               {
                   if (in<NN)
                    {
                         FlipSpin(T1);
                    } else {
                        FlipSpin(T2);
                    }

                    //Proizvodim izmereniya


                    if (in>NN)
                    {
                        Dinamic_Value(it,in);
                         //Static_Value();
                    }
               }//Po shagam Monte-Karlo
               //Otvet_Static(im,it,numproc,myid,n_proc,T);
          //}//Po temperaturam
         // Dinamic_Otvet(im,myid);
         for (auto i=0;i<NT;i++)
        {
            for (auto j=0;j<NK  ;j++)
            {
                D1[i][j]+=D[i][j]/NM;
            }
        }

          for (auto i=0;i<NT;i++)
        {
            for (auto j=0;j<NK  ;j++)
            {
                D[i][j]=0;
            }
        }
		  if (myid==0) printf("--------------------%.3lf\n",((im+1)*(myid+1)*100.0)/(n_proc*(myid+1)*1.0));
     }//Po usredneniyam
    Dinamic_Otvet(im,myid);
     //printf("\a");//Podiem, truba zovet! Otstavit bainki

	printf("\a");//Podiem, truba zovet! Otstavit bainki
	if (myid==0) printf("\nFinite l`a comedeya!\n");

	MPI_Finalize();
	//Konec paralelnoy chasti
	/*
	if (myid==0)
	{
		//getch();
	}
	*/
     return 0;//
}
