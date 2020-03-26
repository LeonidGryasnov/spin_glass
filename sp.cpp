#include <stdio.h>
#include <stdlib.h>
//#include <conio.h>
#include <string.h>
#include <math.h>
#include <time.h>

constexpr const int L=6;//Lineyniy razmer sistemi
constexpr const double D0=1;//Koncentraciya spinov orientirovannih "vverh"
constexpr const int spin_quantity=L*L*L;//Kolichestvo spinov v reshetke
constexpr const int monte_karlo_steps=6000;//Kolichestvo shagov Monte-Karlo
constexpr const int skipped_monte_karlo_steps=3000;//Kolichestvo nachalnopropuskaemih shagov
constexpr const int averaging=2;//Kolichestvo usredneniy po odnoy temperature
constexpr const double T_hide = 5;
constexpr const double T = 1.15;
constexpr const double P=0.5;
///-------------------------------------
constexpr const int steps=monte_karlo_steps-skipped_monte_karlo_steps;
double Xsg[steps];
///-------------------------------------
int J[L][L][L][L][L][L];//Matrica vzaimodeystviya
///-------------------------------------
int Spin[2][L][L][L];//Reshetka spinov
//Peremennie dlya generatora sluchaynih chisel
int mult1, mult2;
int modm1, modm2;
double rn[256];
int ibm1,ibm2;
int xx1,xx2;
///-------------------------------------
int Gener_Init()
//Inicializatiya generatora sluchaynih chisel
{
     //Initializiruem zerna generatora
     xx1=1235*time(NULL);
     xx2=6356;
     mult1=16807;
     mult2=65539;
     modm1=2147483647;
     modm2=2147483647;

     ibm1=2*xx1+1;

     for (auto i=0;i<256;i++)
     {
          ibm1*=mult1;
          if (ibm1<0) ibm1+=(modm1+1);
          rn[i]=ibm1/(modm2*1.0);
     }
     return 0;
}
///-------------------------------------
///*
int Init_Dinamic()
//Inicializatiya dinamicheskih peremennih
{
     for (auto i=0;i<steps;i++)
     {
        Xsg[i]=0;
     }

     return 0;
}
///-------------------------------------
int Gran(int i)
//Periodicheskie granichnie ucloviya
{
     if (i==(-1)) i=L-1;
     if (i==L) i=0;
     return i;
};
///-------------------------------------
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
///-------------------------------------
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
///-------------------------------------
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
///-------------------------------------
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
///-------------------------------------
int FlipSpin(double T)
//Perevorot spina
{
     int i,j,k,z,is;
     int dE;//Energiya ot perevorota spina
     double w;//Veroyatnost perevorota

     for (is=0;is<spin_quantity;is++)
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
///-------------------------------------
int Dinamic_Value(int it)
//Podschet dinamicheskih velichin
{
     double dValue=0;

     for (auto i=0;i<L;i++)
     {
          for (auto j=0;j<L;j++)
          {
               for (auto k=0;k<L;k++)
               {
                    dValue+=(Spin[1][i][j][k]*Spin[0][i][j][k])/(2.0*spin_quantity);
               }
          }
     }
     Xsg[it]=spin_quantity*((dValue*dValue)/steps);
     //Xsg[it]=1;
     return 0;
};
///-------------------------------------
int Dinamic_Otvet(int im)
{
     FILE *fp;

     char s[30];

	 sprintf(s,"%d - Xsg(t,L).dat",im);

	 fp=fopen(s,"w");

     //Zapolnyaem tablicu
     for (auto i=0;i<steps;i++)
     {
          fprintf(fp,"%i",i+1);
          fprintf(fp,"\t%lf",Xsg[i]);
          fprintf(fp,"\n");
     }

     fclose(fp);

     for (auto i=0;i<steps;i++)
     {
         Xsg[i]=0;
     }

     return 0;
}
///-------------------------------------

int main()
{
    Gener_Init();//Inicializatiya generatora sluchaynih chisel

    for (auto im=0;im<averaging;im++) {
        Init_Dinamic();//Inicializatiya dinamicheskih peremennih
        Init_J();//Inicializiruem matricu obmennogo vzaimodeystviya
        Init_Spin();//Inicializiruem nachalnoe pologenie spinov

        for (auto in=0;in<monte_karlo_steps;in++) {
            /*if (in<skipped_monte_karlo_steps)
                    {
                         FlipSpin(T_hide);
                         //printf("- \n");
                    } else {
                        FlipSpin(T);
                        Dinamic_Value(in);
                        //printf("%d \n", Xsg[in]);
                    }
                    //Proizvodim izmereniya
                    */
                     FlipSpin(T_hide);
                    if (in>skipped_monte_karlo_steps)
                    {
                        FlipSpin(T);
                        Dinamic_Value(in-skipped_monte_karlo_steps);
                    }
        }
        Dinamic_Otvet(im);
        printf("%d average\n",im);
    }
    return 0;
}
