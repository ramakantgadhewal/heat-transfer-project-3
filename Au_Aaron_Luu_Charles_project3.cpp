/****************************** MECH 8120 Heat Transfer, Fall 2020******************************
Project 3 - Forced Convection Steady State Laminar Flow with Constant Wall Temperature
 Purpose: To determine the local and average Nusselt numbers for hydrodynamically
          fully-developed, thermally developing, steady state, laminar flow in a pipe
 Author(s):      Aaron Au,    Charles Luu
 Student ID(s):  A00935196,   A00903725
Declaration:
We, Aaron Au and Charles Luu, declare that the following program was written by us.
Date Created: November 10, 2020
Revised:
***********************************************************************************************/
//--------------------------------Standard Libraries-------------------------------- 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <float.h>

//--------------------------------Global Constants-------------------------------- 
const double twentyP = 0.20;       //20%
const double fortyP = 0.40;        //40%
const double sixtyP = 0.60;        //60%
const double eightyP = 0.80;       //80%
const double hundyP = 1.0;         //100%
const double thetaWall = 1.0;      //BC at pipe wall
const double dtdr = 0.0;           //BC at pipe center; represents d(theta)/dr
const int axialNodes = 12001;      //predetermined number of nodes in axial direction

//--------------------------------Structure Definitions--------------------------------
typedef struct PIPEDATA
{
   int J;                  //number of nodes in radial direction
   int L;                  //length of the heated section
   double Re, Pr;          //Reynolds number and Prandtl number
   double tMean, NuMean;   //dimensionless mean temperature and mean Nusselt number
}PIPEDATA;

typedef struct PIPEPOINT
{
   double u, r;               //local velocity and radial position
   double a, b, c, y;         //a, b, c are for the coefficient matrix; y is the result
   double theta;              //dimensionless temperature
   double P1, P2, P3, P4, P5; //temperature profiles
}PIPEPOINT;

//--------------------------------Function Prototypes--------------------------------
void flushInputBuffer();                //flushes input buffer after scanf_s() call
void finiteDiff(PIPEDATA *);            //finite difference solver
void Tridiagonal(int, PIPEPOINT *);     //solve tridiagonal system of equations

//--------------------------------Main--------------------------------
int main()
{
   PIPEDATA PD; //declare PD variable for inputs

   //prompt user for input
   printf("How many radial nodes, Gandalf?: ");
   scanf_s("%d", &PD.J);
   flushInputBuffer();

   printf("\nWhat is the length of the heated section (in # of diameters), Gandalf?: ");
   scanf_s("%d", &PD.L);
   flushInputBuffer();

   printf("\nWhat about the Reynolds number, Gandalf?: ");
   scanf_s("%lf", &PD.Re);
   flushInputBuffer();

   printf("\nWhat about the Prandtl number, Gandalf?: ");
   scanf_s("%lf", &PD.Pr);
   flushInputBuffer();

   //function call to finiteDiff solver
   finiteDiff(&PD);

   //prompt to end program
   printf("\n\nPress ENTER to close the program...\n\n");
   (void)getchar();

   return 0;
}
//--------------------------------Function Definitions--------------------------------

//--------------------------------flushInputBuffer-------------------------------- 
// DESCRIPTION: This function flushes the input buffer to avoid scanf issues
// ** CALL THIS FUNCTION AFTER EVERY CALL TO SCANF!!! **
// ARGUMENTS: none
// RETURN VALUE: none

void flushInputBuffer()
{
   char ch; // temp variable
   if(feof(stdin)) return; // don't do a getchar if stdin is empty!
      // exit loop when all characters are flushed
   while((ch = getchar()) != '\n' && ch != EOF) {}
}

//--------------------------------FiniteDiff-------------------------------- 
//DESCRIPTION: This function calculates the dimensionless temperature at each node
//             as a system of equations
// RETURN VALUE: none
void finiteDiff(PIPEDATA *PD)
{
   PIPEPOINT *PP;      //pointer to PIPEPOINT structure
   int I, J;           //number of nodes in axial and radial directions, respectively
   int i, j;           //counter variables

   double Pe;          //Peclet number
   double dx, dr;      //step sizes in axial and radial directions, respectively
   double L;           //length of heated section
   double sum;         //running total for trapezoid rule

   Pe = (PD->Re) * (PD->Pr);  //initialize Peclet number
   L = PD->L;                 //initialize length of heated section
   I = axialNodes;            //number of axial nodes
   J = PD->J;                 //initialize number of radial nodes 
   dx = 2.0 * L / (I - 1.0);  //initialize axial step size
   dr = 1.0 / (J - 1.0);      //initialize radial step size

   FILE *Nusselt;             //pointer to file for Nusselt plot
   FILE *Profiles;            //pointer to file for profile plots
   FILE *Mean;                //pointer to file for Mean temp. plot
   errno_t ERR;               //checkfor error

   //attempt to open file for writing
   ERR = fopen_s(&Nusselt, "Nusselt.dat", "w");
   ERR = fopen_s(&Profiles, "Profile.dat", "w");
   ERR = fopen_s(&Mean, "Mean.dat", "w");


   //print error message if unable to open file
   if(ERR != 0 || Nusselt == NULL || Profiles == NULL || Mean == NULL)
   {
      printf("\nCannot open file for writing...");
      return;
   }

   //Allocating memory for all PIPEPOINT arrays and Nusselt array
   PP = (PIPEPOINT *)malloc(J * sizeof(PIPEPOINT));
   double *Nu = (double *)malloc(I * sizeof(double));

   //check for errors in allocating memory
   if(PP == NULL || Nu == NULL)
   {
      printf("Cannot allocate memory for ... Press ENTER to exit program");
      return;
   }

   //initializing Nusselt numbers in axial direction to zero
   for(i = 0; i < I; i++)
   {
      Nu[i] = 0.0;
   }

   //initializing radial nodes and velocities
   for(j = 0; j < J; j++)
   {
      PP[j].r = j * dr;
      PP[j].u = 2.0 * (1.0 - pow(PP[j].r, 2.0));
   }

   //initialize PP arrays to zero
   for(j = 0; j < J; j++)
   {
      PP[j].a = 0.0;
      PP[j].b = 0.0;
      PP[j].c = 0.0;
      PP[j].y = 0.0;
      PP[j].theta = 0.0;
      PP[j].P1 = 0.0;
      PP[j].P2 = 0.0;
      PP[j].P3 = 0.0;
      PP[j].P4 = 0.0;
      PP[j].P5 = 0.0;
   }

   i = 0; //reset counter for do-while loop

   do
   {
      //reset sum and tMean for next axial step
      sum = 0.0;
      PD->tMean = 0.0;


      //assigning values to first elements of arrays
      PP[0].a = -3.0 / 2.0 / dr;   //a coefficient
      PP[0].b = 4.0 / 2.0 / dr;    //b coefficient
      PP[0].c = -1.0 / 2.0 / dr;   //c coefficient
      PP[0].y = dtdr;              //boundary condition at pipe center

      //assigning values to last elements of arrays
      PP[J - 1].a = 1.0;           //a coefficient
      PP[J - 1].b = 0.0;           //b coefficient
      PP[J - 1].c = 0.0;           //c coefficient
      PP[J - 1].y = thetaWall;     //boundary condition at pipe wall

      //calculate interior elements of arrays
      for(j = 1; j < J - 1; j++)
      {
         PP[j].a = (PP[j].u / dx) + (2.0 / Pe / pow(dr, 2.0) / PP[j].r) *
            ((PP[j].r + dr / 2.0) + (PP[j].r - dr / 2.0));

         PP[j].b = -2.0 * (PP[j].r + dr / 2.0) / Pe / pow(dr, 2.0) / PP[j].r;
         PP[j].c = -2.0 * (PP[j].r - dr / 2.0) / Pe / pow(dr, 2.0) / PP[j].r;
         PP[j].y = PP[j].u * PP[j].theta / dx;
      }

      //Function call to Tridiagonal
      Tridiagonal(J, PP);


      //Store theta values from y array
      for(j = 0; j < J; j++)
      {
         PP[j].theta = PP[j].y;
      }

      //storing temperature profiles at 5 equally spaced points
      if((double)i / (axialNodes - 1.0) == twentyP)   //20% of pipe length
      {
         for(j = 0; j < J; j++)
         {
            PP[j].P1 = PP[j].y;
         }
      }
      else if((double)i / (axialNodes - 1.0) == fortyP)  //40% of pipe length
      {
         for(j = 0; j < J; j++)
         {
            PP[j].P2 = PP[j].y;
         }
      }

      else if((double)i / (axialNodes - 1.0) == sixtyP)  //60% of pipe length
      {
         for(j = 0; j < J; j++)
         {
            PP[j].P3 = PP[j].y;
         }
      }
      else if((double)i / (axialNodes - 1.0) == eightyP) //80% of pipe length
      {
         for(j = 0; j < J; j++)
         {
            PP[j].P4 = PP[j].y;
         }
      }
      else if((double)i / (axialNodes - 1.0) == hundyP) //100% of pipe length
      {
         for(j = 0; j < J; j++)
         {
            PP[j].P5 = PP[j].y;
         }
      }

      //calculate mean dimensionless temperature from trapezoid rule
      for(j = 0; j < J; j++)
      {
         sum += PP[j].theta * PP[j].u * PP[j].r + PP[j + 1].theta * PP[j + 1].u * PP[j + 1].r;
      }

      PD->tMean = 2.0 * dr / 2.0 * sum;

      //calculate local Nusselt number
      Nu[i] = 2.0 / (PP[J - 1].theta - PD->tMean) / dr * ((3.0 / 2.0 * PP[J - 1].theta) -
         (4.0 / 2.0 * PP[J - 2].theta) + (1.0 / 2.0 * PP[J - 3].theta));

      //write Nusselt numbers and mean temps. to .dat files
      fprintf(Mean, "%lf, %lf\n", i * dx / 2.0, PD->tMean);
      fprintf(Nusselt, "%lf, %lf\n", i * dx / 2.0, Nu[i]);

      i++; //increment counter by 1
   }
   while(i < I);

   //write temperature profiles to .dat file
   for(j = 0; j < J; j++)
   {
      fprintf(Profiles, "%lf, %lf, %lf, %lf, %lf\n", PP[j].P1, PP[j].P2, PP[j].P3, PP[j].P4, PP[j].P5);
   }

   //reset sum to calculate mean Nusselt number
   sum = 0.0;

   //calculate mean Nusselt number
   for(i = 0; i < I; i++)
   {
      sum += Nu[i];
   }

   PD->NuMean = sum / (double)I;

   printf("\nNu mean: %lf", PD->NuMean);

   //close file after writing
   fclose(Nusselt);
   fclose(Profiles);
   fclose(Mean);

   //Free memory
   free(PP);
   free(Nu);
}

//--------------------------------Tridiagonal-------------------------------- 
//DESCRIPTION: Modified Thomas algorithm to solve tridiagonal 
//             system of equations
// RETURN VALUE: none
//
// The non zero terms are:
//    a  b  c
//    c  a  b
//       c  a  b
//          c  a  b
//          b  c  a
void Tridiagonal(int I, PIPEPOINT *PP)
{
   int i;

   // Process first row
   PP[0].b /= PP[0].a;
   PP[0].c /= PP[0].a;

   // Adjust b[1]
   PP[1].b -= PP[1].c * PP[0].c;

   // Process rows 2 to n-1
   for(i = 1; i < (I - 1); i++)
   {
      PP[i].a -= PP[i].c * PP[i - 1].b;
      PP[i].b /= PP[i].a;
   }

   // Last row
   PP[I - 1].c -= PP[I - 1].b * PP[I - 3].b;
   PP[I - 1].a -= PP[I - 1].c * PP[I - 2].b;

   // Forward solution step
   PP[0].y /= PP[0].a;

   for(i = 1; i < (I - 1); i++)
   {
      PP[i].y = (PP[i].y - (PP[i].c * PP[i - 1].y)) / PP[i].a;
   }

   PP[I - 1].y = (PP[I - 1].y - (PP[I - 1].c * PP[I - 2].y) - (PP[I - 1].b * PP[I - 3].y)) / PP[I - 1].a;

   // Backward solution step
   for(i = I - 2; i > 0; i--)
   {
      PP[i].y -= PP[i].b * PP[i + 1].y;
   }

   PP[0].y -= (PP[0].b * PP[1].y) + (PP[0].c * PP[2].y);

   return;
}