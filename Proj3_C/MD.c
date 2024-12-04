// THIS IS FOR TESTING THE FUNCTIONS ONE BY ONE BEFORE PUTTING ALL TOGETHER

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "2DMatAlloc.h" // import the function header that allocates 2D mat in a not painful way
#include "V_LJ.h"
#include "T_energy.h"
#include "xlr8.h"
#include "Verlet.h"

#define Sigma 0.3345
#define Epsilon 0.0661




int main(){

double t_step = 0.1;
double n_step = 1000;
double M = 10;



FILE* input_file = fopen("inp.txt", "r"); // assigns the label 'input_file' to the now open file of input
size_t N; // first line of the "inp.txt" tells the program the num of atoms (thus the matrix dimension)



fscanf(input_file,"%zu\n",&N); // reads N in the file and assigns to the memory address of N
printf("N = %zu \n", N); // is everything ok?

// now we move on





double** P = malloc_2d(N,3); // I allocate a matrix with N rows and 3 columns (xyz)
double *mass = malloc(N * sizeof(double));
double** vel = malloc_2d(N,3);


for (int i=0; i<N;i++){
 for (int j=0;j<N;j++){
  vel[i][j]=0.;
 } 
}

for (int i=0; i<N; i++){
 for (int j=0; j<3;j++){
  fscanf(input_file,"%lf", &P[i][j]); // read the values into the matrix Position (P)
  if (j==2){
  fscanf(input_file, "%lf", &mass[i]);
  }
 }
}

fclose(input_file); // No need to keep the file open... so i close it!



// ---------------------------------------------------
// . . . End of the reading input part . . . 
// ---------------------------------------------------

/*
printf("Coordinates: \n"); // Pretty output 
for (int i=0; i<N; i++){
 for (int j=0; j<3;j++){
  printf("%lf\t",P[i][j]); // Now it displays for debugging the values by row
 }
 printf("\n"); // new line! (to make things prettier)
}
printf("Masses: \n");
for (int i=0; i<N;i++){
 printf("%lf\n",mass[i]);
}
*/


double ** Dist = malloc_2d(N,N); // Now I'm allocating some space for the matric 
double r ; // this is a temporary variable in which I put a pairwise distance





// I exploit the symmetry of the distance matrix with respect to the diagonal
// Thus I run the loop only on the lower triangle of this square Matrix.
// Also, obviously, the diagonal will be 0

for (int i=0; i<N; i++){
 Dist[i][i]=0.; // set the diagonal as 0
 for (int j=0; j<i; j++ ){
  double dx = (P[i][0]-P[j][0]) ; // distance on the 3 directions 
  double dy = (P[i][1]-P[j][1]) ; // this goes for every component
  double dz = (P[i][2]-P[j][2]) ; // row = atom, col = component
  
  r = sqrt(dx*dx + dy*dy + dz*dz); // now compute the square root tnx to function
  
  Dist[i][j] = r ;
  Dist[j][i] = r ; // Let's save some time by taking advantage of symmetries!
 }
}
/*
// Debug time! Let's see if everything is fine...
printf("Distance matrix:\n");
for (int i=0; i<N; i++){
 for (int j=0; j<N; j++){
  printf("%lf \t", Dist[i][j]);
 }
 printf("\n");
}
*/


// Computing the total potential in this conformation

double V_lj; // total potential variable
double T_lj; // total kinetic variable
double E_lj; // total energy varible


V_lj=V(Epsilon,Sigma,N, Dist);
// Hope this is the right number...
//printf("Potential V total is: %lf\n", V_lj);

T_lj=T(N,vel,mass);
// Hope this is the right number...
//printf("Kinetic T total is: %lf\n", T_lj);

E_lj = V_lj + T_lj;
// Hope this is the right number...
//printf("E total is: %lf\n", E_lj);


// COMPUTE THE XLR-ATION
double** xlr8tion = malloc_2d(N,3);
double** next_xlr8tion = malloc_2d(N,3); // prepare another 2D array to store the next accel, necessary to compute the updated vel.

compute_acc(N, P, mass,Dist,xlr8tion,Epsilon,Sigma);

/*
for (int i=0; i<N;i++){
 for (int j=0; j<3;j++){
  printf("%lf\t", xlr8tion[i][j]);
 }
 printf("\n");
}
*/


// VERLET PART:
pos_update(P, vel,xlr8tion, t_step, N);
compute_acc(N, P, mass,Dist,next_xlr8tion,Epsilon,Sigma);
vel_update(P,vel,xlr8tion,next_xlr8tion,t_step,N);










return 0;
}
