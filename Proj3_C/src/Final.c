/// THIS IS FOR TESTING THE FUNCTIONS ONE BY ONE BEFORE PUTTING ALL TOGETHER

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "2DMatAlloc.h" // import the function header that allocates 2D mat in a not painful way
#include "V_LJ.h"
#include "T_energy.h"
#include "xlr8.h"
#include "Verlet.h"

#define Sigma 0.3345
#define Epsilon 0.0661

//#define Sigma 1
//#define Epsilon 1


// THIS IS THE SAME AS LAST ITERATION, I'M PREPARING THE GLOBAL VARIABLES, MOVING THEM IN THE UPPER PART SO THAT DURING THE LOOP
// THEY WILL NOT BE CREATED AGAIN AND AGAIN (OR WORST, ALLOCATED.








int main(){
// Time related simulation variables
double t_step = 0.2;
int n_step = 1000;
int M = 10;


char file_name[32];

printf("Available .txt files:\n");
system("ls *.txt");


printf("\n");
printf("Input file name: ");
scanf("%99s", file_name);


FILE* input_file = fopen(file_name, "r"); // assigns the label 'input_file' to the now open file of input
FILE* traj = fopen("output.xyz", "w");
FILE* Energies = fopen("TV_E.txt","w");
FILE* speed = fopen("vel.txt", "w");
FILE* ac = fopen("accel.txt","w");
FILE* distances = fopen("r_ij.txt", "w");

size_t N; // first line of the "inp.txt" tells the program the num of atoms (thus the matrix dimension)


// read the number of particles in the input
fscanf(input_file,"%zu\n",&N); // reads N in the file and assigns to the memory address of N
//printf("N = %zu \n", N); // is everything ok?




// Now that we know with how many part. we are working with: allocations of all the upcoming useful arrays
double** P = malloc_2d(N,3); // I allocate a matrix with N rows and 3 columns (xyz)
double *mass = malloc(N * sizeof(double));
double** Dist = malloc_2d(N,N); // Now I'm allocating some space for the matric 
double** vel = malloc_2d(N,3);
double** xlr8tion = malloc_2d(N,3);
double** next_xlr8tion = malloc_2d(N,3); // prepare another 2D array to store the next accel, necessary to compute the updated vel.

// initiate the vel vector as 0 for all the particles
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



double r ; // this is a temporary variable in which I put a pairwise distance



for (int step = 0; step<n_step; step++){

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


// Computing the total potential in this conformation

double V_lj; // total potential variable
double T_lj; // total kinetic variable
double E_lj; // total energy varible


V_lj = V(Epsilon,Sigma,N, Dist);
T_lj = T(N,vel,mass);
E_lj = V_lj + T_lj;


// COMPUTE THE XLR-ATION
compute_acc(N, P, mass,Dist,xlr8tion,Epsilon,Sigma); // if this is the first step, then i compute the initial accel


// compute the remainder of the number of current stev and M (I want to save only every M steps)
int modulo = step % M;

if (modulo == 0){
fprintf(Energies,"%d\t%lf\t %lf\t %lf\n",step,T_lj, V_lj, E_lj );
fprintf(traj,"\t%zu\n",N);
fprintf(traj,"\n");

fprintf(speed,"\t%zu\n",N);
fprintf(speed,"\n");

fprintf(distances,"\t%zu\n",N);
fprintf(distances,"\n");

fprintf(ac,"\t%zu\n",N);
fprintf(ac,"\n");

for (int i=0; i<N; i++){
fprintf(traj,"Ar\t");
fprintf(ac,"%lf\t", mass[i]);
for (int j=0; j<3; j++){
fprintf(traj,"%lf\t",P[i][j]);
fprintf(speed,"%lf\t",vel[i][j]);
fprintf(ac,"%lf\t",xlr8tion[i][j]);
}
for (int j=0;j<N;j++){
fprintf(distances,"%lf\t",Dist[i][j]);
}
fprintf(traj,"\n");
fprintf(speed,"\n");
fprintf(distances,"\n");
fprintf(ac,"\n");



}


}



// VERLET PART: first update position, compute the accelleration at the next step (with updated pos), use the previous and current accel to fint the new vel
pos_update(P, vel,xlr8tion, t_step, N);
compute_acc(N, P, mass,Dist,next_xlr8tion,Epsilon,Sigma);
vel_update(P,vel,xlr8tion,next_xlr8tion,t_step,N);



}






// Free the memory (like a chill guy)
free_2d(P);
free(mass); 
free_2d(Dist); 
free_2d(vel); 
free_2d(xlr8tion);
free_2d(next_xlr8tion);




printf("==========================\n");
printf("Job Done!");

return 0;
}
