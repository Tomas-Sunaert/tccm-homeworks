void compute_acc(size_t   Natoms,double** coord,double*  mass,double** distance,double** acceleration, double epsilon, double sigma){
// this subroutine will calculate the acceleration for each particle


for (int i=0;i<Natoms;i++){// cycles for all the N-particles
  double m = mass[i];




 for (int k=0;k<3;k++){// cycles for all 3 components of the xlr-ation
   double a = 0.; // preparest for the accumulation


  for (int j=0; j<Natoms; j++){// sums over all the N-particles contrib. 
   if (j==i){continue;}

   double lj = pow((sigma/distance[i][j]), 6);
   double U = 24*epsilon*(lj-2*lj*lj)/(distance[i][j]); 
   a = a + U * (coord[i][k]-coord[j][k])/distance[i][j];

  }
   acceleration[i][k]=-a/m;  
 } 
}











}
