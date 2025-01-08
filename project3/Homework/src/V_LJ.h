// this is a prototype of function for Lennard Johnes potential


// epsilon and sigma are given from the main program


double V(double epsilon,double sigma,size_t Natoms,double** distance){

double result=0.; // Just initializing the variable as zero, so that i can later accumulate the potential in it

for (int i=0; i<Natoms; i++){
 for (int j=0; j<i; j++){
  double lj = pow((sigma/distance[i][j]),6); // i call only once the function per loop, the other term i just compute it as the square
  result = result + 4 * epsilon * (lj*lj-lj); // in this way the result will be stored and accumulated here
 }
}



return result;


}










