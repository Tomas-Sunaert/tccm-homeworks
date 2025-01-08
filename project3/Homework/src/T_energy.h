double T(size_t Natoms, double** velocity, double* mass){

double KE=0.; // same trick! Initialize this to 0, so that I can accumulate the total KE inside it

for (int i = 0; i<Natoms; i++){
 for (int j = 0; j<3; j++){
 KE+= 0.5 * mass[i] * velocity[i][j] * velocity[i][j];
}}

return KE;



}
