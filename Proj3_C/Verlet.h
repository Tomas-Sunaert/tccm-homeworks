void pos_update(double** coord,double** velox,double** acceleration, double t_step, double Natoms){

for (int i = 0; i<Natoms; i++){
 for (int j = 0; j<3; j++){
  coord[i][j]=coord[i][j] + velox[i][j] * t_step + 0.5 * acceleration[i][j] * t_step * t_step;
 }
}
}

void vel_update(double** coord,double** velox,double** acceleration, double** next_acceleration, double t_step, double Natoms){
for (int i = 0; i<Natoms; i++){
 for (int j = 0; j<3; j++){
  velox[i][j]=velox[i][j] + 0.5 * (acceleration[i][j] + next_acceleration[i][j]) * t_step;
 }
}

}
