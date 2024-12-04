// THIS IS FOR TESTING THE FUNCTIONS ONE BY ONE BEFORE PUTTING ALL TOGETHER

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "2DMatAlloc.h"
//#include <V_LJ.h>




int main(){
FILE* test_input = fopen("input_mat.txt", "r");

double** A= malloc_2d(5,3);

for (int i=0; i<5; i++){
for (int j=0; j<3;j++){
fscanf(test_input,"%lf", &A[i][j]);
}}



printf("Matrix A:\n");
for (int i=0; i<5; i++){
for (int j=0; j<3;j++){
printf("%lf\t",A[i][j]);
}
printf("\n");
}





return 0;
}
