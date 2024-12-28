Here you have 2 .f90 files.



The first is Heatmap.f90.

This program takes as input the input files in the data folder 
(that has fixed path with respect to the .f90 file, so don't move it!) and creates a "Dense" file with the numbers in a 2D array, 
so that gnuplot can read that and create a heat map showing the matrix elements in a nicer way.

Also, there is a .gnuplot file that has all the commands and features of the interactive heatmap, don't delete this.

The second one is Sparse_MatMul.f90.

This reads from name file 2 matrices, checks for errors and then computes the product using the 3 1D array.
At the end it also creates a "dense" matrix rapresentation saved in "My_MaMul". 
In the"MaMul" file you can find the result from the standart MatMul() intrinsic function of F90 to compare.

