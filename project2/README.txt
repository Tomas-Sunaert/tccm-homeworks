So, you know the drill.

You have the data folder, that NEEDS TO STAY THERE (I configured the path with respect to the src).

The code, when it's running, will display the current "data" folder that i just mentioned, so you can choose the matrix you want. 

Then it will generate a dense matrix representation (2D matrix) in the file "Dense_matrix". 

You need to have on your local gnuplot because from the code i call the bash terminal and execute the file .gnupot that generates an interactive interface to visualize the data in a heat map format.



HOW DOES THE CODE OPERATE?

1) Lists from bash shell the content of ../data 
2) After choosing the file, il counts the number of lines in total and allocates the 3 arrays for columns, rows and values with the correct size
3) Now finally it read all the file from the top filling the 3 arrays element by element
4) Since it could happen that the matrix could be a block matrix, it read the max vale between the column and row indices to then allocate the smalles matrix possible to contain the data (not completely necessary)
5) With the i-th index of row in "Row" array and the i-th in "Col" array it fills the 2D matrix with the i-th "Val" element, taking care of the symmetry since Dense(i,j) = Dense(j,i)
6) Calls the bash shell and executes the gnuplot commands.







Have fun!
