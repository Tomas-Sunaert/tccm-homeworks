The file provided cimputes some classical MD simulation given an input coordinate for N atoms. 
Parameters and potential are optimized for Ar atoms with a shallow LJ potential.
The initial velocity is set to 0.

In the src folder you'll find all the source codes to check out, as well as the headers for the program, thus breaking up the problem in smaller, more digestable, functions.

The test run is in the "test" folder, where you'll find a 2 input files to try "inp.txt" and "inp2.txt", as well as some outputs of a previous run. This includes some files of the energies, accelerations, velocities or relative distances, that can be used both for debugging if you want to implement new features or to have a deeper look inside the MD.

In the file acceleration.txt you can find in a xyz-like format the mass of each atom followed by the vector components.

In the vel.txt file the velocities are shown for each step (no mass included).

In the r_ij.txt file the relative distances between different atoms is shown for each step.

In TV_E.txt you'll find 4 columns, the first is the current step, the second is the total kinetic energy, the third is the total potential and the last on is the sum of the previous two. Use this to check the validity of the time step.





- Jacopo
