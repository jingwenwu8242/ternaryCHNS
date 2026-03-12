This project provides the source code for the paper, Axisymmetric simulation of the two-phase fluids in contact with solids. 
The code is used to generate the results shown in Sections 4.2 (Fig. 9) of the paper.
Author
Jingwen Wu, Junxiang Yang*   
*Corresponding Author
File Descriptions
In the root folder
`makefile`:  For macOS, Linux, WSL (Windows) systems, automates the compilation and cleaning process.
`CITATION.cff`: Standard citation format file for this program.
In 4_3_4, 4_3_5 (corresponding to the section number in the paper respectively)
`*.c` & `*.h`: C++ source and header files for the calculation.
`show_figure.m`: MATLAB script to visualize the final 2D x-y state of two-phase flows in contact with solids under the axisymmetric assumption.

Compilation and Execution
MacOS/Linux/WSL (Windows)
Compiling and running the code
Open a terminal and run the following commands:
cd /path/to/the/code
make all
./chnsangle.out

File descriptions
chnsangle.c      Main C source file for the computation of axisymmetric two-phase flows in contact with solids.
main1.h          Header file containing the main variables, parameters, and function declarations.
mainutil1.h      Header file containing utility functions and auxiliary numerical routines.
show_figure1.m   MATLAB script for visualizing the 2D x-y results based on the axisymmetric assumption.

Deleting all the results
Input `make delete` in your terminal.
Warning: Executing this command will permanently remove all generated data files (e.g., `.m` files containing calculated results) and figure files (e.g., `.eps` files) from the output directory. This operation is irreversible. Please ensure you have backed up any critical results before proceeding.
