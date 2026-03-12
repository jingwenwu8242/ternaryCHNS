This project provides the source code for the paper “Axisymmetric simulation of two-phase fluids in contact with solids.” The code is used to generate the results reported in Section 4.2 (Fig. 9) of the paper.

Authors
Jingwen Wu, Junxiang Yang*
*Corresponding author

File descriptions

In the root folder
makefile       For macOS, Linux, and WSL (Windows) systems; automates the compilation and cleaning processes.
CITATION.cff   Standard citation file for this program.

In the code folder
chnsangle.c    Main C source file for simulating axisymmetric two-phase flows in contact with solids.
main1.h        Header file containing the main variables, parameters, and function declarations.
mainutil1.h    Header file containing auxiliary utility functions and numerical routines.
show_figure1.m MATLAB script for visualizing the final 2D x-y state of two-phase flows in contact with solids under the axisymmetric assumption.

Compilation and execution

Open a terminal and run:
cd /path/to/the/code
make all
./chnsangle.out

Deleting all results

Run:
make delete

Warning: This command permanently removes all generated data files (e.g., .m files containing computed results) and figure files (e.g., .eps files) from the output directory. This operation is irreversible. Please back up any important results before proceeding.
