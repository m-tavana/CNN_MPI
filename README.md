# HPC-Course
MPI implementation

Note that implementations 1 and 3 need at least two processing units to execute.

Use ./cnn_impelementationX for program execution #X is 1,2 or 3

Use the following to compile all the three impelementation.

mpiCC -O2  -I.   -c -o conv_layer4_implementation1.o conv_layer4_implementation1.c

mpiCC -O2  -I.   -c -o conv_layer4_implementation2.o conv_layer4_implementation2.c

mpiCC -O2  -I.   -c -o conv_layer4_implementation3.o conv_layer4_implementation3.c

mpiCC -O2  -I.   -c -o utils.o utils.c

mpiCC -O2  -I. -o cnn_implementation1 conv_layer4_implementation1.o utils.o -lm

mpiCC -O2  -I. -o cnn_implementation2 conv_layer4_implementation2.o utils.o -lm

mpiCC -O2  -I. -o cnn_implementation3 conv_layer4_implementation3.o utils.o -lm





