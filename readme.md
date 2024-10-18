Implementation of the minimal split checkerboard method adapted to unsymmetric matrices.

The main point of this Python implementation is to do edge colouring in Python where it is just a single commmand, as opposed to in another language (my use is Fortran) where I would have to manually create a graph data structure and edge colouring algorithm: the results are ready to be printed to a file and read in. See my dqmc repository for an implementation.

This Python implementation also allows for easy testing (seeing how well the approximation holds).

Note: this implementation assumes the matrix to be checkerboarded does not have a diagonal (ie, it has already been Trotterized out: exp(A) ~ exp(diag(A)) * exp(offdiag(A)))

TODO:

Implement diagonal.