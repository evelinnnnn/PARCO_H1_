# **Homework 1: Exploring Implicit and ExplicitParallelism with OpenMP** <br>
### Evelin Begher - 235188 - evelin.begher@studenti.unitn.it <br>

I wrote the code in C, and to simplify the execution process, I created a [PBS script](code.pbs) that allows me to run all the necessary instructions in a single step.
<br>
I can compile my code in another way by directly using the processor and executing the following command: <br>
*$ gcc -o omp H1.c -fopenmp <br>
$ export OMP_NUM_THREADS=**n_thread**; ./omp **n*** <br>
where ***n_thread*** represents the number of threads to be used (in our case 2, 4, 8, or 16) and ***n*** represents the matrix dimension (from $$2^4$$ to $$2^{12}$$)<br>
This procedure takes significantly more time, especially when using multiple threads and changing the matrix dimensions for various configurations. <br>
<br>
These two modes are implemented in my [code](code.c) and executed accordingly. In this code we Initialize a random matrix *M* *n×n*. Then, we implemented a function to check if the matrix is symmetric (**checkSym**) . Additionally, we created a function to compute the transpose of *M* and store the result in a new matrix *T* (**matTranspose**). <br>
The assignment requires us to implement the function **checkSym** and the function **matTranspose** using three different methods: sequential implementation, implicit parallelization, and explicit parallelization using OpenMP. The goal of this assignment is to measure the time taken to perform both the symmetry check and the transpose operation using these three approaches, and to explore how we can optimize performance, especially when working with large datasets or complex calculations.

