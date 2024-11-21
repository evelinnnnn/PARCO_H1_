**Homework 1: Exploring Implicit and ExplicitParallelism with OpenMP** <br>
Evelin Begher - 235188 - evelin.begher@studenti.unitn.it <br>

I wrote the code in C, and to simplify the execution process, I created a PBS script that allows me to run all the necessary instructions in a single step.[[PBS code](code.pbs)]
<br>
I can compile my code in another way by directly using the processor and executing the following command: <br>
*$ gcc -o omp H1.c -fopenmp <br>
$ export OMP_NUM_THREADS=**n_thread**; ./omp **n*** <br>
where ***n_thread*** represents the number of threads to be used (in our case 2, 4, 8, or 16) and ***n*** represents the matrix dimension (from $$2^4$$ to $$2^{12}$$)<br>
This procedure takes significantly more time, especially when using multiple threads and changing the matrix dimensions for various configurations. <br>
This procedure enables the execution of the program in parallel, leveraging multiple threads to optimize performance for larger matrix sizes.
