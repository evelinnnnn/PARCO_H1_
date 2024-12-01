# **Homework 1: Exploring Implicit and Explicit Parallelism with OpenMP** <br>
### Evelin Begher - 235188 - evelin.begher@studenti.unitn.it <br>

I wrote the code in C, and to simplify the execution process, I created a [PBS script](matrix_transpose_H1.pbs) that allows me to run the instructions in a single step.<br>
<br>
I can compile my code in another way by directly using the processor and executing the following command: <br>
*$ gcc -o omp H1.c -fopenmp <br>
$ ./omp **n*** <br>
where ***n*** represents the matrix dimension (from $$2^4$$ to $$2^{12}$$)<br>
This procedure takes more time, especially when changing the matrix dimensions for various configurations. <br>
<br>

These two modes are implemented in my [code](matrix_transpose_H1.c) and executed accordingly. In this code we Initialize a random matrix *M* *n×n*. Then, we implemented a function to check if the matrix is symmetric (**checkSym**) . Additionally, we created a function to compute the transpose of *M* and store the result in a new matrix *T* (**matTranspose**). <br>
The assignment requires us to implement the function **checkSym** and the function **matTranspose** using three different methods: sequential implementation, implicit parallelization, and explicit parallelization using OpenMP. The goal of this assignment is to measure the time taken to perform both the symmetry check and the transpose operation using these three approaches, and to explore how we can optimize performance, especially when working with large datasets or complex calculations.

### Here there are the steps to follow to run my code: 

* Obtain the [MATRIX_TRANSPOSE_ H1_PBS](matrix_transpose_H1.pbs) and [MATRIX_TRANSPOSE_ H1.C](matrix_transpose_H1.c) files. These will be used for the job submission and running the program.

* Connect to the HPC Cluster:

    * Establish a secure VPN connection to the University of Trento network.
    * Open your SSH client and connect to the cluster

*  Enter your university login credentials 
      * Upload the Files: Once connected, navigate to your desired directory on the cluster or create a new one. Use the file transfer feature in your SSH client to upload the [MATRIX_TRANSPOSE_ H1_PBS](matrix_transpose_H1.pbs) and [MATRIX_TRANSPOSE_ H1.C](matrix_transpose_H1.c) to that directory.

* Request an Interactive Session and Reserve a Node:
     * Navigate to the directory where the files are located: cd ./directory
     * Request an interactive session on a compute node with the desired specifications (ncpus=64:ompthreads=64:mem=1mb) using the following command: qsub -I matrix_transpose_H1.pbs
     * Once the job completes, you can find the results in a file named job.o in the same directory where the job was submitted.



