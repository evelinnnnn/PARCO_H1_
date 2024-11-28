#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#define MIN 0.0
#define MAX 10.0
#define TOT_RUN 10

#ifdef _OPENMP
   double wt1,wt2;
#endif

// Variables for storing time
double wstart,wend;

//global variable for the speedup and the efficency (glocal because I use it in two different function)
double avg_time_sequential;
double avg_time_sequential_cs; //for the checksym

int num_threads = 0; 

double bandwidth = 0.0;

// all the function I use for my code
void initializeMatrix(float **mat, int _n);
void print(float **mat, int _n);
int checkSym(float **mat, int _n) ; 
float **matTranspose(float **mat, int _n);
int  checkSymImp(float **mat, int _n);
float **matTransposeImp(float **mat, int _n); 
int  checkSymOMP(float **mat, int _n) ;
float **matTransposeOMP(float **mat, int _n);
int checksquare(int num); 
void checktraspose(float **SEQ, float **B, int _n);
  
int main(int argc, char** argv){
  
  /*Initialize the random number generator with the current time to 
  have different random sequences on each program execution */
  srand(time(NULL)); 
   
  //take a value as input
  int n = atoi(argv[1]); 
   	
  /* Input and check for matrix dimension (square of 2)
  do {
  printf("insert number: ");
  scanf("%d", &n);
  }while (checksquare(n) == 0);*/
  
  // Allocate memory for matrix M
  float **M = (float **)malloc(n * sizeof(float *)); 
  float **TSEQ = (float **)malloc(n * sizeof(float *));
  float **TIMP = (float **)malloc(n * sizeof(float *));
  float **TOMP = (float **)malloc(n * sizeof(float *));
  for (int i = 0; i < n; i++) {
        M[i] = (float *)malloc(n * sizeof(float));
        TSEQ[i] = (float *)malloc(n * sizeof(float));
        TIMP[i] = (float *)malloc(n * sizeof(float));
        TOMP[i] = (float *)malloc(n * sizeof(float));
  } 
   
  //check for allocation success
  if ( M == NULL) {
    printf (" Memory allocation failed \n") ;
    return 1;
  } 
  
  initializeMatrix(M, n);  
  
  printf("\n%d\n", n); 
  
  /* Open CSV file for writing results
    FILE *file = fopen("H1_results.csv", "w");
    if (file == NULL) {
        perror("Failed to open file");
        return 1;
    }

   // Write header for CSV
   fprintf(file, "Num_Threads;Avg_Parallel_Time;Avg_Speedup;Avg_Efficiency\n");*/
  
//sequential
printf("\n_______________________________sequential____________________________\n"); 
 
 //check for the allocation success
  if ( TSEQ == NULL) {
    printf (" Memory allocation failed \n") ;
    return 1;
  } 
  
  if (checkSym(M, n) == 1 ){   
  
    printf("the matrix M is symmetric \n");   
  
  } else {  
    
    TSEQ = matTranspose(M, n);      
    //checktraspose(M, TSEQ, n);
      
  }   
  
//with simd  
printf("\n___________________implicit parallelization (SIMD)___________________\n"); 
   
  if ( TIMP == NULL) {
    printf (" Memory allocation failed \n") ;
    return 1;
  } 
          
  if (checkSymImp(M, n) == 1 ){
   
    printf("the matrix M is symmetric \n"); 
  
  } else { 
  
    TIMP = matTransposeImp(M, n); 
    //checktraspose(TSEQ, TIMP, n); 
  
  }
  
//with OPENMP
printf("\n___________________explicit parallelizarion (OPENMP)__________________\n"); 
    
  if ( TOMP == NULL) {
    printf (" Memory allocation failed \n") ;
    return 1;
  } 
  
  if (checkSymOMP(M, n) == 1) {
  
    printf("the matrix M is symmetric \n"); 
  
  } else {
    TOMP = matTransposeOMP(M, n); 
    //checktraspose(TSEQ, TOMP, n); 
  }
   
   //__________________________________________________________________________________
  // Free the allocated memory
  for (int i = 0; i < n; i++) {
    free(M[i]);
    free(TSEQ[i]);
    free(TIMP[i]);
    free(TOMP[i]);
  }
  free(M); 
  free(TSEQ);
  free(TIMP);
  free(TOMP);
  
  //fclose(file);
  
  return 0;
}

//inizializate the matrix
void initializeMatrix(float **mat, int _n) {  
  for (int i = 0; i < _n; i++) {
  	for(int j = 0; j < _n; j++ ) {
     //insert randomic number in the matrix
  		mat[i][j] = MIN + (float)rand()/(float) (RAND_MAX / (MAX-MIN)); 
  	}
  }
}

//print the matrix
void print(float **mat, int _n) {
  for (int i = 0; i < _n; i++) {
  	for(int j = 0; j < _n; j++ ) {
      //.1 to specify only one digit after the decimal point 
 		  printf("%.1f  ", mat[i][j]);
  	}
  	printf("\n"); 
  }
}

//control if the matrix is symmetric 
int checkSym(float **mat, int _n) {
  
  double total_time_sequential_cs= 0.0;
  double time_sequential_cs = 0.0; 
  int is_symmetric= 1; 
  
  for (int run = 0; run < TOT_RUN; run++) {
  
  //Start time
    wstart = omp_get_wtime();
    
    for (int i = 0; i < _n; i++) {
      for(int j = 0; j < _n; j++ ) {
        if (mat[i][j] != mat[j][i]){
          is_symmetric = 0; 
        } 
      }
    }
    
    //End time 
    wend = omp_get_wtime();
    
    time_sequential_cs = wend - wstart;
    total_time_sequential_cs += time_sequential_cs;
  }
  
  avg_time_sequential_cs = total_time_sequential_cs / TOT_RUN;
  printf("Execution times of the checksym routine: %12.4g seconds\n", avg_time_sequential_cs);
  
  bandwidth = ((2 * (_n *_n) * sizeof(float))/avg_time_sequential_cs)/ 1e9;
  printf("Execution times of the bandwidth: %.3g seconds\n", bandwidth);
  
  return is_symmetric; 		
}

//calculate the traspose of the matrix 
float **matTranspose(float **mat, int _n) {

  double total_time_sequential= 0.0; 
  double time_sequential = 0.0; 
   
  //create a temporary matrix 
  float **temp = (float **)malloc(_n * sizeof(float *)); 
  for (int i = 0; i < _n; i++) {
        temp[i] = (float *)malloc(_n * sizeof(float));
  }; 
  
  for (int run = 0; run < TOT_RUN; run++) {
  
    // Start time
    wstart = omp_get_wtime();
  
    //calculate the transpose    	
    for (int i = 0; i < _n; i++) {
    	for(int j = 0; j < _n; j++ ) {
    	  temp[i][j] = mat[j][i];	
    	}
    }
  
    //End time
    wend = omp_get_wtime(); 
  
    time_sequential = wend - wstart;
    total_time_sequential += time_sequential;
  }
  
  //caclcualate the average of the execution times of the transposition routine check 
  //to have more relevant value
  avg_time_sequential = total_time_sequential / TOT_RUN;
  printf("Execution times of the transposition routine check: %f seconds\n", avg_time_sequential);
  
  bandwidth = ((2 * (_n *_n) * sizeof(float))/avg_time_sequential)/ 1e9;
  printf("Execution times of the bandwidth: %.3g seconds\n", bandwidth);

  return temp;
  
  //free the allocated memory
  for (int i = 0; i < _n; i++) {
        free(temp[i]);
  }
  free(temp); 
}

//control if the matrix is symmetric
#pragma GCC optimized ("O2")
int  checkSymImp(float **mat, int _n) { 
  
  double total_time_implicit_cs= 0.0;
  double time_implicit_cs = 0.0; 
  int is_symmetric = 1; 
  
  for (int run = 0; run < TOT_RUN; run++) {
  
    //Start time
    wstart = omp_get_wtime();
    
    #pragma simd reduction(&:is_symmetric)
    for (int i = 0; i < _n; i++) {
        for (int j = 0; j < i; j++) {
            if (mat[i][j] != mat[j][i]) {
                is_symmetric = 0;
            }
        }
    }
    
    //End time 
    wend = omp_get_wtime();
    
    time_implicit_cs = wend - wstart;
    total_time_implicit_cs += time_implicit_cs;
  }
  
  double avg_time_implicit_cs = total_time_implicit_cs / TOT_RUN;
  printf("Execution times of the checksym routine: %12.4g seconds\n", avg_time_implicit_cs);
  
  bandwidth = ((2 * (_n *_n) * sizeof(float))/avg_time_implicit_cs)/ 1e9;
  printf("Execution times of the bandwidth: %.3g seconds\n", bandwidth);
  
  return is_symmetric;		
}

//calculate the traspose of the matrix
#pragma GCC optimize ("O2")
float **matTransposeImp(float **mat, int _n) {

  double total_time_implicit= 0.0;
  double time_implicit = 0.0; 
   
  //create a temporary matrix 
  float **temp = (float **)malloc(_n * sizeof(float *)); 
  for (int i = 0; i < _n; i++) {
        temp[i] = (float *)malloc(_n * sizeof(float));
    };
    
    
    for (int run = 0; run < TOT_RUN; run++) {
          
      // Start time
      wt1 = omp_get_wtime(); 
          
       #pragma simd 	
      for (int i = 0; i < _n; i++) {
      	for(int j = 0; j < _n; j++ ) {
      	  temp[i][j] = mat[j][i];	
      	}
      }      
  
      // end time
      wt2 = omp_get_wtime();
      
      time_implicit =  wt2 - wt1;
      total_time_implicit += time_implicit;
          
    }
              
  double avg_time_implicit = total_time_implicit / TOT_RUN;
  printf("Execution times of the transposition routine check: %f seconds\n", avg_time_implicit);
  
  bandwidth = ((2 * (_n *_n) * sizeof(float))/avg_time_implicit)/ 1e9;
  printf("Execution times of the bandwidth: %.3g seconds\n", bandwidth);
   
  return temp; 
  
  for (int i = 0; i < _n; i++) {
        free(temp[i]);
  }
  free(temp); 
}

//control if the matrix is symmetric with OPENMP 
int  checkSymOMP(float **mat, int _n) {

  // variable for sppedup and efficency
  double avg_time_parallel_cs=0.0;
  double total_time_parallel_cs; 
  double avg_speedup_cs=0.0;
  double avg_efficiency_cs=0.0;
  double time_parallel_cs = 0.0; 

  int is_symmetric = 1;
  int i, j;  
  
  printf("\n CHECKSYM OMP \n"); 
  
  printf(" Num_Threads | Avg_Parallel_Time | Avg_Speedup | Avg_Efficiency\n"); 
    
   for (num_threads = 1; num_threads <= 64; num_threads *= 2) {
    omp_set_num_threads(num_threads);
    double total_time_parallel_cs = 0.0;
    
    for (int run = 0; run < TOT_RUN; run++) {
    
      #ifdef _OPENMP
      // start time
      wt1 = omp_get_wtime();
           
      #pragma omp parallel for shared(mat, is_symmetric) private(i, j)
      for (i = 0; i < _n; i++) {
          for (j = 0; j < i; j++) {
              if (mat[i][j] != mat[j][i]) {    
                  is_symmetric = 0;
              }
          }
      }     
    
      //end time
      wt2 = omp_get_wtime();
      
      #endif 
  
      time_parallel_cs =  wt2 - wt1;
      total_time_parallel_cs += time_parallel_cs;
      
    }
    
    avg_time_parallel_cs = total_time_parallel_cs / TOT_RUN;
    avg_speedup_cs = avg_time_sequential_cs / avg_time_parallel_cs;
    avg_efficiency_cs = avg_speedup_cs / num_threads;
    bandwidth = ((2 * (_n *_n) * sizeof(float))/avg_time_parallel_cs)/ 1e9;
  
    #ifdef _OPENMP
    printf("%11d  | %17f | %11f | %13.2f%% | %14.3f\n", num_threads, avg_time_parallel_cs, avg_speedup_cs, avg_efficiency_cs * 100, bandwidth); 
    #endif
  
  }

  return is_symmetric; 		
}

//calculate the traspose of the matrix using OPENMP
float **matTransposeOMP(float **mat, int _n) {
  
  // variable for sppedup and efficency
  double avg_time_parallel=0.0;
  double total_time_parallel; 
  double avg_speedup=0.0;
  double avg_efficiency=0.0;
  double time_parallel = 0.0;  
  
  // Allocate memory for the transposed matrix
  float **temp = (float **)malloc(_n * sizeof(float *)); 
  for (int i = 0; i < _n; i++) {
    temp[i] = (float *)malloc(_n * sizeof(float));
  }
  
  printf("\n MATTRANSPOSE OMP \n"); 
  printf(" Num_Threads | Avg_Parallel_Time | Avg_Speedup | Avg_Efficiency | Bandwidth (GB/s)\n");
  
    // Test performance with thread counts of 1, 2, 4, 8, 16, 32, and 64
    
    for (num_threads = 1; num_threads <= 64; num_threads *= 2) {
      omp_set_num_threads(num_threads);
      double total_time_parallel = 0.0;
      
      for (int run = 0; run < TOT_RUN; run++) {
      
        #ifdef _OPENMP
        // start time
        wt1 = omp_get_wtime();
        #endif
         
        #ifdef _OPENMP
        // Parallelize the matrix transpose
        #pragma omp parallel for collapse(2)  
        for (int i = 0; i < _n; i++) {
          for (int j = 0; j < _n; j++) {
            temp[i][j] = mat[j][i];
          }
        } 
        #endif 
        
        #ifdef _OPENMP
        //end time
        wt2 = omp_get_wtime();
        #endif
    
        time_parallel =  wt2 - wt1;
        total_time_parallel += time_parallel;
      
    }
  
    avg_time_parallel = total_time_parallel / TOT_RUN;
    avg_speedup = avg_time_sequential / avg_time_parallel;
    avg_efficiency = avg_speedup / num_threads;
    bandwidth = ((2 * (_n *_n) * sizeof(float))/avg_time_parallel)/ 1e9;
    
    //printf("Execution times of the transposition routine check: %f seconds\n", avg_time_parallel );
    #ifdef _OPENMP
    printf("%11d  | %17f | %11f | %13.2f%% | %14.3f\n", num_threads, avg_time_parallel, avg_speedup, avg_efficiency * 100, bandwidth);
    #endif
    
  
  }
  
  /* Write results to CSV
  fprintf(file, "%d; %f; %f; %f\n", num_threads, avg_time_parallel, avg_speedup, avg_efficiency * 100);
  }*/
  
   

  
  return temp; 
  
  for (int i = 0; i < _n; i++) {
    free(temp[i]);
  }
  free(temp); 
}

//check if dimention of the matrix is a power of two
int checksquare(int num) {

  if (num == 0) {
  return 0; }
  
  while (num%2 == 0) {
  num= num/2;} 
  
  return (num==1); 
}

//check if the transpose of the sequential code is equal to the other 
void checktranspose(float **SEQ, float **B, int _n) {

  int check = 1; 
  
  for (int i = 0; i < _n; i++) {
 	  for(int j = 0; j < _n; j++ ) {
  		if (SEQ[i][j] != B[i][j]){
  			check = 0;  
  		} 
  	}
  }
  
  if (check == 1) {
    printf("\nequal\n"); 
  } else {
    printf("\nnot equal\n"); 
  }

}
