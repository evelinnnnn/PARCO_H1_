#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#define MIN 0.0
#define MAX 10.0

#ifdef _OPENMP
   double wt1,wt2;
#endif

// Variables for storing time
  double wstart,wend;

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
  		printf("%.1f  ", mat[i][j]);
  	}
  	printf("\n"); 
  }
}

//control if the matrix is symmetric 
int checkSym(float **mat, int _n) {

  for (int i = 0; i < _n; i++) {
  	for(int j = 0; j < _n; j++ ) {
  		if (mat[i][j] != mat[j][i]){
  			return 0; 
  		} 
  	}
  }
  
  return 1; 		
}

//calculate the traspose of the matrix 
float **matTranspose(float **mat, int _n) {
   
  //create a temporary matrix 
  float **temp = (float **)malloc(_n * sizeof(float *)); 
  for (int i = 0; i < _n; i++) {
        temp[i] = (float *)malloc(_n * sizeof(float));
    };
    
    // Start time
  wstart = omp_get_wtime(); 
     	
  for (int i = 0; i < _n; i++) {
  	for(int j = 0; j < _n; j++ ) {
  	temp[i][j] = mat[j][i];	
  	}
  }
  
  // End time
  wend = omp_get_wtime();
  
  printf("Execution times of the transposition routine check: %12.4g seconds\n", wend - wstart);

  return temp;
  
  for (int i = 0; i < _n; i++) {
        free(temp[i]);
  }
  free(temp); 
}

//control if the matrix is symmetric
#pragma GCC optimized ("O2")
int  checkSymImp(float **mat, int _n) {

  #pragma simd
  int is_symmetric = 1;
    for (int i = 0; i < _n; i++) {
        for (int j = 0; j < i; j++) {
            if (mat[i][j] != mat[j][i]) {
                is_symmetric = 0;
            }
        }
    }
  
  return is_symmetric; 		
}

//calculate the traspose of the matrix
#pragma GCC optimize ("O2")
float **matTransposeImp(float **mat, int _n) {

  // Start time
  wstart = omp_get_wtime(); 
   
  //create a temporary matrix 
  float **temp = (float **)malloc(_n * sizeof(float *)); 
  for (int i = 0; i < _n; i++) {
        temp[i] = (float *)malloc(_n * sizeof(float));
    };
   
  #pragma simd 	
  for (int i = 0; i < _n; i++) {
  	for(int j = 0; j < _n; j++ ) {
  	temp[i][j] = mat[j][i];	
  	}
  }
  
  // End time
  wend = omp_get_wtime();

  printf("Execution times of the transposition routine check: %12.4g seconds\n", wend - wstart);

  return temp; 
  
  for (int i = 0; i < _n; i++) {
        free(temp[i]);
  }
  free(temp); 
}

//control if the matrix is symmetric with OPENMP 
int  checkSymOMP(float **mat, int _n) {
  int is_symmetric = 1;
  int i, j; 
  
  #ifdef _OPENMP
  #pragma omp parallel for shared(mat, is_symmetric) private(i, j)
    for (i = 0; i < _n; i++) {
        for (j = 0; j < i; j++) {
            if (mat[i][j] != mat[j][i]) {    
                is_symmetric = 0;
            }
        }
    }    
    #endif

  return is_symmetric; 		
}

//calculate the traspose of the matrix using OPENMP
float **matTransposeOMP(float **mat, int _n) {

  #ifdef _OPENMP
    wt1 = omp_get_wtime();
  #endif 
  
  // Allocate memory for the transposed matrix
  float **temp = (float **)malloc(_n * sizeof(float *)); 
  for (int i = 0; i < _n; i++) {
    temp[i] = (float *)malloc(_n * sizeof(float));
  }
        
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
    wt2 = omp_get_wtime();
  #endif
  
  #ifdef _OPENMP
  printf("Execution times of the transposition routine check: %f seconds\n", wt2 - wt1);
  #endif
  
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
 //if n = 1 return 1 (true), in the other case (n != 1) return 0 (false) 
  return (num==1); 
}

//check if the traspose of the sequential code is equal to the other 
void checktraspose(float **SEQ, float **B, int _n) {

  int check = 1; 
  
  for (int i = 0; i < _n; i++) {
 	  for(int j = 0; j < _n; j++ ) {
  		if (SEQ[i][j] != B[i][j]){
                        //the other code is not equal to the sequential one, so the comparation is false
  			check = 0;  
  		} 
  	}
  }
  
  if (check == 1) {
          //the other code is equal to the sequential one (true = 1) 
    printf("\nequal\n"); 
  } else {
    printf("\nnot equal\n"); 
  }

}
  
int main(int argc, char** argv){
  
  /*Initialize the random number generator with the current time to 
  have different random sequences on each program execution */
  srand(time(NULL)); 
  
  int n = atoi(argv[1]); 
//create 3 different matrix for the traspose because I have to check if the traspose is equal in every method
  float **TSEQ, **TIMP, **TOMP; 
  	
  /* Input and check for matrix dimension
  do {
  printf("insert number: ");
  scanf("%d", &n);
  }while (checksquare(n) == 0);*/
  
   // Allocate memory for matrix M and check for allocation success
  float **M = (float **)malloc(n * sizeof(float *)); 
  for (int i = 0; i < n; i++) {
        M[i] = (float *)malloc(n * sizeof(float));
  } 
   
  if ( M == NULL ) {
    printf (" Memory allocation failed \n") ;
    return 1;
  } 
  
  initializeMatrix(M, n);  
  
  printf("\n%d\n", n); 

  
//sequential
printf("\n_______________________________sequential____________________________\n"); 
  
   // Start time
  wstart = omp_get_wtime();
  
  if (checkSym(M, n) == 1 ){
   
   // End time 
  wend = omp_get_wtime();
  
  printf("the matrix M is symmetric \n");
    
  printf("Execution times of the checksym routine: %12.4g seconds\n", wend - wstart);    
    
  } else { 
  
  // End time
  wend = omp_get_wtime();
      
  //printf("Execution times of the checksym routine: %12.4g seconds\n", wend - wstart); 
  
  TSEQ = matTranspose(M, n);
 
  //checktraspose(M, TSEQ, n); 
     
  }
  
//with simd  
printf("\n___________________implicit parallelization (SIMD)___________________\n"); 
  
  // Start time
   wstart = omp_get_wtime();
      
  if (checkSymImp(M, n) == 1 ){
  
  	// End time 
    wend = omp_get_wtime();
    
    printf("the matrix M is symmetric \n"); 
     
    printf("Execution times of the checksym routine: %12.4g seconds\n", wend - wstart); 
    
  } else { 

  // End time
  wend = omp_get_wtime(); 
    
  //printf("Execution times of the checksym routine: %12.4g seconds\n", wend - wstart);
  
  TIMP = matTransposeImp(M, n);
  
  //checktraspose(TSEQ, TIMP, n); 
  
  }
  
//with OPENMP
printf("\n___________________explicit parallelizarion (OPENMP)__________________\n"); 

  
  #ifdef _OPENMP
    wt1 = omp_get_wtime();
    #endif
    
  if (checkSymOMP(M, n) == 1) {
  
    #ifdef _OPENMP
    wt2 = omp_get_wtime();
    #endif
      
    printf("the matrix M is symmetric \n");
    
    #ifdef _OPENMP
    printf("Execution times of the checksym routine: %f seconds\n", wt2 - wt1);
    #endif
    
    } else {
    
    #ifdef _OPENMP
    wt2 = omp_get_wtime();
    #endif
    
    /*#ifdef _OPENMP
    printf("Execution times of the checksym routine: %f seconds\n", wt2 - wt1);
    #endif */
  
    TOMP = matTransposeOMP(M, n);  
    
    //checktraspose(TSEQ, TOMP, n); 
  
    }
   
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
  
  return 0;
}
