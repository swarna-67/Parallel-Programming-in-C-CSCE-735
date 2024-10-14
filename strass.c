#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<time.h>
#include<omp.h>
//#include<ctime.h>
                
#define MAX_THREADS     2000
// Global variables
int num_threads;
int k_limit;
	  
   
       
void split_matrix(int** M, int** R, int nl, int x, int y){
    
        for(int z=0;z<nl;z++)
        R[z] = (int*)malloc(nl*sizeof(int ));
       // omp_set_num_threads(16);       //parllelize the for loop in split matrix function
       #pragma omp parallel for
        for(int i1 = 0, i2 = x; i1 < nl; i1++, i2++)
            #pragma omp parallel for
            for(int j1 = 0, j2 = y; j1 < nl; j1++, j2++)
                R[i1][j1] = M[i2][j2];
    
  
} 
//function to perform addition of 2 matrices
int** add_matrix(int** X, int** Y,int m){
  
   int** res = (int**)malloc(m*sizeof(int*));

  for(int i3 =0;i3<m;i3++){
       res[i3] = (int*)malloc(m*sizeof(int));
     }
      //omp_set_num_threads(16);       //parllelize the for loop in add matrix function
   #pragma omp parallel for
    for(int g=0;g<m;g++)
        #pragma omp parallel for
		for(int h=0;h<m;h++)
			res[g][h]= X[g][h]+ Y[g][h];

   return res;
}  
 
//function to perform subtraction of 2 matrices
int** sub_matrix(int** X1, int** Y1, int m1){
   
  int** res1 = (int**)malloc(m1*sizeof(int*));
 
for(int i4 =0;i4<m1;i4++){
    res1[i4] = (int*)malloc(m1*sizeof(int));
   }
      //omp_set_num_threads(16);       //parllelize the for loop in add matrix function
      #pragma omp parallel for
      for(int g1=0;g1<m1;g1++)
         #pragma omp parallel for
		  for(int h1=0;h1<m1;h1++)
		    	res1[g1][h1]= X1[g1][h1] - Y1[g1][h1];

return res1;
}
   
//function to get add back the C11,C12,C21,C22  to the product matrix that contains the result of multiplication
void combine_to_result(int** child, int **parent, int rw,  int  cl, int n   ){
 
//omp_set_num_threads(16);//parllelize the for loop in add matrix function
//printf("Hello world from threadcombine = %d \n",omp_get_num_threads());
#pragma omp parallel for
for(int z1 =0, z2 = rw; z1<n; z1++, z2++){
     
     #pragma omp parallel for
     for(int y1=0, y2 = cl; y1<n; y1++, y2++ ){
         parent[z2][y2] = child[z1][y1];
     }
}
         

 
}  
   
int** get_product_strassen(int** MATA, int** MATB,int n,int k_dash){
 
  //int tid = omp_get_num_threads();
  
 
int** product = (int**)malloc(n*sizeof(int*));

for(int i1 =0;i1<n;i1++){
    product[i1] = (int*)malloc(n*sizeof(int));
   }
   

   
 
if(n==1 || n<= k_limit ){ //if n<=k_dash
  //perfrom standard multiplication after some k' recursive process
 
   int n_new = n;
  // printf("Value of q here = %d \n",q);
  //omp_set_num_threads(8);
 // #pragma omp parallel 
 // {
 //   int tic = omp_get_num_threads();
  //  printf("Hello world from thread2.1 = %d \n",tic);  
  #pragma omp parallel for   
  for(int p1 =0;p1<n_new;p1++){
    #pragma omp parallel for
    for(int p2=0;p2<n_new;p2++){
      product[p1][p2]=0;
       for(int p3 =0;p3<n_new;p3++){
         int multi_term = MATA[p1][p3]* MATB[p3][p2];
         product[p1][p2] = product[p1][p2] + multi_term;
       }

    }
  } 
  //}
   
}
//do strassens recursive process
else {
 
  int n_nw= n/2;
 
 //declare M1
 int** M1 = (int**)malloc(n_nw*sizeof(int*));

  //declare M2
  int** M2 = (int**)malloc(n_nw*sizeof(int*));
 
   //declare M3
  int** M3 = (int**)malloc(n_nw*sizeof(int *));
 
   //declare M4
   int**M4 = (int**)malloc(n_nw*sizeof(int *));
 
   //declare M5
   int**M5 = (int**)malloc(n_nw*sizeof(int *));
  
   //declare M6
   int**M6 = (int**)malloc(n_nw*sizeof(int *));
 
   //declare M7
  int**M7 = (int**)malloc(n_nw*sizeof(int *));
 //declare A11,A12,A21,A22, B12, B21,B22 ,B11
   

 
 //omp_set_num_threads(16);
 // #pragma omp parallel sections
 // {
 //   #pragma omp section
//   {
   int **A11 = (int**)malloc(n_nw*sizeof(int *));
  split_matrix(MATA, A11,n_nw, 0 , 0);
//     }
//  #pragma omp section
//   {
  int **A12 =  (int**)malloc(n_nw*sizeof(int *)); 
  split_matrix(MATA, A12, n_nw, 0 , n_nw);
//   }
//pragma omp section
 // {
 int **A21 = (int**)malloc(n_nw*sizeof(int *));
  split_matrix(MATA, A21, n_nw, n_nw, 0);
 // }
// #pragma omp section
//  {
  int **A22 =  (int**)malloc(n_nw*sizeof(int *));
  split_matrix(MATA, A22, n_nw, n_nw, n_nw);
// } 
     
  //splitting the same way for B matrix
 // #pragma omp section
 // {
  int **B11 = (int**)malloc(n_nw*sizeof(int *));
  split_matrix(MATB, B11,n_nw, 0 , 0);
 // }
 // #pragma omp section
 // {
   int **B12 =  (int**)malloc(n_nw*sizeof(int *));
  split_matrix(MATB, B12, n_nw, 0, n_nw);
//  }
 // #pragma omp section
//{
  int **B21 =  (int**)malloc(n_nw*sizeof(int *));
  split_matrix(MATB, B21,n_nw, n_nw, 0);
//  }
 // #pragma omp section
// {
  int **B22 =  (int**)malloc(n_nw*sizeof(int *));
  split_matrix(MATB, B22, n_nw, n_nw, n_nw);
 // }
     
// }   
  //creating pointers to M1,M2,M3,M4 ,M5,M6 AND M7n where 7 multiplications are required
  //parallelize the recursive calls using task openmp
 //omp_set_num_threads(8);
 #pragma omp parallel 
   {   
     //int tid = omp_get_num_threads();
     // printf("Hello world from thread2 = %d \n",tid);  
 #pragma omp single
  {
   #pragma omp task 
  {     
        
        M1 = get_product_strassen(add_matrix(A11, A22,n_nw), add_matrix(B11, B22,n_nw), n_nw,k_dash);
 }
  
  #pragma omp task 
  {  
         
         M2 = get_product_strassen(add_matrix(A21, A22,n_nw), B11,n_nw,k_dash);
  } 

 #pragma omp task 
  {
        
        M3 = get_product_strassen(A11, sub_matrix(B12, B22,n_nw),n_nw, k_dash);
  }
 #pragma omp task 
  {
        
        M4 = get_product_strassen(A22, sub_matrix(B21, B11,n_nw),n_nw,k_dash);
  }
 #pragma omp task 
   {   
        
        M5 = get_product_strassen(add_matrix(A11, A12,n_nw), B22,n_nw, k_dash);
  }
#pragma omp task 
  {
       
        M6 = get_product_strassen(sub_matrix(A21, A11,n_nw), add_matrix(B11, B12,n_nw), n_nw, k_dash);
  }
 #pragma omp task 
  {
       
       M7 = get_product_strassen(sub_matrix(A12, A22,n_nw), add_matrix(B21, B22,n_nw),n_nw, k_dash);
  }
 #pragma omp taskwait
    
  }
   } 
     
   
   
  //C11 = M1 +M4 - M5 + M7 
  int** C11 = add_matrix(sub_matrix(add_matrix(M1,M4,n_nw),M5,n_nw),M7,n_nw);
  //join C11 to the product matrix
 // #pragma omp task
 // {
  combine_to_result(C11,product, 0, 0, n_nw);
 // }
  
  //C12 = M3 + M5
  int** C12 = add_matrix(M3, M5,n_nw);
  
  //join C12 to the product matrix
//  #pragma omp task
 // {
  combine_to_result(C12,product, 0, n_nw, n_nw);
 // }
  //C21 = M2 + M4
  int** C21 = add_matrix(M2, M4, n_nw);
  //join C21 to the product matrix
 // #pragma omp task
//  {
  combine_to_result(C21, product, n_nw, 0, n_nw);
 // }
  
  //C22 = M1 -M2 + M3 + M6
  int** C22 = add_matrix(sub_matrix(add_matrix(M1, M3,n_nw),M2, n_nw), M6, n_nw);
   //join C22 to the product matrix
 // #pragma omp task
 // {
  combine_to_result(C22, product, n_nw, n_nw, n_nw);
//  }
//  #pragma omp taskwait
   
        
		   
   //deallocate the memory for all the pointers previously allocated using malloc 

	free(M1); free(M2); free(M3); free(M4); free(M5); free(M6); free(M7);
	free(A11); free(A12); free(A21); free(A22);
	free(B11); free(B12); free(B21); free(B22);
		    
	 	 
	 	  
		 
  } 



 
return product;
 


}



int main(int argc, char* argv[]) {
	  
	int i=0,j=0;
	int k,n,k_dash,q;
	
	if(argc !=4) {
	
		printf("\n Invalid number of arguments!\n\n");
		return 0;
	}
     
  double total_time, time_res;
   k=atoi(argv[1]);
   n= pow(2, k);
  int n1 =n;
  k_dash = atoi(argv[2]);
  k_limit = pow(2, k_dash);
  q = atoi(argv[3]);
  num_threads = q;
  //setting the number of threads to be used in all the parallel regions
  
  if (num_threads  > MAX_THREADS) {
  printf("Maximum number of threads allowed: %d.\n", MAX_THREADS);
  exit(0);
  }
  
  
  if(n<1 || n%2!=0){
		
		printf("\nPlease Enter a  matrix dimension that is power of 2 only!\n\n");
		return 0;
	
	}
 
    //create matrix A
    int** MATA = (int**)malloc(n*sizeof(int *));

    for(int i=0;i<n;i++){
        MATA[i] = (int*)malloc(n*sizeof(int));
        for(int j=0;j<n;j++){
            MATA[i][j] = (rand() % 500 );
        }
    }  
    
   //create matrix b 
   int** MATB = (int**)malloc(n*sizeof(int *));
   
   
   for(int k=0;k<n;k++){
        MATB[k] = (int*)malloc(n*sizeof(int));
        for(int l=0;l<n;l++){
            MATB[k][l] = (rand() % 500 );
        }
    }
     
          
   
   //printf(" %d", k_dash);
   //MATRIX C TO store the product from strassens multiplication computation
   int** MATC = (int**)malloc(n*sizeof(int*));
   
   for(int i6 =0;i6<n;i6++){
       MATC[i6] = (int*)malloc(n*sizeof(int));
   }   
   //clock_gettime(CLOCK_REALTIME, &start);//start time
   //omp_set_nested(0);  //tired to set this nested parallelizm to parallelize the successive recursive calls but gives a high overhead of threads when used for more levels
   //omp_set_max_active_levels(1);
   omp_set_dynamic(0);   
   omp_set_num_threads(q);
   double start_time = omp_get_wtime();
    
       
   MATC = get_product_strassen(MATA,MATB,n,k_dash);
 
    
   //to find the execution time
  // clock_gettime(CLOCK_REALTIME, &stop);
  double end_time = omp_get_wtime();

  total_time = end_time - start_time; //(stop.tv_sec-start.tv_sec);
   
 
 // regular matrix mulitplication to check the correctness  of strassens algo
  int** standardprod  = (int**)malloc(n1*sizeof(int *));
    for(int z4=0;z4<n1;z4++){
         standardprod[z4] = (int*)malloc(n1*sizeof(int));
   } 
  for(int z1 =0;z1<n1;z1++){
    for(int z2=0;z2<n1;z2++){
        standardprod[z1][z2]=0;
         for(int z3 =0;z3<n1;z3++){
          int mul_term = MATA[z1][z3]* MATB[z3][z2];
          standardprod[z1][z2] = standardprod[z1][z2] + mul_term;
       }
 
    }
  }    
   
       
   //print the result matrix
   //printf("\nMatrix C\n");
   for(int t=0;t<n1;t++){
     for(int u =0;u<n1;u++){
            //printf("   %d   ", MATC[t][u]);
            if((int)standardprod[t][u]!=(int)MATC[t][u]){
             printf("\n Algo incorrect\n");
            return 0;
            }
       }
   //  printf("\n");
   }   
   
   printf("\nStrassen's multiplication Is correct\n");
   printf("\nMatrix Size  = %d, k'_size = %d,  Threads = %d,  Total time taken- time (sec) = %8.7f",n , k_limit , q , total_time);
   
   //deallocate the memory previously allocted via malloc   
   free(MATA); free(MATB);
   free(MATC); free(standardprod);
   return 0; 
    
} 

  

