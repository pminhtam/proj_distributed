#include "mpi.h"              /* MPI header file */
#include <stdio.h>
#include<string.h>

#include <stdlib.h>
#include <sys/time.h>
#define MAX_LEN 100000   /* Max array size */

/* Usual search for largest function  */
int find_max(int a[],int len)
{
   int i;   
   int max; /* Current max */
   max = a[0];
   for (i=1;i<len;++i)
		    if (a[i] > max) max = a[i];
   return max;
}
int tag = 27;
int arrDist(int *data,int len,int localLen,int num_procs,MPI_Comm comm){
	int num_per_procs = len/num_procs;

	int **subArr = (int **) malloc(sizeof(int*)*(num_procs-1));
	for(int i=0;i<num_procs-1;i++){
		subArr[i] = (int *) malloc(sizeof(int)*(num_per_procs));
	}
	for (int j=0;j<num_procs-1;j++){
		for (int i=0;i<num_per_procs;i++){
			subArr[j][i] = data[localLen + j*num_per_procs+i];
		}	
	}

	for(int i=0; i<num_procs-1;i++){
		MPI_Send(&num_per_procs,1,MPI_INT,i+1,tag,comm);
	}

	for(int i=0; i<num_procs-1;i++){
		MPI_Send(subArr[i],num_per_procs,MPI_INT,i+1,tag+1,comm);
	}
	printf("Send xong\n");

} 
static int bufsize=1<<20; //256MB line buffer
static char line[1<<20];


int main(int argc, char *argv[]){

	MPI_Init(&argc, &argv);
	MPI_Comm comm = MPI_COMM_WORLD;
	int numParts;
	MPI_Comm_size(comm, &numParts);
	int myRank;
	MPI_Comm_rank(comm, &myRank);
	int *localArr;
	int localLen;
	if (myRank == 0){
		char *fname = "data.txt";
		FILE *fp = fopen(fname, "r");
		char *c;
		int *data;
		c = fgets(line, bufsize, fp);    // read a line 

		int len = atoi(strtok(c, " "));
		data = (int *) malloc(sizeof(int)*len);
		for (int i=0;i<len;i++){
			data[i] = atoi(strtok(NULL, " "));
		}

		// for (int i=0;i<len;i++){
		// 	printf("%d\n",data[i] );
		// }

		//====================
		localLen = len/numParts + len%numParts;
		localArr = (int *) malloc(sizeof(int)*(localLen));
		for(int i=0;i<localLen;i++){
			localArr[i] = data[i];
		}
		arrDist(data,len,localLen,numParts,comm);
		// printf("\nlocalLen 0 ====  %d \n", localLen);
		// for (int i=0;i<localLen;i++){
		// 	printf("%d\n",localArr[i] );
		// }
	}
	else{
		MPI_Recv(&localLen, 1, MPI_INT, 0, tag, comm, MPI_STATUS_IGNORE);
		// printf("\n %d \n", localLen);
		localArr = (int *) malloc(sizeof(int)*(localLen));
		MPI_Recv(localArr, localLen, MPI_INT, 0, tag+1, comm, MPI_STATUS_IGNORE);
		// for (int i=0;i<localLen;i++){
		// 	printf("%d\n",localArr[i] );
		// }
	}

	int max = find_max(localArr,localLen);
	printf("\n max proc %d ==== %d \n",myRank,max );
	free(localArr);
	int global_max;
	MPI_Reduce(&max, &global_max, 1, MPI_INT, MPI_MAX,0, comm);
	if(myRank == 0)
		printf("\n %d \n",global_max );


	MPI_Finalize();

}