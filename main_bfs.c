#include "graph_bfs.c"
#include <mpi.h>
#include <time.h>


void partGraph(GraphStruct *G_graph,int num_procs){
	int i;
	// int num_procs = 3;
	// MPI_Comm_size(comm, &num_procs);
	for(i=0; i<G_graph->numVtxs; i++){	// chia cac canh vao cac process 
		int part = G_graph->vtxGIDs[i] % num_procs;
		G_graph->part[i] = part;
	}
}


void graphSend(MPI_Comm comm, GraphStruct sub_graph, int dest_proc) {

	int info_graph[3];
	int ack1 = 100, ack_tag = 5, count_tag = 10, id_tag = 15;
	int ack2 = 0;
	info_graph[0] = sub_graph.numVtxs;
	info_graph[1] = sub_graph.numNbors;

	MPI_Send(info_graph, 3, MPI_INT, dest_proc, count_tag, comm);
	MPI_Recv(&ack2, 1, MPI_INT, dest_proc, ack_tag, comm, MPI_STATUS_IGNORE);

	if (ack2 != 200) {
		printf("Rank 0 send info_graph from rank %d failed!\n", dest_proc);
	}
	if (info_graph[0] > 0) {
		MPI_Send(sub_graph.vtxGIDs, info_graph[0], MPI_INT, dest_proc, id_tag, comm);
		MPI_Send(sub_graph.nborIndex, info_graph[0] + 1, MPI_INT, dest_proc, id_tag + 1, comm);
		MPI_Send(sub_graph.part, info_graph[0], MPI_INT, dest_proc, id_tag + 2, comm);

		if (info_graph[1] > 0) {
			MPI_Send(sub_graph.nborGIDs, info_graph[1], MPI_INT, dest_proc, id_tag + 3, comm);
			MPI_Send(sub_graph.nborProcs, info_graph[1], MPI_INT, dest_proc, id_tag + 4, comm);
		}
	}
	// MPI_Send(&ack1, 1, MPI_INT, dest_proc, ack_tag, comm);

}
void graphRecv(MPI_Comm comm, GraphStruct * local_graph){
	printf("graphRecv  44 \n ");
	int info_graph[3];
	int ack2 = 200, ack_tag = 5, count_tag = 10, id_tag = 15;
	int ack1 = 0;
	MPI_Recv(info_graph, 3, MPI_INT, 0, count_tag, comm, MPI_STATUS_IGNORE);
	MPI_Send(&ack2, 1, MPI_INT, 0, ack_tag, comm);

	if (info_graph[0] > 0) {
		graphInit(local_graph, info_graph[0], info_graph[1]);

		MPI_Recv(local_graph->vtxGIDs, info_graph[0], MPI_INT, 0, id_tag, comm, MPI_STATUS_IGNORE);
		MPI_Recv(local_graph->nborIndex, info_graph[0] + 1, MPI_INT, 0, id_tag + 1, comm, MPI_STATUS_IGNORE);
		MPI_Recv(local_graph->part, info_graph[0], MPI_INT, 0, id_tag + 2, comm, MPI_STATUS_IGNORE);
		if (info_graph[1] > 0) {
			MPI_Recv(local_graph->nborGIDs, info_graph[1], MPI_INT, 0, id_tag + 3, comm, MPI_STATUS_IGNORE);
			MPI_Recv(local_graph->nborProcs, info_graph[1], MPI_INT, 0, id_tag + 4, comm, MPI_STATUS_IGNORE);
		}
	}

	MPI_Recv(&ack1, 1, MPI_INT, 0, 0, comm, MPI_STATUS_IGNORE);
	if (ack1 != 100 ) {
		int my_rank;
		MPI_Comm_rank(comm, &my_rank);
		printf("Rank %d receive subgraph from rank 0 failed!\n", my_rank);
	}
}

void graphDistribute(MPI_Comm comm, GraphStruct G_graph, GraphStruct * L_graph,int num_procs){
// void graphDistribute(GraphStruct G_graph,int num_procs){

	int i, j, proc_id, vtx_lid, nnbors, nbor_index, nbor_proc;
	GraphStruct *sub_graphs = (GraphStruct *) malloc(sizeof(GraphStruct) * num_procs);
	for (i = 0; i < num_procs; i++) {
		sub_graphs[i].numVtxs = 0;
		sub_graphs[i].numNbors = 0;
	}
printf(" send 80 ");

	int *global_part_vector = malloc(sizeof(int) * G_graph.numVtxs);
	for (i = 0; i < G_graph.numVtxs; i++) {
		proc_id = G_graph.part[i];
		sub_graphs[proc_id].numVtxs++;
		sub_graphs[proc_id].numNbors += (G_graph.nborIndex[i + 1] - G_graph.nborIndex[i]);
		global_part_vector[G_graph.vtxGIDs[i]] = G_graph.part[i];
	}
printf(" send 88 ");

	for (i = 0; i < num_procs; i++) {
		graphInit(&sub_graphs[i], sub_graphs[i].numVtxs, sub_graphs[i].numNbors);
		sub_graphs[i].numVtxs = 0;
		sub_graphs[i].numNbors = 0;
	}
printf(" send 94 ");

	for (i = 0; i < G_graph.numVtxs; i++) {
		proc_id = G_graph.part[i];
		vtx_lid = sub_graphs[proc_id].numVtxs++;
		sub_graphs[proc_id].vtxGIDs[vtx_lid] = G_graph.vtxGIDs[i];
		sub_graphs[proc_id].part[vtx_lid] = G_graph.part[i];

		nnbors = G_graph.nborIndex[i + 1] - G_graph.nborIndex[i];
		sub_graphs[proc_id].nborIndex[vtx_lid + 1] = sub_graphs[proc_id].nborIndex[vtx_lid] + nnbors;
		for (j = G_graph.nborIndex[i]; j < G_graph.nborIndex[i + 1];j++) {
			nbor_index = sub_graphs[proc_id].numNbors++;
			nbor_proc = global_part_vector[G_graph.nborGIDs[j]];

			sub_graphs[proc_id].nborProcs[nbor_index] = nbor_proc;
			sub_graphs[proc_id].nborGIDs[nbor_index] = G_graph.nborGIDs[j];
		}
	}

// graphPrint(sub_graphs[0]);
// graphPrint(sub_graphs[1]);
// graphPrint(sub_graphs[2]);

	free(global_part_vector);
	*L_graph = sub_graphs[0];
printf(" send 119 ");
	for (i = 1; i < num_procs; i++) {
		printf(" send %d ",i);
		graphSend(comm, sub_graphs[i], i);
		graphDeinit(&sub_graphs[i]);
	}
	free(sub_graphs);
	int ack1 = 100;
	for (i = 1; i < num_procs; i++) {
		MPI_Send(&ack1, 1, MPI_INT, i, 0, comm);
	}
}

typedef struct ArrayList {
	int* data;
	int length;
	int size;
} ArrayList_t;
#define DEFAULT_SIZE 64

ArrayList_t* listCreate() {
	ArrayList_t *list = (ArrayList_t*) calloc(1, sizeof(ArrayList_t));
	if (!list)
		return NULL;
	list->size = DEFAULT_SIZE;
	list->length = 0;
	if (!(list->data = (int*) calloc(sizeof(int), list->size))) {
		free(list);
		return NULL;
	}
	return list;
}

void listDestroy(ArrayList_t *list) {
	free(list->data);
	free(list);
}

static int _listExpand(ArrayList_t *list, int max) {
	if (max < list->size)
		return 0;
	int new_size = list->size  * 1.5;
	if (new_size < max)
		new_size = max;

	int *t;
	if (!(t = realloc(list->data, new_size * sizeof(int))))
		return -1;
	list->data = t;
	(void) memset(list->data + list->size, 0, (new_size - list->size) * sizeof(int));
	list->size = new_size;
	return 0;
}
int listPutIdx(ArrayList_t *list, int idx, int data) {
	if (_listExpand(list, idx + 1))
		return -1;
	list->data[idx] = data;
	if (list->length <= idx)
		list->length = idx + 1;
	return 0;
}
int listAppend(ArrayList_t *arr, int data) {
	return listPutIdx(arr, arr->length, data);
}

void listClear(ArrayList_t *list){
	list->length = 0;
}



int* BFS(MPI_Comm comm, GraphStruct L_graph, int srcLid, int srcRank,int num_procs){
	int i, j;
	int myRank;
	MPI_Comm_rank(comm, &myRank);
	 
	int totalVtx;
	MPI_Allreduce(&(L_graph.numVtxs), &totalVtx, 1, MPI_INT, MPI_SUM, comm);
	int *gid2lid = (int *) calloc(sizeof(int), totalVtx);
	for(i=0; i<L_graph.numVtxs; i++){
		int gid = L_graph.vtxGIDs[i];
		gid2lid[gid] = i;
	}

	
	int *recvCount = (int *) malloc(sizeof(int) * num_procs);
	int **recvBuf = (int **) malloc(sizeof(int *) * num_procs);
	ArrayList_t **sendBuf = (ArrayList_t **) malloc(sizeof(ArrayList_t *) * num_procs);
	for(i=0; i<num_procs; i++){
		sendBuf[i] = listCreate();
	}
	/*****************************************************************/
	int *d = (int *) malloc(sizeof(int)* L_graph.numVtxs);
	memset(d, -1, sizeof(int) * L_graph.numVtxs);
	ArrayList_t * FS = listCreate();	// cac canh dang o vi tri duyet
	if(srcRank == myRank){
		d[srcLid] = 0;
		listAppend(FS, srcLid);
	}
	int level = 1;
	long ETPS = 0;
	int numActiveVertices = 0;
	do{
		//duyet
		ArrayList_t * NS = listCreate();	//cac canh o buoc tiep theo
		for(i=0; i<FS->length; i++){
			int lid = FS->data[i];
			for(j=L_graph.nborIndex[lid]; j<L_graph.nborIndex[lid + 1]; j++){
				int nborGID = L_graph.nborGIDs[j];
				int owner = L_graph.nborProcs[j];
				if(owner == myRank){
					int lid = gid2lid[nborGID];
					if(d[lid] == -1){
						listAppend(NS, lid);
						d[lid] = level;
					}
				}else{
					listAppend(sendBuf[owner], nborGID);
				}
			}

			int numNbors = L_graph.nborIndex[lid + 1] - L_graph.nborIndex[lid];
			ETPS += numNbors;
		}
		listDestroy(FS);
		FS = NS;

		MPI_Request request;
		//gui cac dinh duyet toi process chua dinh ay
		for(i=0; i<num_procs; i++){
			if(sendBuf[i]->length){
				MPI_Isend(sendBuf[i]->data, sendBuf[i]->length, MPI_INT, i, 1, comm, &request);
				MPI_Request_free(&request);
			}
		}
		for(i=0; i<num_procs; i++){	// 
			MPI_Gather(&(sendBuf[i]->length), 1, MPI_INT, recvCount, 1, MPI_INT, i, comm);
		}
		for(i=0; i<num_procs; i++){
			recvBuf[i] = (int *) malloc(sizeof(int) * recvCount[i]);
			if(recvCount[i]){
				MPI_Recv(recvBuf[i], recvCount[i], MPI_INT, i, 1, comm, MPI_STATUS_IGNORE);
			}
		}

		// gan gia tri cho cai dinh da duyet
		for(i=0; i<num_procs; i++){
			for(j=0; j<recvCount[i]; j++){
				int gid = recvBuf[i][j];
				int lid = gid2lid[gid];
				if(d[lid] == -1){
					d[lid] = level;
					listAppend(FS, lid);
				}
			}
			free(recvBuf[i]);
		}
		numActiveVertices = FS->length;
		MPI_Allreduce(MPI_IN_PLACE, &numActiveVertices, 1, MPI_INT, MPI_SUM, comm);
		for(i=0; i<num_procs; i++){
			listClear(sendBuf[i]);
		}
//		if(myRank == 0) printf("step-%d\n", level);
		level ++;
	}while(numActiveVertices > 0);

	free(gid2lid);

	for(i=0; i<num_procs; i++){
		listDestroy(sendBuf[i]);
	}
	free(sendBuf);
	free(recvBuf);
	free(recvCount);

	return d;
}

int main(int argc, char *argv[]){
	MPI_Init(&argc, &argv);
	MPI_Comm comm = MPI_COMM_WORLD;
	// char *fname = "graph.txt";
	char *fname;
	fname = argv[1];
	int numParts;
	MPI_Comm_size(comm, &numParts);
	int myRank;
	MPI_Comm_rank(comm, &myRank);

	double time_distributed,time_bfs;


	GraphStruct L_graph;
printf("%d\n",myRank);

 	if (myRank == 0){
		FILE *fp = fopen(fname, "r");
		if (fp == NULL) {
			printf("File %s not found\n", fname);
		}
		GraphStruct G_graph;
		graphLoad(&G_graph, fp);
		printf(" graphLoad done 149\n");
		// graphPrint(G_graph);

		partGraph(&G_graph,numParts);
		printf(" partGraph done 153\n");
		time_distributed = MPI_Wtime();
		graphDistribute(comm,G_graph,&L_graph,numParts);
		time_distributed = MPI_Wtime() - time_distributed;
		printf(" graphDistribute done 155\n");

		graphDeinit(&G_graph);
		fclose(fp);
 	}
 	else{
 		printf("%d\n",myRank);
 		graphRecv(comm, &L_graph);
 	}

 	srand(time(NULL));
	int srcRank = rand() % numParts;
	int srcLid = rand() % L_graph.numVtxs;
	time_bfs = MPI_Wtime();
	int *d = BFS(comm, L_graph, srcLid, srcRank,numParts);
	time_bfs = MPI_Wtime() -time_bfs;
	if (myRank == 0){
	printf("time_distributed  %10.3lf milliseconds \n", time_distributed*1000);

	printf("time_bfs  %10.3lf milliseconds \n", time_bfs*1000);
	}
	for(int i=0; i<L_graph.numVtxs; i++){
		 printf("myRank = %d  d[%d] =  %d \n",myRank,i,d[i]);
	}
 	MPI_Finalize();
}
