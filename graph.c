// #include<iostream>
#include<stdio.h>
#include<string.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>

typedef struct {
  int numVtxs;
  int numNbors;
  int *vtxGIDs;
  int *part;
  int *nborIndex;
  int *nborGIDs;
  int *nborProcs;
} GraphStruct;
void graphInit(GraphStruct *graph,int numVtxs,int numNbors){
  graph->numVtxs = numVtxs;
  graph->numNbors = numNbors;
  if(numVtxs>0){
    graph->vtxGIDs = (int *) malloc(sizeof(int)*numVtxs);
    graph->part = (int *) malloc(sizeof(int) * (numVtxs));
    graph->nborIndex = (int *) calloc(sizeof(int), (numVtxs + 1));
    graph->nborGIDs = (int *) malloc(sizeof(int) * numNbors);
    graph->nborProcs = (int *) malloc(sizeof(int) * numNbors);

  }
}

void graphDeinit(GraphStruct * graph) {
    if(graph->numVtxs > 0){
      free(graph->vtxGIDs);
      free(graph->part);
      free(graph->nborIndex);
      free(graph->nborGIDs);
      free(graph->nborProcs);
    }

    graph->numVtxs = 0;
    graph->numNbors = 0;
}
void graphPrint(GraphStruct graph){
  printf("\n num_vtx=%d, num_nbors=%d \n", graph.numVtxs, graph.numNbors);

  int i, j, num_nbors;
  // for(i=0; i<graph.numVtxs; i++){
  //   num_nbors = graph.nborIndex[i + 1] - graph.nborIndex[i];

  //   printf("%d %d ", graph.vtxGIDs[i], num_nbors);
  //   for(j=graph.nborIndex[i]; j<graph.nborIndex[i + 1]; j++){
  //     printf("cc %d ", graph.nborGIDs[j]);
  //   }
  //   printf("\n");
  // }
  printf("\n graph->vtxGIDs ");
  for (i=0;i<graph.numVtxs;i++){
    
    printf("%d ",graph.vtxGIDs[i] );
  }
  printf("\n graph->part ");
  for (i=0;i<graph.numVtxs;i++){
    
    printf("%d ",graph.part[i] );
  }
  printf("\n graph->nborIndex ");
  for (i=0;i<graph.numVtxs + 1;i++){
    
    printf("%d ",graph.nborIndex[i] );
  }
  printf("\n graph->nborGIDs ");
  for (i=0;i<graph.numNbors;i++){
    
    printf("%d ",graph.nborGIDs[i] );
  }
  printf("\n graph->nborProcs ");
  for (i=0;i<graph.numNbors;i++){
    
    printf("%d ",graph.nborProcs[i] );
  }
}




static int bufsize=1<<28; //256MB line buffer
static char line[1<<28];
/* Function to find next line of information in input file */
int getNextLine(FILE *fp, char *buf, int bufsize) {
  int i, cval, len;
  char *c;
  while (1) {
    c = fgets(buf, bufsize, fp);    // read a line 
    if (c == NULL)
      return 0; /* end of file */
    len = strlen(c);
    // printf(" %d  c:  %s   buf :%s ",len,c,buf );
    for (i = 0, c = buf; i < len; i++, c++) {   // skip space at head line
      cval = (int) *c;
      if (isspace(cval) == 0)
        break;
    }
    if (i == len)
      continue; /* blank line */
    if (*c == '#')
      continue; /* comment */
    if (c != buf) {
      strcpy(buf, c);
    }
    break;
  }
  return strlen(buf); /* number of characters */
}
int graphLoad(GraphStruct * graph, FILE * fp) {
  int numGlobalVertices, numGlobalEdges, numParts;
  int i, j, nnbors;

  /* Get the number of vertices */
  getNextLine(fp, line, bufsize);
  sscanf(line, "%d", &numGlobalVertices);

  /* Get the number of edges  */
  getNextLine(fp, line, bufsize);
  sscanf(line, "%d", &numGlobalEdges);

  /* Allocate arrays to read in entire graph */
  graphInit(graph, numGlobalVertices, numGlobalEdges << 1);

  char * token;

  for (i = 0; i < numGlobalVertices; i++) {
    getNextLine(fp, line, bufsize);

    token = strtok(line, " ");
    graph->vtxGIDs[i] = atoi(token);

    token = strtok(NULL, " ");
    nnbors = atoi(token);

    graph->nborIndex[i + 1] = graph->nborIndex[i] + nnbors;
    for (j = graph->nborIndex[i]; j<graph->nborIndex[i + 1]; j++) {
      token = strtok(NULL, " ");
      graph->nborGIDs[j] = atoi(token);
    }
  }
  return 0;
}



// int main(){
//   char *fname = "graph.txt";
//   FILE *fp = fopen(fname, "r");
//   if (fp == NULL) {
//     printf("File %s not found\n", fname);
//   }
//   GraphStruct G_graph;
//   graphLoad(&G_graph, fp);

//   graphPrint(G_graph);
// }
