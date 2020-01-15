#include <assert.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>



typedef struct {
    int rank;
    int row, col;
    int p, q;
    MPI_Comm comm;
    MPI_Comm row_comm;
    MPI_Comm col_comm;
} MPI_GRID_INFO;

void setup_grid(MPI_GRID_INFO *grid);
int check_can_devide(int p, int n);
int *read_mtrx(FILE * fp,int *n);
void send_sub_mtrx(int *_mtrx, int n, int q);
void *process_mtrx(MPI_GRID_INFO *grid, int *mtrx_A, int n);
void floyd_warshall(int *_A, int *_B, int *_C, int n);
// void fix_final_mtrx(int *_mtrx, int *_mtrx_F, int n, int q);
void graphPrint(int *_mtrx, int n);

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    MPI_GRID_INFO grid;
    setup_grid(&grid);

    int n, *mtrx;

    char *fname;
    fname = argv[1];

    if (grid.rank == 0) {
        FILE *fp = fopen(fname, "r");
        mtrx = read_mtrx(fp,&n);
        // doc so n: la so chieu cua ma tran ke
        // if (!scanf("%d", &n)) {
        //     fprintf(stderr, "Loi 44 \n");
        //     MPI_Abort(MPI_COMM_WORLD, 0);
        //     exit(1);
        // }
        // if (!check_can_devide(grid.p, n)) {
        //     fprintf(stderr, "Loi 49 %d , %d \n", n, grid.p);
        //     MPI_Abort(MPI_COMM_WORLD, 0);
        //     exit(1);
        // }
        // // doc ma tran ke
        // mtrx = read_mtrx(n);
            printf("main 56 \n");

// graphPrint(mtrx, n);

    }

    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    // Neu la process goc va co nhieu process thi gui submatrix
    double time_send;
    if (grid.rank == 0) {
        time_send = MPI_Wtime();
    }

    if (grid.rank == 0 && grid.p > 1) send_sub_mtrx(mtrx, n, grid.q);
    if (grid.rank == 0) {
        time_send = MPI_Wtime() - time_send;
        printf("Thoi gian khoi tao %10.3lf  milliseconds \n", time_send* 1000);
    }
    const int m = n / grid.q;   // kich thuoc cua submatrix
    int *mtrx_A;
    if (grid.p > 1) {       // process recieve ma tran
        mtrx_A = (int *) malloc(m * m * sizeof(int));
        assert(mtrx_A != NULL);
        MPI_Recv(mtrx_A, m * m, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else {
        mtrx_A = mtrx;
    }

    double time;
    // graphPrint(mtrx, n);

    // graphPrint(mtrx_A, n);
    time = MPI_Wtime();
    int *mtrx_C = process_mtrx(&grid, mtrx_A, n);
    time = MPI_Wtime() - time;
    // graphPrint(mtrx_C, n);
    if (grid.p > 1) free(mtrx_A);       // neu co nhieu prcess thi xa luon ma tran ban dau

    int *mtrx_F = malloc(n * n * sizeof(int));
    assert(mtrx_F != NULL);
    MPI_Gather(mtrx_C, m * m, MPI_INT, mtrx_F, m * m, MPI_INT, 0, MPI_COMM_WORLD);
    free(mtrx_C);

    if (grid.rank == 0) {
        // fix_final_mtrx(mtrx, mtrx_F, n, grid.q);
        fprintf(stderr, "\nThoi gian chay: %10.3lf milliseconds.\n\n", time * 1000);
        // graphPrint(mtrx, n);
        graphPrint(mtrx_F, n);

    }
    free(mtrx_F);
    MPI_Comm_free(&grid.comm);
    MPI_Comm_free(&grid.row_comm);
    MPI_Comm_free(&grid.col_comm);
    MPI_Finalize();
    return 0;
}

void setup_grid(MPI_GRID_INFO *grid) {
    MPI_Comm_size(MPI_COMM_WORLD, &(grid->p));
    MPI_Comm_rank(MPI_COMM_WORLD, &(grid->rank));

    grid->q = sqrt(grid->p);
    // printf("%d\n",grid->q );
    int dims[2] = { grid->q, grid->q };
    int periods[2] = { 1, 1 };
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 1, &(grid->comm));

    int coords[2];
    MPI_Comm_rank(grid->comm, &(grid->rank));
    MPI_Cart_coords(grid->comm, grid->rank, 2, coords);
    grid->row = coords[0];
    grid->col = coords[1];
    printf(" 120  rr  %d \n ", grid->row);
    printf(" 121 rr %d \n ", grid->col);

    coords[0] = 0; coords[1] = 1;
    MPI_Cart_sub(grid->comm, coords, &(grid->row_comm));
    coords[0] = 1; coords[1] = 0;
    MPI_Cart_sub(grid->comm, coords, &(grid->col_comm));
}

int check_can_devide(int p, int n) {
    int q = sqrt(p);
    if (q * q == p && n % q == 0) return 1;
    return 0;
}

int *read_mtrx(FILE * fp,int *n) {
    fscanf(fp, "%d", n);
    printf("%d\n", *n);   
    // mtrx = (int *) malloc(*n * *n * sizeof(int));
    // assert(mtrx != NULL);
    int (*_mtrx)[*n] = (int (*)[*n]) malloc(*n * *n * sizeof(int));
    int i, j;
    for(i = 0; i < *n; i++) {
        for(j = 0; j < *n; j++) {
            // if (!scanf("%d", &mtrx[i][j])) {
            if (!fscanf(fp, "%d", &_mtrx[i][j])) {

                printf("Loi 143 \n");
                MPI_Abort(MPI_COMM_WORLD, 0);
                exit(1);
            }
            // if (mtrx[i][j] == 0 && i != j) mtrx[i][j] = INF;
        }
    }
    printf("read_mtrx 161 \n");
    return (int *)_mtrx;
}

void send_sub_mtrx(int *_mtrx, int n, int q) {
    // n: kich thuoc ma tran
    // q can bac 2 cua comm size
    int m = n / q;
    int (*mtrx)[n] = (int (*)[n]) _mtrx;
    int (*sub_mtrx)[m] = (int(*)[m]) malloc(m * m * sizeof(int));
    assert(sub_mtrx != NULL);

    int dst = 0;
    int i, j, k, l;
    for (k = 0; k < q; k++) {
        for (l = 0; l < q; l++) {
            for (i = 0; i < m; i++) {
                for (j = 0; j < m; j++) {
                    sub_mtrx[i][j] = mtrx[i + (k * m)][j + (l * m)];
                }
            }
            MPI_Send(sub_mtrx, m * m, MPI_INT, dst++, 1, MPI_COMM_WORLD);
        }
    }
    free(sub_mtrx);
}

void *process_mtrx(MPI_GRID_INFO *grid, int *mtrx_A, int n) {
    int m = n / grid->q;
    int src = (grid->row + 1) % grid->q;
    int dst = (grid->row - 1 + grid->q) % grid->q;

    int *temp_A = (int *) malloc(m * m * sizeof(int));
    assert(temp_A != NULL);
    int *mtrx_B = (int *) malloc(m * m * sizeof(int));
    assert(mtrx_B != NULL);
    int *mtrx_C = (int *) malloc(m * m * sizeof(int));
    assert(mtrx_C != NULL);
    memcpy(mtrx_C, mtrx_A, m * m * sizeof(int));
// graphPrint(temp_A, m);
    double time_process_floyd;
    double sum_time_process_floyd = 0; 
    int iter;
    for (iter = 1; iter < n; iter <<= 1) {
        memcpy(mtrx_B, mtrx_C, m * m * sizeof(int));
        int stage;
        for (stage = 0; stage < grid->q; stage++) {
            int bcast_root = (grid->row + stage) % grid->q;
            if (bcast_root == grid->col) {
                MPI_Bcast(mtrx_A, m * m, MPI_INT, bcast_root, grid->row_comm);
                time_process_floyd = MPI_Wtime();

                floyd_warshall(mtrx_A, mtrx_B, mtrx_C, m);
                time_process_floyd = MPI_Wtime() - time_process_floyd;
                sum_time_process_floyd += time_process_floyd;
                // printf(" dong 194 bcast_root %d   grid->col %d\n",bcast_root, grid->col);
                // graphPrint(temp_A, m);
            } else {
                // printf(" dong 194 bcast_root %d   grid->col  %d  grid->rank   %d \n",bcast_root, grid->col,grid->rank);
                // graphPrint(temp_A, m);
                MPI_Bcast(temp_A, m * m, MPI_INT, bcast_root, grid->row_comm);
                time_process_floyd = MPI_Wtime();

                floyd_warshall(temp_A, mtrx_B, mtrx_C, m);
                time_process_floyd = MPI_Wtime() - time_process_floyd;
                sum_time_process_floyd += time_process_floyd;
                // printf(" dong 199 bcast_root %d   grid->col %d\n",bcast_root, grid->col);
                // graphPrint(temp_A, m);
            }
            MPI_Sendrecv_replace(mtrx_B, m * m, MPI_INT, dst, 1, src, 1, grid->col_comm, MPI_STATUS_IGNORE);
            // printf(" 211 stage  %d \n",stage);
        }

        // printf(" 214 iter  %d \n",iter);
    }
    printf("Prpcess %d Thoi gian chay floyd  %10.3lf  \n",grid->rank, sum_time_process_floyd*1000);
    free(temp_A);
    free(mtrx_B);
    return mtrx_C;
}

inline void floyd_warshall(int *_A, int *_B, int *_C, int n) {
    int (*A)[n] = (int (*)[n]) _A;
    int (*B)[n] = (int (*)[n]) _B;
    int (*C)[n] = (int (*)[n]) _C;

    int i, j, k;
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            for (k = 0; k < n; k++)
                C[i][j] = C[i][j] < A[i][k] + B[k][j] ? C[i][j] : A[i][k] + B[k][j] ;
}

// void fix_final_mtrx(int *_mtrx, int *_mtrx_F, int n, int q) {
//     int m = n / q;
//     int (*mtrx)[n] = (int (*)[n]) _mtrx;
//     int (*mtrx_F)[n] = (int (*)[n]) _mtrx_F;

//     int count = 0;
//     int i, j, k, l;
//     int a = 0, b = 0;
//     for (k = 0; k < q; k++) {
//         for (l = 0; l < q; l++) {
//             for (i = k * m; i < (k + 1) * m; i++) {
//                 for (j = l * m; j < (l + 1) * m; j++) {
//                     mtrx[i][j] = mtrx_F[a][b] == INF ? 0 : mtrx_F[a][b];
//                     b = ++b == n ? (++a && 0) : b;
//                 }
//             }
//         }
//     }
// }

void graphPrint(int *_mtrx, int n) {
    int (*mtrx)[n] = (int (*)[n]) _mtrx;

    int i, j;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n - 1; j++) {
            printf("%d ", mtrx[i][j]);
        }
        printf("%d\n", mtrx[i][j]);
    }
    fflush(stdout);
}