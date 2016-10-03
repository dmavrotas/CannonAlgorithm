#include <stdio.h>
#include <math.h>
#include "common.h"
#include "omp.h"

#define PROCESSORS 2
#define THREADS 2
#define VIRTUALM 32
#define VIRTUALN 32
#define VIRTUALK 32

#define A_RN 400
#define A_CN 400
#define B_RN 400
#define B_CN 400

#define BLOCK_SIZE (A_RN / VIRTUALM)

int lcm(a, b) {
    int n;
    for(n = 1;; n++) {
        if(n % a == 0 && n % b == 0)
            return n;
    }
}

void RSync_Process_Blocks(Matrix* matrix, Matrix* subMatrix, int id, int push) {
    int r_pq = 0;
    int c_pq = 0;
    int r_mn = 0;
    int c_mn = 0;
    int ri = 0;
    int rj = 0;
    int roffset = A_RN / VIRTUALM;
    int coffset = A_CN / VIRTUALN;
    int rbegin = 0;
    int rend = 0;
    int cbegin = 0;
    int cend = 0;
    int r = 0;
    int s = 0;

    r_pq = (id / PROCESSORS);
    c_pq = (id % THREADS);

    int step = 0;
    int subr = 0;
    int subc = 0;

    for(r_mn = r_pq; r_mn < VIRTUALM; r_mn += PROCESSORS) {
        rbegin = r_mn * roffset;
        rend = rbegin + roffset;

        for(c_mn = c_pq; c_mn < N; c_mn += Q, step++){
            
            cbegin = c_mn * coffset;
            cend = cbegin + coffset;

            subr = (step / (M/P));
            subc = (step % (N/P));

            for(ri = rbegin, r = 0; ri < rend; ri++, r++){
                for(rj = cbegin, s = 0; rj < cend; rj++, s++){
                    
                    if (push) matrix->items[ri][rj].data = subMatrix->items[r + (subr * roffset)][s + (subc * coffset)].data;
                    else subMatrix->items[r + (subr * roffset)][s + (subc * coffset)].data = matrix->items[ri][rj].data;
                }
            }
        }
    }
}

void RSync_Process_SubMatrix(Matrix *local, Matrix* subMatrix, int row, int col, int push) {
    int i = 0;
    int j = 0; 
    int r = 0; 
    int c = 0;
    int offset = subMatrix->sizex;
    int ibegin = row * offset;
    int jbegin = col * offset;

    for(i = ibegin, r = 0; i < (ibegin + offset); i++, r++){
        for(j = jbegin, c = 0; j < (jbegin + offset); j++, c++){
            if (push) local->items[i][j].data = subMatrix->items[r][c].data;
            else subMatrix->items[r][c].data = local->items[i][j].data;
        }
    }
}

int main(int argc, char*argv[]) {
    int ID = 0;
    int sarn = 0;
    int sacn = 0;
    int sbcn = 0;
    int LCM = 0;
    int t = 0;
    int i = 0;
    int j = 0;
    int l = 0;
    int jp = 0;
    int ip = 0;
    int local_block = 0;
    double t1 = 0;
    double t2 = 0;

    Matrix A = NULL;
    Matrix B = NULL;
    Matrix C = NULL;
    Matrix sa = NULL;
    Matrix sb = NULL;
    Matrix sc = NULL;
    Matrix sunbsa = NULL;    
    Matrix sunbsb = NULL;  
    Matrix sunbsc = NULL;

    A = CreateMatrix(A_RN, A_CN, CreateRandomDataSet(A_RN, A_CN));
    B = CreateMatrix(B_RN, B_CN, CreateRandomDataSet(B_RN, B_CN));
    C = CreateMatrix(A_RN, B_CN, CreateEmptyDataSet(A_RN, B_CN));

    sarn = (A_RN / VIRTUALM) * (VIRTUALM / PROCESSORS);
    sacn = (A_CN / VIRTUALK) * (VIRTUALK / THREADS);
    sbrn = (B_RN / VIRTUALK) * (VIRTUALK / PROCESSORS);
    sbcn = (B_CN / VIRTUALN) * (VIRTUALN / THREADS);

    local_block = BLOCK_SIZE;

    LCM = lcm(PROCESSORS, THREADS);
    LeftShiftMatrix(&A, BLOCK_SIZE, 1);
    UpShiftMatrix(&B, BLOCK_SIZE, 1);

    t1 = omp_get_wtime();

    #pragma omp parallel default(none) shared(A, B, C, sarn, sacn, sbrn, sbcn, LCM, local_block) \
                                       private(sa, sb, sc, id, t, i, j, l, jp, ip, subsa, subsb, subsc) num_threads(PROCESSORS * THREADS)
    {
        ID = omp_get_thread_num();

        sa = CreateMatrix(sarn, sacn, CreateEmptyDataSet(sarn, sacn));
        sb = CreateMatrix(sbrn, sbcn, CreateEmptyDataSet(sbrn, sbcn));
        sc = CreateMatrix(sarn, sbcn, CreateEmptyDataSet(sarn, sbrn));

        for(t = 0; t < LCM; t++) {
            subsa = CreateMatrix(local_block, local_block, CreateEmptyDataSet(local_block, local_block));
            subsb = CreateMatrix(local_block, local_block, CreateEmptyDataSet(local_block, local_block));
            subsc = CreateMatrix(local_block, local_block, CreateEmptyDataSet(local_block, local_block));

            RSync_Process_Blocks(&A, &sa, id, 0);
            RSync_Process_Blocks(&B, &sb, id, 0);

            for(i = 0; i < (VIRTUALM / PROCESSORS); i++){
                for(j = 0; j < (VIRTUALN / THREADS); j++){
                    for(l = 0; l < (VIRTUALK / LCM); l++){
                        jp = (j % VIRTUALK + l * LCM / THREADS) % (VIRTUALK / THREADS);
                        ip = (i % VIRTUALK + l * LCM / PROCESSORS) % (VIRTUALK / PROCESSORS);

                        RSync_Process_SubMatrix(&sc, &subsc, i, j, 0);
                        RSync_Process_SubMatrix(&sa, &subsa, i, jp, 0);
                        RSync_Process_SubMatrix(&sb, &subsb, ip, j, 0);
                        
                        MultiplyMatrix(&subsc, &subsa, &subsb);

                        RSync_Process_Blocks(&sc, &subsc, i, j, 1);
                    }
                }
            }

            RSync_Process_Blocks(&c, &sc, id, 1);

            #pragma omp barrier

            #pragma omp single
            {
                LeftShiftMatrix(&A, BLOCK_SZ, 0);
                LeftShiftMatrix(&B, BLOCK_SZ, 0);
            }
        }
    }

    t2 = omp_get_wtime();

    return 0;
}