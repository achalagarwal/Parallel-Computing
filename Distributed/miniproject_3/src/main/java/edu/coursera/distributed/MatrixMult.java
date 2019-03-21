package edu.coursera.distributed;

import edu.coursera.distributed.util.MPI;
import edu.coursera.distributed.util.MPI.MPIException;

/**
 * A wrapper class for a parallel, MPI-based matrix multiply implementation.
 */
public class MatrixMult {
    /**
     * A parallel implementation of matrix multiply using MPI to express SPMD
     * parallelism. In particular, this method should store the output of
     * multiplying the matrices a and b into the matrix c.
     *
     * This method is called simultaneously by all MPI ranks in a running MPI
     * program. For simplicity MPI_Init has already been called, and
     * MPI_Finalize should not be called in parallelMatrixMultiply.
     *
     * On entry to parallelMatrixMultiply, the following will be true of a, b,
     * and c:
     *
     *   1) The matrix a will only be filled with the input values on MPI rank
     *      zero. Matrix a on all other ranks will be empty (initialized to all
     *      zeros).
     *   2) Likewise, the matrix b will only be filled with input values on MPI
     *      rank zero. Matrix b on all other ranks will be empty (initialized to
     *      all zeros).
     *   3) Matrix c will be initialized to all zeros on all ranks.
     *
     * Upon returning from parallelMatrixMultiply, the following must be true:
     *
     *   1) On rank zero, matrix c must be filled with the final output of the
     *      full matrix multiplication. The contents of matrix c on all other
     *      ranks are ignored.
     *
     * Therefore, it is the responsibility of this method to distribute the
     * input data in a and b across all MPI ranks for maximal parallelism,
     * perform the matrix multiply in parallel, and finally collect the output
     * data in c from all ranks back to the zeroth rank. You may use any of the
     * MPI APIs provided in the mpi object to accomplish this.
     *
     * A reference sequential implementation is provided below, demonstrating
     * the use of the Matrix class's APIs.
     *
     * @param a Input matrix
     * @param b Input matrix
     * @param c Output matrix
     * @param mpi MPI object supporting MPI APIs
     * @throws MPIException On MPI error. It is not expected that your
     *                      implementation should throw any MPI errors during
     *                      normal operation.
     */
    public static void parallelMatrixMultiply(Matrix a, Matrix b, Matrix c,
            final MPI mpi) throws MPIException {


        double[] al;
        double[] bl;
        double[][] cl;
        int size = mpi.MPI_Comm_size(mpi.MPI_COMM_WORLD);
        int rows = c.getNRows();
        int cols = c.getNCols();
        int rank = mpi.MPI_Comm_rank(mpi.MPI_COMM_WORLD);
        //two ways
        //method1: send full matrix but assign specific indices
        int factor = (rows+size -1)/size;
        int lb = factor*(rank-1);
        int ub = Math.min(factor*(rank) -1, rows-1);
        int i,j,k;
        if(rank==0) {
            mpi.MPI_Bcast(a.getValues(), 0, a.getNCols()*a.getNRows(), 0, mpi.MPI_COMM_WORLD);
            mpi.MPI_Bcast(b.getValues(), 0,b.getNCols()*b.getNRows(), 0, mpi.MPI_COMM_WORLD);
            cl = new double[0][0];
            al = new double[0];
            bl = new double[0];
            for (i = lb; i <= ub; i++) {
                for (j = 0; j < cols; j++) {
                    c.set(i, j, 0.0);
//                    mpi.MPI_Barrier(mpi.MPI_COMM_WORLD);
                }
            }

            for (i = lb; i < ub; i++) {
                for (j = 0; j < cols; j++) {
                    for (k = 0; k < a.getNCols(); k++)
                        c.incr(i - lb, j, a.get(i - lb, k) * b.get(j, k));
                }
            }

            double crecv[][] = new double[size - 1][(factor + 1) * cols];

            for (i = 0; i < size; i++) {
                mpi.MPI_Recv(crecv[i], 0, (factor + 1) * cols, i + 1, 0, mpi.MPI_COMM_WORLD);
                int factr = (int) crecv[i][0];
                int index = 0;
                for (j = 0; j < factor; j++) {
                    for (k = 0; k < cols; k++) {
                        c.set(i * factr + j - 1, k, crecv[i][index++]);
                    }
                }
            }
        }




//            mpi.MPI_Recv();


        else {
            al = new double[a.getNCols()*a.getNRows()];
            bl = new double[b.getNCols()*b.getNRows()];
            cl = new double[factor][cols];
            mpi.MPI_Bcast(al,0,al.length,0,mpi.MPI_COMM_WORLD);
            mpi.MPI_Bcast(bl,0,bl.length,0,mpi.MPI_COMM_WORLD);

            double am[][] = new double[a.getNRows()][a.getNCols()];
            double bm[][] = new double[b.getNRows()][b.getNCols()];

            for (i = 0; i < am.length; i++) {
                for (j = 0; j < am[0].length; j++) {
                    am[i][j] = al[i * am[0].length + j];
                }
            }

            for (i = 0; i < bm.length; i++) {
                for (j = 0; j < bm[0].length; j++) {
                    bm[i][j] = bl[i * bm[0].length + j];
                }
            }


            for (i = lb; i <= ub; i++) {
                for (j = 0; j < cols; j++) {
                    for (k = 0; k < a.getNCols(); k++) {
                        // i+1 row -> 5*6 : 21%6 = 3
                        cl[i - lb][j] += am[i][k] * bm[k][j];
                    }
                }
            }

            double[] linc = new double[cl.length*cl[0].length+1];
            int index = 1;
            linc[0] = factor;
            for(i = 0;i<cl.length;i++){
                for(j =0;j<cl[0].length;j++){
                    linc[index++] = cl[i][j];
                }
            }

            mpi.MPI_Send(linc,0,linc.length,0,0,mpi.MPI_COMM_WORLD);
        }
//                mpi.MPI_Barrier(mpi.MPI_COMM_WORLD);

//                for (i = lb; i <= ub; i++) {
//                    for(j =0;j<c.getNCols();j++){
//                        c.incr(i, j, al(i*, j) * bl(i, j));
//                    }
//
//                }
//            }


//        for ( i = 0; i < c.getNRows(); i++) {
//            for ( j = 0; j < c.getNCols(); j++) {
//                c.set(i, j, 0.0);
//
//                for (k = 0; k < b.getNRows(); k++) {
//                    c.incr(i, j, a.get(i, k) * b.get(k, j));
//                }
//            }
//        }
    }
}
