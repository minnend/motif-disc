package mdisc.kdd2007;

import no.uib.cipr.matrix.nni.*;
import no.uib.cipr.matrix.*;
import kdm.util.*;

/** test lineary algebra speed via mtj */
public class LinAlgTest
{

   public static void main(String args[])
   {
      // test for NNI BLAS
      {
         double[] x = { 1.1, 2.2, 3.3, 4.4 }, y = { 1.1, 2.2, 3.3, 4.4 };
         int n = x.length;
         double dot = BLAS.dot(n, x, 1, y, 1);
         System.out.println("Answer = " + dot);
      }

      // test for NNI LAPACK
      {
         int[] iseed = { 1998, 1999, 2000, 2001 };
         double[] x = new double[10];
         int[] n = { x.length }; // All arguments are arrays due to the CLAPACK calling convention

         LAPACK.laruv(iseed, n, x);

         System.out.println("Answer = ");
         for(int i = 0; i < x.length; i++)
            System.out.print(x[i] + " ");
         System.out.println();
      }

      TimerMS timer = new TimerMS();
      int N = 100;

      System.err.printf("Creating raw data (N=%d)... ", N);
      timer.reset();
      double[][] raw1 = new double[N][N];
      double[][] raw2 = new double[N][N];
      for(int i = 0; i < N; i++)
         for(int j = 0; j < N; j++){
            raw1[i][j] = Library.random();
            raw2[i][j] = Library.random();
         }
      System.err.printf(" done(%dms).\n", timer.time());

      int nIters = 100;
      System.err.printf("MTJ matrix mult (A*B x%d)... ", nIters);
      timer.reset();
      DenseMatrix A = new DenseMatrix(raw1);
      DenseMatrix B = new DenseMatrix(raw2);
      DenseMatrix C = new DenseMatrix(N, N);
      for(int i = 0; i < nIters; i++)
         A.mult(1.0, B, C);
      System.err.printf("done (%dms)\n", timer.time());
   }

}
