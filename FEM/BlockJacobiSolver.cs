using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Diagnostics;

namespace FEM
{
    class BlockJacobiSolver
    {

        public Vector solve(Matrix A, Vector b, Vector x, int n, int p, int maxIterations, double error) {

            /*
            for (int i = 0; i < n; i++) {
                b.set(i, b.at(i) / A.at(i, i));
                for (int j = 0; j < n; j++) {
                    A.set(i, j, A.at(i,j)/A.at(i,i));
                }
            }
            */


            int m = n / p;

            for (int iteration = 0; iteration < maxIterations; iteration++)
            {
                
                for (int i = 0; i < m; i++)
                {
                    //Stopwatch stopwatch = Stopwatch.StartNew();
                    TriDiagonalMatrix D_part = new TriDiagonalMatrix(A.submatrix(i * p, (i + 1) * p, i * p, (i + 1) * p));

                    Vector b_part = b.subVector(i * p, (i + 1) * p)
                                     .subtract(A.submatrix(i * p, (i + 1) * p, 0, A.getColumns())
                                                .setZeros(0, p, i * p, (i + 1) * p)
                                                .multiply(x));

                    //stopwatch.Stop();
                    //Console.WriteLine("Iteration setup" + i+ " :" + stopwatch.ElapsedMilliseconds);
                    //stopwatch = Stopwatch.StartNew();
                    Vector x_part = jacobiSolver(D_part, b_part, x.subVector(i * p, (i + 1) * p), p, 50, 0.0000001);
                    //stopwatch.Stop();
                    //Console.WriteLine("Iteration solve" + i + " :" + stopwatch.ElapsedMilliseconds);
                                    //solve(D_part, b_part, p);
                    //stopwatch = Stopwatch.StartNew();
                    for (int j = 0; j < p; j++)
                    {
                        x.set(i * p + j, x_part.at(j));
                    }
                    //stopwatch.Stop();
                    //Console.WriteLine("Iteration postetup" + i + " :" + stopwatch.ElapsedMilliseconds);

                }
                

            }

            return x;
        }


        private Vector jacobiSolver(TriDiagonalMatrix A, Vector b, Vector x, int n, int maxIterations, double error)
        {
            int i = 0;

            for (i = 0; i < maxIterations; i++)
            {

                x = jacobiIteration(A, b, x);
                double e = b.subtract(A.multiplyRight(x)).norm(); 
                if (e < error)
                {
                    break;
                }
            }

            if (i == maxIterations) {
                Console.WriteLine("Reached max iterations: " + i);
            }
            
            return x;
        }

        private Vector jacobiIteration(TriDiagonalMatrix A, Vector b, Vector x)
        {
            DiagonalMatrix D = new DiagonalMatrix(A);
            TriDiagonalMatrix N = A.subtract(D).negate();

            DiagonalMatrix Dinv = D.inverse();
            Vector tmp = b.add(N.multiplyRight(x));

            return D.inverse().multiplyRight(b.add(N.multiplyRight(x)));
        }

        static void Mrain(string[] args)
        {
            BlockJacobiSolver solver = new BlockJacobiSolver();

            TriDiagonalMatrix A = new TriDiagonalMatrix(3);
            A.set(0, 1, 1); A.set(1, 2, 1);
            A.set(0, 0, 3); A.set(1, 1, 4); A.set(2, 2, 3);
            A.set(1, 0, 1); A.set(2, 1, 1);
            Vector b = new Vector(3); b.set(0, 4.9); b.set(1, 11.5); b.set(2, 10.8);
            Vector x = new Vector(3); x.set(0, 1); x.set(1, 2); x.set(2, 3);

            for (int i = 0; i < 5; i++)
            {
                x = solver.jacobiIteration(A, b, x);
                Console.WriteLine(b.subtract(A.multiplyRight(x)).norm());
            }

            Console.ReadKey();
        }


        private Vector solve(Matrix M, Vector V, int n)
        {

            for (int i = 0; i < n - 1; i++) {

                for (int j = i+1; j < n; j++) {

                    double coeff = M.at(j,i) / M.at(i,i);

                    for (int k = i; k < n; k++) {
                        M.set(j,k,M.at(j,k) - coeff * M.at(i,k)); 
                    }
                    V.set(j, V.at(j) - coeff * V.at(i));

                }    
            }

            for (int i = n-1; i > 0; i--)
            {
                for (int j = i-1; j >= 0; j--)
                {

                    double coeff = M.at(j, i) / M.at(i, i);

                    for (int k = i; k < n; k++)
                    {
                        M.set(j,k, M.at(j,k) - coeff * M.at(i,k));
                    }
                    V.set(j, V.at(j) - coeff * V.at(i));

                }
            }

            double[] R = new double[n];

            for (int i = 0; i < n; i++) {
                R[i] = V.at(i) / M.at(i,i);   
            }

            return new Vector(R);
        }

    }
}
