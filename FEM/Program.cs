using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using System.Diagnostics;

namespace FEM
{
    class Program
    {
        const int subd = 5;


        const int NNO = subd*subd; //number of nodes
        const int NEL = (subd-1)*(subd-1)*2; //number of elements
        const int NCON = (subd-1)*2; //number of boundary lines on which potentials are known
        const int NMAT = 3; //number of dielectric materials

        int[,] KTRI = new int[NEL,3];// array that indicates the node numbers of each element  +
        int[] MAT = new int[NEL]; // indicates the material number of each element +
        int[] RO = new int[NEL]; // indicates the material static change of each element  (+)
        double[] PERM = new double[NMAT]; // permittivities of the NMAT materials +
        double[] X = new double[NNO]; // x coordinates of the NNO node numbers  +
        double[] Y = new double[NNO]; // y coordinates of the NNO node numbers  +
        double[] VI = new double[NCON]; //imposed potentials on the NCON boundary lines
        //int[,] NOCC = new int[NCON,20]; //node numbers at which the potential VI is imposed (maximum 20 nodes per equipotential line)
        double[,] SS = new double[NNO,NNO]; // global matrix of coefficients of the system of equations
        double[] W = new double[NNO]; // vector of node potentials
        double[] VDR = new double[NNO]; // vector of the right-hand side of the matrix

        public void coords(int m, int n)
        {

            int l = 0;
            double x = 0;
            double y = 10;
            double shift = 10.0 / (subd - 1);//greška

            for (int j = 0; j <= n; j++)
            {
                x = 0;
                for (int i = 1; i <= m+1; i++)
                {
                    X[l] = x;
                    Y[l]= y;
                    l++;
                    x += shift;
                }
                y -= shift;
            }
        }

        public void trianglesIndices(int m, int n) {

             int l = 0;

            for (int j = 0; j <= n-1 ; j++ )
            {
                for (int i = 0; i < m; i++)
                {
                    int k = i + j * (m + 1);

                    KTRI[l, 0] = k;
                    KTRI[l, 1] = k + 1 + m + 1;
                    KTRI[l, 2] = k + 1;
                    l++;

                    KTRI[l, 0] = k;
                    KTRI[l, 1] = k + m + 1;
                    KTRI[l, 2] = k + 1 + m + 1;
                    l++;
                }

            }
        }

        public void materials() { 
        
            for (int i = 0; i < NEL; i++) { MAT[i] = 0; }

            MAT[0] = 1;
            //MAT[19]=1;
            //MAT[20]=1;
            //MAT[21]=1;
            //MAT[22]=1;
            /*
            for (int i = 20; i < 30; i++) {

                for (int j = 40; j < 80; j++)
                {
                    MAT[i * 2 * 50 + j] = 1;
                }            
            }

            for (int i = 20; i < 25; i++)
            {

                for (int j = 40; j < 60; j++)
                {
                    MAT[i * 2 * 50 + j] = 0;
                }

            }*/

            PERM[0] = 1;
            PERM[1] = 10  * PERM[0];
            PERM[2] = 0.01 * PERM[0];
        }

        public void form() {

            double[,] S = new double[3,3];
            int[] NAUX = new int[3];

            //DO FOR NEL ELEMENTS
            for (int i = 0; i < NEL; i++)
            {

                int N1 = KTRI[i,0];
                int N2 = KTRI[i,1];
                int N3 = KTRI[i,2];
                int NM = MAT[i];

                //CALCULATE Ql, Q2, Q3, Rl, R2, R3

                double Q1=Y[N2]-Y[N3];
                double Q2=Y[N3]-Y[N1];
                double Q3=Y[N1]-Y[N2];
                double R1=X[N3]-X[N2];
                double R2=X[N1]-X[N3];
                double R3=X[N2]-X[N1];
                double XPERM = PERM[NM];

                //CALCULATE DETERMINANT, TWICE THE AREA OF TRIANGLE

                double DET = X[N2] * Y[N3] + X[N1] * Y[N2] + X[N3] * Y[N1] - X[N1] * Y[N3] - X[N3] * Y[N2] - X[N2] * Y[N1];

                double COEFF = XPERM / DET / 2;
                double ROEL=-RO[i]*DET/6;

                //CALCULATE THE TERMS S(3,3)

                S[0,0]=COEFF*(Q1*Q1+R1*R1);
                S[0,1]=COEFF*(Q1*Q2+R1*R2);
                S[0,2]=COEFF*(Q1*Q3+R1*R3);
                S[1,0]=S[0,1];
                S[1,1]=COEFF*(Q2*Q2+R2*R2);
                S[1,2]=COEFF*(Q2*Q3+R2*R3);
                S[2,0]=S[0,2];
                S[2,1]=S[1,2];
                S[2,2]=COEFF*(Q3*Q3+R3*R3);

                //ASSEMBLE THE S(3,3) INTO THE MATRIX SS(NNO,NNO) AND SOURCE TERM

                NAUX[0]=N1;
                NAUX[1]=N2;
                NAUX[2]=N3;
                for(int k=0;k<3;k++){
                    int kk=NAUX[k];
                    VDR[kk]=VDR[kk]+ROEL;
                    for(int j=0;j<3;j++){
                    
                        int jj=NAUX[j];

                        SS[kk,jj]=SS[kk,jj]+S[k,j];
                    }
                }

            }
        }

        public void boundary() {

            for (int i = 0; i < subd; i++) {

                for (int j = 0; j < NNO; j++) {

                    SS[i, j] = 0;
                    SS[i + NNO - subd, j] = 0;

                }

                SS[i, i] = 1;

                VDR[i] = 100;
                SS[i + NNO-subd, i + NNO-subd] = 1;
                VDR[i + NNO-subd] = 0;
            }
                
        }

        public double[] jacobiSolver(double[,] A, double[] b, int n, int maxIterations, double error) {

            double[] x = initializeStartingSolution(100.0, 0.0, n*n, n);

            int i = 0;

            for (i = 0; i < maxIterations; i++) {
                if (jacobiIteration(A, b, x) < error)
                {
                    break;
                }        
            }

            Console.WriteLine("Iterations: " + i);            
            return x;
        }

        private double jacobiIteration(double[,] A, double[] b, double[] x) { 
            
            int size = x.Length;

            double[] x_next = new double[size];

            for (int i = 0; i < size; i++) {
                x_next[i] = b[i];
                for (int j = 0; j < size; j++)
                {
                    if(i!=j){
                        x_next[i] += -x[j] * A[i, j]; 
                    }
                }
                x_next[i] = x_next[i] / A[i, i];
            }

            double e_norm = 0;

            for (int i = 0; i < size; i++)
            {
                double e = b[i];
                x[i] = x_next[i];
                for (int j = 0; j < size; j++)
                {
                    e -= x_next[j] * A[i, j];
                }
                e_norm += e * e;
            }

            return Math.Sqrt(e_norm);
        } 

        private double[] initializeStartingSolution(double maxVoltage, double minVoltage, int size, int subdivision) {

            double[] x = new double[size];

            double step = (maxVoltage - minVoltage) / (subdivision - 1);
            for (int i = 0; i < size; i++)
            {
                x[i] = 100 - ( i / subdivision) * step;
            }


            return x;
        }


        public double[] solve(double[,] M, double[] V, int n) {

            for (int i = 0; i < n - 1; i++) {

                for (int j = i+1; j < n; j++) {

                    double coeff = M[j, i] / M[i, i];

                    for (int k = i; k < n; k++) {
                        M[j, k] = M[j, k] - coeff * M[i, k]; 
                    }
                    V[j] = V[j] - coeff * V[i];

                }    
            }

            FileSerializer.saveMatrix(NNO, NNO, M, "L");

            for (int i = n-1; i > 0; i--)
            {
                for (int j = i-1; j >= 0; j--)
                {

                    double coeff = M[j, i] / M[i, i];

                    for (int k = i; k < n; k++)
                    {
                        M[j, k] = M[j, k] - coeff * M[i, k];
                    }
                    V[j] = V[j] - coeff * V[i];

                }
            }

            double[] R = new double[n];

            for (int i = 0; i < n; i++) {
                R[i] = V[i] / M[i, i];   
            }

            return R;
        }

        public double analysis(double[,] A) {

            MatrixSerializer.saveMatrix(new Matrix(A), "A");

            double[,] Dinv = new double[NNO, NNO];
            double[,] N    = new double[NNO, NNO];

            for (int i = 0; i < NNO; i++) {

                Dinv[i, i] = 1/A[i, i];
                for (int j = 0; j < NNO; j++)
                {
                    if (i != j)
                        N[i, j] = -A[i, j];
                }
            }

            MatrixSerializer.saveMatrix(new Matrix(Dinv), "Dinv");
            MatrixSerializer.saveMatrix(new Matrix(N), "N");

            Matrix tmp = new Matrix(Dinv).multiplyRight(new Matrix(N));
            MatrixSerializer.saveMatrix(tmp, "tmp");
            return (tmp.frobeniousNorm());
        }

        public double error(Vector R)
        {
            double[] e = new double[NNO];
            double e_norm = 0;
            for (int i = 0; i < NNO; i++)
            {
                e[i] = VDR[i];
                for (int l = 0; l < NNO; l++)
                {
                    e[i] -= R.at(l) * SS[i, l];
                }
                e_norm += e[i] * e[i];
            }

            return Math.Sqrt(e_norm);
        }


        public double error(double[] R) {
            double[] e = new double[NNO];
            double e_norm = 0;
            for (int i = 0; i < NNO; i++)
            {
                e[i] = VDR[i];
                for (int l = 0; l < NNO; l++)
                {
                    e[i] -= R[l] * SS[i, l];
                }
                e_norm += e[i] * e[i];
            }

            return Math.Sqrt(e_norm);
        }

        static void Main(string[] args)
        {
            Program program = new Program();

            Stopwatch stopwatch_before = Stopwatch.StartNew();
            program.trianglesIndices(subd-1, subd-1);
            program.coords(subd-1,subd-1);
            program.materials();
            program.form();
            program.boundary();

            double[,] SS_tmp = new double[NNO, NNO];
            double[] VDR_tmp = new double[NNO];
            for (int i = 0; i < NNO; i++) {
                VDR_tmp[i] = program.VDR[i];
                for (int l = 0; l < NNO; l++)
                {
                    SS_tmp[i, l] = program.SS[i, l];
                }
            }
            stopwatch_before.Stop();
            Console.WriteLine("Startup:" + stopwatch_before.ElapsedMilliseconds);
            FileSerializer.saveMatrix(NNO,NNO,SS_tmp, "A");
            Stopwatch stopwatch = Stopwatch.StartNew();
            double[] R = program.solve(SS_tmp, VDR_tmp, NNO);
            stopwatch.Stop();
            Console.WriteLine("Gauss elimination time:" + stopwatch.ElapsedMilliseconds);
            Console.WriteLine("Error gauss elim: " + program.error(R));

            
            Stopwatch stopwatch3 = Stopwatch.StartNew();
            Vector initial = new Vector(program.initializeStartingSolution(100.0, 0.0, NNO, subd));
            Vector RBJ = new BlockJacobiSolver().solve(new Matrix(program.SS), new Vector(program.VDR), initial, NNO, subd, 1000, 100);
            stopwatch3.Stop();
            Console.WriteLine("Block jacobi time:" + stopwatch3.ElapsedMilliseconds);
            Console.WriteLine("Error block jac: " + program.error(RBJ));            
            
            /*
            Stopwatch stopwatch2 = Stopwatch.StartNew();
            double[] RJ = program.jacobiSolver(program.SS, program.VDR, subd, 4000, 1E-12);
            stopwatch2.Stop();
            Console.WriteLine("Jacobi time: " + stopwatch2.ElapsedMilliseconds);
            Console.WriteLine("Error jac: " + program.error(RJ));
            */

            Console.ReadKey();
        }
    }
}


/*
 int j = 0;
            int k = 0;
            string s = ""; string si = "";
            for (int i = 0; i < NNO; i++)
            {

                s += "(" + k + "," + R[i].ToString("0.00000").Replace(",", ".") + ")";
                si += "(" + j + "," + R[i].ToString("0.00000").Replace(",", ".") + ")";
                k++;
                if ((i + 1) % subd == 0)
                {
                    j++;
                    k = 0;
                }
            }
            
            using (StreamWriter outfile = new StreamWriter(@"C:\Users\tdubravcevic\Downloads\points.txt")) {
                outfile.WriteLine(s);
            }

            using (StreamWriter outfile = new StreamWriter(@"C:\Users\tdubravcevic\Downloads\pointsInv.txt"))
            {
                outfile.WriteLine(si);
            }
 */
