using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FEM
{
    class MatrixTest
    {
       
        public static void Mrain(string[] args)
        {

            Matrix A = new Matrix(3, 2, new double[] { 1, 2, 3, 4, 5, 6 });
            Matrix B = new Matrix(2, 3, new double[] { 6, 5, 4, 3, 2, 1 });

            Console.WriteLine(A.multiplyLeft(B).frobeniousNorm());

            Console.ReadKey();
        }

    }
}
