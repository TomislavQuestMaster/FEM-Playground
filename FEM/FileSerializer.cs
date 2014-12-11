using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

namespace FEM
{
    class FileSerializer
    {
        public static void saveMatrix(int m, int n, double[,] matrix, string name)
        {
            double[,] data = new double[m, n];
            using (StreamWriter outfile = new StreamWriter(@"C:\Users\tdubravcevic\Downloads\" + name + ".csv"))
            {
                for (int x = 0; x < m; x++)
                {
                    string content = "";
                    for (int y = 0; y < n; y++)
                    {
                        content += matrix[x, y].ToString("0.00000") + ";";
                    }
                    outfile.WriteLine(content);
                }
            }
        }

        public static void saveVector(int m, double[] vector, string name)
        {
            using (StreamWriter outfile = new StreamWriter(@"C:\Users\tdubravcevic\Downloads\" + name + ".csv"))
            {

                for (int x = 0; x < m; x++)
                {
                    string content = "";
                    content += vector[x].ToString("0.00") + ";";
                    outfile.WriteLine(content);
                }
            }
        }

        public static void saveVectorAsMatrix(int n, double[] vector, string name)
        {
            double[,] matrix = new double[n,n];
            int j = 0;
            int k = 0;
            for (int i = 0; i < vector.Length; i++)
            {
                matrix[j, k] = vector[i];

                k++;
                if ((i + 1) % n == 0)
                {
                    j++;
                    k = 0;
                }
            }
            saveMatrix(n, n, matrix, name);
        }
    }
}
