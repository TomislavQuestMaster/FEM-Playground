using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

namespace FEM
{
    class MatrixSerializer
    {
        public static void saveMatrix(Matrix matrix, string name)
        {
            using (StreamWriter outfile = new StreamWriter(@"C:\Users\tdubravcevic\Downloads\" + name + ".csv"))
            {
                for (int x = 0; x < matrix.getRows(); x++)
                {
                    string content = "";
                    for (int y = 0; y < matrix.getColumns(); y++)
                    {
                        content += matrix.at(x, y).ToString("0.00000") + ";";
                    }
                    outfile.WriteLine(content);
                }
            }
        }

        public static void saveVector(Vector vector, string name)
        {
            using (StreamWriter outfile = new StreamWriter(@"C:\Users\tdubravcevic\Downloads\" + name + ".csv"))
            {

                for (int x = 0; x < vector.getSize(); x++)
                {
                    string content = "";
                    content += vector.at(x).ToString("0.00") + ";";
                    outfile.WriteLine(content);
                }
            }
        }

    }
}
