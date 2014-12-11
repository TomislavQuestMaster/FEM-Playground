using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FEM
{
    class DiagonalMatrix
    {
        int size;
        double[] diagonal;

        public DiagonalMatrix(int size) 
        {
            this.size = size;
            this.diagonal = new double[size];
        }

        public DiagonalMatrix(Matrix matrix) 
        {
            this.size = matrix.getRows();
            this.diagonal = new double[this.size];

            for (int i = 0; i < this.size; i++) 
            {
                diagonal[i] = matrix.at(i,i);
            }
        }
        public DiagonalMatrix(TriDiagonalMatrix matrix)
        {
            this.size = matrix.getSize();
            this.diagonal = new double[this.size];

            for (int i = 0; i < this.size; i++)
            {
                diagonal[i] = matrix.at(i, i);
            }
        }
        public DiagonalMatrix(DiagonalMatrix matrix)
        {
            this.size = matrix.getSize();
            this.diagonal = new double[this.size];

            for (int i = 0; i < this.size; i++)
            {
                diagonal[i] = matrix.at(i, i);
            }
        }


        public int getSize()
        {
            return size;
        }

        public double at(int i, int j)
        {

            if (i != j) 
            {
                return 0.0;
            }
            return diagonal[i];
        }

        public void set(int i, int j, double value) 
        {
            if (i == j) 
            {
                diagonal[i] = value;
            }
        }

        public Vector multiplyRight(Vector vector) 
        {
            Vector result = new Vector(this.getSize());

            for(int i=0;i<this.getSize();i++)
            {
                result.set(i, vector.at(i) * diagonal[i]);
            }

            return result;
        }

        public Matrix multiplyRight(Matrix matrix)
        {
            Matrix result = new Matrix(this.getSize(), this.getSize());

            for (int i = 0; i < this.getSize(); i++)
            {
                for (int j = 0; j < this.getSize(); j++)
                {
                    result.set(i, j, matrix.at(i,j) * diagonal[i]);
                }
            }

            return result;
        }

        public DiagonalMatrix inverse()
        {
            DiagonalMatrix result = new DiagonalMatrix(this);

            for (int i = 0; i < this.getSize(); i++)
            {
                result.diagonal[i] = 1 / result.diagonal[i];
            }

            return result;
        }

    }
}
