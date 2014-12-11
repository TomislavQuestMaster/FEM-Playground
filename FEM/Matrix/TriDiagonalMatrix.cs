using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FEM
{
    class TriDiagonalMatrix
    {
        int size;
        double[] upperDiagonal;
        double[] diagonal;
        double[] lowerDiagonal;

        public TriDiagonalMatrix(int size)
        {
            this.size = size;
            this.diagonal = new double[size];
            this.upperDiagonal = new double[size-1];
            this.lowerDiagonal = new double[size-1];
        }

        public TriDiagonalMatrix(Matrix matrix)
        {
            this.size = matrix.getRows();
            this.diagonal = new double[this.size];
            this.upperDiagonal = new double[size - 1];
            this.lowerDiagonal = new double[size - 1];

            for (int i = 0; i < this.size; i++)
            {
                diagonal[i] = matrix.at(i, i);
                if (i != 0) {
                    lowerDiagonal[i - 1] = matrix.at(i, i-1);
                    upperDiagonal[i - 1] = matrix.at(i-1, i);
                }

            }
        }

        public TriDiagonalMatrix(TriDiagonalMatrix matrix)
        {
            this.size = matrix.getSize();
            this.diagonal = new double[this.size];
            this.upperDiagonal = new double[size-1];
            this.lowerDiagonal = new double[size-1];

            for (int i = 0; i < this.size; i++)
            {
                diagonal[i] = matrix.at(i, i);
                if (i != 0)
                {
                    lowerDiagonal[i - 1] = matrix.at(i, i - 1);
                    upperDiagonal[i - 1] = matrix.at(i - 1, i);
                }

            }
        }

        public int getSize()
        {
            return size;
        }

        public double at(int i, int j)
        {

            if (i == j)
            {
                return diagonal[i];
            }
            if (i + 1 == j)
            {
                return upperDiagonal[i];
            }
            if (i == j + 1)
            {
                return lowerDiagonal[j];
            }

            return diagonal[i];
        }

        public void set(int i, int j, double value)
        {
            if (i == j)
            {
                diagonal[i] = value;
            }
            if (i + 1 == j)
            {
                upperDiagonal[i] = value;
            }
            if (i == j + 1)
            {
                lowerDiagonal[j] = value;
            }
        }

        public Vector multiplyRight(Vector vector)
        {
            Vector result = new Vector(this.getSize());

            for (int i = 0; i < this.getSize(); i++)
            {
                if (i == 0) {
                    result.set(i, vector.at(i) * diagonal[i] + vector.at(i + 1) * upperDiagonal[i]);
                    continue;
                }
                if (i == this.getSize()-1)
                {
                    result.set(i, vector.at(i) * diagonal[i] + vector.at(i - 1) * lowerDiagonal[i - 1]);
                    continue;
                }
                result.set(i, vector.at(i) * diagonal[i] + vector.at(i-1) * lowerDiagonal[i-1] + vector.at(i+1) * upperDiagonal[i] );
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
                    if (j == 0)
                    {
                        result.set(i, j, matrix.at(i, j) * diagonal[j] + matrix.at(i, j + 1) * upperDiagonal[j]);
                        continue;
                    }
                    if (j == this.getSize() - 1)
                    {
                        result.set(i, j, matrix.at(i, j) * diagonal[j] + matrix.at(i, j - 1) * lowerDiagonal[j - 1]);
                        continue;
                    }
                    result.set(i, j, matrix.at(i, j) * diagonal[j] + matrix.at(i, j - 1) * lowerDiagonal[j - 1] + matrix.at(i, j + 1) * upperDiagonal[j]);
                }
            }

            return result;
        }

        public TriDiagonalMatrix subtract(DiagonalMatrix matrix)
        {
            TriDiagonalMatrix result = new TriDiagonalMatrix(this);

            for (int i = 0; i < this.getSize(); i++)
            {
                result.set(i, i, diagonal[i] - matrix.at(i, i));
            }

            return result;
        }

        public TriDiagonalMatrix negate() {

            TriDiagonalMatrix result = new TriDiagonalMatrix(this);

            for (int i = 0; i < this.size; i++)
            {
                result.diagonal[i] = -diagonal[i];
                if (i != 0)
                {
                    result.lowerDiagonal[i - 1] = -lowerDiagonal[i - 1];
                    result.upperDiagonal[i - 1] = -lowerDiagonal[i - 1];
                }

            }

            return result;
        }

    }
}
