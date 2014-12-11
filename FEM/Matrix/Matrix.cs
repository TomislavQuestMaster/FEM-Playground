using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FEM
{
    class Matrix
    {
        int rows;
        int columns;
        double[,] matrix;

        public int getRows() {
            return rows;
        }

        public int getColumns()
        {
            return columns;
        }

        public double at(int i, int j) {

            if (i < 0 || i >= rows) {
                throw new Exception("Row index out of bounds " + i);
            }
            if (j < 0 || j > columns)
            {
                throw new Exception("Column index out of bounds " + i);
            }

            return matrix[i, j];
        }

        public void set(int i, int j, double value) {
            matrix[i, j] = value;
        }

        public Matrix(int rows, int columns, double[] values)
        {
            this.rows = rows;
            this.columns = columns;
            matrix = new double[rows, columns];
            for (int i = 0; i < this.rows; i++) {
                for (int j = 0; j < this.columns; j++) {
                    matrix[i,j]=values[i*columns+j];
                }
            }

        }


        public Matrix(int rows, int columns) {
 
            this.rows = rows;
            this.columns = columns;
            matrix = new double[rows, columns];
        }

        public Matrix(double[,] matrix)
        {
            this.rows = matrix.GetLength(0);
            this.columns = matrix.GetLength(1);
            this.matrix = matrix;
        }

        public Matrix setZeros(int startRow, int endRow, int startColumn, int endColumn)
        {
            for (int i = startRow; i < endRow; i++)
            {

                for (int j = startColumn; j < endColumn; j++)
                {
                    this.matrix[i, j] = 0;
                }
            }

            return this;
        }

        public Matrix submatrix(int startRow, int endRow, int startColumn, int endColumn) {

            double[,] result = new double[endRow - startRow, endColumn - startColumn];

            for (int i = 0; i < endRow - startRow; i++) {

                for (int j = 0; j < endColumn - startColumn; j++)
                {
                    result[i, j] = this.at(startRow + i, startColumn + j);
                }
            }

            return new Matrix(result);
        }

        public Vector multiply(Vector vector) {

            double[] result = new double[rows];

            for (int i = 0; i < rows; i++) {

                result[i] = 0;

                for (int j = 0; j < columns; j++) {
                
                    result[i] += this.at(i,j)*vector.at(j);
                }

            }

            return new Vector(result);
        }

        public Matrix multiplyRight(Matrix right) {

            if (this.columns != right.rows) {
                throw new Exception(
                    "Matrix not chained, left matrix has " + this.columns + 
                    " columns and right matrix has " + right.rows +" rows!");
            }

            Matrix result = new Matrix(this.rows, right.columns);

            for (int i = 0; i < this.rows; i++) {

                for (int j = 0; j < right.columns; j++) {

                    result.matrix[i, j] = 0;
                    for (int k = 0; k < this.columns; k++) {
                        result.matrix[i,j]+=this.matrix[i,k]*right.matrix[k,j];
                    }
                }
            }

            return result;
        }

        public Matrix multiplyLeft(Matrix left)
        {

            if (left.columns != this.rows)
            {
                throw new Exception(
                    "Matrix not chained, left matrix has " + left.columns +
                    " columns and right matrix has " + this.rows + " rows!");
            }

            Matrix result = new Matrix(left.rows, this.columns);

            for (int i = 0; i < left.rows; i++)
            {

                for (int j = 0; j < this.columns; j++)
                {

                    result.matrix[i, j] = 0;
                    for (int k = 0; k < left.columns; k++)
                    {
                        result.matrix[i, j] += left.matrix[i, k] * this.matrix[k, j];
                    }
                }
            }

            return result;
        }

        public double frobeniousNorm(){
        
            double norm = 0;

            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < columns; j++)
                {
                    norm += matrix[i,j]*matrix[i,j];
                }
            }       
        
            return Math.Sqrt(norm);
        }
    }
}
