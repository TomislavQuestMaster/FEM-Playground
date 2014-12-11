using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FEM
{
    class Vector
    {
        int size;
        double[] vector;

        public Vector(double[] vector) {
            this.size = vector.Length;
            this.vector = vector;
        }

        public Vector(int size) {
            this.size = size;
            this.vector = new double[size];
        }

        public Vector subVector(int start, int end) { 
        
            double[] result = new double[end-start];
            for (int i = 0; i < end - start; i++) {
                result[i] = this.at(start+i);
            }

            return new Vector(result);
        }

        public Vector subtract(Vector rigth) { 
        
            double[] result = new double[size];

            for(int i=0;i<size;i++){
                result[i] = this.at(i) - rigth.at(i);     
            }

            return new Vector(result);
        }

        public Vector add(Vector rigth)
        {

            double[] result = new double[size];

            for (int i = 0; i < size; i++)
            {
                result[i] = this.at(i) + rigth.at(i);
            }

            return new Vector(result);
        }

        public double norm()
        {
            double result = 0.0;
            for (int i = 0; i < size; i++) {
                result += vector[i] * vector[i];
            }
            return Math.Sqrt(result);
        }

        public int getSize()
        {
        
            return size;
        }

        public double at(int i) {

            return vector[i];
        }

        public void set(int i, double value) {

            vector[i] = value;
        }
    }
}
