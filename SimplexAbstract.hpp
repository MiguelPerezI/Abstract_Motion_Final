
#include <stdarg.h>

#ifndef SIMPLEXABSTRACT
#define SIMPLEXABSTRACT

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <list>
#include <iterator>
#include <sstream>
#include <string>
#include <fstream>

using namespace std;

int factorial(int n);

//Clase para almacenar enteros
class MatrixInt {
    public:
    	int m, n;
    	int * * A;

        void initMatrixInt(int, int);
		void destroyMatrixInt();
        void zeroMatrixInt(int, int);
        int getM();
        int getN();
        int getA(int, int);
        void updateA(int, int, int);
        void makeEqualMatrixInt(MatrixInt);
        void escMatrixInt();
        void escMatrixInt0();

        string stringMatrixInt() {

        	string ret = "";
        	for (int i = 0; i < m; i++) {
        		for (int j = 0; j < n; j++)
        			ret += " " + to_string(A[i][j]);
        		ret += "\n";
        	}

        	return ret;
        }

};

class MatrixIntList {
	public:
		int m;
		MatrixInt * matrixInt;

		void initMatrixIntList(int);
		void initA(int, MatrixInt);
		void push(MatrixInt);
		void escMatrixIntList();
		void escMatrixIntList0();
};

//Clase para coordenadas (i, j)
class VectorInt {
	private:
		int x, y;
	public:
		void initVectorInt(int, int);
		int getX();
		int getY();
		void updateVectorInt(int, int);
		void updateX(int);
		void updateY(int);
		void escVectorInt();
		void makeEqual(VectorInt);
		void zeroVectorInt();
		void upOneX();
		void upOneY();
		void upOneXY();
		void downOneX();
		void downOneY();
		void downOneXY();
		void upCero();
		int areEqual(VectorInt);

		string stringVectorInt() {

			string ret = "";
			ret += "("+ to_string(x) + ", " + to_string(y) + ")";
			return ret;
		}
};

//Clase para almacenar elementos (i, j)
class MatrixVectorInt {
	public:
		int m, n;
		VectorInt * * A;

		void updateMatrixVectorIntSize(int, int);
		void updateRowSize(int, int);
		void initMatrixVectorInt(int, int);
        void zeroMatrixVectorInt(int, int);
        int getM();
        int getN();
        VectorInt getA(int, int);
        void updateA(int, int, VectorInt);
        void makeEqualMatrixVectorInt(MatrixVectorInt);
        void escMatrixVectorInt();
        void updateAX(int, int, int);
        void updateAY(int, int, int);
        void destroyMatrixVectorInt();
        void update(MatrixVectorInt);


        //Only for 1xn matrices
        	VectorInt getRandomA();
        	void push(VectorInt);

        void escMatrixVectorInt0();

        string stringMatrixVectorInt() {

        	string ret = "";
        	for (int i = 0; i < m; i++) {
        		for (int j = 0; j < n; j++)
        			ret += " " + A[i][j].stringVectorInt();
        		ret += "\n";
        	}

        	return ret;
        }
};

//Clase que almacena un simplejo abstracto
//Por el momento todo será publico
class Simplex {
	public:
		int dimension, numVertex, id;
		MatrixInt A;
		void initSimplex(int, ...);
		void initSimplexOne(int);
		void escSimplex();

		string stringSimplex() {

			string ret = "";
			ret += A.stringMatrixInt();
			return ret;
		}
};

//Clase para almacenar Simplejos
class MatrixSimplex {
	public:
		int m, n;
		Simplex * * A;
		void initMatrixSimplex(int, int);
		int getM();
		int getN();
		Simplex getA(int, int);
		void escMatrixSimplex();

		string stringMatrixSimplex() {

			string ret = "";
			for (int i = 0; i < m; i++) {
				for (int j = 0; j < n; j++)
					ret += " " + A[i][j].stringSimplex();
				ret += "\n";
			}

			return ret;
		}
};

#endif
