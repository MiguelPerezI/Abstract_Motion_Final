#include <stdarg.h>

#ifndef MATRIXR
#define MATRIXR

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>

using namespace std;

class Matrix {

	public:
		int m, n;
		double * * A;

		void initMatrix(int m, int n) {

			this->m = m;
			this->n = n;
			this->A = (double **) malloc (m * sizeof(double*));
			for (int i = 0; i < m; i++)
				this->A[i] = (double *) malloc (n * sizeof(double));
		}

		double * operator[] (int k) { 
			return A[k];
		 };

		void updateA(int i, int j, double x) {
			A[i][j] = x;
		};

		void addWeight(int i, int j, double W) {
			updateA(i, j, W);
			updateA(j, i, W);
		};

		void zeroMatrix(int M);
		void escMatrix();
		int getM() {return m;};
		int getN() {return n;};
};

typedef vector<double> doubleList;
typedef vector<int> integerList;
typedef vector<int> :: iterator itInt;

class Dijkstra {

	private:

		doubleList dist;
		integerList Q, S;
		itInt IT;
		integerList prev;
		double aux;
		int minQ, w, s, u, v, alt, ss;

	public:

		double initDijkstra(Matrix graph, int begin, int end);

		void updateAux(double min) {aux = min;};
		double getAux() {return aux;};

		//pushing methods
		void pushIntoDist(double x) { dist.push_back(x);};
		void pushIntoPrev(int x) { prev.push_back(x);};
		void pushIntoPrevDijkstra(int x) { 
			if (x != prev.back())
				prev.push_back(x);
		};
		void pushIntoQ(int x) { Q.push_back(x);};

		//getting methods
		double getDist(int i) {return dist[i];};
		int getPrev(int i) {return prev[i];};
		int getQ(int i) {return Q[i];};
		int getQSize() {return Q.size();};

		//update methods
		void updateDist(int i, double x) { if (i < dist.size()) dist[i] = x;};
		void updatePrev(int i, int x) { if (i < prev.size()) prev[i] = x;};
		void updateQ(int i, int x) { if (i < Q.size()) Q[i] = x;};

		//erase methods
		void eraseInQ(int i) {
			//if (i == 0)
				Q.erase (Q.begin()+i);
		};

		void escListDist() {
			if (dist.size() != 0) 

			for (int i = 0; i < dist.size(); i++)
				std::cout << " " << dist[i];

			 else std::cout << "SANTO VACIO, ARRODILLATE" << std::endl;
			std::cout << "\n";
		};

		void escListPrev() {
			if (prev.size() != 0) 

			for (int i = 0; i < prev.size(); i++)
				std::cout << "  ->  " << prev[i];
			
			 else std::cout << "SANTO VACIO, ARRODILLATE" << std::endl;
			std::cout << "\n";
		};

		void escListS() {
			if (S.size() != 0) 

			for (int i = 0; i < S.size(); i++)
				std::cout << "  ->  " << S[i];
			
			 else std::cout << "SANTO VACIO, ARRODILLATE" << std::endl;
			std::cout << "\n";
		};

		void escListQ() {

			std::cout << Q.size() << "  --> Q := {";
			if (Q.size() != 0) {

			for (int i = 0; i < Q.size(); i++)
				std::cout << " " << Q[i];
				std::cout << "}";
			}
			
			 else std::cout << "SANTO VACIO, ARRODILLATE}" << std::endl;
			std::cout << "\n";
		};

		int minVertexQ() {

			minQ = -1;
			aux = 5e+10;
			for (int i = 0; i < getQSize(); i++) {
				//std::cout << "\ncomparing " << getDist(i) << "  "<<getAux()<<std::endl;
				if (getDist(Q[i]) < getAux()) {
					updateAux(getDist(Q[i]));
					minQ = i;
				}
			}
			return minQ;
		};

		int eraseMinQ() {
			//std::cout << "\nWe find minimum index of Q " << minVertexQ() <<std::endl;
			minVertexQ();
			int ret = Q[minQ];
			eraseInQ(minQ);

			return ret;
		};

		int compareRealNumbers(double x, double y) {
			
			if (sqrt((x-y) * (x-y)) < 1e-22) return 1;
			else return 0;

		};

};

#endif
