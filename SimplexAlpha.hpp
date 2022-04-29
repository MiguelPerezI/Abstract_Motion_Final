#include <stdarg.h>

#ifndef SIMPLEXALPHA
#define SIMPLEXALPHA

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <list>
#include <iterator>
#include "SimplexAbstract.hpp"
#include "Matrix.hpp"
#include <vector>
#include <sstream>
#include <string>
#include <fstream>

using namespace std;

class SimplexAlpha {
	public:
		int dimension, numVertex;
		MatrixVectorInt A;

		void initSimplexAlpha(int);
		void initSimplexAlphaOne(int, ...);
		void updateSimplexSize(int);
		void initVertex(int, int, int);
		void updateVertex(int, int, int);
		void escSimplexAlpha();
		int compareSimplexAlpha(SimplexAlpha);
		void makeCopySimplexAlpha(SimplexAlpha simplex);
		void updateSimplexAlpha(SimplexAlpha simplex);
		void updateNumVertex(int x);
		void destroySimplexAlpha();

		string stringSimplexAlpha() {

			string ret = "";
			ret += this->A.stringMatrixVectorInt();
			return ret;
		}
};

//Clase para almacenar Simplejos
class MatrixSimplexAlpha {
	public:
		int m, n;
		MatrixInt rowLength;
		SimplexAlpha * * A;

		void initMatrixSimplexAlpha(int, int);
		void updateMatrixSimplexAlphaSize(int newM, int newN);
		void initMatrixSimplexAlphaRows(int);
		void initMatrixSimplexAlphaRowsLength(int, int);

		string stringMatrixSimplexAlpha() {

			string ret = "\n";
			for (int i = 0; i < this->m; i++) {
				//ret += " row := " + to_string(i) + "\n";
				for (int j = 0; j < this->n; j++) {
					ret += "[" + to_string(i) + "]" + "[" + to_string(j) + "]" + "           " + this->A[i][j].stringSimplexAlpha();
				}	
			}

			return ret;
		}

		string stringMatrixSimplexAlpha2() {

			string ret = "\n";
			for (int i = 0; i < this->m; i++) {
				ret += " row := " + to_string(i) + "\n";
				for (int j = 0; j < this->rowLength.getA(0, i); j++)
					ret += " " + this->A[i][j].stringSimplexAlpha(); 
			}

			return ret;
		}

		int getM();
		int getN();
		void initA(int, int, int);
		void updateSizeA(int, int, int);
		SimplexAlpha getA(int, int);
		void escMatrixSimplexAlpha();
		void escMatrixSimplexAlphaOne();
		void escMatrixSimplexAlphaOneP();
		void initAByCopyOf(int i, int j, SimplexAlpha);
		void update(MatrixSimplexAlpha);

		//one dimensional case only
		void push(SimplexAlpha);
		SimplexAlpha getRandomA();
		void popAddFacet(SimplexAlpha); // Este método es solamente para AddFacet
		int verifyEmptyAddFacet(); // Este método es solamente para AddFacet
		void destroyMatrixSimplexAlpha();
		void pushZero(SimplexAlpha simplex);
		void pushZeroTypeInit(SimplexAlpha simplex) {

			if (this->n == 0) {
				this->initAByCopyOf(0, 0, simplex);
				this->n += 1;
			}
			else {
				this->A[0] = (SimplexAlpha *) realloc (this->A[0], sizeof(SimplexAlpha) * (this->n + 1));
				this->n += 1;
				this->initAByCopyOf(0, this->n - 1, simplex);
			}
		
		}

		//for m == 1 only
		void updateN(int newN) {

			//printf("-->Pushing Simplex --< sizeBefore := %d\n", this->n);
			
				this->A[0] = (SimplexAlpha *) realloc (this->A[0], sizeof(SimplexAlpha) * newN);
				this->n = newN;
			
		}

		void resetPCovering();

		void replaceWith(MatrixSimplexAlpha M) {

			if (this->n > M.n) {
				this->updateMatrixSimplexAlphaSize(M.m, M.n);
				for (int i = 0; i < M.n; i++)
					this->A[0][i].updateSimplexAlpha(M.A[0][i]);
			}

			if (this->n == M.n)
				for (int i = 0; i < M.n; i++)
					this->A[0][i].updateSimplexAlpha(M.A[0][i]);

			if (this->n < M.n) {
				
			for (int i = 0; i < this->n; i++)
					this->A[0][i].updateSimplexAlpha(M.A[0][i]);	
			for (int i = this->n; i < M.n; i++)
					this->pushZero(M.A[0][i]);
			}

		}
};



// ostream& operator << (ostream& out, const TensorSimplexAlpha& A){
// 	for(int i = 0; i < A.getM(); i++)
// 		for(int j = 0; j < A.getN(); j++ )
// 			out << A[i][j].escMatrixSimplexAlpha();
// 			out << "\n" << endl;

// 	return out;
// }

//Ya estamos listos para definir un Complejo Simplicial
class Complex {
	public:
		int n;
		int numSimplex;
		MatrixSimplex K;
		Matrix graph;
		Dijkstra dijkstra;

		
		void initComplex(int, int);
		void initAdjMat();
		void initAdjMatA(int, int, int);
		int getAdjMatA(int, int);
		MatrixInt getAdjMat();
		void escAdjMat();
		Simplex getSimplex(int, int);
		Simplex getSimplex(int);
		void escComplex();

		string stringComplex() {

			string ret = "\n 					<<COMPLEX K>>\n";
			ret += K.stringMatrixSimplex();
			return ret;
		}
};

class SimplexProd {
	public:
		int maxSimplex;
		int maxOld;
		MatrixInt path;
		MatrixSimplexAlpha maximalSimplices;

		void multiplySimplices(Simplex, Simplex);
		void multiplySimplicesUpdate(Simplex s0, Simplex s1);
		list <SimplexAlpha> getMaximalSimplices();
		void escListSimplexAlpha(list <SimplexAlpha> facets);
};

class ComplexProduct {
	public:
		Complex K;
		SimplexProd aux;
		MatrixVectorInt zero_skeleton;
		MatrixSimplexAlpha listOfFacets;

		void initComplexProduct(Complex K);
		void escMaximalSimplices();

		string stringComplexProduct() {

			string ret = "\n 					<<COMPLEX PRODUCT KxK>>\n";

			ret += "KxK Zero Skeleton \n";
			ret += zero_skeleton.stringMatrixVectorInt();
			ret += "\n\nKxK List Of Maximal Facets\n";
			ret += listOfFacets.stringMatrixSimplexAlpha2() + "\n";

			return ret;
		}
};

class SubComplexJ {
	public:
		int numSimplex;
		MatrixVectorInt zero_skeleton;
		MatrixSimplexAlpha listOfFacets;

		void initSubComplexJ(int numSimplex);
		void initA(int, SimplexAlpha);
		void escSubComplexJ();
		SimplexAlpha getA(int);
		void initZero_Skeleton();
		void updateZero_Skeleton(SimplexAlpha);

		string stringSubComplexJ() {
			string ret = "\nSubComplex::numSimplex = " + to_string(listOfFacets.n) + "\n";
			ret += listOfFacets.stringMatrixSimplexAlpha();

			return ret;
		}

		void escZero_Skeleton();
		VectorInt getRandomZeroSk();
		void destroySubComplexJ();
		SimplexAlpha getRandomJ();

		int isInJ(int i, int j);
		void pushSimplexAlpha(SimplexAlpha);

		void pushSimplexAlphaZero(SimplexAlpha simplex) {

			//cout << "-->Pushing Simplex := "; simplex.escSimplexAlpha(); cout << "\n";
			this->listOfFacets.pushZero(simplex);
			//cout << "--><<<<<<Push successful >> updating zero_skeleton >>\n";
			//cout << "\n-->update_Process :<<>\n";
			this->updateZero_Skeleton(simplex);
			//cout << "\n\n<>>\n\n-->Update SUCSSFUL\n";

			this->numSimplex *= 0;
			this->numSimplex += this->listOfFacets.n;
		}

		void popSimplexAlpha() {
			zero_skeleton.updateMatrixVectorIntSize(1, 1);
			int newN = listOfFacets.n - 1;
			listOfFacets.updateMatrixSimplexAlphaSize(1, newN);
			for (int i = 0; i < newN; i++)
				updateZero_Skeleton(listOfFacets.A[0][i]);
			numSimplex = newN;
		}

		void resetSubComplexJ(SimplexAlpha simplex) {
			zero_skeleton.updateMatrixVectorIntSize(1, 1);
			listOfFacets.updateMatrixSimplexAlphaSize(1, 1);
			listOfFacets.A[0][0].updateSimplexAlpha(simplex);
			//for (int i = 0; i < 1; i++)
			updateZero_Skeleton(listOfFacets.A[0][0]);
			numSimplex = 1;	
		}
		void popSimplexAlpha(SimplexAlpha);
		void copySubComplexJ(SubComplexJ);
		void copySubComplexJ1(SubComplexJ);

};

class SimplicialMap {
	public:
		MatrixInt image;

		void initSimplicialMap();
		void initSimplicialMapZero() {
			image.initMatrixInt(1, 1);
			image.A[0][0] = 0;
		}
		void destroySimplicialMap();
		void projection1(ComplexProduct);
		void projection2(ComplexProduct);
		void escSimplicialMap();
		void escSimplicialMap0();
		void escSimplicialMapArray();
		void evaluateSimplex(SimplexAlpha);
		void evaluateJ(SubComplexJ);
		int contiguous(SimplicialMap b, SubComplexJ J, Complex K);
		void escListInt(list <int> evalMap);
		void initSimplicialMapCopy(SimplicialMap);
		SimplicialMap retSimplicialMapCopy(SimplicialMap toCopy);
		void updateSimplicialMapCopy(SimplicialMap);
		void updateSimplicialMapImageA(int i, int j, int newInt);

		int evaluateSkeleton_J(int j, SubComplexJ J);
		int d(SimplicialMap b);
		int evalOnVertex(VectorInt v);
		int d_k(int f, int g, Complex K);
		int d(SimplicialMap b, SubComplexJ, Complex K);
		int equal(SimplicialMap b);
		int equalJ(SimplicialMap b, SubComplexJ J);

		void escSimplicialMapJ(SubComplexJ J) {
			
			for (int i = 0; i < this->image.getM(); i++) {
				for (int j = 0; j < this->image.getN(); j++) {
					if ( J.isInJ(i, j) == 1)
						printf("   %d", this->image.getA(i, j));
					else 
						printf("   @");
				}
			}
		}

		void escMatrixIntJ(SubComplexJ J);

		string stringSimplicialMap(SubComplexJ J) {

			string ret = "";
			for (int i = 0; i < image.getM(); i++)
				for (int j = 0; j < image.getN(); j++)
					if (J.isInJ(i, j) == 1) {
						ret += "   " + to_string(image.getA(i, j));
					} else {
						ret += "   @";
					}
			return ret;
		}
};


class MatrixSimplicialMap {
	public:
		int m;
		SimplicialMap * map;

		void initMatrixSimplicialMap(int);
		void initA(int, SimplicialMap);
		void updateA(int, SimplicialMap);
		void push(SimplicialMap);
		void escMatrixSimplicialMap();
		void destroyMatrixSimplicialMap();
		void bigPop();
		void resetMatrixSimplicialMap(SimplicialMap);
		void copyNotByInit(MatrixSimplicialMap M);

		void escMatrixSimplicialMapJ(SubComplexJ J) {
		
			printf("\n\n");
			for(int i = 0; i < this->m; i++){
				this->map[i].escSimplicialMapJ(J);
				printf("\n\n");
			}
		}

		void escMatrixSimplicialMapArray() {
		
			printf("\n");
			for(int i = 0; i < this->m; i++){
				this->map[i].escSimplicialMapArray();
				printf("\n\n");
			}
		}


		string stringMatrixSimplicialMap(SubComplexJ J) {

			string ret;
			ret += "\n\nList of SIMPLICIAL MAPS\n";
			for (int k = 0; k < this->m; k++)
				ret += "  " + map[k].stringSimplicialMap(J) + "\n";

			return ret;
		}

		void resetMatrixSimplicialMapZero(SimplicialMap a) {

			this->map = (SimplicialMap *) realloc (this->map, sizeof(SimplicialMap) * (1));
			this->updateA(0, a);
			this->m = 0;
		}

		//careful, this->n must begin a zero, else it creates an empty simplex
		void pushZero(SimplicialMap psy) {
		
			if (this->m == 0) {
				this->updateA(0, psy);
				this->m += 1;
			}
			else {
				this->map = (SimplicialMap *) realloc (this->map, sizeof(SimplicialMap) * (this->m + 1));
				this->m += 1;
				this->initA(this->m - 1, psy);
			}
		}

		void replaceMatrixSimplicialMap(MatrixSimplicialMap psy) {

			this->resetMatrixSimplicialMap(psy.map[0]);
			this->m = 0;

			for (int k = 0; k < psy.m; k++)
				this->pushZero(psy.map[k]);
		}

		void fillEmptyMapList(MatrixSimplicialMap matMap) {

        for (int i = 0; i < matMap.m; i++) {
                this->pushZero(matMap.map[i]);
            }
        }

};

class TensorSimplicialMap {

	public:

		int m;
		MatrixSimplicialMap * mapList;
		SimplicialMap phi;

		MatrixSimplicialMap getMatMap(int i) {return mapList[i];}

		void fillEmptyMapListTensor(int i, MatrixSimplicialMap matMap) {

			this->mapList[i].fillEmptyMapList(matMap);
		}

		void initTensorSimplicialMap() {

			this->m = 0;
			mapList = (MatrixSimplicialMap * ) malloc (1 * sizeof(MatrixSimplicialMap));
			//mapList[0].initMatrixSimplicialMap(1);
			//phi.initSimplicialMapZero();
			//mapList[0].initA(0, phi);
			
		}

		int getM() {return this->m;}

		void escTensorSimplicialMap() {

			cout << "\n";
			for (int k = 0; k < m; k++) {
				cout << "-------------------->List (" << k << ")\n";
				mapList[k].escMatrixSimplicialMapArray();
			}
		}

		void updateA(int i, MatrixSimplicialMap psy) {
			
			this->mapList[i].replaceMatrixSimplicialMap(psy);
		}

		void initAZero(int i, MatrixSimplicialMap psy) {

			this->mapList[i].initMatrixSimplicialMap(1);
			this->mapList[i].initA(0, psy.map[i]);
			this->mapList[i].m = 0;

			for (int j = 0; j < psy.m; j++)
				this->mapList[i].pushZero(psy.map[j]);
		}

		void initA(int i, MatrixSimplicialMap psy) {


			this->mapList[i].resetMatrixSimplicialMap(psy.map[0]);
			this->mapList[i].m = 0;

			for (int i = 0; i < psy.m; i++)
				this->mapList[i].pushZero(psy.map[i]);
		}

		void pushZero(MatrixSimplicialMap psy) {
		
			if (this->m == 0) {
				this->initAZero(0, psy);
				this->m += 1;
			}

			else {
				this->mapList = (MatrixSimplicialMap *) realloc (this->mapList, sizeof(MatrixSimplicialMap) * (this->m + 1));
				this->m += 1;
				this->initAZero(this->m - 1, psy);
			}
		}

		void resetTensorSimplicialMap(MatrixSimplicialMap psy) {
			
			this->mapList = (MatrixSimplicialMap *) realloc (this->mapList, sizeof(MatrixSimplicialMap) * (1));
			this->m = 0;
			this->updateA(0, psy);
		}

		void updateM(int new_m) {
			this->mapList = (MatrixSimplicialMap *) realloc (this->mapList, sizeof(MatrixSimplicialMap) * (new_m));
		}

		void resetTensorSimplicialMapItem(int i, SimplicialMap a) {
			this->mapList[i].resetMatrixSimplicialMapZero(a);
		}
};


class TensorSimplexAlpha {
	
	public:
		int n, m;
		MatrixSimplexAlpha **A;
		MatrixSimplicialMap maps;


		void initTensorSimplexAlpha(int m, int n) {		
			this->m = m;
			this->n = n;
			this->A = (MatrixSimplexAlpha **) malloc (m * sizeof(MatrixSimplexAlpha*));
			
			for (int i = 0; i < m; i++)
				this->A[i] = (MatrixSimplexAlpha *) malloc (n * sizeof(MatrixSimplexAlpha));
		}

		void initTensorSimplexAlpha0(int m, int n, SimplicialMap a) {		
			this->m = m;
			this->n = n;
			this->A = (MatrixSimplexAlpha **) malloc (m * sizeof(MatrixSimplexAlpha*));
			
			for (int i = 0; i < m; i++)
				this->A[i] = (MatrixSimplexAlpha *) malloc (n * sizeof(MatrixSimplexAlpha));

			maps.initMatrixSimplicialMap(1);
			maps.map[0].initSimplicialMapCopy(a);
		}

		void updateN(int newN) {

			this->A[0] = (MatrixSimplexAlpha *) realloc (this->A[0], sizeof(MatrixSimplexAlpha) * newN);
			this->n = newN;
			
		}

		int getM(){ return this->m;}
		int getN() {return this->n;}
		
		MatrixSimplexAlpha* operator [] (int m) {
			return this->A[m%this->m];
		}

		void initA(int i, int j, MatrixSimplexAlpha L) {
			
			this->A[i][j].initMatrixSimplexAlpha(L.m, L.n);
			for (int i0 = 0; i0 < L.m; i0++)
				for (int j0 = 0; j0 < L.n; j0++)
					this->A[i][j].initAByCopyOf(i0, j0, L.A[i0][j0]);
		}

		void initAList(int i, int j, MatrixSimplexAlpha L, MatrixSimplicialMap rccJ, SimplicialMap a) {

			this->A[i][j].initMatrixSimplexAlpha(L.m, L.n);
			for (int i0 = 0; i0 < L.m; i0++)
				for (int j0 = 0; j0 < L.n; j0++)
					this->A[i][j].initAByCopyOf(i0, j0, L.A[i0][j0]);

			//maps.resetMatrixSimplicialMap(a);
			//maps.copyNotByInit(rccJ);
		}

		// Para el caso unidimensional
		void push(MatrixSimplexAlpha A) {

			//printf("\n--> size := %d", getN());

			this->A[0] = (MatrixSimplexAlpha*) realloc(this->A[0], (this->n+1)*sizeof(MatrixSimplexAlpha));
			initA(0, this->n, A);

			//this->n += 1;
			//printf("\n--> size := %d", getN());
		}

	//	void pushSimplexAlpha(int listOfBags_i, bagNum_j, SimplexAlpha sigma) {
			
		//	this->A[listOfBags_i][bagNum_j].pushZero(sigma);
	//	}		

		void escTensorSimplexAlpha() {
			for (int i = 0; i < this->m; i++) {
				for (int j = 0; j < this->n; j++) {
					printf("\n-->List %d\n", j);
					this->A[i][j].escMatrixSimplexAlpha();
				}
			
			printf("\n");	
			}
		}

		//Here we can acces the size of any partition in Tensor
		//one dimensional case only
		int sizeOfPartition(int i) {
			return this->A[0][i].getN();
		}

		//In the case of Covering each partition consists of a
		//one dimensional array of simplicies, this methods gets
		//into a partition and retrives a any simplex needed
		//i stands for the partition number and j for the simplex at
		//that particuliar position

		SimplexAlpha getPartitionSimplex(int i, int j) {

			//A[0][n] is a list of simplicies
			//so we can acces each simplex by
			//A[0][i].get(0, j) which is a simplex
			return this->A[0][i].getA(0, j);
		}

		//one dimensional Case Only
		void popBack() {
			this->A[0][n-1].destroyMatrixSimplexAlpha();

		}

		//one dimensional case only
		void replaceAt(int i, MatrixSimplexAlpha M) {
			this->A[0][i].replaceWith(M);
		}

		string stringTensorSimplexAlpha(){
			string ret = "";

			for (int i = 0; i < this->m; i++) {
				for (int j = 0; j < this->n; j++){
					ret += "\n-->List := " + to_string(j) + "\n\n";
					ret += this->A[i][j].stringMatrixSimplexAlpha();
				}
			}

			return ret;
		}
};


#endif
