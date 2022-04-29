#include "SimplexAlpha.hpp"
#include <iostream>
#include <cstdlib>
#include <math.h>
#include <string>
#include <list>
#include <iterator> 
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include "SimplexAbstract.hpp"
#include <vector>
#include <time.h>
#include "Matrix.hpp"
#include <random>
using namespace std;

/*dK*/

int numOfMaxSimplex(int p, int q) {
	return factorial(p + q) / (factorial(p) * factorial(q));
}

void SimplexAlpha::initSimplexAlpha(int n) {

	this->numVertex = n;
	this->A.initMatrixVectorInt(1, n);
}

void SimplexAlpha::initSimplexAlphaOne(int n, ...) {

	this->numVertex = n;
	this->dimension = n - 1;
	this->A.initMatrixVectorInt(1, n);

	va_list list;
	va_start(list, n);
		for (int i = 0; i < n; i++)
			this->A.updateA(0, i, va_arg(list, VectorInt));
	va_end(list);
}

void SimplexAlpha::destroySimplexAlpha() {

	this->dimension = 0;
	this->numVertex = 0;
	this->A.destroyMatrixVectorInt();
	delete this;
}

void SimplexAlpha::escSimplexAlpha() {
	this->A.escMatrixVectorInt();
}

void SimplexAlpha::initVertex(int i, int x, int y) {
	this->A.A[0][i].initVectorInt(x, y);
}

void SimplexAlpha::updateVertex(int i, int x, int y) {
	this->A.A[0][i].updateX(x);
	this->A.A[0][i].updateY(y);
}

void SimplexAlpha::updateSimplexSize(int newSize) {

	this->dimension = newSize - 1;
	this->numVertex = newSize;
	this->A.updateMatrixVectorIntSize(1, newSize);
}

void SimplexAlpha::updateSimplexAlpha(SimplexAlpha simplex) {

	if (this->A.n > simplex.A.n) {
		this->A.updateMatrixVectorIntSize(simplex.A.m, simplex.A.n);
		for (int i = 0; i < simplex.A.n; i++) {
			this->A.updateA(0, i, simplex.A.A[0][i]);
		}
	}

	if (this->A.n == simplex.A.n) {
		for (int i = 0; i < simplex.A.n; i++) {
			this->A.updateA(0, i, simplex.A.A[0][i]);
		}
	}

	if (this->A.n < simplex.A.n) {
		this->A.updateMatrixVectorIntSize(simplex.A.m, simplex.A.n);
		for (int i = 0; i < this->A.n; i++) {
			this->A.updateA(0, i, simplex.A.A[0][i]);
		}

		for (int i = this->A.n; i < simplex.A.n; i++) {
			this->A.A[0][i].initVectorInt(simplex.A.A[0][i].getX(), simplex.A.A[0][i].getY());
		}
	}
	
	this->numVertex = this->A.n;
	this->dimension = this->numVertex - 1;
	
}

int SimplexAlpha::compareSimplexAlpha(SimplexAlpha beta) {
	
	int i = 0;
	int ret = 1;
	//printf("\n\n-->Checking if Simplicies are equal  ret := %d", ret);
	if (this->numVertex == beta.numVertex) {
		//printf("\n-->They have at least the same number of vertices");
		while (i < this->numVertex) {
			if (this->A.A[0][i].areEqual(beta.A.A[0][i]) == 0) {
				
				ret = 0;
				//printf("\n-->They arent equal ret := %d", ret);
				break;
			}
			i += 1;
		}
	} else {
		ret = 0;
		//printf("\n-->They dont have the same number of vertices ret := %d", ret);
	}

	return ret;
}

void SimplexAlpha::makeCopySimplexAlpha(SimplexAlpha simplex) {

	this->initSimplexAlpha(simplex.numVertex);
	for (int i = 0; i < this->numVertex; i++) 
		this->A.A[0][i].initVectorInt(simplex.A.A[0][i].getX(), simplex.A.A[0][i].getY());
	
}

void SimplexAlpha::updateNumVertex(int x) {
	this-> numVertex = x;
}

void MatrixSimplexAlpha::initMatrixSimplexAlpha(int m, int n) {

	this->m = m;
	this->n = n;
	this->A = (SimplexAlpha **) malloc (m * sizeof(SimplexAlpha*));
	for (int i = 0; i < m; i++)
		this->A[i] = (SimplexAlpha *) malloc (n * sizeof(SimplexAlpha));

}

void MatrixSimplexAlpha::initMatrixSimplexAlphaRows(int m) {

	this->m = m;
	this->A = (SimplexAlpha **) malloc (m * sizeof(SimplexAlpha*));
	this->rowLength.initMatrixInt(1, m);
}

void MatrixSimplexAlpha::initMatrixSimplexAlphaRowsLength(int i, int n) {

	this->rowLength.updateA(0, i, n);
	this->A[i] = (SimplexAlpha *) malloc (n * sizeof(SimplexAlpha));
}

void MatrixSimplexAlpha::updateMatrixSimplexAlphaSize(int newM, int newN) {
	
	this->m = newM;
	this->n = newN;
	this->A = (SimplexAlpha **) realloc (this->A, sizeof(SimplexAlpha*) * newM);
	for (int i = 0; i < newM; i++)
		this->A[i] = (SimplexAlpha *) realloc (this->A[i], sizeof(SimplexAlpha) * newN);
}

//for this and M Matrices unidimensional only
void MatrixSimplexAlpha::update(MatrixSimplexAlpha M) {
	
	if (this->n > M.n) {
			
		this->updateMatrixSimplexAlphaSize(M.m, M.n);
		for (int i = 0; i < M.n; i++)
			this->A[0][i].updateSimplexAlpha(M.A[0][i]);
	} 
	else {

			if (this->n != M.n)
				this->updateMatrixSimplexAlphaSize(M.m, M.n);
			for (int i = 0; i < this->n; i++)
				this->A[0][i].updateSimplexAlpha(M.A[0][i]);
			for (int i = this->n; i < M.n; i++)
				this->initAByCopyOf(0, i, M.A[0][i]);
		}
}

int MatrixSimplexAlpha::getM() {
	return this->m;
}

int MatrixSimplexAlpha::getN() {
	return this->n;
}

SimplexAlpha MatrixSimplexAlpha::getA(int i, int j) {
	return this->A[i][j];
}


void MatrixSimplexAlpha::escMatrixSimplexAlpha() {

	for (int i = 0; i < this->m; i++) {
		for (int j = 0; j < this->n; j++){
			//printf("\n");
			this->A[i][j].escSimplexAlpha();
		}
		printf("\n");
	}
}

void MatrixSimplexAlpha::escMatrixSimplexAlphaOne() {

	for (int i = 0; i < this->m; i++) {
		for (int j = 0; j < this->rowLength.getA(0, i); j++){
			//printf("\n");
			this->A[i][j].escSimplexAlpha();
		}
		printf("\n");
	}
}

void MatrixSimplexAlpha::escMatrixSimplexAlphaOneP() {

	for (int i = 0; i < this->m; i++) {
		for (int j = 0; j < this->rowLength.getA(0, i); j++){
			printf("-->(%d, %d) := ", i, j);
			this->A[i][j].escSimplexAlpha();
		}
		printf("\n");
	}
}

void MatrixSimplexAlpha::initA(int i, int j, int n) {

	this->A[i][j].initSimplexAlpha(n);
}

void MatrixSimplexAlpha::updateSizeA(int i, int j, int newN) {

	this->A[i][j].updateSimplexSize(newN);
}

void MatrixSimplexAlpha::initAByCopyOf(int i, int j, SimplexAlpha simplex) {
	this->A[i][j].makeCopySimplexAlpha(simplex);
}

//one-dimensional case only, with Matrix completely empty
void MatrixSimplexAlpha::push(SimplexAlpha simplex) {

	//printf("-->Pushing Simplex --< sizeBefore := %d\n", this->n);
	
		this->A[0] = (SimplexAlpha *) realloc (this->A[0], sizeof(SimplexAlpha) * (this->n + 1));
		this->n += 1;
		this->initAByCopyOf(0, this->n - 1, simplex);
	
}


//careful, this->n must begin a zero, else it creates an empty simplex
void MatrixSimplexAlpha::pushZero(SimplexAlpha simplex) {

	if (this->n == 0) {
		this->A[0][0].updateSimplexAlpha(simplex);
		this->n += 1;
	}
	else {
		this->A[0] = (SimplexAlpha *) realloc (this->A[0], sizeof(SimplexAlpha) * (this->n + 1));
		this->n += 1;
		this->initAByCopyOf(0, this->n - 1, simplex);
	}
}

void MatrixSimplexAlpha::resetPCovering() {

	this->updateMatrixSimplexAlphaSize(1, 1);
	this->n = 0;

}

//one-dimensional case only
//void MatrixSimplexAlpha::reset() {
//
//	
//}

// Solamente para AddFacet con MatrixSimplexAlpha.m = 1
void MatrixSimplexAlpha::popAddFacet(SimplexAlpha S) {
	for(int i = 0; i < this->n;  i++)
		if(this->A[0][i].compareSimplexAlpha(S) == 1)
			this->A[0][i].updateNumVertex(0);
}

// Este método es solamente para AddFacet
int MatrixSimplexAlpha::verifyEmptyAddFacet() {
	int  v = 0;
	for(int i = 0; i < this->n; i++)
		if(this->A[0][i].numVertex != 0)
			v += 1;

	if(v == 0) return 1;
	else return 0;
}

//one-dimensional case only
SimplexAlpha MatrixSimplexAlpha::getRandomA() {
	
	////srand(time(NULL));
	//srand((unsigned) time(0));
	int i = rand() % (this->n);
	//printf("-->ranInteger := %d\n", i);
	if (this->A[0][i].numVertex == 0) {
		i = 0;
		
		while (this->A[0][i].numVertex == 0) {
			i += 1;
		}
		 
	}

	return this->A[0][i];
}

void MatrixSimplexAlpha::destroyMatrixSimplexAlpha() {

	for (int i = 0; i < this->m; i++) {
		//for (int j = 0; j < this->n; j++)
			//this->A[i][j].destroySimplexAlpha();
		delete this->A[i];
		//this->A[i] = NULL;
	}

	delete this->A;
	//this->A = NULL;

}



void Complex::initComplex(int n, int num) {

	this->n = num;
	this->numSimplex = n;
	this->K.initMatrixSimplex(1, n);
}

void Complex::initAdjMat() {
	graph.initMatrix(n, n);
}

void Complex::escComplex() {
	printf("\nComplex {\n");
	this->K.escMatrixSimplex();
	printf("}\n");
}

Simplex Complex::getSimplex(int i, int j) {
	return K.A[i][j];
}

Simplex Complex::getSimplex(int j) {
	return K.A[0][j];
}







void SimplexProd::multiplySimplices(Simplex s0, Simplex s1) {

	this->maxSimplex = numOfMaxSimplex(s0.dimension, s1.dimension);
	this->maxOld = numOfMaxSimplex(s0.dimension, s1.dimension);
	this->path.initMatrixInt(this->maxSimplex, s0.dimension + s1.dimension);

	this->path.updateA(0, 0, 0); this->path.updateA(0, 1, 0); this->path.updateA(0, 2, 1); this->path.updateA(0, 3, 1);
	this->path.updateA(1, 0, 0); this->path.updateA(1, 1, 1); this->path.updateA(1, 2, 0); this->path.updateA(1, 3, 1);
	this->path.updateA(2, 0, 0); this->path.updateA(2, 1, 1); this->path.updateA(2, 2, 1); this->path.updateA(2, 3, 0);
	this->path.updateA(3, 0, 1); this->path.updateA(3, 1, 0); this->path.updateA(3, 2, 0); this->path.updateA(3, 3, 1);
	this->path.updateA(4, 0, 1); this->path.updateA(4, 1, 0); this->path.updateA(4, 2, 1); this->path.updateA(4, 3, 0);
	this->path.updateA(5, 0, 1); this->path.updateA(5, 1, 1); this->path.updateA(5, 2, 0); this->path.updateA(5, 3, 0);
	
	//this->path.updateA(0, 0, 0); this->path.updateA(0, 1, 1);
	//this->path.updateA(1, 0, 1); this->path.updateA(1, 1, 0);

	this->maximalSimplices.initMatrixSimplexAlpha(1, this->maxSimplex);
	int simplexNum = 1;
	int pathI;
	int pathJ;

	for (int i = 0; i < this->maxSimplex; i++)
		this->maximalSimplices.initA(0, i, s0.dimension + s1.dimension + 1);

	for (int simplexNum = 0; simplexNum < this->maxSimplex; simplexNum++) {

		pathI = 0;
		pathJ = 0;

		this->maximalSimplices.A[0][simplexNum].A.A[0][0].initVectorInt(s0.A.getA(0, pathI), s1.A.getA(0, pathJ));
		for (int i = 1; i < this->maximalSimplices.A[0][simplexNum].numVertex; i++) {
	
			if (this->path.getA(simplexNum, i-1) == 0) {
				pathI +=1;
			} else pathJ += 1;
			this->maximalSimplices.A[0][simplexNum].A.A[0][i].initVectorInt(s0.A.getA(0, pathI), s1.A.getA(0, pathJ));
		}
	}
	//this->maximalSimplices.escMatrixSimplexAlpha();
}

void SimplexProd::multiplySimplicesUpdate(Simplex s0, Simplex s1) {

	this->maxSimplex = numOfMaxSimplex(s0.dimension, s1.dimension);
	
	this->path.updateA(0, 0, 0); this->path.updateA(0, 1, 0); this->path.updateA(0, 2, 1); this->path.updateA(0, 3, 1);
	this->path.updateA(1, 0, 0); this->path.updateA(1, 1, 1); this->path.updateA(1, 2, 0); this->path.updateA(1, 3, 1);
	this->path.updateA(2, 0, 0); this->path.updateA(2, 1, 1); this->path.updateA(2, 2, 1); this->path.updateA(2, 3, 0);
	this->path.updateA(3, 0, 1); this->path.updateA(3, 1, 0); this->path.updateA(3, 2, 0); this->path.updateA(3, 3, 1);
	this->path.updateA(4, 0, 1); this->path.updateA(4, 1, 0); this->path.updateA(4, 2, 1); this->path.updateA(4, 3, 0);
	this->path.updateA(5, 0, 1); this->path.updateA(5, 1, 1); this->path.updateA(5, 2, 0); this->path.updateA(5, 3, 0);
	
	//this->path.updateA(0, 0, 0); this->path.updateA(0, 1, 1);
	//this->path.updateA(1, 0, 1); this->path.updateA(1, 1, 0);
	
	this->maximalSimplices.updateMatrixSimplexAlphaSize(1, this->maxSimplex);
	int simplexNum = 1;
	int pathI;
	int pathJ;

	for (int i = 0; i < this->maxSimplex; i++)
		this->maximalSimplices.updateSizeA(0, i, s0.dimension + s1.dimension + 1);

	if (this->maxOld < this->maxSimplex) {
		for (int i = this->maxOld; i < this->maxSimplex; i++)
			this->maximalSimplices.initA(0, i, s0.dimension + s1.dimension + 1);		
	}

	for (int simplexNum = 0; simplexNum < this->maxSimplex; simplexNum++) {

		pathI = 0;
		pathJ = 0;

		this->maximalSimplices.A[0][simplexNum].A.A[0][0].updateX(s0.A.getA(0, pathI));
		this->maximalSimplices.A[0][simplexNum].A.A[0][0].updateY(s1.A.getA(0, pathJ));
		for (int i = 1; i < this->maximalSimplices.A[0][simplexNum].numVertex; i++) {
	
			if (this->path.getA(simplexNum, i-1) == 0) {
				pathI +=1;
			} else pathJ += 1;
			this->maximalSimplices.A[0][simplexNum].A.A[0][i].updateX(s0.A.getA(0, pathI));
			this->maximalSimplices.A[0][simplexNum].A.A[0][i].updateY(s1.A.getA(0, pathJ));
		}
	}
	//this->maximalSimplices.escMatrixSimplexAlpha();
	this->maxOld = numOfMaxSimplex(s0.dimension, s1.dimension);
}

void ComplexProduct::initComplexProduct(Complex K) {

	this->zero_skeleton.initMatrixVectorInt(K.n, K.n);
	for (int i = 0; i < K.n; i++)
		for (int j = 0; j < K.n; j++)
			this->zero_skeleton.A[i][j].initVectorInt(i, j);

	this->listOfFacets.initMatrixSimplexAlphaRows(K.numSimplex * K.numSimplex);

	int count = 0;
	for (int i = 0; i < K.numSimplex; i++) {
		for (int j = 0; j < K.numSimplex; j++) {

			if (i == 0 && j == 0) {
				this->aux.multiplySimplices(K.getSimplex(i), K.getSimplex(j));
				this->listOfFacets.initMatrixSimplexAlphaRowsLength(count, this->aux.maximalSimplices.n);

				for (int k = 0; k < this->aux.maximalSimplices.n; k++) {
					this->listOfFacets.A[count][k].initSimplexAlpha(this->aux.maximalSimplices.getA(0, k).numVertex);// = this->aux.maximalSimplices.getA(0, k);				
				
					for (int l = 0; l < this->aux.maximalSimplices.getA(0, k).numVertex; l++) 
						this->listOfFacets.A[count][k].initVertex(l, this->aux.maximalSimplices.getA(0, k).A.A[0][l].getX(), this->aux.maximalSimplices.getA(0, k).A.A[0][l].getY());		 
				}
				count += 1;
			} else {

				this->aux.multiplySimplicesUpdate(K.getSimplex(i), K.getSimplex(j));
				this->listOfFacets.initMatrixSimplexAlphaRowsLength(count, this->aux.maximalSimplices.n);

				for (int k = 0; k < this->aux.maximalSimplices.n; k++) {
					this->listOfFacets.A[count][k].initSimplexAlpha(this->aux.maximalSimplices.getA(0, k).numVertex);// = this->aux.maximalSimplices.getA(0, k);				
				
					for (int l = 0; l < this->aux.maximalSimplices.getA(0, k).numVertex; l++) 
						this->listOfFacets.A[count][k].initVertex(l, this->aux.maximalSimplices.getA(0, k).A.A[0][l].getX(), this->aux.maximalSimplices.getA(0, k).A.A[0][l].getY());		 
				}
				count += 1;

			}

		}
	}
}


void SubComplexJ::initSubComplexJ(int numSimplex) {

	this->numSimplex = numSimplex;
	this->listOfFacets.initMatrixSimplexAlpha(1, numSimplex);
}

void SubComplexJ::destroySubComplexJ() {

	this->zero_skeleton.destroyMatrixVectorInt();
	this->listOfFacets.destroyMatrixSimplexAlpha();
	this->numSimplex = 0;
}

void SubComplexJ::initA(int i, SimplexAlpha simplex) {
	
	this->listOfFacets.A[0][i].initSimplexAlpha(simplex.numVertex);
	for (int k = 0; k < simplex.numVertex; k++)
		this->listOfFacets.A[0][i].initVertex(k, simplex.A.A[0][k].getX(), simplex.A.A[0][k].getY());
}


SimplexAlpha SubComplexJ::getA(int i) {
	return this->listOfFacets.A[0][i];
}

SimplexAlpha SubComplexJ::getRandomJ() {
	return this->listOfFacets.getRandomA();
}

void SubComplexJ::copySubComplexJ(SubComplexJ J) {
	
	this->initSubComplexJ(J.numSimplex);
	for(int i = 0; i < J.numSimplex; i++)
		this->listOfFacets.A[0][i].makeCopySimplexAlpha(J.listOfFacets.A[0][i]);
	this->initZero_Skeleton();

	//printf("\n-->done Copy");
}

void SubComplexJ::copySubComplexJ1(SubComplexJ J) {
	
	printf("\n-----(1)");
	this->listOfFacets.update(J.listOfFacets);
	printf("\n-----(2)");
	this->zero_skeleton.update(J.zero_skeleton);
	printf("\n-----(3)");
	this->numSimplex = J.listOfFacets.n;
	printf("\n-----(4)");
}

void SubComplexJ::escSubComplexJ() {

	if (this->listOfFacets.n == 0)
		cout << "-->SubComplexJ info :: { } Empty Set";
	if (this->listOfFacets.n > 0) {
		cout << "-->SubComplexJ info :: {\n";
		this->listOfFacets.escMatrixSimplexAlpha();
		cout << "}\n\n";
	}
}

void SubComplexJ::initZero_Skeleton() {

	//printf("\n      -->initZero_Skeleton memory again");
	this->zero_skeleton.initMatrixVectorInt(1, this->listOfFacets.A[0][0].numVertex);

	//printf("\n      -->Adding beginning facet again");
	for (int i = 0; i < this->listOfFacets.A[0][0].numVertex; i++) {
		this->zero_skeleton.A[0][i].initVectorInt(this->listOfFacets.A[0][0].A.A[0][i].getX(), this->listOfFacets.A[0][0].A.A[0][i].getY());	
	}

	//printf("\n      -->Adding the other Facets");
	for (int i = 1; i < this->listOfFacets.n; i++) { //Here we iterate through all the simplicies s of J
		for (int j = 0; j < this->listOfFacets.A[0][i].numVertex; j++) {//we are in a s_i simplex of J
			int oldN = this->zero_skeleton.n;
			int count = oldN;
			int check = 0;
			for (int k = 0; k < oldN; k++) {  
				if (this->zero_skeleton.A[0][k].areEqual(this->listOfFacets.A[0][i].A.A[0][j]) == 0) {
					check += 0;
				} else {check += 1;}
			}
			if (check == 0) {
				count += 1;
				this->zero_skeleton.updateRowSize(0, count);
				this->zero_skeleton.A[0][count-1].initVectorInt(this->listOfFacets.A[0][i].A.A[0][j].getX(), this->listOfFacets.A[0][i].A.A[0][j].getY());		
			}
		}
	}



}

void SubComplexJ::updateZero_Skeleton(SimplexAlpha simplex) {

	for (int i = 0; i < simplex.numVertex; i++) {

		
					int j = 0;
					int count = 0;
					//printf("\n\n 			-->NumberOfSimplices : = %d", this->zero_skeleton.n);
					while (j < this->zero_skeleton.n) {
						//printf("\n 			-->(%d)Comparing ", j);simplex.A.A[0][i].escVectorInt();printf("  vs. "); this->zero_skeleton.A[0][j].escVectorInt(); 
						if (simplex.A.A[0][i].areEqual(this->zero_skeleton.A[0][j]) == 1) {
							count += 1;
							//printf("\n 			-->FOUND EQUAL PAIR");
						}

						j ++;
					}
					//printf("\n 			-->CYCLE := %d\n", count);
					if (count == 0) {
						//printf("-->pushing");
						this->zero_skeleton.push(simplex.A.A[0][i]);
						//printf("\n-->done pushing");
					}
	}


}

void SubComplexJ::escZero_Skeleton() {

	printf("-->SubComplexJ Zero Skeleton ::\n");
	this->zero_skeleton.escMatrixVectorInt();
}

VectorInt SubComplexJ::getRandomZeroSk() {

	//srand(time(NULL));
	return this->zero_skeleton.getRandomA();
}

int SubComplexJ::isInJ(int i, int j) {

	int ret;

	int count = 0;
	int a , b;

	//Vamos con la idea de que (i, j) no este en J entonces in == 0
	int in = 0;
	while (count < this->zero_skeleton.n) {

		a = this->zero_skeleton.getA(0, count).getX();
		b = this->zero_skeleton.getA(0, count).getY();
		if (i == a && j == b) {
			in = 1;
			break;
		}

		count += 1;
	}

	return in;

}

void SubComplexJ::pushSimplexAlpha(SimplexAlpha simplex) {

	this->listOfFacets.push(simplex);
	this->updateZero_Skeleton(simplex);
	this->numSimplex *= 0;
	this->numSimplex += this->listOfFacets.n;
	//printf("\n-->pushed Sim\n");
}

void ComplexProduct::escMaximalSimplices() {

	int sum = 0;
	for (int i = 0; i < this->listOfFacets.m; i++)
		for (int j = 0; j < this->listOfFacets.rowLength.getA(0, i); j++)
			sum += 1;

	this->listOfFacets.escMatrixSimplexAlphaOneP();	
	printf("\n-->This KxK is defined by %d maximal facets\n\n-->Skeleton\n", sum);
	this->zero_skeleton.escMatrixVectorInt();
}

void SimplicialMap::destroySimplicialMap() {
	this->image.destroyMatrixInt();
}

void SimplicialMap::projection1(ComplexProduct KxK) {

	this->image.initMatrixInt(KxK.zero_skeleton.getM(), KxK.zero_skeleton.getN());
	for (int i = 0; i < KxK.zero_skeleton.getM(); i++)
		for (int j = 0; j < KxK.zero_skeleton.getN(); j++)
			this->image.updateA(i, j, i);

		//this->Equal.initMatrixInt(KxK.zero_skeleton.getM(), KxK.zero_skeleton.getN());
}

void SimplicialMap::projection2(ComplexProduct KxK) {

	this->image.initMatrixInt(KxK.zero_skeleton.getM(), KxK.zero_skeleton.getN());
	for (int i = 0; i < KxK.zero_skeleton.getM(); i++)
		for (int j = 0; j < KxK.zero_skeleton.getN(); j++)
			this->image.updateA(i, j, j);

		//this->Equal.initMatrixInt(KxK.zero_skeleton.getM(), KxK.zero_skeleton.getN());
}

void SimplicialMap::escSimplicialMap() {
	this->image.escMatrixInt();
}

void SimplicialMap::escSimplicialMapArray() {

	for (int i = 0; i < this->image.m; i++)
		for (int j = 0; j < this->image.n; j++)
			cout << " " << this->image.A[i][j];
}

void SimplicialMap::escSimplicialMap0() {
	this->image.escMatrixInt0();
}

//void SimplicialMap::escSimplicialMapJ(SubComplexJ J) {
//	this->image.escMatrixIntJ(J);
//}

void SimplicialMap::evaluateSimplex(SimplexAlpha simplex) {

	//printf("\n\n-->Simplex given ");
	//simplex.escSimplexAlpha();
	//printf("\n-->phi1(simplex) :=  ");
	for (int i = 0; i < simplex.numVertex; i++) {
		printf("%d  ", this->image.getA(simplex.A.A[0][i].getX(), simplex.A.A[0][i].getY()));
	}
}

void escMatrixIntJ(MatrixInt image, SubComplexJ J) {
	printf("\n");
	for (int i = 0; i < image.getM(); i++) {
		for (int j = 0; j < image.getN(); j++)
			if ( J.isInJ(i, j) == 1)
				printf("   %d", image.getA(i, j));
	}
	printf("\n");
}

void escMatrixIntJ0(MatrixInt image, SubComplexJ J) {
	printf("\n");
	for (int i = 0; i < image.getM(); i++) {
		for (int j = 0; j < image.getN(); j++) {
			if ( J.isInJ(i, j) == 1)
				printf("   %d", image.getA(i, j));
			else 
				printf("   @");
		}
		printf("\n");
	}
	printf("\n");
}

int SimplicialMap::evaluateSkeleton_J(int j, SubComplexJ J) {

	int i0 = J.zero_skeleton.A[0][j%(J.zero_skeleton.n)].getX();
	int j0 = J.zero_skeleton.A[0][j%(J.zero_skeleton.n)].getY();

	return this->image.getA(i0, j0);
}


int SimplicialMap::contiguous(SimplicialMap b, SubComplexJ J, Complex K) {

	list <int> evalMap;
	//printf("\n\n-->Contiguous:\n");
	//printf("-->Scanning the subComplex J \n\n");
	
	int i = 0;
	int result = 0;
	while (i <= J.numSimplex - 1){
		//printf("-->phi1("); J.getA(i).escSimplexAlpha(); printf(") := ");
		//printf("projection_1 := "); this->evaluateSimplex(J.getA(i));printf("\n");
		for (int k = 0; k < J.getA(i).numVertex; k++)
			evalMap.push_back(this->image.getA(J.getA(i).A.A[0][k].getX(), J.getA(i).A.A[0][k].getY()));
		
		//printf("\n-->phi2("); J.getA(i).escSimplexAlpha(); printf(") := ");
		//printf("projection_2 := "); b.evaluateSimplex(J.getA(i));printf("\n");
		for (int k = 0; k < J.getA(i).numVertex; k++)
			evalMap.push_back(b.image.getA(J.getA(i).A.A[0][k].getX(), J.getA(i).A.A[0][k].getY()));
		
		evalMap.sort();
		evalMap.unique();

		int i0 = 0;
		int contiguous = 0;
		while (i0 < K.numSimplex) {

					if (evalMap.size() == K.getSimplex(i0).numVertex) {
						int check = 0;
		
		
								list <int> :: iterator it;
								it = evalMap.begin();
		
								for (int l = 0; l < evalMap.size(); l++) {
				
									if (*it == K.getSimplex(i0).A.getA(0, l)) check += 0;
									else check += 1;
				
									advance(it, 1);
								}
		
								if (check == 0) {
									contiguous = 111;
									//printf("-->Contiguous face\n");
									break;
								}
					}

			i0 += 1;
		}

		//printf("Union(projection_1, projection_2) := "); this->escListInt(evalMap);
		if (contiguous == 111) {
			//printf("\n-->EvalSimplex is in K");
		} else {
			//printf("-->EvalSimplex not in K\n");
			result += 1;
		}
		//printf("\n\n\n");

		evalMap.clear();
		i += 1;
	} 

	if (result == 0) return 1;
	else return 0;

}


list <SimplexAlpha> SimplexProd::getMaximalSimplices() {

	list <SimplexAlpha> facets;
	for (int i = 0; i < this->maximalSimplices.n; i++) {
		facets.push_back(this->maximalSimplices.getA(0, i));
	}

	return facets;
}

void SimplexProd::escListSimplexAlpha(list <SimplexAlpha> facets) {
	
	list <SimplexAlpha> :: iterator it;
	it = facets.begin();

	for (int i = 0; i < facets.size(); i++) {
		it->escSimplexAlpha();
		advance(it, 1);
	}
}

void SimplicialMap::escListInt(list <int> evalMap) {
	
	list <int> :: iterator it;
	it = evalMap.begin();

	for (int i = 0; i < evalMap.size(); i++) {
		//printf(" %d", it);
		cout << "  " << '\t' << *it; 
		advance(it, 1);
	}
}

void SimplicialMap::initSimplicialMapCopy(SimplicialMap toCopy) {

	this->image.initMatrixInt(toCopy.image.getM(), toCopy.image.getN());
	this->image.makeEqualMatrixInt(toCopy.image);	
}

SimplicialMap SimplicialMap::retSimplicialMapCopy(SimplicialMap toCopy) {

	SimplicialMap ret;
	ret.image.initMatrixInt(toCopy.image.getM(), toCopy.image.getN());
	ret.image.makeEqualMatrixInt(toCopy.image);	
	return ret;
}

void SimplicialMap::updateSimplicialMapCopy(SimplicialMap toCopy) {
	//this->image.makeEqualMatrixInt(toCopy.image);

	for (int i = 0; i < this->image.m; i++)
		for (int j = 0; j < this->image.n; j++) {
			this->image.A[i][j] *= 0;
			this->image.A[i][j] += toCopy.image.A[i][j];
		}
}

void SimplicialMap::updateSimplicialMapImageA(int i, int j, int newInt) {
	this->image.updateA(i, j, newInt);
}

int SimplicialMap::evalOnVertex(VectorInt v) {
	return this->image.getA(v.getX(), v.getY());
}

int SimplicialMap::d_k(int f, int g, Complex K) {
	
	//return dk(K.B, f, g);
	return K.dijkstra.initDijkstra(K.graph, f, g);
}

int SimplicialMap::d(SimplicialMap b, SubComplexJ J, Complex K) {

	int f, g;
	int sum = 0;

	//printf("-->>>>>>>\n");
	for (int i = 0; i < J.zero_skeleton.n; i++) {
		//printf(">>>>>>>>>>>>>>>>>>>"); J.zero_skeleton.getA(0, i).escVectorInt();printf("\n");
		f = this->evalOnVertex(J.zero_skeleton.getA(0, i));
		g = b.evalOnVertex(J.zero_skeleton.getA(0, i));
		sum += this->d_k(f, g, K);
		//printf("\n");

	}

	return sum;
}

int SimplicialMap::equal(SimplicialMap b) {

	
	int equal = 0;
			for (int i = 0; i < b.image.getM(); i++)
				for (int j = 0; j < b.image.getN(); j++)
					if (this->image.getA(i, j) == b.image.getA(i, j))
						equal += 2;
					else 
						equal = 1;

	if (equal %2 == 1)
		return 0;
	else return 1;

}

int SimplicialMap::equalJ(SimplicialMap b, SubComplexJ J) {

	int equal = 0;
			for (int i = 0; i < b.image.getM(); i++)
				for (int j = 0; j < b.image.getN(); j++)
					if (J.isInJ(i, j)) {
						if (this->image.getA(i, j) == b.image.getA(i, j))
							equal += 2;
						else 
							equal = 1;
					}

	if (equal %2 == 1)
		return 0;
	else return 1;

}

void MatrixSimplicialMap::initMatrixSimplicialMap(int m) {

	this->m = m;
	this->map = (SimplicialMap *) malloc (m * sizeof(SimplicialMap));

}

void MatrixSimplicialMap::initA(int i, SimplicialMap matI) {

	this->map[i].initSimplicialMapCopy(matI);	
}

void MatrixSimplicialMap::updateA(int i, SimplicialMap matI) {
	map[i].image.makeEqualMatrixInt(matI.image);
}


void MatrixSimplicialMap::push(SimplicialMap M) {

	this->map = (SimplicialMap *) realloc (this->map, sizeof(SimplicialMap) * (this->m + 1));
	this->m += 1;
	
	this->initA(this->m-1, M);
}

void MatrixSimplicialMap::bigPop() {

	this->map = (SimplicialMap *) realloc (this->map, sizeof(SimplicialMap) * (0));
	this->m = 0;
	
	//this->initA(this->m-1, matInt);	
}

void MatrixSimplicialMap::copyNotByInit(MatrixSimplicialMap M) {

	if (this->m >= M.m) {
		for (int i = 0; i < M.m; i++)
			this->map[i].updateSimplicialMapCopy(M.map[i]);

		this->map = (SimplicialMap *) realloc (this->map, sizeof(SimplicialMap) * (M.m));
		this->m = M.m;
	} else { //if this->m < M

		this->map = (SimplicialMap *) realloc (this->map, sizeof(SimplicialMap) * (M.m));
		for (int i = 0; i < this->m; i++)
			this->map[i].updateSimplicialMapCopy(M.map[i]);

		for (int i = this->m; i < M.m; i++)
			this->map[i].initSimplicialMapCopy(M.map[i]);

		this->m = M.m;
	}

}

void MatrixSimplicialMap::resetMatrixSimplicialMap(SimplicialMap a) {

	this->map = (SimplicialMap *) realloc (this->map, sizeof(SimplicialMap) * (1));
	this->m = 1;
	this->updateA(this->m-1, a);
}

void MatrixSimplicialMap::destroyMatrixSimplicialMap(){
	
	for(int i = 0; i < this->m; i++)
		this->map[i].destroySimplicialMap();

	this->bigPop();
}
