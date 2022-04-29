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
#include <random>
using namespace std;

template<typename T>
T random(std::vector<T> const &v)
{
	auto it = v.cbegin();
	int random = rand() % v.size();
	std::advance(it, random);

	return *it;
}

int factorial(int n) {

	int ret = 1;
	for (int i = 1; i < n + 1; i++)
		ret *= i;
	return ret;
}

void MatrixInt::initMatrixInt(int m, int n) {

	this->m = m;
	this->n = n;
	this->A = (int **) malloc (m * sizeof(int*));
	for (int i = 0; i < m; i++)
		this->A[i] = (int *) malloc (n * sizeof(int));

}

void MatrixInt::destroyMatrixInt() {

	for (int i = 0; i < this->m; i++) {
		delete this->A[i];
		this->A[i] = NULL;
	}

	delete this->A;
	this->A = NULL;

	this->m = 0;
	this->n = 0;
}

int MatrixInt::getM() {
	return this->m;
}

int MatrixInt::getN() {
	return this->n;
}

int MatrixInt::getA(int i, int j) {
	return this->A[i][j];
}

void MatrixInt::updateA(int i, int j, int x) {
	this->A[i][j] = x;
}

void MatrixInt::makeEqualMatrixInt(MatrixInt mat) {

	if (this->m == mat.getM() && this->n == mat.getN()) {

		for (int i = 0; i < this->m; i++)
			for (int j = 0; j < this->n; j++)
				this->updateA(i, j, mat.getA(i, j));
	} else {
		printf("\n\n--> In makeEqualMatrixInt :: unable to make equal\n--> Please check dimensions\n\n");
	}
}

void MatrixInt::zeroMatrixInt(int m, int n) {

	this->initMatrixInt(m, n);
	for (int i = 0; i < this->m; ++i)
		for (int j = 0; j < this->n; ++j)
			this->A[i][j] = 0;
}

void MatrixInt::escMatrixInt() {

	for (int i = 0; i < this->m; i++) {
		for (int j = 0; j < this->n; j++)
			printf("   %d", this->getA(i, j));
		printf("\n\n");
	}
}

void MatrixInt::escMatrixInt0() {
	printf("\n");
	for (int i = 0; i < this->m; i++) {
		for (int j = 0; j < this->n; j++)
			printf("   %d", this->getA(i, j));
	}
	printf("\n");
}

void MatrixIntList::initMatrixIntList(int m) {

	this->m = m;
	this->matrixInt = (MatrixInt *) malloc (m * sizeof(MatrixInt));

}

void MatrixIntList::initA(int i, MatrixInt matI) {

	this->matrixInt[i].initMatrixInt(matI.getM(), matI.getN());
	
	for (int i0 = 0; i0 < matI.getM(); i0++)
		for (int j0 = 0; j0 < matI.getN(); j0++)
			this->matrixInt[i].A[i0][j0] = matI.getA(i0, j0);		
}

void MatrixIntList::push(MatrixInt matInt) {

	this->matrixInt = (MatrixInt *) realloc (this->matrixInt, sizeof(MatrixInt) * (this->m + 1));
	this->m += 1;
	
	this->initA(this->m-1, matInt);
}

void MatrixIntList::escMatrixIntList() {

	for (int i = 0; i < this->m; i++) {
		printf("%d-->\n", i);
		this->matrixInt[i].escMatrixInt();
	}
}

void MatrixIntList::escMatrixIntList0() {

	for (int i = 0; i < this->m; i++) {
		this->matrixInt[i].escMatrixInt0();
	}
}

//método get para consegir el elemento i-ésimo?


void VectorInt::initVectorInt(int x, int y) {
	this->x = x;
	this->y = y;
}

int VectorInt::getX() {
	return this->x;
}

int VectorInt::getY() {
	return this->y;
}

void VectorInt::updateVectorInt(int newX, int newY) {

	this->x = newX;
	this->y = newY;
}

void VectorInt::updateX(int newX) {
	this->x = newX;
}

void VectorInt::updateY(int newY) {
	this->y = newY;
}

void VectorInt::escVectorInt() {

	printf("(%d, %d)", this->x, this->y);
}

void VectorInt::makeEqual(VectorInt v) {
	
	this->x = v.getX();
	this->y = v.getY();
}

void VectorInt::zeroVectorInt() {
	this->initVectorInt(0, 0);
}

void VectorInt::upOneX() {
	this->updateX(this->getX() + 1);
}

void VectorInt::upOneY() {
	this->updateY(this->getY() + 1);
}

void VectorInt::upOneXY() {
	this->updateX(this->getX() + 1);
	this->updateY(this->getY() + 1);
}

void VectorInt::downOneX() {
	this->updateX(this->getX() - 1);
}

void VectorInt::downOneY() {
	this->updateY(this->getY() - 1);
}

void VectorInt::downOneXY() {
	this->updateX(this->getX() - 1);
	this->updateY(this->getY() - 1);
}

int VectorInt::areEqual(VectorInt v) {
	if (this->x == v.getX() && this->y == v.getY())
		return 1;
	else
		return 0;
}



void MatrixVectorInt::initMatrixVectorInt(int m, int n) {

	this->m = m;
	this->n = n;
	this->A = (VectorInt **) malloc (m * sizeof(VectorInt*));
	for (int i = 0; i < m; i++)
		this->A[i] = (VectorInt *) malloc (n * sizeof(VectorInt));

}

//WARNING ONLY FOR THE SPECIAL CASE WHEN this->m == 1
void MatrixVectorInt::updateRowSize(int i, int newN) {
	this->n = newN;
	this->A[i] = (VectorInt *) realloc (this->A[i], sizeof(VectorInt) * newN);
}

void MatrixVectorInt::updateMatrixVectorIntSize(int newM, int newN) {
	
	this->m = newM;
	this->n = newN;
	this->A = (VectorInt **) realloc (this->A, sizeof(VectorInt*) * newM);
	for (int i = 0; i < m; i++)
		this->A[i] = (VectorInt *) realloc (this->A[i], sizeof(VectorInt) * newN);
}

void MatrixVectorInt::update(MatrixVectorInt M) {
	if (this->n > M.n) {
			
		this->updateMatrixVectorIntSize(M.m, M.n);
		for (int i = 0; i < M.n; i++)
			this->updateA(0, i, M.A[0][i]);

		this->m = M.m;
		this->n = M.n;
	} else {


		if (this->n == M.n) {
				
			for (int i = 0; i < M.n; i++)
				this->updateA(0, i, M.A[0][i]);
		}
	
	
		else {
			this->updateMatrixVectorIntSize(M.m, M.n);
			for (int i = 0; i < this->n; i++)
				this->updateA(0, i, M.A[0][i]);
	
			for (int i = this->n; i < M.n; i++)
				this->A[0][i].initVectorInt(M.A[0][i].getX(), M.A[0][i].getY());
			this->m = M.m;
			this->n = M.n;
		}
	}
}

int MatrixVectorInt::getM() {
	return this->m;
}

int MatrixVectorInt::getN() {
	return this->n;
}

VectorInt MatrixVectorInt::getA(int i, int j) {
	return this->A[i][j];
}

void MatrixVectorInt::updateA(int i, int j, VectorInt x) {
	this->A[i][j].makeEqual(x);
}

void MatrixVectorInt::makeEqualMatrixVectorInt(MatrixVectorInt mat) {

	if (this->m == mat.getM() && this->n == mat.getN()) {

		for (int i = 0; i < this->m; i++)
			for (int j = 0; j < this->n; j++)
				this->updateA(i, j, mat.getA(i, j));
	} else {
		printf("\n\n--> In makeEqualMatrixVectorInt :: unable to make equal\n--> Please check dimensions\n\n");
	}
}

void MatrixVectorInt::zeroMatrixVectorInt(int m, int n) {

	this->initMatrixVectorInt(m, n);
	for (int i = 0; i < this->m; ++i)
		for (int j = 0; j < this->n; ++j)
			this->A[i][j].zeroVectorInt();
}

void MatrixVectorInt::escMatrixVectorInt() {

	//printf("\n");
	for (int i = 0; i < this->m; i++) {
		for (int j = 0; j < this->n; j++) {
			printf("     ");
			this->A[i][j].escVectorInt();
		}
		printf("\n");
	}
}

void MatrixVectorInt::escMatrixVectorInt0() {

	//printf("\n");
	for (int i = 0; i < this->m; i++) {
		for (int j = 0; j < this->n; j++) {
			printf("     ");
			this->A[i][j].escVectorInt();
		}
	}
}

void MatrixVectorInt::updateAX(int i, int j, int X) {
	this->A[i][j].updateX(X);
}

void MatrixVectorInt::updateAY(int i, int j, int X) {
	this->A[i][j].updateY(X);
}

VectorInt MatrixVectorInt::getRandomA() {
	
	////srand(time(NULL));
	int i = rand() % (this->n); 
	return this->A[0][i];
}

void MatrixVectorInt::destroyMatrixVectorInt() {

	for (int i = 0; i < this->m; i++) {
		delete this->A[i];
	}

	delete this->A;

	this->m = 0;
	this->n = 0;
}

//one dimensional case
void MatrixVectorInt::push(VectorInt v) {

	int oldN = this->n;
	this->updateMatrixVectorIntSize(1, oldN + 1);
	this->A[0][oldN].initVectorInt(v.getX(), v.getY());
}


void Simplex::initSimplex(int n, ...) {

	va_list list;
	this->id = 0;
	this->numVertex = n;
	this->dimension = n - 1;
	this->A.initMatrixInt(1, this->numVertex);

	va_start(list, n);

		for (int i = 0; i < n; i++)
			this->A.updateA(0, i, va_arg(list, int));
	va_end(list);
}

void Simplex::initSimplexOne(int n) {

	this->numVertex = n;
	this->dimension = n - 1;
	this->A.initMatrixInt(1, this->numVertex);
}

void Simplex::escSimplex() {

	printf("Simplex(");
	for (int i = 0; i < this->numVertex; i++) {
		if (i < this->numVertex - 1)
			printf("%d, ", this->A.getA(0, i));
		else
			printf("%d", this->A.getA(0, i));
	}
	printf(")");
}



void MatrixSimplex::initMatrixSimplex(int m, int n) {

	this->m = m;
	this->n = n;
	this->A = (Simplex **) malloc (m * sizeof(Simplex*));
	for (int i = 0; i < m; i++)
		this->A[i] = (Simplex *) malloc (n * sizeof(Simplex));

}

int MatrixSimplex::getM() {
	return this->m;
}

int MatrixSimplex::getN() {
	return this->n;
}

Simplex MatrixSimplex::getA(int i, int j) {
	return this->A[i][j];
}

void MatrixSimplex::escMatrixSimplex() {

	for (int i = 0; i < this->m; i++) {
		for (int j = 0; j < this->n; j++){
			printf("\n");
			this->A[i][j].escSimplex();
		}
		printf("\n\n");
	}
}
