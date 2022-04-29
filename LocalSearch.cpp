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
#include "SimplexAlpha.hpp"
#include "LocalSearch.hpp"
#include <vector>
#include <time.h>
#include <random>

using namespace std;

void escMatrixIntJ11(MatrixInt image, SubComplexJ J) {
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

int LocalSearch::randomPhi_vi(SubComplexJ J, Complex K, SimplicialMap a) {

	int ranIndex = rand() % (J.zero_skeleton.n);
	//printf("\n-->i := %d\n", ranIndex);
	this->vi.makeEqual(J.zero_skeleton.getA(0, ranIndex));
	int f_vi = a.evalOnVertex(this->vi);
	//printf("-->vi_TRUE := \n"); this->vi.escVectorInt();printf("\n");

	return f_vi;
}

void LocalSearch::setVkBuild(Complex K) {

	int i = 0;
	while (i < K.n) {

		if (i == this->f_vi) {
			i += 1;
		} else {
			if (i < this->f_vi)
				this->Vk_minus_vert.updateA(0, i, i);
			if (this->f_vi < i)
				this->Vk_minus_vert.updateA(0, i-1, i);
			i += 1;
		}
	}
}

void LocalSearch::ranVkVertex(Complex K) {

	int ranIndex = rand() % (K.n - 1);
	this->randVK = this->Vk_minus_vert.getA(0, ranIndex);
	//printf("-->random ii := %d\n\n", this->randVK);
}


void LocalSearch::initLocalSearch(SubComplexJ J, Complex K, SimplicialMap a, SimplicialMap b, int M, double r) {

	this->PSY.initMatrixSimplicialMap(1);
	this->PSY.initA(0, a);
	this->PSY.m = 0;
	this->Vk_minus_vert.initMatrixInt(1, K.n - 1);
	this->vi.initVectorInt(0, 0);
	int ii = 0;

	if (a.contiguous(b, J, K) == 1) {
		this->PSY.push(b);
	}
	else {

		this->Faux.initSimplicialMapCopy(a);
		f.initSimplicialMapCopy(a);
		int i = 0;
	}

	this->psy_reduced.initMatrixSimplicialMap(1);
	this->psy_reduced.initA(0, a);
}

string LocalSearch::updateLocalSearch(SubComplexJ J, Complex K, SimplicialMap a, SimplicialMap b, int M, double r) {

	printf("\n\n\n-->>LocalSearch start >>>>>>>>>\n");
	//printf("-->We start by defining a list Psy = {a}.\n");
	
	this->destroyPsy_Reduced(a);
	this->vi.updateVectorInt(0, 0);
	int ii = 0;
	string ret = "";

	//printf("-->Our parameters are given by M := %d and r := %lg\n", M, r);
	this->PSY.pushZero(a);
	if (a.contiguous(b, J, K) == 1) {
		cout << "\n-->TRIVIAL CASE FOR LOCALSEARCH\n\n";
		this->PSY.push(b);
	}
	else {
 
		//printf("\n-------> a and b not Contiguous\n");
		this->Faux.updateSimplicialMapCopy(a);
		f.updateSimplicialMapCopy(a);
		int i = 0;
			int ranIndex111 = rand() % 100;
			int ii = 0;
			while (i < M) {
				
				f.updateSimplicialMapCopy(Faux);
				this->f_vi = this->randomPhi_vi(J, K, Faux);
				this->setVkBuild(K);
				this->ranVkVertex(K);
				f.image.A[this->vi.getX()][this->vi.getY()] *= 0;
				f.image.A[this->vi.getX()][this->vi.getY()] += this->randVK;
				this->p = (double) rand()/RAND_MAX;

				if ( (f.contiguous(Faux, J, K) == 1)  && (this->p < r || f.d(b, J, K) < Faux.d(b, J, K))) {
					this->PSY.push(f);
					ii += 1;
					Faux.updateSimplicialMapCopy(f);
					if (this->PSY.map[PSY.m-1].equalJ(b, J) == 1) {

						printf("\n\n\n\n-->//////////////////////////////\n\n//////////////////////////////\n\nContinuity %d-Chain Found\n\n//////////////////////////////\n\n//////////////////////////////\n\n\n\n", PSY.m);
						printf("\n-->>>>updateLocalSearch found a chain");
							this->reduceChain(J, K);
							ret += this->stringReducedChainResults(J, a);
						break;
					}
				}
				i += 1;
			}

			if (this->PSY.map[PSY.m-1].equalJ(b, J) == 0) {
				//cout << "\n-->LocalSearch failed to find a map in subcomplex J\n";
				this->destroyPsy_Reduced(a);
			}

	}


	this->destroyPSY(a);
	printf("-->>LocalSearch finish >>>>>>>>>>>>>\n");
	return ret;
}
 
void LocalSearch::destroyPSY(SimplicialMap a) {
	this->PSY.resetMatrixSimplicialMapZero(a);
}

void LocalSearch::destroyPsy_Reduced(SimplicialMap a) {
	this->psy_reduced.resetMatrixSimplicialMapZero(a);
}

MatrixInt LocalSearch::getPsy(int i){
	return this->PSY.map[i].image;
}
