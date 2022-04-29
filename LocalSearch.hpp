#include <stdarg.h>

#ifndef LOCALSEARCH
#define LOCALSEARCH

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <list>
#include <iterator>
#include "SimplexAbstract.hpp"
#include "SimplexAlpha.hpp"
using namespace std;

class LocalSearch {
	public:
		SimplicialMap f, Faux;
		MatrixInt Vk_minus_vert;
		MatrixSimplicialMap PSY, psy_reduced;
		VectorInt vi;

		int vA, vB;

		integerList VK_m_vert;
		int f_vi;	
		int randVK;
		double p;
		int suc;

		int ghost;

		void start(SimplicialMap);
		void initLocalSearch(SubComplexJ, Complex, SimplicialMap a, SimplicialMap b, int M, double r);
		string updateLocalSearch(SubComplexJ, Complex, SimplicialMap a, SimplicialMap b, int M, double r);
		int randomPhi_vi(SubComplexJ, Complex K, SimplicialMap a);
		void setVkBuild(Complex);
		void ranVkVertex(Complex);
		void escLocalSearchResults(SubComplexJ);
		MatrixInt getPsy(int i);
		void destroyPsy_Reduced(SimplicialMap a);

		void escPsy_Reduced(SubComplexJ J, SimplicialMap a) {
			this->escReducedChainResults(J, a);
		}

		void reduceChain(SubComplexJ J, Complex K) {
			
			int j = 0;
			this->psy_reduced.pushZero(PSY.map[0]);
		
			while (j != PSY.m-1){
				
				int i;
				i = PSY.m-1;
		
				while (PSY.map[i].contiguous(PSY.map[j], J, K) != 1) 
					i -= 1;
		
				this->psy_reduced.pushZero(PSY.map[i]);
				j = i;
			}
			
		}

	void escReducedChainResults(SubComplexJ J, SimplicialMap a) {
		
		//this->psy_reduced.resetMatrixSimplicialMap(a);
		//this->psy_reduced.m = 0;
		this->psy_reduced.escMatrixSimplicialMapJ(J);
		cout << "\n\n-->> Reduced found " << this->psy_reduced.m << " contiguity chains -->>\n"; 

		

	}

	string stringReducedChainResults(SubComplexJ J, SimplicialMap a) {

                //this->psy_reduced.resetMatrixSimplicialMap(a);
                //this->psy_reduced.m = 0;
		string ret = "\n\nMap\n";
                ret += this->psy_reduced.stringMatrixSimplicialMap(J);
                cout << "\n\n-->> Reduced found " << this->psy_reduced.m << " contiguity chains -->>\n";
		ret += "\n";
		return ret;
        }
	
	void destroyPSY(SimplicialMap);

	void ghostLocalSearch(SubComplexJ J, Complex K, SimplicialMap a, SimplicialMap b, int M, double r) {

	this->ghost = 0;																						
	//cout << "\n\n\n-->>GhostLocalSearch start >>>>>>>>> CODE := " << this->ghost << endl;
	//printf("-->We start by defining a list Psy = {a}.\n");
	
	this->vi.updateVectorInt(0, 0);
	int ii = 0;

	//printf("-->Our parameters are given by M := %d and r := %lg\n", M, r);
	if (a.contiguous(b, J, K) == 1)
		this->ghost = 420;
	else {
 
		//printf("\n-------> a and b not Contiguous\n");
		this->Faux.updateSimplicialMapCopy(a);
		f.updateSimplicialMapCopy(a);
		int i = 0;

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
					ii += 1;
					Faux.updateSimplicialMapCopy(f);
					if (Faux.equalJ(b, J) == 1) {
						this->ghost = 420;
						//printf("\n\n\n\n-->//////////////////////////////\n\n//////////////////////////////\n\nContinuity %d-Chain Found\n\n//////////////////////////////\n\n//////////////////////////////\n\n\n\n", ii);
						//cout << "\n-->>>>GhostLocalSearch found a chain | CODE := " << this->ghost;
						break;
					}
				}
				i += 1;
			}

		}
	//cout << "-->>LocalSearch finish >>>>>>>>>>>>> CODE := " << this->ghost;
	}
};

#endif
