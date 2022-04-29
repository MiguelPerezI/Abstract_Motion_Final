#include <stdarg.h>

#include <mpi.h>
#ifndef RCCADD
#define RCCADD
/*for i in {1..5}; do ./a.out  && sleep 1; done*/
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <list>
#include <iterator>
#include "SimplexAbstract.hpp"
#include "SimplexAlpha.hpp"
#include "LocalSearch.hpp"

using namespace std;

int get(list<int> _list, int _i){
    list<int>::iterator it = _list.begin();
    for(int i=0; i<_i; i++){
        ++it;
    }
    return *it;
}

class RCC {

private:
	
	MatrixSimplicialMap mapeosContiguos;
	LocalSearch suchen;
	SimplexAlpha sigma;
	int message[1000], message2[1000];
	int bigOsize, bigOsize2; 
	list <int> O;
	int success;

	int thread, core, coreS;
	MPI_Status status;

	int sendList[10000];
	int sendState;
	int receiveState;
	int receiveList[10000];
	int sendSimplex[10000];
	int sendRandom[10000];
	int receiveSimplex[10000];

public:

	SubComplexJ J, Jprime;

	int successRet() {return success;}

	void destroyMapeos(SimplicialMap a) {
		this->mapeosContiguos.resetMatrixSimplicialMapZero(a);
	}

	void fillMapeos() {

		//cout << "\n-->LocalSearch found " << suchen.psy_reduced.m << " contiguous maps";
		for (int i = 0; i < suchen.psy_reduced.m; i++) {
			this->mapeosContiguos.pushZero(suchen.psy_reduced.map[i]);
		}
	}

	void escMapeos(SubComplexJ J) {
		this->mapeosContiguos.escMatrixSimplicialMapJ(J);
	}

	void initRCC_CORE(SubComplexJ L, Complex K, SimplicialMap a, SimplicialMap b, int M, double r, int argc, char **argv) {

                thread = MPI_Init(&argc, &argv);
                thread = MPI_Comm_rank(MPI_COMM_WORLD, &core);
                thread = MPI_Comm_size(MPI_COMM_WORLD, &coreS);
                srand(time(NULL));

                //printf("\n\n<<<<< RCC >>>>>>\n\n");
                J.initSubComplexJ(1);
                Jprime.initSubComplexJ(1);
                //printf("-->Fetching random SimplexAlpha from L and saving to J.\n");
                J.initA(0, L.getRandomJ());
                J.initZero_Skeleton();
                Jprime.initA(0, L.getRandomJ());
                Jprime.initZero_Skeleton();
                //J.escSubComplexJ();
                sigma.makeCopySimplexAlpha(L.listOfFacets.A[0][0]);

		srand(time(NULL) + core);
		for (int i = 0; i < core; i++) {
			sendList[i] = 0;
			receiveList[i] = 0;
		}

		mapeosContiguos.initMatrixSimplicialMap(1);
                mapeosContiguos.map[0].initSimplicialMapCopy(a);
                suchen.initLocalSearch(L, K, a, b, M, r);
        };

	void initAddFacets_CORE(SubComplexJ JsubL, SubComplexJ L, Complex K, SimplicialMap a, SimplicialMap b, int M, double r,  int argc, char **argv) {

		thread = MPI_Init(&argc, &argv);
                thread = MPI_Comm_rank(MPI_COMM_WORLD, &core);
                thread = MPI_Comm_size(MPI_COMM_WORLD, &coreS);
                srand(time(NULL));

                //cout << "\n\n<<<<< ADDFACETS >>>>>>\n\n";
                J.initSubComplexJ(JsubL.numSimplex);
                //cout << "\n numSimplex := " << J.listOfFacets.n << endl;

                for (int i = 0; i < JsubL.numSimplex; i++)
                        J.initA(i, JsubL.listOfFacets.A[0][i]);

                J.initZero_Skeleton();
                //cout << "\nAddFacets given SubComplexJ:\n";
                //J.escSubComplexJ();
                sigma.makeCopySimplexAlpha(L.listOfFacets.A[0][0]);
		
		srand(time(NULL) + core);
                for (int i = 0; i < core; i++) {
                        sendList[i] = 0;
                        receiveList[i] = 0;
                }

                mapeosContiguos.initMatrixSimplicialMap(1);
                mapeosContiguos.map[0].initSimplicialMapCopy(a);
                suchen.initLocalSearch(L, K, a, b, M, r);
        };

	void resetRCCSubComplexJ(SubComplexJ L) {
		
		int RR = rand()%(L.listOfFacets.n);
		sendState = RR;
                sendRandom[core] = RR;
                for (int an_id = 0; an_id < coreS; an_id++)//////////TELL THE OVERS
                        if (an_id != core)
                                thread = MPI_Send(&sendState, 1, MPI_INT, an_id, 2001, MPI_COMM_WORLD);

                for (int iter = 0; iter < coreS; iter++)/////////////////////RECEIVING NEWS FROM THE REST AND UPDATING
                                if (core != iter) {
                                        thread = MPI_Recv(&receiveState, 1, MPI_INT, iter, 2001, MPI_COMM_WORLD, &status);
                                        sendRandom[iter] = receiveState;
                                }
		J.resetSubComplexJ(L.listOfFacets.A[0][sendRandom[0]]);
	}


	void resetAddFacetsSubComplexJ(SubComplexJ JsubL) {

		J.resetSubComplexJ(JsubL.listOfFacets.A[0][0]);
		for (int i = 1; i < JsubL.numSimplex; i++)
			J.pushSimplexAlpha(JsubL.listOfFacets.A[0][i]);
	}


	void resetSubcomplex_CORE(list <int> aaa0, MatrixSimplexAlpha AA0, int ccc) {

                int ran = rand()%(aaa0.size());
                int RR = get(aaa0, ran);
                sendState = RR;
                sendRandom[core] = RR;
                for (int an_id = 0; an_id < coreS; an_id++)//////////TELL THE OVERS
                        if (an_id != core)
                                thread = MPI_Send(&sendState, 1, MPI_INT, an_id, 2001, MPI_COMM_WORLD);

                for (int iter = 0; iter < coreS; iter++)/////////////////////RECEIVING NEWS FROM THE REST AND UPDATING
                                if (core != iter) {
                                        thread = MPI_Recv(&receiveState, 1, MPI_INT, iter, 2001, MPI_COMM_WORLD, &status);
                                        sendRandom[iter] = receiveState;
                                }

                J.resetSubComplexJ(AA0.A[0][sendRandom[0]]);
                cout << "\n\n\n                         XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX CORE  := " << core  << "    cycle :=  " << ccc << "      " << sendRandom[0] << "\n\n";
        }


	void runRCC_CORE(SubComplexJ L, Complex K, SimplicialMap a, SimplicialMap b, int M, double r, int argc, char **argv) {

                int count = 0;

                cout << "\n<<<<<<<<<< RCC RUNNING >>>>>>>>>>>\n";

		int b1 = 0;
                int nn1 = 0;
                for (int i = 0; i < L.numSimplex; i++) {
                                for (int j = 0; j < J.numSimplex; j++)
                                        b1 += L.listOfFacets.A[0][i].compareSimplexAlpha(J.listOfFacets.A[0][j]);
                                if (b1 == 0) {
                                        if (nn1 == 0)
                                                O.push_back(i);
                                        else
                                                O.push_back(i);
                                        nn1 += 1;
                                        }
                                b1 = 0;
                }
		int cycleRCC = 0;
	
		while (  O.size() > 0) {

			
			//////////////////////////////////////////////////////////////Choosing random O element	
			int i0 = rand()%(O.size());
			int RR0 = get(O, i0);
			
			//////////////////////////////////////////////////////////////Unión of J and {random(O)}
			sigma.updateSimplexAlpha(L.listOfFacets.A[0][RR0]);
			J.pushSimplexAlpha(sigma);

			
			//////////////////////////////////////////////////////////////FINDING A MAP OVER JU{random(O)}
			suchen.ghostLocalSearch(J, K, a, b, M, r);
	
			int callack = 0;

			if (suchen.ghost == 420) {////////////////////////////////////MAP FOUND!! DO:
				sendState = RR0;

				for (int an_id = 0; an_id < coreS; an_id++)///////////TELL THE REST
					if (an_id != core)
					thread = MPI_Send(&sendState, 1, MPI_INT,
						        an_id, 2001, MPI_COMM_WORLD	);
			} else {//////////////////////////////////////////////////////NO MAP FOUND
				J.popSimplexAlpha();

				sendState = -RR0 - 1;
				for (int an_id = 0; an_id < coreS; an_id++)//////////TELL THE OVERS
					if (an_id != core)
                                        thread = MPI_Send(&sendState, 1, MPI_INT,
                                                        an_id, 2001, MPI_COMM_WORLD     );
			}

			sendList[core] = sendState;
			
			for (int iter = 0; iter < coreS; iter++)/////////////////////RECEIVING NEWS FROM THE REST AND UPDATING
				if (core != iter) {
					thread = MPI_Recv(&receiveState, 1, MPI_INT, iter, 2001, MPI_COMM_WORLD, &status);
					sendList[iter] = receiveState;
				}
			
			//////////////////////////////////////////////SHOW CURRENT STATE OF AFFAIRS
			for (int iter = 0; iter < coreS; iter++)
				cout << sendList[iter] << ", ";

			int check = 0;
			for (int iter = 0; iter < coreS; iter++){

				if (sendList[iter] <= -1) {check += 1;}//////////////WE MUST CHOOSE ONLY ONE SIMPLEX
				else {
                                        if (sendList[core] > -1) J.popSimplexAlpha();//SO WE REMOVE WHAT WE UNITED

                                        sigma.updateSimplexAlpha(L.listOfFacets.A[0][sendList[iter]]);//AND WE KEEP ONLY ONE
                                        J.pushSimplexAlpha(sigma);
                                        O.remove(sendList[iter]);//////////////////////MAKING SURE THAT THE O-LIST's REMAIN PARALLEL
                                        iter += 10000000000;
				}


			}

				for (int iter = 0; iter < coreS; iter++) {////////////REMOVING all FAILED ATEMPS
					if (sendList[iter] <= -1)
						O.remove(-1 * (sendList[iter] + 1));
				}



			cycleRCC++;
		}



		O.clear();

        };

	void endRCC_CORE() {
		thread = MPI_Finalize();
		cout << "<<<<<<<<<< RCC END >>>>>>>>>>>\n\n";
	}


	MatrixSimplicialMap getMapeosC() {
		return mapeosContiguos;
	}


	void escRCCResults() {

		cout << "\n\n-->RCC found \n";
		//escMapeos(J);
		cout << "\n Defined over the subcomplex J :\n";
		J.escSubComplexJ();
	}

	SubComplexJ getJ() {
		return J;
	}

	MatrixSimplexAlpha getJListOfFacets() {
		return J.listOfFacets;
	}

	MatrixSimplexAlpha getJ_Facets() {
		return J.listOfFacets;
	}

	//for listOfFacets oneDimensional
	int getJSize() {
		return J.listOfFacets.getN();
	}

	//for listOfFacets oneDimensional
	SimplexAlpha getJ_Facet(int k) {
		return J.listOfFacets.A[0][k];
	}


};


#endif
