#include <stdarg.h>

#ifndef COVERING
#define COVERING

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <list>
#include <iterator>
#include "SimplexAbstract.hpp"
#include "SimplexAlpha.hpp"
#include "RCC.hpp"

using namespace std;

int getAcomp(list<int> _list, int _i){
    list<int>::iterator it = _list.begin();
    for(int i=0; i<_i; i++){
        ++it;
    }
    return *it;
}

class Covering {

    	public:

        RCC rcc;

        MatrixSimplexAlpha AA;
        MatrixSimplexAlpha P;
        SubComplexJ I;
        list <int> Acomp;

        TensorSimplexAlpha p, pOrder;
        TensorSimplicialMap pMap, pOrderMap;

        void initAcomp() {
            for (int i = 0; i < AA.getN(); i++) 
                Acomp.push_back(i);
        }

        TensorSimplicialMap getPOrderMap() {return pOrderMap;}
        TensorSimplicialMap getPMap() {return pMap;}



        void initCovering(SubComplexJ L, Complex K, SimplicialMap a, SimplicialMap b, int M, double r, int argc, char **argv) {
            cout << "\n\n\n<<<<<<<< Begining Covering >>>>>>>>>> " << endl;
            
            //L.escSubComplexJ();
            P.initMatrixSimplexAlpha(1, 1);
            P.initAByCopyOf(0, 0, L.listOfFacets.A[0][0]);
            P.n = 0;

            AA.initMatrixSimplexAlpha(1, 1);
            AA.initAByCopyOf(0, 0, L.listOfFacets.A[0][0]);

            //std :: cout << "-->We begin by copying the facets of L into A:";
            for(int i = 1; i < L.numSimplex; i++){
                AA.push(L.listOfFacets.A[0][i]);
            }

            
            //for (int i = 0; i < AA.getN(); i++) 
            //    Acomp.push_back(i);

            //cout << "\n-->Initializing auxilary list";
            //a.escListInt(Acomp);
            //cout << endl;


            cout << "\n--> L List Given By :\n";
            AA.escMatrixSimplexAlpha();

            //std :: cout << "\n-->Now we define an empty P = {}";
            p.initTensorSimplexAlpha(1, 1);
            pOrder.initTensorSimplexAlpha(1, 1);
            pMap.initTensorSimplicialMap();
            pOrderMap.initTensorSimplicialMap();
            //p.push(L.listOfFacets);
            //p.push(L.listOfFacets);
            //p.A[0][1].A[0][1].updateSimplexAlpha(L.listOfFacets.A[0][0]);
            //p.escTensorSimplexAlpha();

            //std :: cout << "\n-->Initializing a subComplex I that inherits all of A's facets";
            I.initSubComplexJ(1);
            I.listOfFacets.initAByCopyOf(0, 0, AA.A[0][0]);

            for(int i = 1; i < AA.getN(); i++)
                I.listOfFacets.push(AA.A[0][i]);

            I.initZero_Skeleton();
            
            //I.escSubComplexJ();
	    //
	    
	    

            rcc.initRCC_CORE(L, K, a, b, M, r, argc, argv);
            //cout << "\n\n-->RCC ready for use";
        }


        void setPAuxilary(int ccc) {
            
            //cout << "\n-->Check for facets in J that intersect with facets in L";
            //compareSimplexAlpha(SimplexAlpha)
            if (ccc == 1) {
                //p.push(rcc.getJ_Facets());
                for (int iterJ = 0; iterJ < rcc.getJSize(); iterJ++)
                    P.pushZero(rcc.getJ_Facet(iterJ));
                p.initA(0, 0, P);
                pOrder.initA(0, 0, P);
                //p.escTensorSimplexAlpha();
                //cout << "\n-->Tensor of size 1, pushing J directly";
            }
            else {
                //cout << "\nccc = " << ccc << endl;
                //cout << "\n-->Tensor of size "<< p.n <<"       --> Building P";
                //cout << "\n-->Number of Partitions " << p.n << endl;
                //p.escTensorSimplexAlpha();
                //cout << "\nnewJFacets := \n";
                //rcc.getJ_Facets().escMatrixSimplexAlpha();



                for (int iterJ = 0; iterJ < rcc.getJSize(); iterJ += 1) {
                    int bit = 0;
                    for (int i = 0; i < p.getN(); i++) {
                        //printf("\n-->Partition := %d of size %d", i, p.sizeOfPartition(i)); //could be zero
                        

                        int bb = 0;
                        while (bb < p.sizeOfPartition(i)) {

                                    if (rcc.getJ_Facet(iterJ).compareSimplexAlpha(p.getPartitionSimplex(i, bb)) == 1) {
                                        bit = 1;
                                        i = p.getN() + 1;
                                        break;
                                    }

                            bb += 1;
                        }
                    }

                    if (bit == 0) {
                        //cout << "\n-->No Match, Adding Facet";
                        P.pushZero(rcc.getJ_Facet(iterJ));
                    }
                }


                //cout << "\n-->newP := ";
                //P.escMatrixSimplexAlpha();
    
                //cout << "\n-->Pushing P into Tensor";
                p.push(P);
                p.n += 1;

                pOrder.push(P);
                pOrder.n += 1;

                //cout << "size >>>>>>>>>:= " << p.n;
                //cout << "\n-->Tensor newSize " << p.getN();
            }


            
            //cout << "\n\n\n";
        }

        void resetSubComplexI_FromA() {

            I.resetSubComplexJ(AA.A[0][getAcomp(Acomp, 0)]);
            if (AA.getN() > 0)
                for (int i = 1; i < Acomp.size(); i++)
                    I.pushSimplexAlpha(AA.A[0][getAcomp(Acomp, i)]);
        }

        void eliminateP_From_A() {

            for (int i = 0; i < P.n; i++) {
                //cout << "-----------------------------------> " << i << endl;
                int count = 0;
                while (count < Acomp.size()) {

                    int RR = getAcomp(Acomp, count);
                    //cout << "\n-->Comparing \n";
                    //printf("-->A(%d) := ", RR); AA.getA(0, RR).escSimplexAlpha(); cout << "\n";
                    //printf("-->P(%d) := ", i);P.getA(0, i).escSimplexAlpha(); cout << "\n";
                    if (AA.getA(0, RR).compareSimplexAlpha(P.getA(0, i)) == 1) {
                        //cout << ">>>>>Difference FOUND<<<<< Removing From A:\n";
                        //printf("-->A(%d) := ", RR); AA.getA(0, RR).escSimplexAlpha(); cout << "\n";
                        Acomp.remove(RR);
                        break;
                    }

                    count += 1;
                }
            }

            //cout << "\n>>>>>>>>Elimination Finish<<<<<<<<\n";
        }


        void runCovering(SubComplexJ L, Complex K, SimplicialMap a, SimplicialMap b, int M, double r, int argc, char **argv) {

            cout << "\n\n--> Covering Running";
            cout << "\n-->Running Main Loop";
            int cycle = 1;

            initAcomp();

            while (Acomp.size() != 0) {

                //cout << "\n\n\n -->Cycle " << cycle;

                cout << "\n debug (0)---------------------------------------\n";
                resetSubComplexI_FromA();

                //cout << "\n\n--> reset of Subcomplex I := \n";
                //I.escSubComplexJ(); cout << "\n\n";
                cout << "\n debug (1)---------------------------------------\n";
                if (Acomp.size() > 1) { 

                            cout << "\n debug (2)---------------------------------------\n";
                            rcc.runRCC_CORE(I, K, a, b, M, r, argc, argv);
			    cout << "\n debug (3)---------------------------------------\n";
                            rcc.escRCCResults();
                            cout << "\n debug (4)---------------------------------------\n";
                            //pMap.pushZero(rcc.getMapeosC());
                            cout << "\n debug (5)---------------------------------------\n";
                            //pOrderMap.pushZero(rcc.getMapeosC());
                            cout << "\n debug (6)---------------------------------------\n";
                            //cout << "\n-->RCC Result J\n";
                            //rcc.getJ().escSubComplexJ();
                            setPAuxilary(cycle);
                            cout << "\n debug (7)---------------------------------------\n";
                            //cout << "\n\n J Facets are : \n";
                            //rcc.getJ().escSubComplexJ();
                            //cout << "\n\n";
                            cout << "\n debug (8)---------------------------------------\n";
                            eliminateP_From_A();
                            cout << "\n debug (9)---------------------------------------\n";
                            //cout << "\n-->new A list := ";
                            //a.escListInt(Acomp);
                            cout << "\n debug (10)---------------------------------------\n";
                            P.resetPCovering();
                            cout << "\n debug (11)---------------------------------------\n";
                            if (Acomp.size() > 1)
                                rcc.resetSubcomplex_CORE(Acomp, AA, cycle);
                            cout << "\n debug (12)---------------------------------------\n";
                            
                } else {
                            //for (int iterJ = 0; iterJ < rcc.getJSize(); iterJ++)
                            cout << "\n debug (13)---------------------------------------\n";
                            P.pushZero(AA.getA(0, getAcomp(Acomp, 0)));
                            cout << "\n debug (14)---------------------------------------\n";
                            p.push(P);
                            cout << "\n debug (15)---------------------------------------\n";
                            p.n += 1;

                            cout << "\n debug (16)---------------------------------------\n";
                            pOrder.push(P);
                            cout << "\n debug (17)---------------------------------------\n";
                            pOrder.n += 1;

                            //break;
                            cout << "\n debug (18)---------------------------------------\n";
                            eliminateP_From_A();
                            cout << "\n debug (19)---------------------------------------\n";
                            //cout << "\n-->new A list := ";
                            //a.escListInt(Acomp);
                            P.resetPCovering();
                            cout << "\n debug (20)---------------------------------------\n";
                }

                //Acomp.pop_back();
                //cout << "\n\n";
                //I.escSubComplexJ();
                //I.escZero_Skeleton();

                cycle += 1;
                //cout << "\n---------------------------------------------------------------------------------------------------------->\n";

            }


            //I.escSubComplexJ();

            cout << "\n\n<<<<<<<<<<<End of Covering>>>>>>>>>>>>>>>\n\n";
            
            p.escTensorSimplexAlpha();


            int si = 0;
            for (int i = 0; i < p.n; i++) {
                si += p.sizeOfPartition(i);
            }

            if (si != L.numSimplex) {
                throw runtime_error("WROOOOOOOOOOOOONG");
		cout << "\n==================>> " << si << "\n\n";
	    }

	    return;
        }

	void endCovering() {
		cout << "\n\n===================>>Covering Finalized\n\n";
		rcc.endRCC_CORE();
	}


    void orderPartition(SimplicialMap a) {

//        list <int> order;
//        for (int i = 0; i < p.n; i++) 
//                order.push_back(p.sizeOfPartition(i));
//
//        order.sort();
//        cout << "\n-->Prepararing order algorithm\n";
//        a.escListInt(order);
//
//        int count = 0;
//        while (order.size() != 0) {
//            int lastS = order.size() - 1;
//            for (int i = 0; i < p.n; i++) {
//                if (p.sizeOfPartition(i) == getAcomp(order, lastS)) {
//                    //cout << "\n-->Found Partition in P of size := <<" << p.sizeOfPartition(i) << ">>  comparble with last <<" << lastS << ">>";
//                    if (count == 0)
//                        pOrder.initA(0, 0, p.A[0][i]);
//                    else {
//                        pOrder.push(p.A[0][i]);
//                        pOrder.n += 1;
//
//                    }        
//                    order.pop_back();
//                    count += 1;
//                    i = p.n+1;
//                }
//            }
//        }
//
//        order.clear();

        int count = 0;
        for (int i = 0; i < p.n; i++) {
                if (count == 0)
                    pOrder.initA(0, 0, p.A[0][i]);
                else {
                    pOrder.push(p.A[0][i]);
                    pOrder.n += 1;
                }
                count += 1;
            }

        for (int i = 0; i < p.n; i++) {
            pOrder.replaceAt(i, p.A[0][i]);
        }

        cout << "\n\n-->ordered tensor should be:\n";
        pOrder.escTensorSimplexAlpha();
    }

    int getPM() {
        return pOrder.getM();
    }

    int getPN() {
        return pOrder.getN();
    }



    int getPartitionSize(int i) {
        return  pOrder.sizeOfPartition(i);
    }



    // { [m0, n0], [m1, n1],  [m2, n2],  [m3, n3]}
    //[0, ni] = sigma(0, ni)

    SimplexAlpha getPartitionSimplex(int i, int j) {

    //    if (pOrder.n > i && pOrder.A[0][j].n > j) {
    //        //cout << "\n<<NumTensor Size numT := " << pOrder.n  << " <>> and i = " << i << "\n";
    //        //cout << "\n<<NumTensorIndexMat Size numTIM := " << pOrder.A[0][i].n  << " <>> and j = " << j << "\n\n";
    //        pOrder.A[0][i].escMatrixSimplexAlpha();
    //    }

        //cout << "\npetition(i, j) := "; this->pOrder.getPartitionSimplex(i, j).escSimplexAlpha();
        //cout << "\n";
        return pOrder.getPartitionSimplex(i, j);
    }






    void replaceAt_i(int i, MatrixSimplexAlpha M) {
        pOrder.replaceAt(i, M);
    }
	
    void escPOrder() {
        pOrder.escTensorSimplexAlpha();
    }
};

#endif
