#include <stdarg.h>

#ifndef OPTIMIZEDCOVERING
#define OPTIMIZEDCOVERING

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <list>
#include <iterator>
#include "SimplexAbstract.hpp"
#include "SimplexAlpha.hpp"
//#include "RCC.hpp"
#include <unistd.h>
#include "Covering_CORE.hpp"

using namespace std;

class OptimizedCovering {
    private:
        //RCC rcc, addFacets;
        Covering c;
        MatrixSimplexAlpha D, AUX;
        SubComplexJ JD, JQ, JP, JQs, joker;
        LocalSearch suchen, suchen0;
        
	public:

        void initOptimizedCovering_CORE(SubComplexJ L, Complex K, SimplicialMap a, SimplicialMap b, int M, double r,  int argc, char **argv) {

            cout << "\n\n<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>\n";
            cout << "<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>\n";
            cout << "<<<<<<<OPTIMIZEDCOVERING>>>>>>>\n";
            cout << "<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>\n";
            cout << "<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>\n";

            cout << "M := " << M << "\n";
            cout << "r := " << r << "\n";

            //Memory ALLOCATION SETUP

            c.initCovering(L, K, a, b, M, r, argc, argv);
            suchen.initLocalSearch(L, K, a, b, M, r);
            suchen0.initLocalSearch(L, K, a, b, M, r);
            
            D.initMatrixSimplexAlpha(1, 1);
            
            AUX.initMatrixSimplexAlpha(1, 1);
            AUX.initAByCopyOf(0, 0, L.listOfFacets.A[0][0]);
            AUX.n = 0;

            JD.initSubComplexJ(1);
            JD.initA(0, L.listOfFacets.A[0][0]);
            JD.initZero_Skeleton();
            JD.listOfFacets.n = 0;

            JQ.initSubComplexJ(1);
            JQ.initA(0, L.listOfFacets.A[0][0]);
            JQ.initZero_Skeleton();
            JQ.listOfFacets.n = 0;
            
            JP.initSubComplexJ(1);
            JP.initA(0, L.listOfFacets.A[0][0]);
            JP.initZero_Skeleton();
            JP.listOfFacets.n = 0;

            JQs.initSubComplexJ(1);
            JQs.initA(0, L.listOfFacets.A[0][0]);
            JQs.initZero_Skeleton();
            JQs.listOfFacets.n = 0;

        }

        
                                            //////////////////////////
                                            //////////////////////////
                                            // runOptimizedCovering //
                                            //////////////////////////
                                            //////////////////////////


        void runOptimizedCovering_CORE(SubComplexJ L, Complex K, SimplicialMap a, SimplicialMap b, int M, int Ns, double r, int argc, char **argv) {

            cout << "\n\n               <<RUNNING COVERING FOR OPTIMIZED>>\n";
            c.runCovering(L, K, a, b, M, r, argc, argv);////////////////////////////MPI STILL OPENED

            //cout << "                 <<ORDER PARTITION>>\n";
            //c.orderPartition(a);
            cout << "                   <<PARTITION P'>>\n";
            	//							c.p.escTensorSimplexAlpha();

            cout << "                   <<REORDERING P>>\n";
                                        				orderP();


            	int j = 0;
            				int N = 3;
									int i = 0;


            while (i < Ns) {
				i += 1;
                            //cout << "                   WHILE   <<MAIN LOOP>>\n\n";
                            //cout << "                   <<      "<< " " << "<>>--------------------------------------------------------------------\n";
                            //cout << "                   <<CYCLE "<< i << "<>>--------------------------------------------------------------------\n";            
                            //cout << "                   <<WORKING WITH PARTITION P >>\n"; 
                                        //c.pOrder.escTensorSimplexAlpha();
                            //cout << "                   <<REORDERED P>>\n";
                                        //c.pOrder.escTensorSimplexAlpha();

                            
                            //cout << "                   <<BIG UNION OF PARTITION INTO JD from j >>  j := " << j << "\n";
                                        unionJD(j);
                            
                            //cout << "                   <<JD AFTER UNION k := j to n >>\n";
                                        //JD.escSubComplexJ(); 
                            
                            //cout << "                   <<RESET THE RCC SubComplexJ>>\n";
                                        c.resetRCCSubComplexJ(JD);
                            //cout << "                   <<Running RCC WITH SUBCOMPLEX JD>>\n";
                                        //rcc.runRCC_CORE(JD, K, a, b, M, r, argc, argv);
                           		//c.simplexScrutiny(JD, K, a, b, M, r, argc, argv);
				        c.runRCC(JD, K, a, b, M, r, argc, argv);	
                            //cout << "                   <<RCCJD := RCC_RESULTS(SUBCOMPLEX JD) := >>\n";
                                        //rcc.getJ().escSubComplexJ();
                            //            cout << "\n\n";
                                        
                            //cout << "                   <<nBuilding JP subcomplex >>\n";
                                        buildP();
                            
                            //cout << "                   Is card(P" << j << ") > card(JP)??\n";
                                        largest(j, L);
                            
                            //cout << "\nNOW WORKING ON Q--------------------------------------------------------------------------------------------------------------\n\n                   <<WORKING WITH list P of partitions>>\n"; 
                                        //c.pOrder.escTensorSimplexAlpha();
                            
                            //cout << "                   <<BIG UNION OF PARTITION INTO JD from j-1 >>  j-1 := " << j-1 << "\n";
                             //          cout << "\n                   <<JQ>>\n";
                                        
					
					int R0 = rand()%(L.listOfFacets.n);
					JD.resetSubComplexJ(L.listOfFacets.A[0][R0]);
                                        JD.listOfFacets.n = 0;
                                        //JD.escSubComplexJ();

					//c.AUX_CORE(L, &JD);

                                        unionJD2(j-1);

                            
                            //cout << "                   <<JD AFTER UNION k:=j-1 to n>>\n";
                                        //JD.escSubComplexJ();
                            
                            //cout << "                   <<RESET THE ADDFACET SubComplexJ J from JP>>\n";
                                        c.resetAddFacetsSubComplexJ(JP);
                            		
                            //cout << "                   <<Running ADDFACETS WITH SUBCOMPLEX JD covering JP>>\n";
                                       	// addFacets.runRCC_CORE(JD, K, a, b, M, r, argc, argv);
                            		//c.simplexScrutinyGeneral(JD, K, a, b, M, r, argc, argv);
					c.runAddFacets(JD, K, a, b, M, r, argc, argv);
                            //cout << "                   <<ADDFACET RESULTS := JP U {s0, s1, ..., sm}>>\n";
                                        //addFacets.getJ().escSubComplexJ();

                            //cout << "                   <<Building SubComplex JQs>>\n";
                                        buildQ();
                            
                            if (j > 0) {
                              //  cout << "                   <<(j > 0) UPDATE JQs WITH P" << j-1 <<" if card(P" << j-1 << ") > card(JQs) >>\n";
                                        largestJQ(j-1, L);
                            }
                            
                                    if (j == 0) cout << "                   <<(j == 0) NO UPDATE to JQs>>\n";
                            


                            deleteVoids(j-1, JQs);

                            if (c.p.n < c.pOrder.n) {
                                //cout << "                   <<card(P') < card(P)>>\n";
                                //cout << "                   << P <- P' >>\n";
                                //we replace P with P'
                                c.pOrder.updateN(c.p.n);
                                for (int l = 0; l < c.p.n; l++) {
                                    c.pOrder.A[0][l].replaceWith(c.p.A[0][l]);
                                }

                                orderP();
                                //cout << "                   <<P is now := >>\n";///////////////////////////////////////////////////////
                                //c.pOrder.escTensorSimplexAlpha();

                            }

                            
                            
                                                        ///////////////////////////////
                                                        ///////////////////////////////
                                                        /////////               ///////
                                                        /////////     REFRESH   ///////
                                                        /////////    FOR NEXT   ///////
                                                        /////////    ITERATION  ///////
                                                        ///////////////////////////////
                                                        ///////////////////////////////
                            
                            
                            cout << "\n\n\n                   <<RESETING JD AND JQ FOR NEXT ITERATION>>\n";
                            cout << "                   <<JD>>\n";
                                        
                                        R0 = rand()%(L.listOfFacets.n);
                                        JD.resetSubComplexJ(L.listOfFacets.A[0][R0]);
                                        JD.listOfFacets.n = 0;
                                        //JD.escSubComplexJ();

                                        
                            cout << "\n                   <<JQ>>\n";
                                        R0 = rand()%(L.listOfFacets.n);
                                        JQ.resetSubComplexJ(L.listOfFacets.A[0][R0]);
                                        JQ.listOfFacets.n = 0;
                                        //JQ.escSubComplexJ();

                            cout << "\n                   <<JP>>\n";
                                        R0 = rand()%(L.listOfFacets.n);
                                        JP.resetSubComplexJ(L.listOfFacets.A[0][R0]);
                                        JP.listOfFacets.n = 0;
                                        //JP.escSubComplexJ();

                                        //matMapJP.resetMatrixSimplicialMapZero(destroyAux);

                            cout << "\n                   <<JQs>>\n";
                                        R0 = rand()%(L.listOfFacets.n);
                                        JQs.resetSubComplexJ(L.listOfFacets.A[0][R0]);
                                        JQs.listOfFacets.n = 0;
                                        //JQs.escSubComplexJ();
                                        //matMapJQs.resetMatrixSimplicialMapZero(destroyAux);
                            
                            
                            j += 1;
                            if (j > c.pOrder.n - 1) j = 0;


                            cout << "                   <<END OF LOOP>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><\n";

            }//End of Loop
		
	    //c.pOrder.escTensorSimplexAlpha();
	    c.WriteTxt(c.pOrder.stringTensorSimplexAlpha());
	
		/////////////////
            //Preparing Buffer SubComplex
            ////////////////
            

            joker.initSubComplexJ(1);
            joker.initA(0, c.pOrder.getPartitionSimplex(0, 0));
            joker.initZero_Skeleton();

            /////////////////
            //Copying 1st subcomplex of pOrder
            ////////////////
	    string ret = "";

            	//for (int I = 0; I < 3; I++){
                //        for (int iii = 1; iii < c.pOrder.sizeOfPartition(I); iii++)
                //            joker.pushSimplexAlpha(c.pOrder.getPartitionSimplex(I, iii));
            
                //        ret += "\n\n------------------------------------<<<For Partition := " + to_string(I) + "s <<RESULTS>>---------------------------------------------------\n\n";
                //        //joker.escSubComplexJ();
            
                //        /////////////////
                //        //Running LocalSearch
                //        ////////////////
                //        if (joker.listOfFacets.n > 1) {
                //            while(suchen.psy_reduced.m < 2) 
                //                ret += suchen.updateLocalSearch(joker, K, a, b, M, 0.1);
                //        }
            
                //        if (joker.listOfFacets.n == 1) ret += "\nSize of complex == 1, trivial case\n\n";
            
                //        if (I < c.pOrder.m - 1)
                //            joker.resetSubComplexJ(c.pOrder.getPartitionSimplex(I+1, 0));



		//}






		c.WriteTxt(ret);
	    c.endCovering();
        }//End of runOptimizedCovering

        void closeWriterO(){
        c.closeWriter();
        }

        //used in optimizedCovering j is considered to come as j-1
        //updateSimplexAlpha(SimplexAlpha simplex)
        void deleteVoids(int j, SubComplexJ L) {

            //cout << "\n                   <<deleteVoids>>\n";
            //cout << "                   <<P list:>>\n";
            //c.escPOrder();

            //cout <<   "            Q is defined as:\n";
            //JQs.escSubComplexJ();

             if (j < 0) j = 0;
            if (JQs.listOfFacets.n == 0) cout << "\nWARNING NO INFO ON J GIVEN\n";

            for (int k = j; k < c.pOrder.n; k++) {
             //   cout << "       -->list(" << k << ")       \n";

                            for (int i = 0; i < c.getPartitionSize(k); i++) {
                           //     cout << " " << i;

                                int bit = 0;
                                for (int l = 0; l < JQs.listOfFacets.n; l++) {
                                    if (c.getPartitionSimplex(k, i).compareSimplexAlpha(JQs.listOfFacets.A[0][l]) == 1)
                                        bit += 1;
                                }

                                cout << "\n";
                                if (bit == 0) {
                                    cout << "--> " << i << "  is not in JQs\n";
                                    AUX.pushZero(c.getPartitionSimplex(k, i));
                                }

                                if (bit > 0) {
                                    //cout << "--> " << i << "  is in JQs\n";
                                }
                            }
                            //cout << "\n";

                            //cout << "           AUX Matrix as:\n";
                            //AUX.escMatrixSimplexAlpha();

                            //cout << "           We empty Pk: k:= j, ..., n\n";
                            c.pOrder.A[0][k].updateMatrixSimplexAlphaSize(1, 1);
                            c.pOrder.A[0][k].A[0][0].updateSimplexAlpha(AUX.A[0][0]);
                            c.pOrder.A[0][k].n = 0;

                            //cout << "           Refill Pk with no JQs simplicies\n";
                            for (int l = 0; l < AUX.getN(); l++)
                                c.pOrder.A[0][k].pushZero(AUX.A[0][l]);



                            //cout << "          P list is now :\n";
                            //c.pOrder.escTensorSimplexAlpha();

                            AUX.updateMatrixSimplexAlphaSize(1, 1);
                            AUX.A[0][0].updateSimplexAlpha(L.listOfFacets.A[0][0]);
                            AUX.n = 0;

                            //cout << "          AUX list is now :\n";
                           // AUX.escMatrixSimplexAlpha();
            }

                            //cout << "          PUSHING Q INTO P :\n";

                            c.pOrder.push(JQs.listOfFacets);
                            c.pOrder.n += 1;

                            //cout << "          P list is now :\n";
                            //c.pOrder.escTensorSimplexAlpha();

                           // cout << "                   <<REORDERING P>>\n";
                                        orderP();

                            //cout << "          P list is now :\n";
                            //c.pOrder.escTensorSimplexAlpha();

                            int voidCount = 0;
                            for (int k0 = 0; k0 < c.pOrder.n; k0++)
                                if (c.getPartitionSize(k0) == 0)
                                    voidCount += 1;
                            //cout << "\n          NUMBER OF VOIDS IN P := " << voidCount << "\n";
                            
                            if (voidCount == 0) cout << "       NO VOIDS in P\n\n";
                            if (voidCount > 0) {
                                //cout << "       DELETING VOIDS in P\n\n";
                                int newN = c.pOrder.n - voidCount;
                                c.pOrder.updateN(newN);
                                //cout << "          P list is now :\n";
                                //c.pOrder.escTensorSimplexAlpha();
                            }
        }

        
        void buildP() {

            //cout << "\nBuilding JP from RCC \n";
            for (int j = 0; j < c.J.listOfFacets.n; j++) 
                    JP.pushSimplexAlphaZero(c.J.listOfFacets.A[0][j]);
            
            //cout << "Our JP subcomplex given by the RCC is:\n";
            //JP.escSubComplexJ();

        }

        void buildQ() {

            //cout << "\nBuilding JQs from ADDFACETS \n";
            for (int j = 0; j < c.Jprime.listOfFacets.n; j++) 
                    JQs.pushSimplexAlphaZero(c.Jprime.listOfFacets.A[0][j]);

            //cout << "Our JQs subcomplex given by the ADDFACETS is:\n";
            //JQs.escSubComplexJ();  
        }      

        void unionJD(int J) {

            //cout << "\nBig Union         numPartitions := " << c.pOrder.n << "\n";
            for (int i = J; i < c.pOrder.n; i++) {
                if (c.getPartitionSize(i) != 0)
                for (int j = 0; j < c.getPartitionSize(i); j++) {
                    JD.pushSimplexAlphaZero(c.getPartitionSimplex(i, j));
                }
            }
        }

        void unionJD2(int J) {

            if (J < 0) J = 0;
            //cout << "\nBig Union of pOrder to JD        numPartitions := " << c.pOrder.n << "\n";
            for (int i = J; i < c.pOrder.n; i++) {
                if (c.getPartitionSize(i) != 0)
                for (int j = 0; j < c.getPartitionSize(i); j++) {
                    JD.pushSimplexAlphaZero(c.getPartitionSimplex(i, j));
                }
            }
        }

        void unionJQ(int J) {

            //cout << "\nBig Union         numPartitions := " << c.pOrder.n << "\n";
            for (int i = J; i < c.pOrder.n; i++) {
                if (c.getPartitionSize(i) != 0)
                for (int j = 0; j < c.getPartitionSize(i); j++) {
                    JQ.pushSimplexAlphaZero(c.getPartitionSimplex(i, j));//c.getPartitionSimplex(i, j));
                }                
            }   
        }

        int findMax(int J) {

            int max = 0;
            int index = 0;

            for (int i = J; i < c.pOrder.n; i++)
                if (max <= c.getPartitionSize(i)) {
                    max = c.getPartitionSize(i);
                    index = i;
                }

            return index;
        }

        void orderP() {
            
            for (int i = 0; i < c.pOrder.n; i++) {
                int maxI = findMax(i);

                if (maxI != i) {
                        //cout << "-->The Max Partition from " << i << " to " << c.pOrder.n - 1 << " is := " << maxI << "\n";
                        for (int l = 0; l < c.getPartitionSize(maxI); l++)
                            AUX.pushZero(c.getPartitionSimplex(maxI, l));
                        
                        c.pOrder.A[0][maxI].updateMatrixSimplexAlphaSize(1, 1);
                        c.pOrder.A[0][maxI].A[0][0].updateSimplexAlpha(AUX.A[0][0]);
                        c.pOrder.A[0][maxI].n = 0;

                        for (int l = 0; l < c.pOrder.A[0][i].getN(); l++)
                            c.pOrder.A[0][maxI].pushZero(c.pOrder.A[0][i].A[0][l]);

                        c.pOrder.A[0][i].updateMatrixSimplexAlphaSize(1, 1);
                        c.pOrder.A[0][i].A[0][0].updateSimplexAlpha(AUX.A[0][0]);
                        c.pOrder.A[0][i].n = 0;

                        for (int l = 0; l < AUX.getN(); l++)
                            c.pOrder.A[0][i].pushZero(AUX.A[0][l]);

                        AUX.updateMatrixSimplexAlphaSize(1, 1);
                        AUX.A[0][0].updateSimplexAlpha(c.getPartitionSimplex(i, 0));
                        AUX.n = 0;
                        //AUX.escMatrixSimplexAlpha();
                }
            }
        }


        //We replace JP with Pj if card(Pj)>card(JP)

        void largest(int J, SubComplexJ L) {

            //cout << "\n--> card(P" << J <<") := " << c.getPartitionSize(J) << "        card(JP) := " << JP.listOfFacets.n << "\n\n"; 

            if (c.getPartitionSize(J) > JP.listOfFacets.n) {
                //cout << "\n--> card(JP) < card(P" << J <<")     UPDATE JP from P" << J << "\n";
                //cout << "\n                   <<deleting JP>>\n";

                int R0 = rand()%(L.listOfFacets.n);
                JP.resetSubComplexJ(L.listOfFacets.A[0][R0]);
                JP.listOfFacets.n = 0;
                //JP.escSubComplexJ();
                
		for (int i = 0; i < c.getPartitionSize(J); i++)
                    JP.pushSimplexAlphaZero(c.pOrder.A[0][J].A[0][i]);

                //cout << "\n\n-->New JP := \n";
                //JP.escSubComplexJ();

            } else {
                //cout << "\n--> card(P" << J <<") <= card(JP)   NO UPDATE to JP\n\n";
            }
            
        }

        void largestJQ(int J, SubComplexJ L) {
            
            //cout << "\n--> card(P" << J <<") := " << c.getPartitionSize(J) << "        card(JQs) := " << JP.listOfFacets.n << "\n\n"; 

            if (c.getPartitionSize(J) > JQs.listOfFacets.n) {
                //cout << "\n--> card(JQs) < card(P" << J <<")     UPDATE JP from P" << J << "\n";
                //cout << "\n                   <<deleting JP>>\n";

                int R0 = rand()%(L.listOfFacets.n);
                JQs.resetSubComplexJ(L.listOfFacets.A[0][R0]);
                JQs.listOfFacets.n = 0;
                //JQs.escSubComplexJ();
                
		for (int i = 0; i < c.getPartitionSize(J); i++)
                    JQs.pushSimplexAlphaZero(c.pOrder.A[0][J].A[0][i]);

                //cout << "\n\n-->New JQs := \n";
                //JQs.escSubComplexJ();

            } else {
                //cout << "\n--> card(P" << J <<") <= card(JQs)   NO UPDATE to JQs\n\n";
            }
            
        }


};

#endif
