#include <iostream>
#include <cstdlib>
#include <math.h>
#include <list>
#include "Matrix.hpp"

using namespace std;

void Matrix::zeroMatrix(int M) {
	
	for (int i = 0; i < M; i++)
		for (int j = 0; j < M; j++) 
			updateA(i, j, 0.0);		
}

void Matrix::escMatrix() {

	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++)
			std::cout << A[i][j] << " ";

		std::cout << std::endl;
	}
}

double Dijkstra::initDijkstra(Matrix graph, int s1, int f) {

	s = f;
	f = s1;

	if (s == f) {
		//std::cout << "\n<<No Action in Dijkstra: s == f>>\n";
		return 0.0;
	}
	
	else {

		//std::cout << "\n<<< Dijkstra >>>\n" <<std::endl;
		//std::cout << "\n\n-->Initializing dynamic auxiliary arrays" <<std::endl;
		
		for (int i = 0; i < graph.getM(); i++) {
			pushIntoQ(i);		//Q must have all the verticies in the graph
			pushIntoPrev(-1); 	//prev must be undefined
			//dist must be set at infinity except at vertex s
			if (i == s) pushIntoDist(0.0);
			else pushIntoDist(10000000000);
		}

		//escListDist();
		//escListPrev();
		//escListQ();

		//int u, v;
		//double alt;

		//std::cout << "\nBegin Main loop\n\n";
		while (getQSize() != 0) { 	// The main loop
			//escListQ();
			u = eraseMinQ();		//delete minimum from Q and asign it to u
				
				//if (u == f) break;
				//else {
								//std::cout << "u = " << u << " minimum element in Q\n";
								for (int i = 0; i < getQSize(); i++) { // only v that are still in Q
									
									v = getQ(i); 
									if (compareRealNumbers(graph[u][v], 0.0) == 0) {
									
											//std::cout << "   v = "<< v << " in Q\n";
											alt = getDist(u) + graph[u][v];
											//std::cout << "   -->alt = "<< getDist(u) << " + " << graph[u][v] << " = " << alt << "\n";
											if (alt < getDist(v)) {
												//std::cout << "         " << alt << " < " << dist[v] << "\n\n";
												updateDist(v, alt);
												updatePrev(v, u);
												//pushIntoPrevDijkstra(u);
											}
									}
					
								}
								//std::cout << "\n\n";
		
				//}
		}

		w = prev[f];
		S.push_back(w);

		while (w != -1) {
			w = prev[w];
			if (w != -1) {
      			S.push_back(w);
    		}
		}

		//std::cout << f;
		//escListS();
		ss = S.size();
		S.clear();
		prev.clear();
		dist.clear();
		return ss;
	}


}
