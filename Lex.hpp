#ifndef LEX_HPP
#define LEX_HPP

#include <iostream>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <complex>
#include <sstream>
#include <string>
#include <fstream>

using namespace std;

//Escribimos una clase que escribe el texto arrojado por Expr.cpp el cual esta en LaTeX

class SimplexAlphaWriter {

      private:
        string filePath;

      public:
        ofstream writer;

        SimplexAlphaWriter(){};

        SimplexAlphaWriter(string path) {

            writer.open(path.c_str());
  
            if (writer.is_open()) {
                
                filePath = path;
                writer << "\n\n                       <<Optimized Covering Layout>>\n\n\n" << endl;
            }

        }

        void write(string inputTexto) {

          writer << inputTexto;
        }

        void closeSimplexAlphaWriter() {
          writer << "\n\n                       <<End Of Document>>\n\n" << endl;
          writer.close();
        }

};
#endif
