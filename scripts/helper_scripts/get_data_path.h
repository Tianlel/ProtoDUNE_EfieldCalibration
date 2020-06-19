#ifndef _GET_DATA_PATH
#define _GET_DATA_PATH

/* Given a txt file and a line number, return the file path (string literal) on that line */

#include <iostream>
#include <fstream>
#include <string>
#include <limits>
#include <stdlib.h>

using namespace std;

string get_data_path(string filename, unsigned int line_number)
{
    ifstream file (filename.c_str(), std::ifstream::in);
    if (file.is_open()) { 
        string outstr;
        file.seekg(ios::beg);
        for (int i=0; i < line_number -1; ++i){
            file.ignore(numeric_limits<streamsize>::max(), '\n');
        }
        getline(file, outstr);
        return outstr;
    }
    else {
        cout << "Error opening file " << filename << endl;
        exit(EXIT_FAILURE);
    } 
}

#endif /* _GET_DATA_PATH */



