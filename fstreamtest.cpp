#include <iostream>
#include <fstream>

using namespace std;

ofstream outstream("myfile.txt");
outstream << "write this text" << endl;
outstream.close();
