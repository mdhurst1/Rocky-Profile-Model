#include <iostream>
#include <fstream>

using namespace std;

ofstream outstream;
outstream.open("myfile.txt");
outstream << "write this text" << endl;
outstream.close();
