#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstring>

using namespace std;

int main()
{
    ofstream WriteFile;
    string Filename = "myfile.txt";

    WriteFile.open(Filename.c_str());
    WriteFile << "write this text" << endl;
    WriteFile.close();
}