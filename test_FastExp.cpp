#include <cmath>
#include <cstdlib>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <unistd.h>
#include "./FastExp.hpp"

using namespace std;

int main()
{
    vector<double> Exponents{-0.000001,-0.00001,-0.0001,-0.001,-0.01,-0.1,-1.,-10.,-100.,-1000.};
    vector<double> Results(Exponents.size());
    vector<double> FastResults(Exponents.size());

    for (int i=0; i<Exponents.size(); ++i)
    {
        Results[i] = exp(Exponents[i]);
        FastResults[i] = fastexp(Exponents[i]);
        
    }

    for (int i=0; i<Exponents.size(); ++i)
    {
        cout << Exponents[i] << ",";
    }
    cout << endl;
        for (int i=0; i<Exponents.size(); ++i)
    {
        cout << Results[i] << ",";
    }
    cout << endl;
        for (int i=0; i<Exponents.size(); ++i)
    {
        cout << FastResults[i] << ",";
    }
    cout << endl;
    
    //cout << Exponents << endl;
    //cout << Results << endl;
    //cout << FastResults << endl;
}
