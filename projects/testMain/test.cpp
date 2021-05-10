#include <vector>
#include <iostream>
#include <fstream>
#include <string>

using namespace std;

int main(int argc, char* argv[])
{
	if(argc != 3)
	{
		cout << "Wrong number of input arguments. Correct is 2." << endl;
		exit(-1);
	}

	double a,b;
	a = stod(argv[1]);
	b = stod(argv[2]);

	ofstream wFile;
	wFile.open("result.txt");
	wFile << a*a + b*b;
	wFile.close();

	return 0;
}