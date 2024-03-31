#include <string>
#include <fstream>
#include <iostream>

using namespace std;

// g++ -I/usr/local/include -o scriptGen ScriptGen.cpp

int main (int argc, char* argv[])
{
	ofstream	writer;
	string		outputFileName;
	double		kAA;
	double		kAB;
	double		kBB;
	double		kBA;
	
	outputFileName = string(argv[1]);
	kAA	= atof(argv[2]);
	kAB	= atof(argv[3]);
	kBB = atof(argv[4]);
	kBA = atof(argv[5]);
	
	cout	<< "Creating execution script for " << outputFileName 
			<< " with rA = " << kAA/kAB << " and rB = " << kBB/kBA << endl;
	
	writer.open(outputFileName.c_str(),ios::out);
	writer	<< "./rr " << outputFileName + "01" << " 0.01 0.01 0.99 " 
			<< kAA << " " << kAB << " " << kBB << " " << kBA 
			<< " " << "4000 0.000000000001" << endl;
	writer	<< "./rr " << outputFileName + "03" << " 0.01 0.03 0.97 " 
			<< kAA << " " << kAB << " " << kBB << " " << kBA 
			<< " " << "4000 0.000000000001" << endl;
	writer	<< "./rr " << outputFileName + "05"<< " 0.01 0.05 0.95 " 
			<< kAA << " " << kAB << " " << kBB << " " << kBA 
			<< " " << "4000 0.000000000001" << endl;
	writer	<< "./rr " << outputFileName + "10" << " 0.01 0.10 0.90 " 
			<< kAA << " " << kAB << " " << kBB << " " << kBA 
			<< " " << "4000 0.000000000001" << endl;
	writer	<< "./rr " << outputFileName + "15" << " 0.01 0.15 0.85 " 
			<< kAA << " " << kAB << " " << kBB << " " << kBA 
			<< " " << "4000 0.000000000001" << endl;
	writer	<< "./rr " << outputFileName + "20" << " 0.01 0.20 0.80 " 
			<< kAA << " " << kAB << " " << kBB << " " << kBA 
			<< " " << "4000 0.000000000001" << endl;
	writer	<< "./rr " << outputFileName + "30" << " 0.01 0.30 0.70 " 
			<< kAA << " " << kAB << " " << kBB << " " << kBA 
			<< " " << "4000 0.000000000001" << endl;
	writer	<< "./rr " << outputFileName + "40" << " 0.01 0.40 0.60 " 
			<< kAA << " " << kAB << " " << kBB << " " << kBA 
			<< " " << "4000 0.000000000001" << endl;
	writer	<< "./rr " << outputFileName + "50" << " 0.01 0.50 0.50 " 
			<< kAA << " " << kAB << " " << kBB << " " << kBA 
			<< " " << "4000 0.000000000001" << endl;
	writer	<< "./rr " << outputFileName + "60" << " 0.01 0.60 0.40 " 
			<< kAA << " " << kAB << " " << kBB << " " << kBA 
			<< " " << "4000 0.000000000001" << endl;
	writer	<< "./rr " << outputFileName + "70" << " 0.01 0.70 0.30 " 
			<< kAA << " " << kAB << " " << kBB << " " << kBA 
			<< " " << "4000 0.000000000001" << endl;
	writer	<< "./rr " << outputFileName + "80" << " 0.01 0.80 0.20 " 
			<< kAA << " " << kAB << " " << kBB << " " << kBA 
			<< " " << "4000 0.000000000001" << endl;
	writer	<< "./rr " << outputFileName + "90" << " 0.01 0.90 0.10 " 
			<< kAA << " " << kAB << " " << kBB << " " << kBA 
			<< " " << "4000 0.000000000001" << endl;
	writer	<< "./rr " << outputFileName + "95" << " 0.01 0.95 0.05 " 
			<< kAA << " " << kAB << " " << kBB << " " << kBA 
			<< " " << "4000 0.000000000001" << endl;
	writer	<< "./rr " << outputFileName + "97" << " 0.01 0.97 0.03 " 
			<< kAA << " " << kAB << " " << kBB << " " << kBA 
			<< " " << "4000 0.000000000001" << endl;
	writer	<< "./rr " << outputFileName + "99" << " 0.01 0.99 0.01 " 
			<< kAA << " " << kAB << " " << kBB << " " << kBA 
			<< " " << "4000 0.000000000001" << endl;
	writer.close();
	return 0;
}