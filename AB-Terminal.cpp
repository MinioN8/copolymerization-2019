/*****************************************************************************************
	
	Written by Nate Lynd
	This program simulates a copolymerization
	g++ -I/usr/local/include/eigen3 -O3 -msse3 -o rr AB-Terminal.cpp

*****************************************************************************************/
#include <cmath>
#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include <iostream>

#include <Eigen/Core>
#include "odeint.h"
#include "stepper.h"

using namespace std;
using namespace Eigen;

/*****************************************************************************************
	g++ -I /Users/nl7597/Documents/Calculations/Eigen/ -O3 -msse3 -o fit fits.cpp

	I/Users/nl7597/Documents/Publications/2020/2020-Imbrogno-Mechanism/

	The following generates the RHS derivatives (f(t,y)) for problems of the form:
	y' = f(x,y), where we seek y(x) as a solution.
	This is then passed to machinery that can integrate it.
*****************************************************************************************/
struct rhs_ter {
	/*
		The constructor of f(x,y):
	
		y_in() = {I/0, A/1, B/2, XA/3, XB/4}
		dydx_out() = {dI/dt, dA/dt, dB/dt, dXA/dt, dXB/dt}
	
		// parameters are double {kAA,kAB,kBA,kBB}
	*/
	double kA,kB,kAA,kBA,kAB,kBB;
	rhs_ter(double kA,double kB,double kAA,double kBA,double kAB,double kBB) 
			: kA(kA),kB(kB),kAA(kAA),kBA(kBA),kAB(kAB),kBB(kBB) {}
	void operator() (const double x, VectorXd &y_in, VectorXd &dydx_out)
	{
	// I - initiator concentration
	dydx_out(0) =	-kA*y_in(0)*y_in(1)-kB*y_in(0)*y_in(2);
	// A - monomer A concentration
	dydx_out(1) =	-kA*y_in(0)*y_in(1)-kAA*y_in(3)*y_in(1)-kBA*y_in(4)*y_in(1);
	// B - monomer B concentration
	dydx_out(2) =	-kB*y_in(0)*y_in(2)-kBA*y_in(3)*y_in(2)-kBB*y_in(4)*y_in(2);
	// XA - concentration of chains ending in A monomer
	dydx_out(3) =	kA*y_in(0)*y_in(1) + kBA*y_in(4)*y_in(1)-kAB*y_in(3)*y_in(2);
	// XB - concentration of chains ending in B monomer
	dydx_out(4) =	kB*y_in(0)*y_in(2) + kAB*y_in(3)*y_in(2)-kBA*y_in(4)*y_in(1);
	}
};

/*****************************************************************************************

*****************************************************************************************/
int main(int argc, char* argv[])
{
	string		kinetics_file;
	string		Skeist_file;
	string		FinemanRoss_file;
	string		KelenTudos_file;
	string		MeyerLowry_file;
	string		BeckinghamLyndA_file;
	string		BeckinghamLyndB_file;
	ofstream	kinetics;
	double		I,A,B;
	double		I0,A0,B0;
	double		kA,kB;
	double		kAA,kAB,kBA,kBB;
	double		termination_time;
	double		accuracy_goal;
	
	const int	iI=0,iA=1,iB=2,iXA=3,iXB=4;
	// to use integrator; number of variables:
	const	int		nvar = 5;

	// read input from command line:
	kinetics_file				=	string(argv[1]);
	I0							=	atof(argv[2]);
	A0							=	atof(argv[3]);
	B0							=	atof(argv[4]);
	kAA							=	atof(argv[5]);
	kAB							=	atof(argv[6]);
	kBB							=	atof(argv[7]);
	kBA							=	atof(argv[8]);
	termination_time			=	atof(argv[9]);
	accuracy_goal				=	atof(argv[10]);
	
	// for simplicity, we set kA = kAA and kB = kBB
	kA = 20.;
	kB = 20.;
	cout << "rA*rB = " << kAA/kAB*kBB/kBA << endl;
	cout << "kA = " << kA << endl;
	cout << "kB = " << kB << endl;
	cout << "kAA = " << kAA << endl;
	cout << "kAB = " << kAB << endl;
	cout << "rA = " << kAA/kAB << endl;
	cout << "kBB = " << kBB << endl;
	cout << "kBA = " << kBA << endl;
	cout << "rB = " << kBB/kBA << endl;
	cout << "I0 = " << I0 << endl;
	cout << "A0 = " << A0 << endl;
	cout << "B0 = " << B0 << endl;

	const	double	atol=accuracy_goal, rtol=atol, h1=0.002, hmin=0.0, x1=0.0, x2=termination_time;
	VectorXd		ystart = VectorXd::Zero(nvar);

	Skeist_file = string(kinetics_file) + ".Copoly.dat";
	FinemanRoss_file=string(kinetics_file) + ".FR.dat";
	KelenTudos_file=string(kinetics_file) + ".KT.dat";
	MeyerLowry_file = string(kinetics_file) + ".ML.dat";
	BeckinghamLyndA_file = string(kinetics_file) + ".BSL.A.dat";
	BeckinghamLyndB_file = string(kinetics_file) + ".BSL.B.dat";
	kinetics_file = kinetics_file + ".dat";

	// set up initial conditions
	ystart[0] = I0; // M(0)
	ystart[1] = A0; // A(0)
	ystart[2] = B0; // B(0)
	ystart[3] = 0.;	// XA(0)
	ystart[4] = 0.; // XB(0)

	Output out(-1);		// Dense output at 20 points plus x1.
	rhs_ter d(kA,kB,kAA,kBA,kAB,kBB);	// give rate constants to functor (above)
	Odeint<StepperDopr5<rhs_ter> > ode(ystart,x1,x2,atol,rtol,h1,hmin,out,d);
	ode.integrate();
	
	// file to store kinetic data:
	kinetics.open(kinetics_file.c_str(),ios::out);
	kinetics.precision(14);
	kinetics << " Full parameters: " << endl;
	kinetics << "rA*rB = " << kAA/kAB*kBB/kBA << endl;
	kinetics << "kA = " << kA << endl;
	kinetics << "kB = " << kB << endl;
	kinetics << "kAA = " << kAA << endl;
	kinetics << "kAB = " << kAB << endl;
	kinetics << "rA = " << kAA/kAB << endl;
	kinetics << "kBB = " << kBB << endl;
	kinetics << "kBA = " << kBA << endl;
	kinetics << "rB = " << kBB/kBA << endl;
	kinetics << "I0 = " << I0 << endl;
	kinetics << "A0 = " << A0 << endl;
	kinetics << "B0 = " << B0 << endl;
	kinetics << "time \t X(t) \t A(t) \t B(t) \t XA(t) \t XB(t) \t conv \t A/A0 \t conv \t B/B0" << endl;
	for (int i=0;i<out.count;i++)
		kinetics 	<< out.xsave(i) << "\t" 
					<< out.ysave(iI,i) << "\t" 
					<< out.ysave(iA,i) << "\t" 
					<< out.ysave(iB,i) << "\t"
					<< out.ysave(iXA,i) << "\t"
					<< out.ysave(iXB,i) << "\t"
					<< (A0-out.ysave(iA,i)+B0-out.ysave(iB,i))/(A0+B0) << "\t" // Conv
					<< (out.ysave(iA,i))/A0 << "\t"  // A/A0
					<< (A0-out.ysave(iA,i)+B0-out.ysave(iB,i))/(A0+B0) << "\t" // Conv
					<< (out.ysave(iB,i))/B0 << endl;  // B/B0
	kinetics.close();
	
	double x, y, F, G, f, f0, conv;
	// Skeist equation:  conversion, feed, polymer composition
	cout << "Writing Copolymer... " << endl;
	kinetics.open(Skeist_file.c_str(),ios::out);
	kinetics.precision(14);
	kinetics << fixed << "conv.\tf\tF" << endl;
	for (int i=0;i<out.count;i++)
	{
		conv = (A0-out.ysave(iA,i)+B0-out.ysave(iB,i))/(A0+B0);
		f = (out.ysave(iA,i)/(out.ysave(iA,i)+out.ysave(iB,i)));
		F = (A0-out.ysave(iA,i))/(A0-out.ysave(iA,i)+B0-out.ysave(iB,i));
		kinetics << fixed << conv << "\t" << f << "\t" << F << endl;
	}
	kinetics.close();
	
	// Fineman-Ross file: FA^2/fA vs FA/fA*(fA-1)
	cout << "Writing Fineman-Ross..." << endl;
	double ffr, Ffr;
	kinetics.open(FinemanRoss_file.c_str(),ios::out);
	kinetics.precision(14);
	kinetics << fixed << "conv\t f(t) \t FA^2/ffrA \t FA/fA(t)*(fA-1)" << endl;
	for (int i=0;i<out.count;i++)
	{
// 		conv = (A0-out.ysave(iA,i)+B0-out.ysave(iB,i))/(A0+B0);
// 		F = (A0-out.ysave(iA,i))/(B0-out.ysave(iB,i));
// 		f = (out.ysave(iA,i)/(out.ysave(iA,i)+out.ysave(iB,i))); 
// 		kinetics << fixed << conv << "\t" << f << "\t"  << F*F/(A0/(B0)) << "\t" << F/(A0/(B0))*((A0/(B0))-1) << endl;
		conv = (A0-out.ysave(iA,i)+B0-out.ysave(iB,i))/(A0+B0);
		F = (A0-out.ysave(iA,i))/(B0-out.ysave(iB,i));
		f = (out.ysave(iA,i)/(out.ysave(iA,i)+out.ysave(iB,i))); 
		ffr = (out.ysave(iA,i)/(out.ysave(iB,i)));
		kinetics << fixed << conv << "\t" << f << "\t"  << ffr*ffr/F << "\t" << ffr*(F-1)/F << endl;
	}
	kinetics.close();
	
	// Kelen-Tudos file: eta = (rA + rB/a) xi - rB/a
	// eta = G / (a + F)
	// xi = F / (a + F)
	// F = x*x/y
	// G = x*(y-1)/y
	// x = MA/MB // monomer concentrations
	// y = mA/mB // repeat unit concentrations
	// a = sqrt(Fmin*Fmax)
		// file.KL
	cout << "Writing Kelen Tudos..." << endl;
	kinetics.open(KelenTudos_file.c_str(),ios::out);
	kinetics.precision(14);
	// This can't be done until a conversion is decided upon. Therefore, for KT, this
	// will only spit out x, y, F, G
	kinetics << fixed << "conv\tx\ty\tF\tG" << endl;
	for (int i=0;i<out.count;i++)
	{
		conv = (A0-out.ysave(iA,i)+B0-out.ysave(iB,i))/(A0+B0);
		x = out.ysave(iA,i)/out.ysave(iB,i);
		y = (A0-out.ysave(iA,i))/(B0-out.ysave(iB,i));
		F = x*x/y;
		G = x*(y-1)/y;
		kinetics << fixed << conv << "\t" << x << "\t" << y << "\t" << F << "\t" << G << endl;
	}
	kinetics.close();
	
	// Meyer-Lowry: conversion vs. 1 - (A(t)/A0)^(a)* ((1-A(t))/(1-A0))^(b) * ((A0 - d)/(A(t)-d))^g
	// a = rB/(1-rB)
	// b = rA/(1-rA)
	// d = (1-rB)/(2-rA-rB)
	// g = (1-rA*rB)/((1-rA)*(1-rB))
		// file.ML
	cout << "Writing Meyer-Lowry..." << endl;
	kinetics.open(MeyerLowry_file.c_str(),ios::out);
	kinetics.precision(14);
	kinetics << fixed << "conv.\tf/f0" << endl;
	for (int i=0;i<out.count;i++)
	{
		conv = (A0-out.ysave(iA,i)+B0-out.ysave(iB,i))/(A0+B0);
		f = out.ysave(iA,i)/(out.ysave(iA,i)+out.ysave(iB,i));
		f0 = A0/(A0+B0);
		kinetics << fixed << f << "\t" << conv << endl;
	}
	kinetics.close();
	// Beckingham-Lynd:
		// file.BSLA
	cout << "Writing Beckingham-Sanoja-Lynd..." << endl;
	kinetics.open(BeckinghamLyndA_file.c_str(),ios::out);
	kinetics.precision(14);
	kinetics << fixed << "(A/A0)\tPAB" << endl;
	for (int i=0;i<out.count;i++)
	{
		kinetics << fixed << (out.ysave(iA,i))/A0 << "\t" 
			<<(A0-out.ysave(iA,i)+B0-out.ysave(iB,i))/(A0+B0) << endl;
	}
	kinetics.close();
	kinetics.open(BeckinghamLyndB_file.c_str(),ios::out);
	kinetics.precision(14);
	kinetics << fixed << "(B/A0)\tPAB" << endl;
	for (int i=0;i<out.count;i++)
	{
		kinetics << fixed << (out.ysave(iB,i))/B0 << "\t" 
			<< (A0-out.ysave(iA,i)+B0-out.ysave(iB,i))/(A0+B0) << endl;
	}
	kinetics.close();

	return 0;
}