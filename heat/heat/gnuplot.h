#pragma once
// gnuplot.h
#ifndef _GNUPLOT_H_
#define _GNUPLOT_H_

#include <cstdio>
#include <string>
#include <iostream>
#include <stdexcept>



#define GNUPLOT_NAME "F:\Programes\gnuplot\bin\wgnuplot.exe"


using std::string;

class Gnuplot
{
public:
	Gnuplot();
	~Gnuplot();
	void operator ()(const string & command); // отправить команду gnuplot

protected:
	FILE *gnuplotpipe;
};

Gnuplot::Gnuplot()
{
	gnuplotpipe = _popen("\"F:\\Programes\\gnuplot\\bin\\gnuplot.exe\"", "w");


	if (!gnuplotpipe)
	{
		std::cout << ("Gnuplot not found !");
	}
}
Gnuplot::~Gnuplot()
{
	fprintf(gnuplotpipe, "exit\n");

	_pclose(gnuplotpipe);

}
void Gnuplot::operator()(const string & command)
{
	fprintf(gnuplotpipe, "%s\n", command.c_str());
	fflush(gnuplotpipe); //без fflush ничего рисоваться не будет
};

#endif // #ifndef _GNUPLOT_H_