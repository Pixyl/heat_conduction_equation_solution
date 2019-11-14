// ConsoleApplication1.cpp : Defines the entry point for the console application.
//



#include <iostream>
#include <math.h>
#include <fstream> 
#include <string> 
#include <sstream>
#include "windows.h"
#include "matrix.h"
#include "matrix.cpp"
#include <chrono>


const double pi = atan(1) * 4;
int N;
int M;


using namespace std;

Vector<double>  Sweeping(TridiagMatrix<double>& matA, Vector<double>& matB)
{

	const int N = matA.size();
	const int N1 = N - 1;
	double y;
	Vector<double> a, B, matRes;
	matRes.assign(N, 0);
	a.assign(N, 0);
	B.assign(N, 0);
	y = matA.get(0, 0);
	a[0] = -matA.get(0, 1) / y;
	B[0] = matB[0] / y;
	for (int i = 1; i < N1; i++) {
		y = matA.get(1, 1) + matA.get(i, i - 1) * a[i - 1];
		a[i] = -matA.get(i, i + 1) / y;
		B[i] = (matB[i] - matA.get(i, i - 1)* B[i - 1]) / y;
	}

	matRes[N1] = (matB[N1] - matA.get(N1, N1 - 1)* B[N1 - 1]) / (matA.get(N1, N1) + matA.get(N1, N1 - 1)* a[N1 - 1]);
	for (int i = N1 - 1; i >= 0; i--) {
		matRes[i] = a[i] * matRes[i + 1] + B[i];
	}
	return matRes;
}

Vector<double> UUoverU(int k, int l, int n, Vector<double> f, double C)
{
	Vector<double> R;
	R.assign(f.size(), 0);
	for (int s = 1; s <= n; s++)
	{
		if ((((k + 1)*s % (n + 1)) == 0) || (((l + 1)*s % (n + 1)) == 0)) continue;
		double as = 2 * pow(-1, s - 1)*sin((k + 1)*pi*s / (n + 1))*sin((l + 1)*pi*s / (n + 1)) / (n + 1);
		double S = C - 2 * cos(pi*s / (n + 1));
		TridiagMatrix<double> I(f.size(), S, -1, -1);
		R = R + Sweeping(I, f*as);
	}
	return R;
}

Vector<double> UoverU(int k, int n, Vector<double> f, double C)
{
	Vector<double> R, F;
	R.assign(f.size(), 0);
	for (int s = 1; s <= n; s++)
	{
		double as = 2 * pow(-1, s - 1)*sin((k + 1)*pi*s / (n + 1))*sin(pi*s / (n + 1)) / (n + 1);
		if ((k + 1)*s % (n + 1) == 0) continue;
		double S = C - 2 * cos(pi*s / (n + 1));
		TridiagMatrix<double> I(f.size(), S, -1, -1);
		R = R + Sweeping(I, f*as);
	}
	return R;
}



void Direct(Vector<Vector<double> >& P, double C)
{
	for (int k = 0; k < (int)log2(N); k++)
	{  
		int step = (int)pow(2, k + 1);
		
		for (int i = step - 1; i < N; i += step)
		{
			int l = i - (int)(step / 2);
			int r = (i + (int)(step / 2) > N) ? (N) : (i + (int)(step / 2));
			int ll = i - step;
			int rr = ((i + step) > N) ? (N) : (i + step);

			if (r < N) P[i] = P[i] + UoverU(rr - r - 1, rr - i - 1, P[r], C);
			if (l >= 0) P[i] = P[i] + UoverU(l - ll - 1, i - ll - 1, P[l], C);
		}

	}
}



void Reverse(Vector<Vector<double> >& P, double C)
{
	for (int k = (int)log2(N); k >= 0; k--)
	{
		int step = (int)pow(2, k + 1);
		

		for (int i = (int)(step/2) - 1; i < N; i += step)
		{
			int l = i - (int)(step / 2);
			int r = (i + (int)(step / 2) > N) ? (N) : (i + (int)(step / 2));
			P[i] = UUoverU(i - l - 1, r - i - 1, r - l - 1, P[i], C);
			if (l >= 0) P[i] = P[i] + UoverU(r - i - 1, r - l - 1, P[l], C);
			if (r < N) P[i] = P[i] + UoverU(i - l - 1, r - l - 1, P[r], C);
		}
	}
}

Vector<Vector<double>> tick(Vector<Vector<double>> P, Vector<double> Tup, Vector<double> Tdown, Vector<double> Tleft, Vector<double> Tright, double C)
{
	Vector<Vector<double>> F(P);
	//////////////////////////////////////

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < M; j++)
		{
			F[i][j] = F[i][j] * C;
		}
	}

	for (int i = 0; i <= N - 1; i++)
	{
		F[i][0] = F[i][0] + Tleft[i];
		F[i][M - 1] = F[i][M - 1] + Tright[i];
	}

	for (int j = 0; j <= M - 1; j++)
	{
		F[0][j] = F[0][j] + Tup[j];
		F[N - 1][j] = F[N - 1][j] + Tdown[j];
	}
	//////////////////////////////////////


	Direct(F, C + 4);
	Reverse(F, C + 4);
	return F;
}

void printMat(Vector<Vector<double>> F)
{
	system("cls");
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < M; j++)
		{
			cout << F[i][j] << " ";
		}
		cout << endl;
	}

	cout << endl;
	Sleep(1000);
}



int main()
{
	double p = 1, c = 1, lambda = 1;
	cout << "Enter N,M" << endl;
	cin >> N >> M;
	cout << "Enter p, c, lambda" << endl;
	cin >> p >> c >> lambda;
	double a = p*c / lambda;
	double h1 = 1 / (double)N, h2 = 1 / (double)M, t = 0.001;
	double C = a*h1*h2 / t;

	cout << "C = " << C << endl;

	Vector<Vector<double> > P, F;
	Vector<double> Pi;
	Pi.assign(M, 1);
	P.assign(N, Pi);



	//printMat(P);

	Vector<double> Tup, Tdown, Tleft, Tright;
	Tup.assign(M, 5);
	Tdown.assign(M, 5);
	Tleft.assign(N, 5);
	Tright.assign(N, 5);

	F = P;
	//while (true)
	//{
	std::chrono::steady_clock::time_point start;
	start = std::chrono::steady_clock::now();
	F = tick(F, Tup, Tdown, Tleft, Tright, C);
	std::chrono::steady_clock::time_point finish;
	finish = std::chrono::steady_clock::now();
	int time = std::chrono::duration_cast<std::chrono::milliseconds>(finish - start).count();
	std::cout << "time = " << time << " ms" << std::endl;
	//printMat(F);
	//}

	char ch;
	cin >> ch;
	return 0;
}