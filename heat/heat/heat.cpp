// ConsoleApplication1.cpp : Defines the entry point for the console application.
//


#include "stdafx.h"
#include <iostream>
#include <math.h>
#include "matrix.h"
#include "matrix.cpp"
const double pi = atan(1) * 4;
int N = 6;
int M = 6;


using namespace std;


/**

*/
Vector<double>  Sweeping(TridiagMatrix<double>& matA, Vector<double>& matB)
{
	//matA.consolePrint();
	//matB.consolePrint();
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
	//matRes.consolePrint();
	return matRes;
}

Vector<double> UUoverU(int k, int l, int n, Vector<double> f, double C)
{
	Vector<double> R;
	R.assign(f.size(), 0);
	for (int s = 1; s <= n; s++)
	{
		if ((((k + 1)*s % (n + 1)) == 0) || (((l + 1)*s % (n + 1))==0)) continue;
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
	//cout << k << "!!" << n << "!" << endl;
	for (int s = 1; s <= n; s++)
	{
		double as = 2 * pow(-1, s - 1)*sin((k + 1)*pi*s / (n + 1))*sin(pi*s / (n + 1)) / (n + 1);
		if ((k + 1)*s % (n + 1) == 0) continue;
		double S = C - 2 * cos(pi*s / (n + 1));
		TridiagMatrix<double> I(f.size(), S, -1,-1);
		R = R + Sweeping(I, f*as);
		//cout << "!\n";
	}
	//R.consolePrint();
	return R;
}

void Direct(Vector<Vector<double> >& P, double C)
{
	//cout << (int)log2(N) << endl;
	for (int k = 0; k < (int)log2(N); k++)
	{
		for (int i = (int)pow(2, k)-1; i < N; i += (int)pow(2, k + 1))
		{
			
			int l = i - (int)pow(2, k) ;
			int r = (i + (int)pow(2, k) > N) ? (N) : (i + (int)pow(2, k));
			
			//cout << l << endl;
			//cout << i << endl;
			//cout << r << endl;
			
			if (l >= 0) P[l] = P[l] + UoverU(r - i - 1, r - l - 1, P[i], C);
			if (r < N) P[r] = P[r] + UoverU(i - l - 1, r - l - 1, P[i], C);
		}
	}
}

void Reverse(Vector<Vector<double> >& P, double C)
{
	for (int k = (int)log2(N); k >= 0; k--)
	{
		//cout << k << endl;
		for (int i = (int)pow(2, k)-1; i < N; i += (int)pow(2, k + 1))
		{
			
			int l = i - (int)pow(2, k);
			int r = (i + (int)pow(2, k) > N) ? (N) : (i + (int)pow(2, k));
			//cout << l << endl;
			//cout << i << endl;
			//cout << r << endl;
			P[i] = UUoverU(i - l - 1, r - i - 1, r - l - 1, P[i], C);
			if (l >= 0) P[i] = P[i] + UoverU(r - i - 1, r - l - 1, P[l], C);
			if (r < N) P[i] = P[i] + UoverU(i - l - 1, r - l - 1, P[r], C);
		}
	}
}



int main()
{
	Vector<Vector<double> > P;
	Vector<double> Pi;
	Pi.assign(M, 16);
	P.assign(N, Pi);

	Vector<double> Tup, Tdown, Tleft, Tright;
	Tup.assign(M, 0);
	Tdown.assign(M, 0);
	Tleft.assign(N, 0);
	Tright.assign(N, 0);

	Vector<Vector<double>> F(P);
	//////////////////////////////////////
	for (int i = 0; i <= N-1; i++)
	{
		F[i][0] = F[i][0] + Tleft[i];
		F[i][M-1] = F[i][M-1] + Tright[i];
	}

	for (int j = 0; j <= M-1; j++)
	{
		F[0][j] = F[0][j] + Tup[j];
		F[N-1][j] = F[N-1][j] + Tdown[j];
	}
	//////////////////////////////////////
	Direct(F, 4);
	Reverse(F, 4);
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < M; j++)
		{
			cout << F[i][j] << " ";
		}
		cout << endl;
	}
	char a;
	cin >> a;
	return 0;
}

