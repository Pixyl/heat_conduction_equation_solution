#pragma once
#include <iostream>
#include <vector>

template<class T>
class Vector
{
public:
	std::vector<T> v;
	Vector<T>& operator= (std::vector<T>& v);
	T& operator [](int idx);
	void assign(int size, T value);
	int size();
	void consolePrint();
};


template<class T>
class TridiagMatrix
{
private:
	Vector<T> c;
	Vector<T> a;
	Vector<T> b;
	int n;

public:
	TridiagMatrix();
	TridiagMatrix(int size);
	TridiagMatrix(int size, T data_c);
	TridiagMatrix(int size, T data_c, T data_a, T data_b);
	TridiagMatrix(Vector<T> data_c, Vector<T> data_a, Vector<T> data_b);
	TridiagMatrix(Vector<T> data_c);
	int size();
	Vector<T> getA();
	Vector<T> getB();
	Vector<T> getC();
	T get(int i, int j);
	void setA(Vector<T> data);
	void setB(Vector<T> data);
	void setC(Vector<T> data);
	void set(int i, int j, T data);
	void consolePrint();
};


