#include "stdafx.h"
#include "matrix.h"


template<class T>
Vector<T> operator+ (Vector<T>& v1, Vector<T>& v2)
{
	Vector<T> res(v1);
	if (v1.v.size() == v2.v.size())
	{
		int n = v1.v.size();
		res.v.assign(n, 0);
		for (int i = 0; i < n; i++)
		{
			res.v[i] = v1.v[i] + v2.v[i];

		}
	}
	return res;
}
template<class T>
void Vector<T>::consolePrint() 
{
	int n = size();
	for (int i = 0; i < n; i++)
	{
		std::cout << (*this)[i] << " ";
	}
	std::cout << std::endl;
}

template<class T>
Vector<T> operator- (Vector<T>& v2, Vector<T>& v1)
{
	Vector<T> res(v1);
	if (v1.v.size() == v2.v.size())
	{
		int n = v1.v.size();
		res.v.assign(n, 0);
		for (int i = 0; i < n; i++)
		{
			res.v[i] = v1.v[i] - v2.v[i];
		}
	}
	return res;
}

template<class T>
Vector<T> operator* (Vector<T>& v2, T a)
{
	Vector<T> res(v2);
	int n = res.size();
	for (int i = 0; i < n; i++)
	{
		res[i] = res[i] * a;
	}
	return res;
}


template<class T>
Vector<T>& Vector<T>::operator= (std::vector<T>& v2)
{
	(*this).v = v2;
}


template<class T>
T& Vector<T>::operator[] (int idx)
{
	return (*this).v[idx];
}

template<class T>
void Vector<T>::assign(int size, T value)
{
	(*this).v.assign(size, value);
}

template<class T>
int Vector<T>::size()
{
	return (*this).v.size();
}




template<class T>
TridiagMatrix<T>::TridiagMatrix(int size, T data_c, T data_a, T data_b)
{
	n = size;
	a.assign(size - 1, data_a);
	b.assign(size - 1, data_b);
	c.assign(size, data_c);

}

template<class T>
TridiagMatrix<T>::TridiagMatrix() : TridiagMatrix(3, 0, 0, 0) {}

template<class T>
TridiagMatrix<T>::TridiagMatrix(int size) : TridiagMatrix(size, 0, 0, 0) {}

template<class T>
TridiagMatrix<T>::TridiagMatrix(int size, T data_c) : TridiagMatrix(size, data_c, 0, 0) {}

template<class T>
TridiagMatrix<T>::TridiagMatrix(Vector<T> data_c, Vector<T> data_a, Vector<T> data_b)
{
	if ((data_a.size() == data_b.size()) && ((data_a.size() + 1) == data_c.size()))
	{
		n = data_c.size();
		a = data_a;
		b = data_b;
		c = data_c;
	}
}

template<class T>
TridiagMatrix<T>::TridiagMatrix(Vector<T> data_c)
{
	n = data_c.size();
	c = data_c;
	a.assign(n - 1, 0);
	b.assign(n - 1, 0);
}

template<class T>
Vector<T> TridiagMatrix<T>::getA()
{
	return a;
}

template<class T>
Vector<T> TridiagMatrix<T>::getB()
{
	return b;
}

template<class T>
Vector<T> TridiagMatrix<T>::getC()
{
	return c;
}

template<class T>
void TridiagMatrix<T>::setA(Vector<T> data)
{
	if (data.size() == n - 1)
		a = data;
}
template<class T>
void TridiagMatrix<T>::setB(Vector<T> data)
{
	if (data.size() == n - 1)
		b = data;
}

template<class T>
void TridiagMatrix<T>::setC(Vector<T> data)
{
	if (data.size() == n)
		c = data;
}
template<class T>
int TridiagMatrix<T>::size()
{
	return n;
}
template<class T>
T TridiagMatrix<T>::get(int i, int j)
{
	if (i < n && j < n)
	{
		if (i > -1 && j > -1)
		{
			if (i == j) return c[i];
			if (i == j - 1) return b[i];
			if (i == j + 1) return a[j];
			return 0;
		}
	}
	return NULL;
}

template<class T>
void TridiagMatrix<T>::set(int i, int j, T data)
{

	if (i < n && j < n)
	{
		if (i > -1 && j > -1)
		{
			if (i == j) c[i] = data;
			if (i == j - 1) b[i] = data;
			if (i == j + 1) a[j] = data;
		}
	}
}

template<class T>
void TridiagMatrix<T>::consolePrint()
{
	for (int i = 0; i < size(); i++)
	{
		for (int j = 0; j < size(); j++)
			std::cout << get(i, j) << " ";
		std::cout << std::endl;
	}
}
template<class T>
TridiagMatrix<T> operator+ (TridiagMatrix<T>& T1, TridiagMatrix<T>& T2)
{
	if (T1.size() == T2.size())
	{
		TridiagMatrix<T> Res(T1.getC() + T2.getC(), T1.getA() + T2.getA(), T1.getB() + T2.getB());
		return Res;
	}
	return T1;
}
template<class T>
TridiagMatrix<T> operator- (TridiagMatrix<T> T1, TridiagMatrix<T> T2)
{
	if (T1.size() == T2.size())
	{
		TridiagMatrix<T> Res(T1.getC() - T2.getC(), T1.getA() - T2.getA(), T1.getB() - T2.getB());
		return Res;
	}
	return T1;
}

template<class T>
TridiagMatrix<T> operator* (TridiagMatrix<T>& T1, T a)
{
	TridiagMatrix<T> Res(T1.getC()*a, T1.getA()*a, T1.getB()*a);

	return Res;
}


template<class T>
Vector<T> operator* (TridiagMatrix<T>& T1, Vector<T>& a)
{
	Vector<T> Rez(a);
	if (T1.size() == a.size())
	{
		int n = T1.size();
		Rez.assign(n, 0);
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				Rez[i] = Rez[i] + a[j] * T1.get(i, j);
			}
		}
	}
	return Rez;
}

