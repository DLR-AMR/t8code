// Matrixklasse mit LR-Zerlegung, entnommen aus http://www.igpm.rwth-aachen.de/Download/ss13/na2/na2-base.h
// Numerische Analysis II, SS 2013, Prof. Dr. Wolfgang Dahmen, Dr. Markus Bachmayr

#ifndef __VECMAT_INCLUDE__
#define __VECMAT_INCLUDE__

#include <t8.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <cassert>
#include <cstring>

////////////////////
/// Vektorklasse ///
////////////////////
//T8_EXTERN_C_BEGIN ();

class vec
{
	double * a;
	int n;

public:
	vec() : a(0), n(0)
	{}

	/// Uninitialisierten Vektor der Laenge n anlegen
	explicit vec(int n) : a(0), n(n)
	{
		if (n > 0) {
			a = new double[n];
		}
	}

	/// Copy-Konstruktor
	vec(const vec& v) : a(0), n(v.n)
	{
		if (n > 0) {
			a = new double[n];
		}
		for (int i = 0; i < n; i++) {
			a[i] = v.a[i];
		}
	}

	~vec()
	{
		delete[] a;
		a = nullptr;
	}

	/// Zugriffsoperator
	double& operator()(int i)
	{
		assert(i >= 0 && i < n);
		return a[i];
	}
	/// Operator fuer const-Zugriff
	double operator()(int i) const
	{
		assert(i >= 0 && i < n);
		return a[i];
	}

	vec& operator=(const vec& v)
	{
		if (v.n != n) {
			delete[] a;
			a = new double[v.n];
		}
		n = v.n;
		for (int i = 0; i < n; i++) {
			a[i] = v.a[i];
		}
		return *this;
	}

	/// Alle Eintraege auf den Wert v setzen
	vec& operator=(double v)
	{
		for (int i = 0; i < n; i++) {
			a[i] = v;
		}
		return *this;
	}

	vec& operator+=(const vec& y)
	{
		assert(y.n >= n);
		for (int i = 0; i < n; i++) {
			a[i] += y(i);
		}
		return *this;
	}

	vec& operator-=(const vec& y)
	{
		assert(y.n >= n);
		for (int i = 0; i < n; i++) {
			a[i] -= y(i);
		}
		return *this;
	}

	vec& operator*=(double v)
	{
		for (int i = 0; i < n; i++) {
			a[i] *= v;
		}
		return *this;
	}

	/// Laenge auf nn aendern, voriger Inhalt geht verloren
	/// falls sich die Groesse aendert
	void resize(int nn)
	{
		if (n != nn) {
			delete[] a;
			a = nullptr;
			a = new double[nn];
			n = nn;
		}
	}

	/// returns the length of the vector
	int size() const
	{
		return n;
	}

	friend std::ostream& operator<<(std::ostream& str, const vec& v)
	{
		for (int i = 0; i < v.n; i++) {
			str << v(i) << "\n";
		}
		return str;
	}
};


/// Inneres Produkt zweier Vektoren
inline double inner(const vec& v1, const vec& v2)
{
	assert(v1.size() == v2.size());
	double s = 0.;
	for (int i = 0; i < v1.size(); i++) {
		s += v1(i) * v2(i);
	}
	return s;
}

/// 1-Norm eines Vektors
inline double l1norm(const vec& v)
{
	double s = 0.;
	for (int i = 0; i < v.size(); i++) {
		s += std::abs(v(i));
	}
	return s;
}

/// Euclidische Norm eines Vektors
inline double l2norm(const vec& v)
{
	return std::sqrt(inner(v,v));
}

/// Maximumnorm eines Vektors
inline double linftynorm(const vec& v)
{
	double m = 0.;
	for (int i = 0; i < v.size(); i++) {
		double av = std::abs(v(i));
		if (av > m) {
			m = av;
		}
	}
	return m;
}


/////////////////////////////////////////
/// Klasse fuer vollbesetzte Matrizen ///
/////////////////////////////////////////

class mat
{
	double * a;
	int m, n;

public:
	mat() : a(0), m(0), n(0)
	{}

	/// Erzeugt uninitialisierte Matrix mit m Zeilen und n Spalten
	mat(int m, int n) : a(0), m(m), n(n)
	{
		if (n*m) {
			a = new double[n*m];
		}
	}

	/// Copy-Konstruktor
	mat(const mat& mat) : a(0), m(mat.m), n(mat.n)
	{
		if (n*m) {
			a = new double[n*m];
		}
		std::memcpy(a, mat.a, n*m*sizeof(double));
	}

	~mat()
	{
		delete[] a;
		a = nullptr;
	}

	/// Zugriffsoperator
	double& operator()(int i, int j);
	/// Operator fuer const-Zugriff
	double operator()(int i, int j) const;

	mat& operator=(const mat& mat);

	/// Alle Eintraege auf Wert v setzen
	mat& operator=(double v);

	/// Groesse auf mm Zeilen, nn Spalten setzen,
	/// voriger Inhalt geht verloren falls sich die Groesse aendert
	void resize(int mm, int nn);

	/// Anzahl Zeilen
	int rows() const;

	/// Anzahl Spalten
	int cols() const;


	friend std::ostream& operator<<(std::ostream& str, const mat& M)
	{
		for (int i=0; i<M.m; i++) {
			for (int j=0; j<M.n; j++) {
				str << M(i,j) << ' ';
			}
			str << std::endl;
		}
		return str;
	}
	/// LR-Zerlegung fuer vollbesetzte Matrizen
	// (Zeile i vertauscht mit Zeile r[i] >= i)
	void lr_factors(mat& A, std::vector<int>& r);

	void lr_solve(const mat& A, const std::vector<int>& r, vec& x);

};

//T8_EXTERN_C_END ();
#endif
