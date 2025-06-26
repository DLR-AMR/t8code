#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <cassert>
#include <cstring>
#include "vecmat.hxx"
#include <t8.h>
/// Zugriffsoperator
double&
mat::operator() (int i, int j)
{
  assert (i >= 0 && i < m && j >= 0 && j < n);
  return a[m * j + i];
}
/// Operator fuer const-Zugriff
double
mat::operator() (int i, int j) const
{
  assert (i >= 0 && i < m && j >= 0 && j < n);
  return a[m * j + i];
}

mat&
mat::operator= (const mat& mat)
{
  if (mat.n * mat.m != m * n) {
    delete[] a;
    a = new double[mat.m * mat.n];
  }
  m = mat.m;
  n = mat.n;
  for (int k = 0; k < m * n; k++) {
    a[k] = mat.a[k];
  }
  return *this;
}

/// Alle Eintraege auf Wert v setzen
mat&
mat::operator= (double v)
{
  for (int k = 0; k < m * n; k++) {
    a[k] = v;
  }
  return *this;
}

/// Groesse auf mm Zeilen, nn Spalten setzen,
/// voriger Inhalt geht verloren falls sich die Groesse aendert
void
mat::resize (int mm, int nn)
{
  if (m != mm || n != nn) {
    delete[] a;
    a = nullptr;
    a = new double[mm * nn];
    m = mm;
    n = nn;
  }
}

/// Anzahl Zeilen
int
mat::rows () const
{
  return m;
}
/// Anzahl Spalten
int
mat::cols () const
{
  return n;
}

/// LR-Zerlegung fuer vollbesetzte Matrizen
// (Zeile i vertauscht mit Zeile r[i] >= i)
void
mat::lr_factors (mat& A, std::vector<int>& r)
{
  const int n = A.rows ();
  r.resize (n);

  for (int j = 0; j < n; j++) {
    int piv = j;
    double Ajmax = std::abs (A (j, j));
    for (int p = j + 1; p < n; p++) {
      double Ap = std::abs (A (p, j));
      if (Ap > Ajmax) {
        piv = p;
        Ajmax = Ap;
      }
    }
    r[j] = piv;
    if (piv != j) {
      for (int k = 0; k < n; k++) {
        std::swap (A (piv, k), A (j, k));
      }
    }
    for (int i = j + 1; i < n; i++) {
      A (i, j) /= A (j, j);
      for (int k = j + 1; k < n; k++) {
        A (i, k) -= A (i, j) * A (j, k);
      }
    }
  }
}

void
mat::lr_solve (const mat& A, const std::vector<int>& r, vec& x)
{
  const int n = A.rows ();
  for (int i = 0; i < n; i++) {
    if (i != r[i]) {
      std::swap (x (i), x (r[i]));
    }
  }
  for (int i = 1; i < n; i++) {
    for (int j = 0; j < i; j++) {
      x (i) -= A (i, j) * x (j);
    }
  }
  for (int i = n - 1; i >= 0; i--) {
    for (int j = i + 1; j < n; j++) {
      x (i) -= A (i, j) * x (j);
    }
    x (i) /= A (i, i);
  }
}
