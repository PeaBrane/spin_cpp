#ifndef HELPER_H_INCLUDED
#define HELPER_H_INCLUDED

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <vector>
#include <stdlib.h>
#include <omp.h>
#include <algorithm>
#include <random>
#include <ctime>
using namespace std;

typedef vector< vector< vector<int> > > Tensor;
typedef vector< vector<int> > Matrix;
typedef vector<int> Row;

static uniform_int_distribution<int> uid(0,1);
static uniform_int_distribution<int> uid2(0,2^32);
static uniform_real_distribution<double> u01(0.0,1.0);

void print_row(Row &row);
void print_mat(Matrix &mat);
void print_mat_order(Matrix &mat);
void print_ten(Tensor &ten);

Row row_empty(int n);
Matrix mat_empty(int n, int m);
Tensor ten_empty(int n, int m, int k);

Row getpart(int nthreads);
Row part(int n, int k);

int randnum(int x, int y);
Matrix randmat(int n, int m);

Row rowdiff(Row &row1, Row &row2);
Row getrow(Matrix &matrix, int dir);
Matrix getmat(Tensor &tensor, int k);
Row subvec(Row &v_in, int i1, int i2);
Matrix submat(Matrix &m_in, int i1, int i2, int j1, int j2);
Tensor subten(Tensor &t_in, int i1, int i2, int j1, int j2);

#endif // HELPER_H_INCLUDED
