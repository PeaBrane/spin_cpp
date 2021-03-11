#include "helper.h"

void print_row(Row &row)
{

    int n = row.size();
    for (int i=0; i<n; i++)
    {
        cout << row[i] << " ";
    }

}

void print_mat(Matrix &mat)
{

    int n = mat.size();
    int m = mat[0].size();

    for (int i=0; i<n; i++)
    {
        for (int j=0; j<m; j++)
        {
            cout << mat[i][j] << " ";
            if (mat[i][j]>=0)
            {
                cout << " ";
            }
        }
        cout << endl;
    }
}

void print_mat_order(Matrix &mat)
{

    int n = mat.size();
    int m = mat[0].size();

    for (int j=m-1; j>=0; j--)
    {
        for (int i=0; i<n; i++)
        {
            cout << mat[i][j] << " ";
            if (mat[i][j]>=0)
            {
                cout << " ";
            }
        }
        cout << endl;
    }
}

void print_ten(Tensor &ten)
{

    int n = ten.size();
    int m = ten[0].size();
    int l = ten[0][0].size();

    for (int k=0; k<l; k++)
    {
        for (int i=0; i<n; i++)
        {
            for (int j=0; j<m; j++)
            {
                cout << ten[i][j][k] << " ";
            }
            cout << endl;
        }
        cout << endl;
    }

}

Row row_empty(int n)
{
    Row temp;
    temp.resize(n);
    fill(temp.begin(),temp.end(),0);
    return temp;
}

Matrix mat_empty(int n, int m)
{
    Matrix temp;
    temp.resize(n);
    fill(temp.begin(),temp.end(),row_empty(m));
    return temp;
}

Tensor ten_empty(int n, int m, int k)
{
    Tensor temp;
    temp.resize(n);
    fill(temp.begin(),temp.end(),mat_empty(m,k));
    return temp;
}

Row getpart(int nthreads)
{

    Row row;

    if (nthreads == 2)
    {
        row = {1,2};
    }
    else if (nthreads == 3)
    {
        row = {1,3};
    }
    else if (nthreads == 4)
    {
        row = {2,2};
    }
    else if (nthreads == 5)
    {
        row = {1,5};
    }
    else if (nthreads == 6)
    {
        row = {2,3};
    }
    else if (nthreads == 7)
    {
        row = {1,7};
    }
    else if (nthreads == 8)
    {
        row = {2,4};
    }
    else if (nthreads == 9)
    {
        row = {3,3};
    }
    else if (nthreads == 10)
    {
        row = {2,5};
    }
    else if (nthreads == 11)
    {
        row = {1,11};
    }
    else if (nthreads == 12)
    {
        row = {3,4};
    }
    else if (nthreads == 13)
    {
        row = {1,13};
    }
    else if (nthreads == 14)
    {
        row = {2,7};
    }
    else if (nthreads == 15)
    {
        row = {3,5};
    }
    else if (nthreads == 16)
    {
        row = {4,4};
    }
    else if (nthreads == 17)
    {
        row = {1,17};
    }
    else if (nthreads == 18)
    {
        row = {3,6};
    }
    else if (nthreads == 19)
    {
        row = {1,19};
    }
    else if (nthreads == 20)
    {
        row = {4,5};
    }
    else if (nthreads == 21)
    {
        row = {3,7};
    }
    else if (nthreads == 22)
    {
        row = {2,11};
    }
    else if (nthreads == 23)
    {
        row = {1,23};
    }
    else if (nthreads == 24)
    {
        row = {4,6};
    }

    return row;

}

Row part(int n, int k)
{

    Row temp;
    double incre = (double)n/k;
    for (int i=1; i<=k; i++)
    {
        int x = (int)round(incre*i);
        temp.push_back(x-1);
    }
    return temp;

}

int randnum(int x, int y)
{
    int diff = y-x;
    diff = rand() % (diff+1);
    return x+diff;
}

Matrix randmat(int n, int m)
{

    Row temp;
    for (int i=0; i<n; i++)
    {
        for (int j=0; j<m; j++)
        {
            temp.push_back(i*m+j);
        }
    }
    random_shuffle(temp.begin(),temp.end());

    Matrix mat;
    for (int i=0; i<n; i++)
    {
        Row row;
        for (int j=0; j<m; j++)
        {
            row.push_back(temp[i*m+j]);
        }
        mat.push_back(row);
    }

    return mat;

}

Row rowdiff(Row &row1, Row &row2)
{

    int n = row1.size();
    Row diff;

    for (int i=0; i<n; i++)
    {
        int x = row1[i]-row2[i];
        diff.push_back(x);
    }

    return diff;

}

Row getrow(Matrix &matrix, int dir)
{

    Row row;
    int n = matrix.size(); int m = matrix[0].size();

    if (dir == 1)
    {
        row = matrix[0];
    }
    else if (dir == 3)
    {
        row = matrix[n-1];
    }
    else if (dir == 2)
    {
        for (int i=0; i<n; i++)
        {
            int ele = matrix[i][0];
            row.push_back(ele);
        }
    }
    else if (dir == 4)
    {
        for (int i=0; i<n; i++)
        {
            int ele = matrix[i][m-1];
            row.push_back(ele);
        }
    }

    return row;

}

Matrix getmat(Tensor &tensor, int k)
{

    Matrix matrix;
    int n = tensor.size(); int m = tensor[0].size();

    for (int i=0; i<n; i++)
    {
        Row row;
        for (int j=0; j<m; j++)
        {
            int ele = tensor[i][j][k];
            row.push_back(ele);
        }
        matrix.push_back(row);
    }

    return matrix;

}

Row subvec(Row &v_in, int i1, int i2)
{
    Row::const_iterator first = v_in.begin() + i1;
    Row::const_iterator last = v_in.begin() + i2;
    Row temp(first, last);
    return temp;
}

Matrix submat(Matrix &m_in, int i1, int i2, int j1, int j2)
{

    Matrix::const_iterator first = m_in.begin() + i1;
    Matrix::const_iterator last = m_in.begin() + i2;
    Matrix m_cut(first,last);
    int l = m_cut.size();

    Row vec;
    Matrix temp;
    for (int i=0; i<l; i++)
    {
        vec = m_cut[i];
        Row::const_iterator first = vec.begin() + j1;
        Row::const_iterator last = vec.begin() + j2;
        Row vtemp(first, last);
        temp.push_back(vtemp);
    }

    return temp;

}

Tensor subten(Tensor &t_in, int i1, int i2, int j1, int j2)
{

    Tensor::const_iterator first = t_in.begin() + i1;
    Tensor::const_iterator last = t_in.begin() + i2;
    Tensor t_cut(first,last);
    int l = t_cut.size();

    Matrix mat;
    Tensor temp;
    for (int i=0; i<l; i++)
    {
        mat = t_cut[i];
        Matrix::const_iterator first = mat.begin() + j1;
        Matrix::const_iterator last = mat.begin() + j2;
        Matrix mtemp(first, last);
        temp.push_back(mtemp);
    }

    return temp;

}
