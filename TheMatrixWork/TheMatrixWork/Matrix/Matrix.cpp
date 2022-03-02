#include "Matrix.h"

std::ofstream& operator<<(std::ofstream& out, Matrix const& matrix){
    int maxLength=0;
    for(int i = 0;i < matrix.m;i++){
        for(int j = 0;j < matrix.n;j++){
            std::stringstream ss;
            std::string count;
            ss << std::setprecision(15) << matrix.matrix[i][j];
            count = ss.str();
            if(count.size()>maxLength){
                maxLength = static_cast<int>(count.size());
            }
        }
    }
    maxLength++;
    if(matrix.n==1&&matrix.m==1){
        out << std::setw(maxLength) << "( " << matrix.matrix[0][0] << " )";
    }
    else if(matrix.n==1&&matrix.m!=1){
        for (int i = 0; i < matrix.m; i++)
        {
            for (int j = 0; j < matrix.n; j++)
            {
                if(j==0&&i==0){
                    out << std::setw(maxLength) << "/ " << matrix.matrix[i][j];
                    std::stringstream ss;
                    std::string count;
                    ss << std::setprecision(15) << matrix.matrix[i][j];
                    count = ss.str();
                    out << std::setw(2*static_cast<int>(maxLength)-static_cast<int>(count.size())) << " \\";
                }
                else if(j==0&&i==matrix.m-1){
                    out << std::setw(maxLength) << "\\ " << matrix.matrix[i][j];
                    std::stringstream ss;
                    std::string count;
                    ss << std::setprecision(15) << matrix.matrix[i][j];
                    count = ss.str();
                    out << std::setw(2*static_cast<int>(maxLength)-static_cast<int>(count.size())) << " /";
                }
                else if(j==0){
                    out << std::setw(maxLength) << "| " << matrix.matrix[i][j];
                    std::stringstream ss;
                    std::string count;
                    ss << std::setprecision(15) << matrix.matrix[i][j];
                    count = ss.str();
                    out << std::setw(2*static_cast<int>(maxLength)-static_cast<int>(count.size())) << " |";
                }
            }
            out << "\n";
        }
    }
    else if(matrix.m==1&&matrix.n!=1){
        for (int i = 0; i < matrix.m; i++)
        {
            for (int j = 0; j < matrix.n; j++)
            {
                if(j==0&&i==0){
                    out << std::setw(maxLength) << "( " << matrix.matrix[i][j];
                }
                else if(j==matrix.n-1&&i==0){
                    out << std::setw(maxLength) << matrix.matrix[i][j] << " )";
                }
                else{
                    std::stringstream ss;
                    std::string count;
                    ss << std::setprecision(15) << matrix.matrix[i][j-1];
                    count = ss.str();
                    if(count.size() == maxLength){
                        out << std::setw(static_cast<int>(maxLength)) << matrix.matrix[i][j];
                    }
                    else{
                        out << std::setw(2*static_cast<int>(maxLength)-static_cast<int>(count.size())) << matrix.matrix[i][j];
                    }
                }
            }
            out << "\n";
        }
    }
    else{
        for (int i = 0; i < matrix.m; i++)
        {
            for (int j = 0; j < matrix.n; j++)
            {
                if(j==0&&i==0){
                    out << std::setw(maxLength) << "/ " << matrix.matrix[i][j];
                }
                else if(j==0&&i==matrix.m-1){
                    out << std::setw(maxLength) << "\\ " << matrix.matrix[i][j];
                }
                else if(j==0){
                    out << std::setw(maxLength) << "| " << matrix.matrix[i][j];
                }
                else if(j==matrix.n-1&&i==0){
                    out << std::setw(maxLength) << matrix.matrix[i][j] << " \\";
                }
                else if(j==matrix.n-1&&i==matrix.m-1){
                    out << std::setw(maxLength) << matrix.matrix[i][j] << " /";
                }
                else if(j==matrix.n-1){
                    out << std::setw(maxLength) << matrix.matrix[i][j] << " |";
                }
                else{
                    std::stringstream ss;
                    std::string count;
                    ss << std::setprecision(15) << matrix.matrix[i][j-1];
                    count = ss.str();
                    if(count.size() == maxLength){
                        out << std::setw(static_cast<int>(maxLength)) << matrix.matrix[i][j];
                    }
                    else{
                        out << std::setw(2*static_cast<int>(maxLength)-static_cast<int>(count.size())) << matrix.matrix[i][j];
                    }
                }
            }
            out << "\n";
        }
    }
    out<<"\n";
    return out;
}


std::ostream& operator<<(std::ostream& out, const Matrix& matrix){
    int maxLength=0;
    for(int i = 0;i < matrix.m;i++){
        for(int j = 0;j < matrix.n;j++){
            std::stringstream ss;
            std::string count;
            ss << std::setprecision(15) << matrix.matrix[i][j];
            count = ss.str();
            if(count.size()>maxLength){
                maxLength = static_cast<int>(count.size());
            }
        }
    }
    maxLength++;
    if(matrix.n==1&&matrix.m==1){
        out << std::setw(maxLength) << "( " << matrix.matrix[0][0] << " )";
    }
    else if(matrix.n==1&&matrix.m!=1){
        for (int i = 0; i < matrix.m; i++)
        {
            for (int j = 0; j < matrix.n; j++)
            {
                if(j==0&&i==0){
                    out << std::setw(maxLength) << "/ " << matrix.matrix[i][j];
                    std::stringstream ss;
                    std::string count;
                    ss << std::setprecision(15) << matrix.matrix[i][j];
                    count = ss.str();
                    out << std::setw(2*static_cast<int>(maxLength)-static_cast<int>(count.size())) << " \\";
                }
                else if(j==0&&i==matrix.m-1){
                    out << std::setw(maxLength) << "\\ " << matrix.matrix[i][j];
                    std::stringstream ss;
                    std::string count;
                    ss << std::setprecision(15) << matrix.matrix[i][j];
                    count = ss.str();
                    out << std::setw(2*static_cast<int>(maxLength)-static_cast<int>(count.size())) << " /";
                }
                else if(j==0){
                    out << std::setw(maxLength) << "| " << matrix.matrix[i][j];
                    std::stringstream ss;
                    std::string count;
                    ss << std::setprecision(15) << matrix.matrix[i][j];
                    count = ss.str();
                    out << std::setw(2*static_cast<int>(maxLength)-static_cast<int>(count.size())) << " |";
                }
            }
            out << std::endl;
        }
    }
    else if(matrix.m==1&&matrix.n!=1){
        for (int i = 0; i < matrix.m; i++)
        {
            for (int j = 0; j < matrix.n; j++)
            {
                if(j==0&&i==0){
                    out << std::setw(maxLength) << "( " << matrix.matrix[i][j];
                }
                else if(j==matrix.n-1&&i==0){
                    out << std::setw(maxLength) << matrix.matrix[i][j] << " )";
                }
                else{
                    std::stringstream ss;
                    std::string count;
                    ss << std::setprecision(15) << matrix.matrix[i][j-1];
                    count = ss.str();
                    if(count.size() == maxLength){
                        out << std::setw(static_cast<int>(maxLength)) << matrix.matrix[i][j];
                    }
                    else{
                        out << std::setw(2*static_cast<int>(maxLength)-static_cast<int>(count.size())) << matrix.matrix[i][j];
                    }
                }
            }
            out << std::endl;
        }
    }
    else{
        for (int i = 0; i < matrix.m; i++)
        {
            for (int j = 0; j < matrix.n; j++)
            {
                if(j==0&&i==0){
                    out << std::setw(maxLength) << "/ " << matrix.matrix[i][j];
                }
                else if(j==0&&i==matrix.m-1){
                    out << std::setw(maxLength) << "\\ " << matrix.matrix[i][j];
                }
                else if(j==0){
                    out << std::setw(maxLength) << "| " << matrix.matrix[i][j];
                }
                else if(j==matrix.n-1&&i==0){
                    out << std::setw(maxLength) << matrix.matrix[i][j] << " \\";
                }
                else if(j==matrix.n-1&&i==matrix.m-1){
                    out << std::setw(maxLength) << matrix.matrix[i][j] << " /";
                }
                else if(j==matrix.n-1){
                    out << std::setw(maxLength) << matrix.matrix[i][j] << " |";
                }
                else{
                    std::stringstream ss;
                    std::string count;
                    ss << std::setprecision(15) << matrix.matrix[i][j-1];
                    count = ss.str();
                    if(count.size() == maxLength){
                        out << std::setw(static_cast<int>(maxLength)) << matrix.matrix[i][j];
                    }
                    else{
                        out << std::setw(2*static_cast<int>(maxLength)-static_cast<int>(count.size())) << matrix.matrix[i][j];
                    }
                }
            }
            out << std::endl;
        }
    }
    out << '\n';
    return out;
}


int MaxLengthOfVector(std::vector<double> vector){
    int maxLength = 0;
    for(int i = 0;i < vector.size();i++){
        std::stringstream ss;
        std::string count;
        ss << std::setprecision(15) << vector[i];
        count = ss.str();
        if(count.size() > maxLength){
            maxLength = static_cast<int>(count.size());
        }
    }
    return maxLength;
}


std::istream& operator>>(std::istream &in, Matrix& matrix){
    matrix.matrix.resize(matrix.m);
    for(int i = 0;i<matrix.m;++i){
        matrix.matrix[i].resize(matrix.n);
        for(int j = 0;j<matrix.n;++j){
            in>>matrix[i][j];
        }
    }
    return in;
}


double SquareMatrix::det() {//Определитель
    long int n = matrix.size();
    double res = 1;
    
    std::vector<std::vector<double>> tmp = matrix;


       for(int col = 0; col < n; ++col) {
          bool found = false;
          for(int row = col; row < n; ++row) {
             if(tmp[row][col]) {
                if ( row != col )
                {
                    tmp[row].swap(tmp[col]);
                }
                found = true;
                break;
             }
          }
          if(!found) {
             return 0;
          }
          for(int row = col + 1; row < n; ++row) {
             while(true) {
                int del = tmp[row][col] / tmp[col][col];
                for (int j = col; j < n; ++j) {
                    tmp[row][j] -= del * tmp[col][j];
                }
                if (tmp[row][col] == 0)
                {
                   break;
                }
                else
                {
                    res*=-1;
                    tmp[row].swap(tmp[col]);
                }
             }
          }
       }
       for(int i = 0; i < n; ++i) {
          res *= tmp[i][i];
       }
       return res;
    }

Matrix operator+(const Matrix A, const Matrix B){
    Matrix c(A.m, B.n);
    assert(A.n==B.n);
    assert(A.m==B.m);
    c.matrix.resize(A.m);
    for(int i = 0;i<A.m;++i){
        c.matrix[i].resize(B.n);
        for(int j = 0;j<B.n;++j){
            c.matrix[i][j] = A.matrix[i][j] + B.matrix[i][j];
        }
    }
    return c;
}


Matrix operator-(const Matrix A, const Matrix B){
    Matrix c(A.m, B.n);
    assert(A.n==B.n);
    assert(A.m==B.m);
    c.matrix.resize(A.m);
    for(int i = 0;i<A.m;++i){
        c.matrix[i].resize(B.n);
        for(int j = 0;j<B.n;++j){
            c.matrix[i][j] = A.matrix[i][j] - B.matrix[i][j];
        }
    }
    return c;
}

Matrix operator*(const Matrix A, const Matrix B){
    Matrix c(A.m, B.n);
    assert(A.n==B.m);
    c.matrix.resize(A.m);
    for(int i = 0;i<A.m;++i){
        c.matrix[i].resize(B.n);
        for(int j = 0;j<B.n;++j){
            c.matrix[i][j] = 0;
            for(int s=0;s<A.n;++s){
                c.matrix[i][j] += A.matrix[i][s]*B.matrix[s][j];
            }
        }
    }
    return c;
}

Matrix& operator*(Matrix& A, const double alpha){
    for(int i = 0;i<A.m;++i){
        for(int j = 0;j<A.n;++j){
            A.matrix[i][j]*=alpha;
        }
    }
    return A;
}

Matrix& operator*(const double alpha, Matrix& A){
    for(int i = 0;i<A.m;++i){
        for(int j = 0;j<A.n;++j){
            A.matrix[i][j]*=alpha;
        }
    }
    return A;
}

Matrix operator%(const Matrix A, const Matrix B){
    Matrix c(A.m, B.n);
    assert(A.n==B.n);
    assert(A.m==B.m);
    c.matrix.resize(A.m);
    for(int i = 0;i<A.m;++i){
        c.matrix[i].resize(B.n);
        for(int j = 0;j<B.n;++j){
            c.matrix[i][j] = A.matrix[i][j] * B.matrix[i][j];
        }
    }
    return c;
}


double ScalarProduct(Matrix vector1, Matrix vector2){
    double scalarProduct = 0;
    assert(vector1.GetN() == vector2.GetM());
    assert(vector1.GetM() == vector2.GetN());
    assert(vector1.GetM()==1);
    for(int i = 0;i<vector1.GetM();++i){
        for(int j = 0;j<vector2.GetN();++j){
            for(int s=0;s<vector1.GetN();++s){
                scalarProduct += vector1.GetMatrixElement(i,s)*vector2.GetMatrixElement(s,j);
            }
        }
    }
    return scalarProduct;
}

double operator*(const Vector a, const Vector b){
    double scalarProduct = 0;
    assert(a.n == b.m);
    assert(a.m==1);
    for(int i = 0;i<a.m;++i){
        for(int j = 0;j<b.n;++j){
            for(int s=0;s<a.n;++s){
                scalarProduct += a.matrix[i][s]*b.matrix[s][j];
            }
        }
    }
    return scalarProduct;
}

double SquareMatrix::tr(){
    double trace=0;
    for (int i=0; i<m; ++i) {
        trace+=matrix[i][i];
    }
    
    return trace;
}

int Matrix::rank(){
    const double EPS = 1E-9;

    int rank = fmax(n,m);
    std::vector<char> line_used (n);
    for (int i=0; i<m; ++i) {
        int j;
        for (j=0; j<n; ++j)
            if (!line_used[j] && fabs(matrix[j][i]) > EPS)
                break;
        if (j == n)
            --rank;
        else {
            line_used[j] = true;
            for (int p=i+1; p<m; ++p)
                matrix[j][p] /= matrix[j][i];
            for (int k=0; k<n; ++k)
                if (k != j && abs (matrix[k][i]) > EPS)
                    for (int p=i+1; p<m; ++p)
                        matrix[k][p] -= matrix[j][p] * matrix[k][i];
        }
    }
    return rank;
}

double Vector::euclideanNorm(){
    double euclideanNorm = 0;
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            euclideanNorm += matrix[i][j]*matrix[i][j];
        }
    }
    euclideanNorm = std::sqrt(euclideanNorm);
    return euclideanNorm;
}

double Vector::maxNorm(){
    double maxNorm = 0;
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if(fabs(matrix[i][j])>maxNorm){
                maxNorm = fabs(matrix[i][j]);
            }
        }
    }
    return maxNorm;
}

double Matrix::frobenNorm(){
    double sum=0;
    for(int i=0;i<n;i++){
        for(int j=0;j<m;j++){
            sum+=matrix[i][j]*matrix[i][j];
        }
    }
    sum = std::sqrt(sum);
    return sum;
}

double AngleBetweenVectors(Vector vector1, Vector vector2){
    double angle;
    angle = ScalarProduct(vector1, vector2)/(vector1.euclideanNorm()*vector2.euclideanNorm());
    angle = acos(angle);
    return angle;
}

SquareMatrix SquareMatrix::inverse(){
    double temp;
    
    double res = this->det();
    
    assert(res!=0);
     
    SquareMatrix E(static_cast<int>(matrix.size()));
     
     
        for (int i = 0; i < matrix.size(); i++)
            for (int j = 0; j < matrix.size(); j++)
            {
                E.matrix[i][j] = 0.0;
     
                if (i == j)
                    E.matrix[i][j] = 1.0;
            }
     
        for (int k = 0; k < matrix.size(); k++)
        {
            temp = matrix[k][k];
     
            for (int j = 0; j < matrix.size(); j++)
            {
                matrix[k][j] /= temp;
                E.matrix[k][j] /= temp;
            }
     
            for (int i = k + 1; i < matrix.size(); i++)
            {
                temp = matrix[i][k];
     
                for (int j = 0; j < matrix.size(); j++)
                {
                    matrix[i][j] -= matrix[k][j] * temp;
                    E.matrix[i][j] -= E.matrix[k][j] * temp;
                }
            }
        }
     
        for (int k = static_cast<int>(matrix.size()) - 1; k > 0; k--)
        {
            for (int i = k - 1; i >= 0; i--)
            {
                temp = matrix[i][k];
     
                for (int j = 0; j < matrix.size(); j++)
                {
                    matrix[i][j] -= matrix[k][j] * temp;
                    E.matrix[i][j] -= E.matrix[k][j] * temp;
                }
            }
        }
     
        for (int i = 0; i < matrix.size(); i++)
            for (int j = 0; j < matrix.size(); j++)
                matrix[i][j] = E.matrix[i][j];
    
    return E;
}

void Matrix::toBinary(std::string bin){
    std::ofstream ss(bin+".bin");
    
    ss << (*this);
}
