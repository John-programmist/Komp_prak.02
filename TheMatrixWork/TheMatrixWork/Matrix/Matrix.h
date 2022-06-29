
#ifndef Matrix_h
#define Matrix_h
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <map>
#include <cmath>
#include <ctime>
#include <string>
#include <cstdlib>
#include <cassert>
#include <initializer_list>


enum vector{
    VERTICAL,
    HORIZONTAL
};


class Matrix{
public:
    
    int m; //Строки
    int n; //Столбцы
    std::vector<std::vector<double>> matrix;
    Matrix(int m=0, int n=0):n(n), m(m){//Конструктор по умолчанию
        matrix.resize(m);
        for(int i = 0;i<m;++i){
            matrix[i].resize(n);
            for(int j = 0;j<n;++j){
                matrix[i][j] = 0;
            }
        }
    }
    
    
    
    ~Matrix() = default;//Деструктор
    
    virtual int GetN(){
        return n;
    }
    
    virtual int GetM(){
        return m;
    }
    
    virtual std::vector<std::vector<double>> GetMatrix(){
        return matrix;
    }
    
    virtual int GetMatrixElement(int i=0, int j=0){
        return matrix[i][j];
    }
    
    virtual void toBinary(std::string name);
    
    Matrix(Matrix& A):n(A.n),m(A.m){//Конструктор копирования
        matrix.resize(m);
        for(int i = 0;i<m;++i){
            matrix[i].resize(n);
            for(int j = 0;j<n;++j){
                matrix[i][j] = A.matrix[i][j];
            }
        }
    }
    
    Matrix(vector vector, int N=0){//Вектор
        if(vector){
            n = N;
            m = 1;
            matrix.resize(1);
            matrix[0].resize(N);
            for(int j = 0;j<N;++j){
                matrix[0][j] = 0;
            }
        }
        else{
            n = 1;
            m = N;
            matrix.resize(N);
            for (int i = 0; i < N; ++i) {
                matrix[i].resize(1);
                matrix[i][0] = 0;
            }
        }
        matrix[0][0] = 1;
    }
    
    
    Matrix(std::initializer_list<double> l){
        n = static_cast<int>(l.size());
        m = 1;
        auto it = l.begin();
        matrix.resize(1);
        matrix[0].resize(n);
        for(int j = 0;j<n;++j){
            matrix[0][j] = *it;
            it++;
        }
    }
    
    friend std::ostream& operator<<(std::ostream& out, const Matrix& matrix);//Вывод
    friend std::istream& operator>>(std::istream& in, Matrix& matrix);//Ввод
    friend std::ofstream& operator<<(std::ofstream& out, Matrix const& matrix);
    
    
    friend Matrix operator+(const Matrix A, const Matrix B);//Сложение матриц
    friend Matrix operator-(const Matrix A, const Matrix B);//Вычитание матриц
    friend Matrix operator*(const Matrix A, const Matrix B);//Умноженеи матриц
    friend Matrix operator%(const Matrix A, const Matrix B);//Произведение Адамара
    friend Matrix& operator*(Matrix& A, const double alpha);//Умножение матриц на число
    friend Matrix& operator*(const double alpha, Matrix& A);//Умножение число на матрицу
    
    Matrix& operator=(Matrix& A){//Оператор присваивания
        n = A.n;
        m = A.m;
        matrix.resize(m);
        for(int i = 0;i<m;++i){
            matrix[i].resize(n);
            for(int j = 0;j<n;++j){
                matrix[i][j] = A.matrix[i][j];
            }
        }
        return *this;
    }
    
    
    std::vector<double>& operator[](int i){//Оператор индексации
        return matrix[i];
    }
    
    Matrix trans(){
        std::swap(m,n);
        matrix.resize(m);
        for(int i = 0;i<m;++i){
            matrix[i].resize(n);
            for(int j = 0;j<i;++j){
                std::swap(matrix[i][j],matrix[j][i]);
            }
        }
        return *this;
    }
    
    virtual int rank();
    
    double frobenNorm();
    
    
    

};


class PCA{
public:
    double centr(Matrix matrix, int j);
    double shkala(Matrix matrix, int i, int j);
};



class Vector: public Matrix{
public:
    Vector(vector vector, int N=0):Matrix(vector,N){}
    ~Vector() = default;
    
    Vector(std::initializer_list<double> l):Matrix(l){}
    
    
    double euclideanNorm();
    double maxNorm();
    

    
    
    friend double operator*(const Vector a, const Vector b);
    
};


class SquareMatrix: public Matrix{
public:
    SquareMatrix(int m=0, int n=0){
        this->m = m;
        this->n = m;
        matrix.resize(m);
        for(int i = 0;i<m;++i){
            matrix[i].resize(m);
            for(int j = 0;j<m;++j){
                matrix[i][j] = 0;
            }
        }
    }
    ~SquareMatrix() = default;
    
    
    virtual SquareMatrix inverse();
    virtual double det();
    virtual double tr();
    
};

class Unitmatrix: public SquareMatrix{
public:
    Unitmatrix(int m=0):SquareMatrix(m){
        matrix.resize(m);
        srand(static_cast<unsigned int>(time(nullptr)));
        for(int i = 0;i<m;++i){
            matrix[i].resize(m);
            for(int j = 0;j<m;++j){
                matrix[i][j] = (i==j);
            }
        }
    }
    ~Unitmatrix() = default;
    
};

class DiagonalMatrix: public SquareMatrix{
    DiagonalMatrix(int m=0):SquareMatrix(m){
        matrix.resize(m);
        srand(static_cast<unsigned int>(time(nullptr)));
        for(int i = 0;i<m;++i){
            matrix[i].resize(m);
            for(int j = 0;j<m;++j){
                matrix[i][j] = (i==j)*std::rand()%10;
            }
        }
    }
    ~DiagonalMatrix() = default;
    
};

class UpperTriangularMatrix: public SquareMatrix{
public:
    UpperTriangularMatrix(int m=0):SquareMatrix(m){
        matrix.resize(m);
        srand(static_cast<unsigned int>(time(nullptr)));
        for(int i = 0;i<m;++i){
            matrix[i].resize(m);
            for(int j = 0;j<m;++j){
                matrix[i][j] = (i<=j)*(std::rand()%10+1);
            }
        }
    }
    ~UpperTriangularMatrix() = default;
    
};

class LowerTriangularMatrix: public SquareMatrix{
public:
    LowerTriangularMatrix(int m=0):SquareMatrix(m){
        matrix.resize(m);
        srand(static_cast<unsigned int>(time(nullptr)));
        for(int i = 0;i<m;++i){
            matrix[i].resize(m);
            for(int j = 0;j<m;++j){
                matrix[i][j] = (i>=j)*(std::rand()%10+1);
            }
        }
    }
    ~LowerTriangularMatrix() = default;
    
    
};

class SymmetricMatrix: public SquareMatrix{
public:
    SymmetricMatrix(int m=0):SquareMatrix(m){
        matrix.resize(m);
        srand(static_cast<unsigned int>(time(nullptr))); 
        for(int i = 0;i<m;++i){
            matrix[i].resize(m);
            for(int j = 0;j<m;++j){
                matrix[i][j] = (i>=j)*(std::rand()%10);
                matrix[j][i] = matrix[i][j];
            }
        }
    }
    ~SymmetricMatrix() = default;
    
};

double ScalarProduct(Matrix vector1, Matrix vector2);
double AngleBetweenVectors(Vector vector1, Vector vector2);


#endif
