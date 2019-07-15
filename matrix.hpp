#ifndef MATRIX_HPP_INCLUDED
#define MATRIX_HPP_INCLUDED

#include <iostream>
#include <stdexcept>
#include <algorithm>

void index(int row, int column, int h, int w) { //check the index
    if (row < 0 || column < 0 || row >= h || column >= w)
        throw std::out_of_range("Indexing fault");
}

class mtx {
    int h;
    int w;
    double *data;
    mtx Gauss_e(double&, mtx&);                //Gauss-elimination
public:
    int r;
    mtx(): h(0), w(0), data(NULL), r(0) {}

    mtx(int h, int w): h(h), w(w), data(new double[h*w]), r(1) {
        for(iterator it = begin(); it != end(); ++it)
            *it = 0;
    }

    mtx(mtx const& m): h(m.h), w(m.w), data(new double[h*w]), r(m.r) {
        iterator it1, it2;
        for(it1 = begin(), it2 = m.begin(); it1 != end(); ++it1, ++it2)
            *it1 = *it2;
    }

    mtx& operator=(mtx const&);

    double& operator()(int row, int column) {  //indexing
        index(row, column, h, w);
        return data[row * w + column];
    }

    double operator()(int row, int column) const {  //indexing
        index(row, column, h, w);
        return data[row * w + column];
    }

    mtx operator*(double) const;              //scalar multiplication

    mtx transp() const;                    //transposition

    mtx operator+(mtx const&) const;       //addition

    mtx operator-() const {                   //negative
        return (*this)*(-1);
    }

    mtx operator-(mtx& rhs) const {        //subtraction
        return ((*this) + (-rhs));
    }

    mtx operator*(mtx const&) const;       //matrix multiplication

    mtx& operator*=(double scalar) {
        return (*this) = ((*this) * scalar);
    }

    mtx& operator+=(mtx& rhs) {
        return (*this) = ((*this) + rhs);
    }

    mtx& operator-=(mtx& rhs) {
        return (*this) = ((*this) - rhs);
    }

    mtx& operator*=(mtx& rhs) {
        if(h != w || rhs.h != rhs.w || h != rhs.w)
            throw std::invalid_argument("The matrices are not square matrices!");
        return (*this) = ((*this) * rhs);
    }

    double det();   //determinant

    mtx inverse();   //inverse

    void Gauss();

    int mrank();

    int geth() const {  //getter
        return h;
    }

    int getw() const {  //getter
        return w;
    }

    class iterator {    //iterator class
        double * p;
    public:
        iterator(double * p = NULL)
            : p(p) {
        }

        iterator& operator++() {
            if(p != NULL)
                ++p;
            return *this;
        }
        iterator operator++(int) {
            iterator masolat = *this;
            ++(*this);
            return masolat;
        }

        double& operator*() const {
            return *p;
        }

        double* operator->() const {
            return p;
        }

        bool operator==(iterator rhs) const {
            return p == rhs.p;
        }
        bool operator!=(iterator rhs) const {
            return p != rhs.p;
        }
    };

    iterator begin() const {
        return iterator(data);
    }

    iterator end() const {
        return iterator(data + h * w);
    }

    ~mtx() {
        delete[] data;
    }
};

mtx operator*(double scalar, mtx& rhs) {
    return (rhs*scalar);
}

std::ostream& operator<<(std::ostream& os, mtx const& m)  {
    typename mtx::iterator it;
    int i = 0;
    for(it = m.begin(); it != m.end(); ++it) {
        if (i == 0)
            os << std::endl;
        os << *it;
        if (i == (m.getw() - 1)) {
            i = 0;
        } else {
            os << '\t';
            ++i;
        }
    }
    return os;
}

std::istream& operator>>(std::istream& is, mtx& m)  {
    for (int i = 0; i < m.geth(); ++i)
        for (int j = 0; j < m.getw(); ++j)
            is >> m(i, j);
    m.r = m.mrank();
    return is;
}

mtx& mtx::operator=(mtx const& m) {
    if (this != &m){
        delete[] data;
        h = m.h;
        w = m.w;
        data = new double[h*w];
        iterator it1, it2;
        for(it1 = begin(), it2 = m.begin(); it1 != end(); ++it1, ++it2)
            *it1 = *it2;
    }
    return *this;
}

mtx mtx::operator*(double scalar) const {
    mtx m(h, w);
    iterator it1, it2;
    for(it1 = (*this).begin(), it2 = m.begin(); it1 != end(); ++it1, ++it2)
        *it2 = *it1 * scalar;
    return m;
}

mtx mtx::transp() const {
    mtx m(w, h);
    for (int i = 0; i < w; ++i)
        for (int j = 0; j < h; ++j)
            m(i, j) = (*this)(j, i);
    return m;
}

mtx mtx::operator+(mtx const& rhs) const {
    if (h != rhs.h || w != rhs.w)
        throw std::invalid_argument("The matrices do not have the same number of rows and columns!");
    mtx m(h, w);
    iterator it1, it2, it3;
    for(it1 = begin(), it2 = rhs.begin(), it3 = m.begin(); it1 != end(); ++it1, ++it2, ++it3)
        *it3 = *it1 + *it2;
    return m;
}

mtx mtx::operator*(mtx const& rhs) const {
    if(w != rhs.h)
        throw std::invalid_argument("The number of columns of the left matrix is not the same as the number of rows of the right matrix!");
    mtx m(h, rhs.w);
    for(int i = 0; i < m.h; ++i)
        for(int j = 0; j < m.w; ++j)
            for(int k = 0; k < w; ++k)
                m(i, j) += (*this)(i, k) * rhs(k, j);
    return m;
}

mtx mtx::Gauss_e(double& D, mtx& e) {
    bool inv = false;
    if (h == w)
        inv = true;
    mtx m = (*this);
    if (inv)
        for (int i = 0; i < e.h; ++i)
            (e(i, i)) = 1;
    int i = 0;
    int j = 0;
    bool exit = false;
    int step = 1;
    double tmp = 0;
    while (!exit) {
        switch (step) {
            case 1:
                if (m(i, j) == 0){
                    step = 2;
                    break;
                }
                if (i == j)
                    D *= m(i,j);
                tmp = m(i, j);
                for (int t = 0; t < w; ++t) {
                    m(i, t) /= tmp;
                    if(inv)
                        e(i, t) /= tmp;
                }
                if (i == h-1) {
                    exit = true;
                    break;
                }
                for (int t = i + 1; t < h; ++t) {
                    tmp = m(t, j);
                    for (int k = 0; k < w; ++k){
                        m(t, k) -= (m(i, k)*tmp);
                        if (inv)
                            e(t, k) -= (m(i, k)*tmp);
                    }
                }
                if (j == w - 1) {
                    exit = true;
                    break;
                }
                ++i;
                ++j;
                break;
            case 2:
                if (i < h-1) {
                    int t;
                    for (t = i + 1; t < h-1 && m(t, j) == 0; ++t)
                        ;
                    if(m(t, j) != 0) {
                        for (int k = 0; k < w; ++k) {
                            std::swap(m(i, k), m(t, k));
                            if (inv)
                                std::swap(e(i, k), e(t, k));
                        }
                        D *= -1;
                        step = 1;
                        break;
                    }

                }
                D = 0;
                if (j == w-1) {
                    --i;
                    exit = true;
                    break;
                }
                ++j;
                step = 1;
                break;
            default:
                exit = true;
                break;
        }
    }
    return m;
}

double mtx::det() {
    if (h != w)
        throw std::invalid_argument("The matrix is not a square matrix!");
    double det = 1;
    mtx e(h, w);
    Gauss_e(det, e);
    return det;
}

mtx mtx::inverse() {
    double det = 1;
    mtx e(h, w);
    Gauss_e(det, e);
    if (det == 0)
        throw std::invalid_argument("The determinant of the matrix is 0, the matrix does not have an inverse!");
    return e;
}

void mtx::Gauss() {
    double det = 1;
    mtx e(h, w);
    *this = Gauss_e(det, e);
    bool zero = true;
    for(int i = h-1; i >= 0; --i) {
        for(int j = 0; j < w; ++j)
            if ((*this)(i, j) != 0)
                zero = false;
        if (zero) {
             mtx m(h - 1, w);
             for (int k = 0; k < h - 1; ++k)
                for (int l = 0; l < w; ++l)
                    m(k, l) = (*this)(k, l);
             *this = m;
            }
        zero = true;
    }
}

int mtx::mrank() {
    mtx m(*this);
    m.Gauss();
    return m.h;
}

#endif // MATRIX_HPP_INCLUDED
