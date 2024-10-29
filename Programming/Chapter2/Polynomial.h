#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

#include <vector>
#include <initializer_list>
#include <iostream>

class Polynomial {
public:
    // Constructors
    Polynomial();
    Polynomial(const std::vector<double>& coef);
    Polynomial(std::initializer_list<double> coef);

    // Methods
    Polynomial diff() const;
    void set_coef(const std::vector<double>& coef);
    std::vector<double> get_coef() const;
    double evaluate(double x) const;

    // Overloaded operators
    Polynomial operator+(const Polynomial& P) const;
    Polynomial operator-(const Polynomial& P) const;
    Polynomial operator*(const Polynomial& P) const;
    Polynomial operator*(double scalar) const;
    Polynomial& operator+=(const Polynomial& P);
    Polynomial& operator-=(const Polynomial& P);
    Polynomial& operator*=(const Polynomial& P);

    // Friend functions
    friend std::ostream& operator<<(std::ostream& os, const Polynomial& P);

private:
    // Coefficients in increasing order of degree
    std::vector<double> coefficient;
};

#endif // POLYNOMIAL_H