#include "Polynomial.h"
#include <cmath>

// Constructors
Polynomial::Polynomial() {}

Polynomial::Polynomial(const std::vector<double>& coef) : coefficient(coef) {}

Polynomial::Polynomial(std::initializer_list<double> coef) : coefficient(coef) {}

// Methods
Polynomial Polynomial::diff() const {
    std::vector<double> deriv_coef;
    int n = coefficient.size();
    for (int i = 1; i < n; ++i) {
        deriv_coef.push_back(coefficient[i] * i);
    }
    return Polynomial(deriv_coef);
}

void Polynomial::set_coef(const std::vector<double>& coef) {
    coefficient = coef;
}

std::vector<double> Polynomial::get_coef() const {
    return coefficient;
}

double Polynomial::evaluate(double x) const {
    double result = 0.0;
    for (int i = coefficient.size() - 1; i >= 0; --i) {
        result = result * x + coefficient[i];
    }
    return result;
}

// Overloaded operators
Polynomial Polynomial::operator+(const Polynomial& P) const {
    size_t n = std::max(coefficient.size(), P.coefficient.size());
    std::vector<double> result_coef(n, 0.0);
    for (size_t i = 0; i < coefficient.size(); ++i) {
        result_coef[i] += coefficient[i];
    }
    for (size_t i = 0; i < P.coefficient.size(); ++i) {
        result_coef[i] += P.coefficient[i];
    }
    return Polynomial(result_coef);
}

Polynomial Polynomial::operator-(const Polynomial& P) const {
    size_t n = std::max(coefficient.size(), P.coefficient.size());
    std::vector<double> result_coef(n, 0.0);
    for (size_t i = 0; i < coefficient.size(); ++i) {
        result_coef[i] += coefficient[i];
    }
    for (size_t i = 0; i < P.coefficient.size(); ++i) {
        result_coef[i] -= P.coefficient[i];
    }
    return Polynomial(result_coef);
}

Polynomial Polynomial::operator*(const Polynomial& P) const {
    size_t n = coefficient.size();
    size_t m = P.coefficient.size();
    std::vector<double> result_coef(n + m - 1, 0.0);
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < m; ++j) {
            result_coef[i + j] += coefficient[i] * P.coefficient[j];
        }
    }
    return Polynomial(result_coef);
}

Polynomial Polynomial::operator*(double scalar) const {
    std::vector<double> result_coef = coefficient;
    for (double& c : result_coef) {
        c *= scalar;
    }
    return Polynomial(result_coef);
}

Polynomial& Polynomial::operator+=(const Polynomial& P) {
    *this = *this + P;
    return *this;
}

Polynomial& Polynomial::operator-=(const Polynomial& P) {
    *this = *this - P;
    return *this;
}

Polynomial& Polynomial::operator*=(const Polynomial& P) {
    *this = *this * P;
    return *this;
}

// Friend functions
std::ostream& operator<<(std::ostream& os, const Polynomial& P) {
    const auto& coef = P.get_coef();
    bool first_term = true;
    for (int i = coef.size() - 1; i >= 0; --i) {
        if (coef[i] != 0.0) {
            if (!first_term) {
                os << (coef[i] > 0 ? " + " : " - ");
            } else if (coef[i] < 0) {
                os << "-";
            }
            first_term = false;
            double abs_coef = std::abs(coef[i]);
            if (abs_coef != 1.0 || i == 0) {
                os << abs_coef;
            }
            if (i > 0) {
                os << "x";
                if (i > 1) {
                    os << "^" << i;
                }
            }
        }
    }
    if (first_term) {
        os << "0";
    }
    return os;
}