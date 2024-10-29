#ifndef INTERP_H
#define INTERP_H

#include "InterpCondition.h"
#include "Polynomial.h"

class Interp {
public:
    Interp(const InterpCondition& con);

    // Methods
    void computeDividedDifferences(); // Compute the table of divided differences
    Polynomial getInterpolationPolynomial(); // Return the interpolation polynomial
    double evaluate(double x); // Evaluate interpolation polynomial at x

private:
    InterpCondition condition_;
    Polynomial interPoly_;
    std::vector<std::vector<double>> tableOfDividedDiffs;
    std::vector<double> coefficients; // Divided difference coefficients
};

#endif // INTERP_H