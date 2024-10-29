#include <iostream>
#include <vector>

// Function to compute divided differences coefficients
void computeDividedDifferences(const std::vector<double>& x, const std::vector<double>& y, std::vector<double>& coef) {
    int n = x.size();
    coef = y; // Initialize coefficients with y values
    for (int j = 1; j < n; ++j) {
        for (int i = n - 1; i >= j; --i) {
            coef[i] = (coef[i] - coef[i - 1]) / (x[i] - x[i - j]);
        }
    }
}

// Function to evaluate Newton's interpolating polynomial at x_val
double evaluateNewtonPolynomial(const std::vector<double>& x, const std::vector<double>& coef, double x_val) {
    int n = x.size();
    double result = coef[n - 1];
    for (int i = n - 2; i >= 0; --i) {
        result = result * (x_val - x[i]) + coef[i];
    }
    return result;
}

int main() {
    // Given data
    std::vector<double> day = {0, 6, 10, 13, 17, 20, 28};
    std::vector<double> sp1 = {6.67, 17.3, 42.7, 37.3, 30.1, 29.3, 28.7};
    std::vector<double> sp2 = {6.67, 16.1, 18.9, 15.0, 10.6, 9.44, 8.89};

    // Compute coefficients for Sample 1
    std::vector<double> coef1;
    computeDividedDifferences(day, sp1, coef1);

    // Compute coefficients for Sample 2
    std::vector<double> coef2;
    computeDividedDifferences(day, sp2, coef2);

    // Predict weights after another 15 days (at day = 28 + 15 = 43)
    double day_future = 43;

    double weight1_future = evaluateNewtonPolynomial(day, coef1, day_future);
    double weight2_future = evaluateNewtonPolynomial(day, coef2, day_future);

    std::cout << "Predicted weight for Sample 1 at day " << day_future << " is " << weight1_future << " mg" << std::endl;
    std::cout << "Predicted weight for Sample 2 at day " << day_future << " is " << weight2_future << " mg" << std::endl;

    // Determine if larvae will die (assuming they die if weight <= 0)
    if (weight1_future <= 0) {
        std::cout << "Sample 1 larvae are predicted to die after another 15 days." << std::endl;
    } else {
        std::cout << "Sample 1 larvae are predicted to survive after another 15 days." << std::endl;
    }

    if (weight2_future <= 0) {
        std::cout << "Sample 2 larvae are predicted to die after another 15 days." << std::endl;
    } else {
        std::cout << "Sample 2 larvae are predicted to survive after another 15 days." << std::endl;
    }

    return 0;
}