#include <iostream>
#include <vector>
#include <cmath>

// Function to compute Hermite coefficients
void computeHermiteCoefficients(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& y_prime, std::vector<double>& coef) {
    int n = x.size();
    std::vector<double> z(2 * n);
    std::vector<double> Q(2 * n);

    // Construct z and Q
    for (int i = 0; i < n; ++i) {
        z[2 * i] = z[2 * i + 1] = x[i];
        Q[2 * i] = y[i];
        Q[2 * i + 1] = y[i];
    }

    // Compute divided differences
    coef.resize(2 * n);
    coef[0] = Q[0];
    std::vector<double> divided_diff(2 * n);
    divided_diff[0] = Q[0];

    for (int i = 1; i < 2 * n; ++i) {
        if (i % 2 == 1 && z[i] == z[i - 1]) {
            divided_diff[i] = y_prime[i / 2];
        } else {
            divided_diff[i] = (Q[i] - Q[i - 1]) / (z[i] - z[i - 1]);
        }
    }

    coef[1] = divided_diff[1];

    for (int i = 2; i < 2 * n; ++i) {
        double numerator = divided_diff[i] - divided_diff[i - 1];
        double denominator = z[i] - z[i - 2];
        divided_diff[i] = numerator / denominator;
        coef[i] = divided_diff[i];
    }
}

// Function to evaluate Hermite polynomial
double evaluateHermitePolynomial(const std::vector<double>& x, const std::vector<double>& coef, double x_val) {
    int n = coef.size();
    double result = coef[n - 1];
    for (int i = n - 2; i >= 0; --i) {
        result = result * (x_val - x[i]) + coef[i];
    }
    return result;
}

// Function to evaluate derivative of Hermite polynomial
double evaluateHermitePolynomialDerivative(const std::vector<double>& x, const std::vector<double>& coef, double x_val) {
    int n = coef.size();
    double result = 0.0;
    for (int i = n - 1; i > 0; --i) {
        double term = coef[i];
        for (int j = i - 1; j >= 0; --j) {
            term *= (x_val - x[j]);
        }
        result += term;
    }
    return result;
}

int main() {
    // Given data
    std::vector<double> t = {0, 3, 5, 8, 13};
    std::vector<double> s = {0, 225, 383, 623, 993};
    std::vector<double> v = {75, 77, 80, 74, 72}; // s'(t)

    int n = t.size();

    // Prepare data for Hermite interpolation
    std::vector<double> z(2 * n);
    std::vector<double> Q(2 * n);
    std::vector<double> s_prime(2 * n);

    for (int i = 0; i < n; ++i) {
        z[2 * i] = z[2 * i + 1] = t[i];
        Q[2 * i] = Q[2 * i + 1] = s[i];
        s_prime[2 * i] = v[i];
        s_prime[2 * i + 1] = v[i];
    }

    // Compute coefficients
    std::vector<double> coef;
    computeHermiteCoefficients(z, Q, s_prime, coef);

    // Predict position and speed at t = 10
    double t_val = 10.0;
    double s_val = evaluateHermitePolynomial(z, coef, t_val);
    double v_val = evaluateHermitePolynomialDerivative(z, coef, t_val);

    std::cout << "At t = " << t_val << " seconds:" << std::endl;
    std::cout << "Predicted position s(t) = " << s_val << " feet" << std::endl;
    std::cout << "Predicted speed v(t) = " << v_val << " feet/second" << std::endl;

    // Check if car ever exceeds 81 feet per second
    int num_points = 1000;
    double t_start = t.front();
    double t_end = t.back();
    double dt = (t_end - t_start) / (num_points - 1);
    bool exceeds_speed_limit = false;

    for (int i = 0; i < num_points; ++i) {
        double t_curr = t_start + i * dt;
        double v_curr = evaluateHermitePolynomialDerivative(z, coef, t_curr);
        if (v_curr > 81.0) {
            exceeds_speed_limit = true;
            break;
        }
    }

    if (exceeds_speed_limit) {
        std::cout << "The car exceeds the speed limit of 81 feet/second." << std::endl;
    } else {
        std::cout << "The car does not exceed the speed limit of 81 feet/second." << std::endl;
    }

    return 0;
}