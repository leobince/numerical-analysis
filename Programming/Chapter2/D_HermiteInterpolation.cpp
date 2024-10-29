#include <iostream>
#include <vector>
#include "Polynomial.h"
#include "InterpCondition.h"
#include "Interp.h"

int main() {
    // Given data
    std::vector<double> t = {0, 3, 5, 8, 13};
    std::vector<double> s = {0, 225, 383, 623, 993};
    std::vector<double> v = {75, 77, 80, 74, 72}; // s'(t)

    size_t n = t.size();

    // Duplicate data for Hermite interpolation
    std::vector<double> z(2 * n);
    std::vector<double> y(2 * n);
    std::vector<double> y_prime(2 * n);

    for (size_t i = 0; i < n; ++i) {
        z[2 * i] = z[2 * i + 1] = t[i];
        y[2 * i] = y[2 * i + 1] = s[i];
        y_prime[2 * i] = y_prime[2 * i + 1] = v[i];
    }

    // Create interpolation condition
    InterpCondition cond(z, y, y_prime);

    // Perform Hermite interpolation
    Interp interp(cond);

    // Predict position and speed at t = 10
    double t_val = 10.0;
    double s_val = interp.evaluate(t_val);

    // Since the interpolation polynomial is generated, we can compute its derivative
    // For simplicity, compute derivative numerically
    double h = 1e-5;
    double v_val = (interp.evaluate(t_val + h) - interp.evaluate(t_val - h)) / (2.0 * h);

    std::cout << "At t = " << t_val << " seconds:" << std::endl;
    std::cout << "Predicted position s(t) = " << s_val << " feet" << std::endl;
    std::cout << "Predicted speed v(t) = " << v_val << " feet/second" << std::endl;

    // Check if car ever exceeds 81 feet per second
    int num_points = 1000;
    double t_start = t.front();
    double t_end = t.back();
    double dt = (t_end - t_start) / (num_points - 1);
    bool exceeds_speed_limit = false;

    double max_speed = 0.0;

    for (int i = 0; i < num_points; ++i) {
        double t_curr = t_start + i * dt;
        double v_curr = (interp.evaluate(t_curr + h) - interp.evaluate(t_curr - h)) / (2.0 * h);
        if (v_curr > 81.0) {
            exceeds_speed_limit = true;
            break;
        }
        if (v_curr > max_speed) {
            max_speed = v_curr;
        }
    }

    if (exceeds_speed_limit) {
        std::cout << "The car exceeds the speed limit of 81 feet/second." << std::endl;
    } else {
        std::cout << "The car does not exceed the speed limit of 81 feet/second." << std::endl;
    }

    return 0;
}