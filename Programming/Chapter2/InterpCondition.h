#ifndef INTERPCONDITION_H
#define INTERPCONDITION_H

#include <vector>

class InterpCondition {
public:
    // Constructors
    InterpCondition();
    InterpCondition(const std::vector<std::pair<double, double>>& points);
    InterpCondition(const std::vector<double>& x_vals, const std::vector<double>& y_vals);
    InterpCondition(const std::vector<double>& x_vals, const std::vector<double>& y_vals, const std::vector<double>& y_derivatives);

    // Methods to add conditions
    void addPoint(double x, double y);
    void addDerivative(double x, double y_prime);

    int getOrder() const;

    // Accessors
    const std::vector<double>& getX() const;
    const std::vector<double>& getY() const;
    const std::vector<double>& getYPrime() const; // For derivatives (Hermite interpolation)

private:
    int max_order_; // The maximum order of derivative information available
    std::vector<double> x_;
    std::vector<double> y_;
    std::vector<double> y_prime_; // First derivatives, used for Hermite interpolation
};

#endif // INTERPCONDITION_H