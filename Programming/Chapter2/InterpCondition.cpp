#include "InterpCondition.h"

InterpCondition::InterpCondition() : max_order_(0) {}

InterpCondition::InterpCondition(const std::vector<std::pair<double, double>>& points) : max_order_(0) {
    for (const auto& p : points) {
        x_.push_back(p.first);
        y_.push_back(p.second);
    }
}

InterpCondition::InterpCondition(const std::vector<double>& x_vals, const std::vector<double>& y_vals)
    : max_order_(0), x_(x_vals), y_(y_vals) {}

InterpCondition::InterpCondition(const std::vector<double>& x_vals, const std::vector<double>& y_vals, const std::vector<double>& y_derivatives)
    : max_order_(1), x_(x_vals), y_(y_vals), y_prime_(y_derivatives) {}

void InterpCondition::addPoint(double x, double y) {
    x_.push_back(x);
    y_.push_back(y);
}

void InterpCondition::addDerivative(double x, double y_prime) {
    x_.push_back(x); // For Hermite, x values may appear multiple times
    y_.push_back(y_[y_.size() - 1]); // Duplicate y value
    y_prime_.push_back(y_prime);
    max_order_ = 1;
}

int InterpCondition::getOrder() const {
    return max_order_;
}

const std::vector<double>& InterpCondition::getX() const {
    return x_;
}

const std::vector<double>& InterpCondition::getY() const {
    return y_;
}

const std::vector<double>& InterpCondition::getYPrime() const {
    return y_prime_;
}