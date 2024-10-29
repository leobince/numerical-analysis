#include <iostream>
#include <vector>
#include "Polynomial.h"
#include "InterpCondition.h"
#include "Interp.h"

int main() {
    // Given data
    std::vector<double> day = {0, 6, 10, 13, 17, 20, 28};
    std::vector<double> sp1 = {6.67, 17.3, 42.7, 37.3, 30.1, 29.3, 28.7};
    std::vector<double> sp2 = {6.67, 16.1, 18.9, 15.0, 10.6, 9.44, 8.89};

    // Create interpolation conditions
    InterpCondition cond1(day, sp1);
    InterpCondition cond2(day, sp2);

    // Perform interpolation
    Interp interp1(cond1);
    Interp interp2(cond2);

    // Predict weights after another 15 days (at day = 28 + 15 = 43)
    double day_future = 43;
    double weight1_future = interp1.evaluate(day_future);
    double weight2_future = interp2.evaluate(day_future);

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