#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Function to generate marker points on the heart curve
void generateMarkerPoints(std::vector<double>& x_vals, std::vector<double>& y_vals,
                          std::vector<double>& t_vals, int m) {
    const double sqrt_3 = std::sqrt(3.0);
    double delta_t = (2 * M_PI) / m;

    for (int i = 0; i <= m; ++i) {
        double t = i * delta_t; // Parameter t
        t_vals.push_back(t);

        double x = sqrt_3 * std::cos(t);
        double y = (2.0 / 3.0) * (sqrt_3 * std::sin(t) + std::sqrt(sqrt_3 * std::fabs(std::cos(t))));

        x_vals.push_back(x);
        y_vals.push_back(y);
    }
}

// Function to compute the tangent vectors at each marker point
void computeTangents(const std::vector<double>& t_vals,
                     std::vector<double>& dx_vals, std::vector<double>& dy_vals, int m) {
    const double sqrt_3 = std::sqrt(3.0);
    double delta_t = (2 * M_PI) / m;

    for (size_t i = 0; i < t_vals.size(); ++i) {
        double t = t_vals[i];

        // Compute dx/dt
        double dx_dt = -sqrt_3 * std::sin(t);

        // Compute dy/dt
        double cos_t = std::cos(t);
        double sin_t = std::sin(t);
        double abs_cos_t = std::fabs(cos_t);

        double dy_dt;

        if (abs_cos_t > 1e-8) { // Avoid division by zero
            double sqrt_term = std::sqrt(sqrt_3 * abs_cos_t);
            double sqrt_derivative = - (sqrt_3 * sin_t * ((cos_t >= 0) ? 1.0 : -1.0)) / (2.0 * sqrt_term);
            dy_dt = (2.0 / 3.0) * (sqrt_3 * cos_t + sqrt_derivative);
        } else {
            // When cos_t is close to zero
            dy_dt = (2.0 / 3.0) * sqrt_3 * cos_t;
        }

        // Scale the tangent vector by delta_t
        dx_dt *= delta_t;
        dy_dt *= delta_t;

        dx_vals.push_back(dx_dt);
        dy_vals.push_back(dy_dt);
    }
}

// Function to compute the control points for cubic Bézier curves
void computeControlPoints(const std::vector<double>& x_vals, const std::vector<double>& y_vals,
                          const std::vector<double>& dx_vals, const std::vector<double>& dy_vals,
                          std::vector<std::vector<double>>& qx_vals,
                          std::vector<std::vector<double>>& qy_vals) {
    size_t m = x_vals.size() - 1; // Number of segments

    for (size_t j = 0; j < m; ++j) {
        double px0 = x_vals[j];
        double py0 = y_vals[j];

        // Scaled tangent vectors
        double T_x_j = dx_vals[j];
        double T_y_j = dy_vals[j];

        double T_x_j1 = dx_vals[j + 1];
        double T_y_j1 = dy_vals[j + 1];

        // Compute control points
        double px1 = px0 + T_x_j / 3.0;
        double py1 = py0 + T_y_j / 3.0;

        double px3 = x_vals[j + 1];
        double py3 = y_vals[j + 1];

        double px2 = px3 - T_x_j1 / 3.0;
        double py2 = py3 - T_y_j1 / 3.0;

        // Store the control points for this segment
        qx_vals.push_back({px0, px1, px2, px3});
        qy_vals.push_back({py0, py1, py2, py3});
    }
}

// Function to evaluate Bézier curve at parameter t
void evaluateBezierCurve(const std::vector<double>& qx, const std::vector<double>& qy,
                         std::vector<double>& x_curve, std::vector<double>& y_curve, int num_points) {
    for (int k = 0; k <= num_points; ++k) {
        double t = static_cast<double>(k) / num_points;
        double one_minus_t = 1.0 - t;

        double b0 = one_minus_t * one_minus_t * one_minus_t;
        double b1 = 3 * one_minus_t * one_minus_t * t;
        double b2 = 3 * one_minus_t * t * t;
        double b3 = t * t * t;

        double x = b0 * qx[0] + b1 * qx[1] + b2 * qx[2] + b3 * qx[3];
        double y = b0 * qy[0] + b1 * qy[1] + b2 * qy[2] + b3 * qy[3];

        x_curve.push_back(x);
        y_curve.push_back(y);
    }
}

// Function to write the approximated curve to file
void writeCurveToFile(const std::vector<double>& x_curve, const std::vector<double>& y_curve, const std::string& filename) {
    std::ofstream file(filename);
    for (size_t i = 0; i < x_curve.size(); ++i) {
        file << x_curve[i] << " " << y_curve[i] << std::endl;
    }
    file.close();
}

int main() {
    std::vector<int> m_values = {10, 40, 160};

    for (int m : m_values) {
        std::vector<double> x_vals;
        std::vector<double> y_vals;
        std::vector<double> t_vals;

        // Generate marker points on the heart curve
        generateMarkerPoints(x_vals, y_vals, t_vals, m);

        // Compute tangent vectors at marker points
        std::vector<double> dx_vals;
        std::vector<double> dy_vals;
        computeTangents(t_vals, dx_vals, dy_vals, m);

        // Compute control points for cubic Bézier curves
        std::vector<std::vector<double>> qx_vals;
        std::vector<std::vector<double>> qy_vals;
        computeControlPoints(x_vals, y_vals, dx_vals, dy_vals, qx_vals, qy_vals);

        // For each segment, evaluate the Bézier curve and collect points
        std::vector<double> x_curve;
        std::vector<double> y_curve;
        int num_points_per_segment = 10; // Adjust for smoother curves

        for (size_t j = 0; j < qx_vals.size(); ++j) {
            std::vector<double> x_segment;
            std::vector<double> y_segment;
            evaluateBezierCurve(qx_vals[j], qy_vals[j], x_segment, y_segment, num_points_per_segment);

            // Append segment points to the overall curve (except the last point to avoid duplicates)
            if (j < qx_vals.size() - 1) {
                x_curve.insert(x_curve.end(), x_segment.begin(), x_segment.end() - 1);
                y_curve.insert(y_curve.end(), y_segment.begin(), y_segment.end() - 1);
            } else {
                // For the last segment, include all points
                x_curve.insert(x_curve.end(), x_segment.begin(), x_segment.end());
                y_curve.insert(y_curve.end(), y_segment.begin(), y_segment.end());
            }
        }

        // Write the approximated curve to file
        std::string filename = "heart_bezier_m" + std::to_string(m) + ".txt";
        writeCurveToFile(x_curve, y_curve, filename);
    }

    std::cout << "Heart curve approximations generated using cubic Bézier curves for m = 10, 40, 160." << std::endl;
    std::cout << "Use a plotting tool to visualize the results." << std::endl;

    return 0;
}