// F_HeartBezier.cpp
#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Function to generate marker points on the heart curve
void generateMarkerPoints(std::vector<double>& x_vals, std::vector<double>& y_vals, std::vector<double>& t_vals, int m) {
    for (int i = 0; i <= m; ++i) {
        double t = (2 * M_PI * i) / m; // Parameter t
        t_vals.push_back(t);

        double x = pow(sin(t), 3);
        double y = (13 * cos(t) - 5 * cos(2 * t) - 2 * cos(3 * t) - cos(4 * t)) / 17.0;

        x_vals.push_back(x);
        y_vals.push_back(y);
    }
}

// Function to compute the tangent vectors at each marker point
void computeTangents(const std::vector<double>& x_vals, const std::vector<double>& y_vals,
                     std::vector<double>& dx_vals, std::vector<double>& dy_vals,
                     const std::vector<double>& t_vals) {
    size_t n = x_vals.size();

    for (size_t i = 0; i < n; ++i) {
        double t = t_vals[i];
        double dx_dt = 3 * pow(sin(t), 2) * cos(t);
        double dy_dt = (1.0 / 17.0) * (-13 * sin(t) + 10 * sin(2 * t) + 6 * sin(3 * t) + 4 * sin(4 * t));

        // Tangent vector components dx/dt and dy/dt
        double dx = dx_dt;
        double dy = dy_dt;

        // Compute dy/dx
        // Not needed since we're using parametric form, but we can normalize the vector
        double length = sqrt(dx * dx + dy * dy);
        dx_vals.push_back(dx / length);
        dy_vals.push_back(dy / length);
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

        double px1 = px0 + dx_vals[j] / 3.0;
        double py1 = py0 + dy_vals[j] / 3.0;

        double px3 = x_vals[j + 1];
        double py3 = y_vals[j + 1];

        double px2 = px3 - dx_vals[j + 1] / 3.0;
        double py2 = py3 - dy_vals[j + 1] / 3.0;

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
        computeTangents(x_vals, y_vals, dx_vals, dy_vals, t_vals);

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
