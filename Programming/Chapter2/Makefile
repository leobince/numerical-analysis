# Makefile for Compiling, Running, and Plotting Refactored Tasks A-F and Building the Report

# Compiler and flags
CXX = g++
CXXFLAGS = -std=c++11 -O2 -Wall

# Source files
SOURCES = Polynomial.cpp InterpCondition.cpp Interp.cpp

# Header files
HEADERS = Polynomial.h InterpCondition.h Interp.h

# Programs
PROGRAMS = A_NewtonInterpolation \
           B_RungeInterpolation \
           C_ChebyshevInterpolation \
           D_HermiteInterpolation \
           E_LarvaeGrowth \
           F_HeartBezier

# Gnuplot scripts
GNUPLOTS = plot_run.plt \
           plot_cheb.plt \
           plot_heart.plt

# Plots
PLOTS = runge_plot.png \
        chebyshev_plot.png \
        heart_plot_m10.png \
        heart_plot_m40.png \
        heart_plot_m160.png

# Report
REPORT = Report.pdf

.PHONY: all run plot report clean

# Default target
all: $(PROGRAMS)

# Compiling
A_NewtonInterpolation: A_NewtonInterpolation.cpp $(SOURCES) $(HEADERS)
	$(CXX) $(CXXFLAGS) -o $@ $< $(SOURCES)

B_RungeInterpolation: B_RungeInterpolation.cpp $(SOURCES) $(HEADERS)
	$(CXX) $(CXXFLAGS) -o $@ $< $(SOURCES)

C_ChebyshevInterpolation: C_ChebyshevInterpolation.cpp $(SOURCES) $(HEADERS)
	$(CXX) $(CXXFLAGS) -o $@ $< $(SOURCES) -lm

D_HermiteInterpolation: D_HermiteInterpolation.cpp $(SOURCES) $(HEADERS)
	$(CXX) $(CXXFLAGS) -o $@ $< $(SOURCES)

E_LarvaeGrowth: E_LarvaeGrowth.cpp $(SOURCES) $(HEADERS)
	$(CXX) $(CXXFLAGS) -o $@ $< $(SOURCES)

F_HeartBezier: F_HeartBezier.cpp
	$(CXX) $(CXXFLAGS) -o $@ $< -lm

# Running programs
run: all
	@echo "Running A_NewtonInterpolation..."
	./A_NewtonInterpolation
	@echo "Running B_RungeInterpolation..."
	./B_RungeInterpolation
	@echo "Running C_ChebyshevInterpolation..."
	./C_ChebyshevInterpolation
	@echo "Running D_HermiteInterpolation..."
	./D_HermiteInterpolation > hermite_output.txt
	@echo "Running E_LarvaeGrowth..."
	./E_LarvaeGrowth > larvae_output.txt
	@echo "Running F_HeartBezier..."
	./F_HeartBezier

# Plotting
plot: run $(PLOTS)

# Explicit plot rules
runge_plot.png: plot_run.plt runge_exact.txt runge_poly_n2.txt runge_poly_n4.txt runge_poly_n6.txt runge_poly_n8.txt
	@echo "Generating runge_plot.png..."
	gnuplot plot_run.plt

chebyshev_plot.png: plot_cheb.plt chebyshev_exact.txt chebyshev_poly_n5.txt chebyshev_poly_n10.txt chebyshev_poly_n15.txt chebyshev_poly_n20.txt
	@echo "Generating chebyshev_plot.png..."
	gnuplot plot_cheb.plt

heart_plot_m10.png heart_plot_m40.png heart_plot_m160.png: plot_heart.plt heart_bezier_m10.txt heart_bezier_m40.txt heart_bezier_m160.txt
	@echo "Generating heart plots..."
	gnuplot plot_heart.plt

# Build the report
report: plot $(REPORT)

$(REPORT): Report.tex $(PLOTS) hermite_output.txt larvae_output.txt
	pdflatex Report.tex
	pdflatex Report.tex

# Clean generated files
clean:
	rm -f $(PROGRAMS)
	rm -f hermite_output.txt
	rm -f larvae_output.txt
	rm -f runge_exact.txt runge_poly_n*.txt
	rm -f chebyshev_exact.txt chebyshev_poly_n*.txt
	rm -f heart_bezier_m*.txt
	rm -f $(PLOTS)
	rm -f *.o
	rm -f Report.aux Report.log Report.pdf
