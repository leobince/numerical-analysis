# Compiler
CXX = g++
CXXFLAGS = -std=c++11 -Wall

# Targets
TARGET = solver
SOURCES = main.cpp
HEADERS = Function.hpp EquationSolver.hpp

# Report
REPORT_TEX = report.tex
REPORT_PDF = report.pdf

.PHONY: all run report clean

all: $(TARGET)

$(TARGET): $(SOURCES) $(HEADERS)
	$(CXX) $(CXXFLAGS) $(SOURCES) -o $(TARGET) -lm

run: $(TARGET)
	./$(TARGET) > output.txt

report: run $(REPORT_TEX)
	pdflatex $(REPORT_TEX)
	pdflatex $(REPORT_TEX)

clean:
	rm -f $(TARGET) *.o *.aux *.log *.out *.synctex.gz *.toc *.bbl *.blg $(REPORT_PDF) output.txt