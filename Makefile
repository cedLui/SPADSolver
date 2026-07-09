CXX      := g++
CXXFLAGS := -std=c++17 -Wall -Wextra -I. -Itests
BUILD    := build

INTERP_SRC     := interpolation.cpp
INTERP_TEST    := tests/test_interpolation.cpp
GAUSS_SRC      := GaussianElimination.cpp
GAUSS_TEST     := tests/test_gaussian_elimination.cpp

INTERP_BIN     := $(BUILD)/test_interpolation
GAUSS_BIN      := $(BUILD)/test_gaussian_elimination
PLOT_BIN       := $(BUILD)/interpolation_plot
TEST_BINS      := $(INTERP_BIN) $(GAUSS_BIN)

.PHONY: all test plot clean

all: $(TEST_BINS)

plot: $(PLOT_BIN)

test: all
	@echo "Running interpolation tests..."
	@$(INTERP_BIN)
	@echo "Running Gaussian elimination tests..."
	@$(GAUSS_BIN)

$(BUILD):
	mkdir -p $(BUILD)

$(INTERP_BIN): $(INTERP_SRC) $(INTERP_TEST) | $(BUILD)
	$(CXX) $(CXXFLAGS) -o $@ $(INTERP_SRC) $(INTERP_TEST)

$(GAUSS_BIN): $(GAUSS_SRC) $(GAUSS_TEST) | $(BUILD)
	$(CXX) $(CXXFLAGS) -o $@ $(GAUSS_SRC) $(GAUSS_TEST)

$(PLOT_BIN): interpolation_plot.cpp $(INTERP_SRC) | $(BUILD)
	$(CXX) $(CXXFLAGS) -o $@ interpolation_plot.cpp $(INTERP_SRC)

clean:
	rm -rf $(BUILD)
