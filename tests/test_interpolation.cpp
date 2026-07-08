#include <algorithm>
#include <cmath>
#include <iostream>

#include "interpolation.h"
#include "interpolation_reference.h"
#include "test_helpers.h"

namespace {

constexpr double k_field = 3.0e6;  // V/cm
constexpr double k_temp = 298.0;   // K
constexpr double k_tol = 1e-3;
constexpr double k_gamma_tol = 1.0;

constexpr double k_mole_fracs[6] = {0.0, 0.11, 0.20, 0.38, 0.5, 0.86};
constexpr double k_sampled_energies[6][5] = {
    {2.25, 3.125, 3.25, 3.5, 3.7},
    {2.625, 3.25, 3.5, 3.7, 3.875},
    {2.25, 3.0, 3.5, 3.875, 4.0},
    {2.5, 3.5, 3.75, 4.125, 4.3},
    {2.25, 3.5, 4.0, 4.625, 4.9},
    {2.875, 4.0, 4.75, 5.5, 5.8}};

double mole_frac_midpoint(int interval_index) {
    return (k_mole_fracs[interval_index] + k_mole_fracs[interval_index + 1]) / 2.0;
}

double energy_below_interval(int mole_index) {
    const double candidate = k_sampled_energies[mole_index][0] - 0.15;
    return std::max(candidate, 0.01);
}

double energy_interval_midpoint(int mole_index, int interval_index) {
    return (k_sampled_energies[mole_index][interval_index] +
            k_sampled_energies[mole_index][interval_index + 1]) /
           2.0;
}

double energy_above_interval(int mole_index) {
    return k_sampled_energies[mole_index][4] + 0.15;
}

void check_gamma(double al_frac, double energy) {
    const double actual = gamma(al_frac, energy);
    const double expected = reference::gamma(al_frac, energy);
    CHECK_NEAR(actual, expected, k_gamma_tol);
    CHECK_TRUE(actual > 0.0);
    CHECK_TRUE(std::isfinite(actual));
}

void test_alpha_endpoints() {
    CHECK_NEAR(alpha(0.0, k_field, k_temp),
               reference::alpha_gan_endpoint(k_field, k_temp), k_tol);
    CHECK_NEAR(alpha(0.65, k_field, k_temp),
               reference::alpha_algan_endpoint(k_field), k_tol);
}

void test_beta_endpoints() {
    CHECK_NEAR(beta(0.0, k_field, k_temp),
               reference::beta_gan_endpoint(k_field, k_temp), k_tol);
    CHECK_NEAR(beta(0.65, k_field, k_temp),
               reference::beta_algan_endpoint(k_field), k_tol);
}

void test_alpha_midpoint_interpolation() {
    const double al_frac = 0.325;
    const double expected = reference::alpha(al_frac, k_field, k_temp);
    CHECK_NEAR(alpha(al_frac, k_field, k_temp), expected, k_tol);
}

void test_beta_midpoint_interpolation() {
    const double al_frac = 0.325;
    const double expected = reference::beta(al_frac, k_field, k_temp);
    CHECK_NEAR(beta(al_frac, k_field, k_temp), expected, k_tol);
}

void test_alpha_decreases_with_al_fraction() {
    const double low = alpha(0.1, k_field, k_temp);
    const double high = alpha(0.5, k_field, k_temp);
    CHECK_TRUE(low > high);
}

void test_beta_is_positive_for_typical_inputs() {
    const double value = beta(0.3, k_field, k_temp);
    CHECK_TRUE(value > 0.0);
    CHECK_TRUE(std::isfinite(value));
}

void test_alpha_and_beta_differ() {
    const double a = alpha(0.3, k_field, k_temp);
    const double b = beta(0.3, k_field, k_temp);
    CHECK_TRUE(std::abs(a - b) > 1.0);
}

void test_temperature_dependence() {
    const double cold = beta(0.2, k_field, 250.0);
    const double hot = beta(0.2, k_field, 350.0);
    CHECK_TRUE(std::abs(cold - hot) > 1.0);
}

void test_gamma_at_mole_fraction_interval_midpoints() {
    for (int i = 0; i < 5; ++i) {
        const double al_frac = mole_frac_midpoint(i);
        check_gamma(al_frac, energy_below_interval(i));
        for (int j = 0; j < 4; ++j) {
            check_gamma(al_frac, energy_interval_midpoint(i, j));
            check_gamma(al_frac, k_sampled_energies[i][j]);
            check_gamma(al_frac, k_sampled_energies[i][j + 1]);
        }
        check_gamma(al_frac, energy_above_interval(i));

        check_gamma(al_frac, energy_below_interval(i + 1));
        for (int j = 0; j < 4; ++j) {
            check_gamma(al_frac, energy_interval_midpoint(i + 1, j));
            check_gamma(al_frac, k_sampled_energies[i + 1][j]);
            check_gamma(al_frac, k_sampled_energies[i + 1][j + 1]);
        }
        check_gamma(al_frac, energy_above_interval(i + 1));
    }
}

void test_gamma_at_mole_fraction_endpoints() {
    for (int i = 0; i < 6; ++i) {
        const double al_frac = k_mole_fracs[i];
        const int mole_interval = std::min(i, 4);
        check_gamma(al_frac, energy_interval_midpoint(mole_interval, 1));
        check_gamma(al_frac, energy_above_interval(mole_interval));
    }
}

void test_gamma_mole_fraction_extrapolation_above_range() {
    const double al_frac = 0.9;
    check_gamma(al_frac, energy_below_interval(4));
    for (int j = 0; j < 4; ++j) {
        check_gamma(al_frac, energy_interval_midpoint(4, j));
        check_gamma(al_frac, energy_interval_midpoint(5, j));
    }
    check_gamma(al_frac, energy_above_interval(5));
}

void test_gamma_energy_intervals_per_mole_fraction_index() {
    for (int mole_index = 0; mole_index < 6; ++mole_index) {
        const int mole_interval = std::min(mole_index, 4);
        const double al_frac = mole_frac_midpoint(mole_interval);

        check_gamma(al_frac, energy_below_interval(mole_index));
        for (int j = 0; j < 4; ++j) {
            check_gamma(al_frac, energy_interval_midpoint(mole_index, j));
            check_gamma(al_frac, k_sampled_energies[mole_index][j]);
            check_gamma(al_frac, k_sampled_energies[mole_index][j + 1]);
        }
        check_gamma(al_frac, energy_above_interval(mole_index));
    }
}

void test_gamma_all_mole_and_energy_interval_combinations() {
    for (int mole_interval = 0; mole_interval < 5; ++mole_interval) {
        const double al_frac = mole_frac_midpoint(mole_interval);
        for (int lower_energy_interval = -1; lower_energy_interval <= 4;
             ++lower_energy_interval) {
            for (int upper_energy_interval = -1; upper_energy_interval <= 4;
                 ++upper_energy_interval) {
                const double lower_energy =
                    lower_energy_interval < 0
                        ? energy_below_interval(mole_interval)
                        : lower_energy_interval == 4
                              ? energy_above_interval(mole_interval)
                              : energy_interval_midpoint(mole_interval,
                                                           lower_energy_interval);
                const double upper_energy =
                    upper_energy_interval < 0
                        ? energy_below_interval(mole_interval + 1)
                        : upper_energy_interval == 4
                              ? energy_above_interval(mole_interval + 1)
                              : energy_interval_midpoint(mole_interval + 1,
                                                           upper_energy_interval);
                check_gamma(al_frac, lower_energy);
                check_gamma(al_frac, upper_energy);
                check_gamma(al_frac, (lower_energy + upper_energy) / 2.0);
            }
        }
    }
}

}  // namespace

int main() {
    test_alpha_endpoints();
    test_beta_endpoints();
    test_alpha_midpoint_interpolation();
    test_beta_midpoint_interpolation();
    test_alpha_decreases_with_al_fraction();
    test_beta_is_positive_for_typical_inputs();
    test_alpha_and_beta_differ();
    test_temperature_dependence();
    test_gamma_at_mole_fraction_interval_midpoints();
    test_gamma_at_mole_fraction_endpoints();
    test_gamma_mole_fraction_extrapolation_above_range();
    test_gamma_energy_intervals_per_mole_fraction_index();
    test_gamma_all_mole_and_energy_interval_combinations();

    std::cout << test_count() << " checks, " << test_failures() << " failures\n";
    return test_failures() == 0 ? 0 : 1;
}
