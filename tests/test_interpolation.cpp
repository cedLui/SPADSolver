#include <cmath>
#include <iomanip>
#include <iostream>

#include "interpolation.h"
#include "interpolation_reference.h"
#include "test_helpers.h"

namespace {

constexpr double k_field = 3.0e6;  // V/cm
constexpr double k_temp = 298.0;   // K
constexpr double k_tol = 1e-3;
constexpr double k_gamma_tol = 1.0;

struct GammaIntervalCase {
    double al_frac_low;
    double al_frac_high;
    double al_frac;
    double energy;
};

constexpr GammaIntervalCase k_gamma_intervals[] = {
    {0.0, 0.11, 0.05, 3.75},   {0.11, 0.20, 0.15, 3.875}, {0.20, 0.38, 0.29, 4.125},
    {0.38, 0.50, 0.44, 4.625}, {0.50, 1.0,  0.65, 5.5},
};
constexpr GammaIntervalCase k_gamma_extrapolation = {0.86, 1.0, 0.9, 5.5};

void print_gamma_difference(double al_frac, double energy, const char* label) {
    const double low = gammaLowEnergy(al_frac, energy);
    const double high = gammaHighEnergy(al_frac, energy);
    const double diff = gammaDifference(al_frac, energy);
    std::cout << std::fixed << std::setprecision(6) << label << '\n'
              << "     al_frac         energy         gamma_low        gamma_high        difference\n"
              << std::setw(12) << al_frac << std::setw(15) << energy << std::setw(18)
              << low << std::setw(18) << high << std::setw(18) << diff << '\n';
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

void test_gamma_low_energy_at_sampled_energy() {
    const double al_frac = 0.0;
    const double energy = 3.125;
    const double expected = reference::gamma_low_energy(al_frac, energy);
    CHECK_NEAR(gammaLowEnergy(al_frac, energy), expected, k_gamma_tol);
}

void test_gamma_low_energy_interpolation_at_gan() {
    const double al_frac = 0.0;
    const double energy = 3.0;
    const double expected = reference::gamma_low_energy(al_frac, energy);
    CHECK_NEAR(gammaLowEnergy(al_frac, energy), expected, k_gamma_tol);
}

void test_gamma_low_energy_mole_fraction_interpolation() {
    const double al_frac = 0.055;
    const double energy = 3.0;
    const double expected = reference::gamma_low_energy(al_frac, energy);
    CHECK_NEAR(gammaLowEnergy(al_frac, energy), expected, k_gamma_tol);
}

void test_gamma_low_energy_extrapolation_above_highest_energy() {
    const double al_frac = 0.0;
    const double energy = 4.0;
    const double expected = reference::gamma_low_energy(al_frac, energy);
    CHECK_NEAR(gammaLowEnergy(al_frac, energy), expected, k_gamma_tol);
}

void test_gamma_low_energy_extrapolation_below_lowest_energy() {
    const double al_frac = 0.0;
    const double energy = 2.0;
    const double expected = reference::gamma_low_energy(al_frac, energy);
    CHECK_NEAR(gammaLowEnergy(al_frac, energy), expected, k_gamma_tol);
}

void test_gamma_low_energy_mid_interval() {
    const double al_frac = 0.29;
    const double energy = 3.6;
    const double expected = reference::gamma_low_energy(al_frac, energy);
    CHECK_NEAR(gammaLowEnergy(al_frac, energy), expected, k_gamma_tol);
}

void test_gamma_low_energy_is_positive() {
    const double value = gammaLowEnergy(0.38, 3.75);
    CHECK_TRUE(value > 0.0);
    CHECK_TRUE(std::isfinite(value));
}

void test_gamma_high_energy_at_gan() {
    const double al_frac = 0.0;
    const double energy = 4.5;
    const double expected = reference::gamma_high_energy(al_frac, energy);
    CHECK_NEAR(gammaHighEnergy(al_frac, energy), expected, k_gamma_tol);
}

void test_gamma_high_energy_mid_al_fraction() {
    const double al_frac = 0.3;
    const double energy = 5.0;
    const double expected = reference::gamma_high_energy(al_frac, energy);
    CHECK_NEAR(gammaHighEnergy(al_frac, energy), expected, k_gamma_tol);
}

void test_gamma_high_energy_extrapolation_high_al_fraction() {
    const double al_frac = 0.9;
    const double energy = 4.5;
    const double expected = reference::gamma_high_energy(al_frac, energy);
    CHECK_NEAR(gammaHighEnergy(al_frac, energy), expected, k_gamma_tol);
}

void test_gamma_high_energy_is_positive() {
    const double value = gammaHighEnergy(0.35, 4.8);
    CHECK_TRUE(value > 0.0);
    CHECK_TRUE(std::isfinite(value));
}

void test_gamma_selects_low_energy_per_interval() {
    for (const auto& interval : k_gamma_intervals) {
        const double energy = interval.energy - 0.1;
        CHECK_NEAR(gamma(interval.al_frac, energy),
                   gammaLowEnergy(interval.al_frac, energy), k_gamma_tol);
    }
    const double energy = k_gamma_extrapolation.energy - 0.1;
    CHECK_NEAR(gamma(k_gamma_extrapolation.al_frac, energy),
               gammaLowEnergy(k_gamma_extrapolation.al_frac, energy), k_gamma_tol);
}

void test_gamma_selects_high_energy_per_interval() {
    for (const auto& interval : k_gamma_intervals) {
        const double energy = interval.energy + 0.1;
        CHECK_NEAR(gamma(interval.al_frac, energy),
                   gammaHighEnergy(interval.al_frac, energy), k_gamma_tol);
    }
    const double energy = k_gamma_extrapolation.energy + 0.1;
    CHECK_NEAR(gamma(k_gamma_extrapolation.al_frac, energy),
               gammaHighEnergy(k_gamma_extrapolation.al_frac, energy), k_gamma_tol);
}

void test_gamma_uses_low_energy_at_threshold() {
    for (const auto& interval : k_gamma_intervals) {
        CHECK_NEAR(gamma(interval.al_frac, interval.energy),
                   gammaLowEnergy(interval.al_frac, interval.energy), k_gamma_tol);
    }
    CHECK_NEAR(gamma(k_gamma_extrapolation.al_frac, k_gamma_extrapolation.energy),
               gammaLowEnergy(k_gamma_extrapolation.al_frac, k_gamma_extrapolation.energy),
               k_gamma_tol);
}

void test_gamma_dispatcher_matches_reference() {
    const double cases[][2] = {{0.05, 3.0},  {0.15, 3.5},  {0.29, 3.6},
                               {0.44, 4.0},  {0.65, 4.5},  {0.9, 5.0},
                               {0.05, 3.85}, {0.15, 4.0},  {0.29, 4.25},
                               {0.44, 4.75}, {0.65, 5.6},  {0.9, 5.6}};
    for (const auto& point : cases) {
        const double al_frac = point[0];
        const double energy = point[1];
        const double expected = reference::gamma(al_frac, energy);
        CHECK_NEAR(gamma(al_frac, energy), expected, k_gamma_tol);
    }
}

void test_gamma_difference_matches_reference_in_intervals() {
    for (const auto& interval : k_gamma_intervals) {
        print_gamma_difference(interval.al_frac, interval.energy,
                               "gammaDifference in interval (reference check):");
        const double expected =
            reference::gamma_difference(interval.al_frac, interval.energy);
        CHECK_NEAR(gammaDifference(interval.al_frac, interval.energy), expected,
                   k_gamma_tol);
    }
    print_gamma_difference(k_gamma_extrapolation.al_frac, k_gamma_extrapolation.energy,
                           "gammaDifference in interval (reference check, extrapolation):");
    const double expected = reference::gamma_difference(k_gamma_extrapolation.al_frac,
                                                        k_gamma_extrapolation.energy);
    CHECK_NEAR(gammaDifference(k_gamma_extrapolation.al_frac, k_gamma_extrapolation.energy),
               expected, k_gamma_tol);
}

void test_gamma_difference_at_energy_thresholds() {
    for (const auto& interval : k_gamma_intervals) {
        print_gamma_difference(interval.al_frac, interval.energy,
                               "gammaDifference at energy threshold:");
        const double expected =
            reference::gamma_difference(interval.al_frac, interval.energy);
        CHECK_NEAR(gammaDifference(interval.al_frac, interval.energy), expected,
                   k_gamma_tol);
    }
    print_gamma_difference(k_gamma_extrapolation.al_frac, k_gamma_extrapolation.energy,
                           "gammaDifference at energy threshold (extrapolation):");
    const double expected = reference::gamma_difference(
        k_gamma_extrapolation.al_frac, k_gamma_extrapolation.energy);
    CHECK_NEAR(gammaDifference(k_gamma_extrapolation.al_frac, k_gamma_extrapolation.energy),
               expected, k_gamma_tol);
}

void test_gamma_difference_equals_high_minus_low() {
    for (const auto& interval : k_gamma_intervals) {
        print_gamma_difference(interval.al_frac, interval.energy,
                               "gammaDifference equals high - low:");
        const double expected = gammaHighEnergy(interval.al_frac, interval.energy) -
                                gammaLowEnergy(interval.al_frac, interval.energy);
        CHECK_NEAR(gammaDifference(interval.al_frac, interval.energy), expected,
                   k_gamma_tol);
    }
    print_gamma_difference(k_gamma_extrapolation.al_frac, k_gamma_extrapolation.energy,
                           "gammaDifference equals high - low (extrapolation):");
    const double expected =
        gammaHighEnergy(k_gamma_extrapolation.al_frac, k_gamma_extrapolation.energy) -
        gammaLowEnergy(k_gamma_extrapolation.al_frac, k_gamma_extrapolation.energy);
    CHECK_NEAR(gammaDifference(k_gamma_extrapolation.al_frac, k_gamma_extrapolation.energy),
               expected, k_gamma_tol);
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
    test_gamma_low_energy_at_sampled_energy();
    test_gamma_low_energy_interpolation_at_gan();
    test_gamma_low_energy_mole_fraction_interpolation();
    test_gamma_low_energy_extrapolation_above_highest_energy();
    test_gamma_low_energy_extrapolation_below_lowest_energy();
    test_gamma_low_energy_mid_interval();
    test_gamma_low_energy_is_positive();
    test_gamma_high_energy_at_gan();
    test_gamma_high_energy_mid_al_fraction();
    test_gamma_high_energy_extrapolation_high_al_fraction();
    test_gamma_high_energy_is_positive();
    test_gamma_selects_low_energy_per_interval();
    test_gamma_selects_high_energy_per_interval();
    test_gamma_uses_low_energy_at_threshold();
    test_gamma_dispatcher_matches_reference();
    test_gamma_difference_matches_reference_in_intervals();
    test_gamma_difference_at_energy_thresholds();
    test_gamma_difference_equals_high_minus_low();

    std::cout << test_count() << " checks, " << test_failures() << " failures\n";
    return test_failures() == 0 ? 0 : 1;
}
