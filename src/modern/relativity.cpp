// src/modern/relativity.cpp
#include "physicsforge/modern/relativity.h"

#include <cmath>

#include "physicsforge/constants.h"
#include "physicsforge/detail/validation.h"

namespace physicsforge::modern {

    using constants::SPEED_OF_LIGHT;
    using detail::require_non_negative;

    namespace {
        // beta^2 = (v/c)^2 after checking |v| < c.
        double beta_squared(const double velocity) {
            detail::require_finite(velocity, "velocity");
            const double beta = velocity / SPEED_OF_LIGHT;
            if (std::abs(beta) >= 1.0) {
                detail::throw_invalid("|velocity|", "less than the speed of light", std::abs(velocity));
            }
            return beta * beta;
        }

        // gamma - 1 computed without cancellation for small beta:
        //   gamma - 1 = beta^2 / (sqrt(1 - beta^2) (1 + sqrt(1 - beta^2))).
        double gamma_minus_one(const double velocity) {
            const double b2 = beta_squared(velocity);
            const double s = std::sqrt(1.0 - b2);
            return b2 / (s * (1.0 + s));
        }
    } // namespace

    double lorentz_factor(const double velocity) { return 1.0 / std::sqrt(1.0 - beta_squared(velocity)); }

    double time_dilation(const double proper_time, const double velocity) {
        require_non_negative(proper_time, "proper_time");
        return proper_time * lorentz_factor(velocity);
    }

    double length_contraction(const double proper_length, const double velocity) {
        require_non_negative(proper_length, "proper_length");
        return proper_length * std::sqrt(1.0 - beta_squared(velocity));
    }

    double rest_energy(const double mass) {
        require_non_negative(mass, "mass");
        return mass * SPEED_OF_LIGHT * SPEED_OF_LIGHT;
    }

    double relativistic_momentum(const double mass, const double velocity) {
        require_non_negative(mass, "mass");
        return lorentz_factor(velocity) * mass * velocity;
    }

    double relativistic_kinetic_energy(const double mass, const double velocity) {
        require_non_negative(mass, "mass");
        return gamma_minus_one(velocity) * mass * SPEED_OF_LIGHT * SPEED_OF_LIGHT;
    }

} // namespace physicsforge::modern
