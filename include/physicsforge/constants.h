/**
 * @file constants.h
 * @brief Physical and mathematical constants in SI units (CODATA 2022).
 *
 * Since the 2019 SI redefinition, c, h, e, k_B and N_A are exact by definition, and constants derived
 * only from them are exact up to double rounding. Measured constants carry the CODATA relative
 * standard uncertainty noted on each.
 */

#ifndef PHYSICSFORGE_CONSTANTS_H
#define PHYSICSFORGE_CONSTANTS_H

#include <numbers>

namespace physicsforge::constants {

    /** @name Mathematical constants */

    /** @brief Archimedes' constant, pi. */
    inline constexpr double PI = std::numbers::pi;
    /** @brief Full turn, 2*pi. */
    inline constexpr double TAU = 2.0 * std::numbers::pi;
    /** @brief Euler's number, e (base of the natural logarithm). */
    inline constexpr double EULER = std::numbers::e;

    /** @name Exact SI defining constants */

    /** @brief Speed of light in vacuum, c [m s^-1]. Exact. */
    inline constexpr double SPEED_OF_LIGHT = 299'792'458.0;
    /** @brief Planck constant, h [J s]. Exact. */
    inline constexpr double PLANCK = 6.626'070'15e-34;
    /** @brief Reduced Planck constant, hbar = h / (2 pi) [J s]. */
    inline constexpr double REDUCED_PLANCK = PLANCK / TAU;
    /** @brief Elementary charge, e [C]. Exact. */
    inline constexpr double ELEMENTARY_CHARGE = 1.602'176'634e-19;
    /** @brief Boltzmann constant, k_B [J K^-1]. Exact. */
    inline constexpr double BOLTZMANN = 1.380'649e-23;
    /** @brief Avogadro constant, N_A [mol^-1]. Exact. */
    inline constexpr double AVOGADRO = 6.022'140'76e23;

    /** @name Derived exact constants */

    /** @brief Molar gas constant, R = N_A k_B [J mol^-1 K^-1]. */
    inline constexpr double GAS_CONSTANT = AVOGADRO * BOLTZMANN;
    /** @brief Electron volt [J]. Numerically equal to the elementary charge. */
    inline constexpr double ELECTRON_VOLT = ELEMENTARY_CHARGE;
    /** @brief Stefan-Boltzmann constant, sigma [W m^-2 K^-4]. */
    inline constexpr double STEFAN_BOLTZMANN = 5.670'374'419e-8;
    /** @brief Wien wavelength displacement law constant, b [m K]. */
    inline constexpr double WIEN_DISPLACEMENT = 2.897'771'955e-3;
    /** @brief Standard acceleration of gravity, g_0 [m s^-2]. Exact by convention (CGPM 1901). */
    inline constexpr double STANDARD_GRAVITY = 9.806'65;
    /** @brief Standard atmosphere [Pa]. Exact by convention. */
    inline constexpr double STANDARD_ATMOSPHERE = 101'325.0;

    /** @name Measured constants (CODATA 2022) */

    /** @brief Newtonian constant of gravitation, G [m^3 kg^-1 s^-2]. Rel. uncertainty 2.2e-5. */
    inline constexpr double GRAVITATIONAL_CONSTANT = 6.674'30e-11;
    /** @brief Vacuum electric permittivity, epsilon_0 [F m^-1]. Rel. uncertainty 1.6e-10. */
    inline constexpr double VACUUM_PERMITTIVITY = 8.854'187'8188e-12;
    /** @brief Vacuum magnetic permeability, mu_0 [N A^-2]. Rel. uncertainty 1.6e-10. */
    inline constexpr double VACUUM_PERMEABILITY = 1.256'637'061'27e-6;
    /** @brief Coulomb constant, k_e = 1 / (4 pi epsilon_0) [N m^2 C^-2]. */
    inline constexpr double COULOMB_CONSTANT = 1.0 / (4.0 * PI * VACUUM_PERMITTIVITY);
    /** @brief Electron rest mass, m_e [kg]. */
    inline constexpr double ELECTRON_MASS = 9.109'383'7139e-31;
    /** @brief Proton rest mass, m_p [kg]. */
    inline constexpr double PROTON_MASS = 1.672'621'925'95e-27;
    /** @brief Neutron rest mass, m_n [kg]. */
    inline constexpr double NEUTRON_MASS = 1.674'927'500'56e-27;

} // namespace physicsforge::constants

#endif // PHYSICSFORGE_CONSTANTS_H
