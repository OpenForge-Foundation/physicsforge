// include/physicsforge/modern/quantum.h
//
// Foundational quantum relations built on Planck's constant: photon energy and
// momentum, matter waves, the photoelectric effect, black-body radiation, the Bohr
// hydrogen atom and the Heisenberg uncertainty bound. SI units (J, m, Hz, K).

#ifndef PHYSICSFORGE_MODERN_QUANTUM_H
#define PHYSICSFORGE_MODERN_QUANTUM_H

namespace physicsforge::modern {

    // Planck-Einstein relation: E = h f. Requires f >= 0.
    double photon_energy_from_frequency(double frequency);

    // Photon energy from wavelength: E = h c / lambda. Requires lambda > 0.
    double photon_energy_from_wavelength(double wavelength);

    // Photon momentum: p = h / lambda. Requires lambda > 0.
    double photon_momentum(double wavelength);

    // Wavelength of light in vacuum from frequency: lambda = c / f. Requires f > 0.
    double wavelength_from_frequency(double frequency);

    // de Broglie wavelength of a (non-relativistic) massive particle: lambda = h / (m v).
    // Requires m > 0 and v != 0; the sign of v is ignored.
    double de_broglie_wavelength(double mass, double velocity);

    // Einstein photoelectric equation: K_max = h f - phi. Returns 0 below the threshold
    // frequency (no emission). Requires f >= 0, phi >= 0.
    double photoelectric_kinetic_energy(double frequency, double work_function);

    // Planck's law, spectral radiance per unit wavelength [W sr^-1 m^-3]:
    //   B(lambda, T) = (2 h c^2 / lambda^5) / (exp(h c / (lambda k_B T)) - 1).
    // Requires lambda > 0, T >= 0 (returns 0 at T = 0).
    double planck_spectral_radiance(double wavelength, double temperature);

    // Wien's displacement law, peak wavelength: lambda_max = b / T. Requires T > 0.
    double wien_peak_wavelength(double temperature);

    // Bohr model energy of hydrogen level n (infinite nuclear mass):
    //   E_n = -m_e e^4 / (8 eps_0^2 h^2 n^2). Requires n >= 1.
    double bohr_energy_level(double principal_quantum_number);

    // Heisenberg bound, minimum momentum uncertainty: dp >= hbar / (2 dx). Requires dx > 0.
    double heisenberg_min_momentum_uncertainty(double position_uncertainty);

} // namespace physicsforge::modern

#endif // PHYSICSFORGE_MODERN_QUANTUM_H
