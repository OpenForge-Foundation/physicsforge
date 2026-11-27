/**
 * @file quantum.h
 * @brief Quantum relations built on Planck's constant.
 *
 * Photon energy and momentum, matter waves, the photoelectric effect, black-body radiation,
 * the Bohr hydrogen atom and the Heisenberg uncertainty bound. SI units (J, m, Hz, K).
 *
 * @ingroup quantum
 */

#ifndef PHYSICSFORGE_QUANTUM_H
#define PHYSICSFORGE_QUANTUM_H

/** @defgroup quantum Quantum
 *  @brief Planck–Einstein relation, matter waves, black-body radiation and the Bohr atom.
 */

namespace physicsforge::quantum {

    /**
     * @brief Planck–Einstein relation, \f$ E = h f \f$.
     * @param frequency Frequency \f$ f \f$ [Hz]; must be non-negative.
     * @return Photon energy [J].
     * @throws std::invalid_argument If @p frequency is negative.
     */
    double photon_energy_from_frequency(double frequency);

    /**
     * @brief Photon energy from wavelength, \f$ E = h c / \lambda \f$.
     * @param wavelength Vacuum wavelength \f$ \lambda \f$ [m]; must be positive.
     * @return Photon energy [J].
     * @throws std::invalid_argument If @p wavelength is not positive.
     */
    double photon_energy_from_wavelength(double wavelength);

    /**
     * @brief Photon momentum, \f$ p = h / \lambda \f$.
     * @param wavelength Vacuum wavelength \f$ \lambda \f$ [m]; must be positive.
     * @return Momentum [kg·m/s].
     * @throws std::invalid_argument If @p wavelength is not positive.
     */
    double photon_momentum(double wavelength);

    /**
     * @brief Vacuum wavelength of light from frequency, \f$ \lambda = c / f \f$.
     * @param frequency Frequency \f$ f \f$ [Hz]; must be positive.
     * @return Wavelength [m].
     * @throws std::invalid_argument If @p frequency is not positive.
     */
    double wavelength_from_frequency(double frequency);

    /**
     * @brief de Broglie wavelength of a non-relativistic massive particle, \f$ \lambda = h / (m |v|) \f$.
     * @param mass Particle mass \f$ m \f$ [kg]; must be positive.
     * @param velocity Particle velocity \f$ v \f$ [m/s]; must be non-zero (its sign is ignored).
     * @return Wavelength [m].
     * @throws std::invalid_argument If @p mass is not positive or @p velocity is zero.
     */
    double de_broglie_wavelength(double mass, double velocity);

    /**
     * @brief Einstein photoelectric equation, \f$ K_{max} = \max(0,\, h f - \phi) \f$.
     * @param frequency Incident photon frequency \f$ f \f$ [Hz]; must be non-negative.
     * @param work_function Work function of the surface \f$ \phi \f$ [J]; must be non-negative.
     * @return Maximum kinetic energy of emitted electrons [J]; zero below the threshold frequency.
     * @throws std::invalid_argument If either input is negative.
     */
    double photoelectric_kinetic_energy(double frequency, double work_function);

    /**
     * @brief Planck's law, spectral radiance per unit wavelength,
     *        \f$ B_\lambda(T) = \dfrac{2 h c^2}{\lambda^5} \dfrac{1}{e^{h c / (\lambda k_B T)} - 1} \f$.
     *
     * Evaluated with @c expm1; the deep Wien tail underflows cleanly to zero.
     *
     * @param wavelength Wavelength \f$ \lambda \f$ [m]; must be positive.
     * @param temperature Absolute temperature \f$ T \f$ [K]; must be non-negative (returns 0 at 0 K).
     * @return Spectral radiance [W/(sr·m³)].
     * @throws std::invalid_argument If @p wavelength is not positive or @p temperature is negative.
     */
    double planck_spectral_radiance(double wavelength, double temperature);

    /**
     * @brief Wien's displacement law, \f$ \lambda_{max} = b / T \f$.
     * @param temperature Absolute temperature \f$ T \f$ [K]; must be positive.
     * @return Peak wavelength [m].
     * @throws std::invalid_argument If @p temperature is not positive.
     */
    double wien_peak_wavelength(double temperature);

    /**
     * @brief Bohr-model energy of hydrogen level @p n (infinite nuclear mass),
     *        \f$ E_n = -\dfrac{m_e e^4}{8 \varepsilon_0^2 h^2 n^2} \f$.
     * @param principal_quantum_number Principal quantum number \f$ n \f$; must be an integer \f$ \ge 1 \f$.
     * @return Level energy [J] (negative; bound state).
     * @throws std::invalid_argument If @p principal_quantum_number is not a positive integer.
     */
    double bohr_energy_level(double principal_quantum_number);

    /**
     * @brief Heisenberg bound on momentum uncertainty, \f$ \Delta p \ge \hbar / (2 \Delta x) \f$.
     * @param position_uncertainty Position uncertainty \f$ \Delta x \f$ [m]; must be positive.
     * @return Minimum momentum uncertainty [kg·m/s].
     * @throws std::invalid_argument If @p position_uncertainty is not positive.
     */
    double heisenberg_min_momentum_uncertainty(double position_uncertainty);

} // namespace physicsforge::quantum

#endif // PHYSICSFORGE_QUANTUM_H
