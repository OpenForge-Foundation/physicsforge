/**
 * @file rocketry.h
 * @brief Rocket propulsion: the Tsiolkovsky rocket equation and related quantities.
 *
 * The ideal rocket equation \f$ \Delta v = v_e \ln(m_0 / m_f) \f$ assumes a constant effective
 * exhaust velocity and neglects gravity and drag losses.
 *
 * @ingroup mechanics
 */

#ifndef PHYSICSFORGE_MECHANICS_ROCKETRY_H
#define PHYSICSFORGE_MECHANICS_ROCKETRY_H

namespace physicsforge::mechanics {

    /**
     * @brief Tsiolkovsky rocket equation, \f$ \Delta v = v_e \ln(m_0 / m_f) \f$.
     * @param exhaust_velocity Effective exhaust velocity \f$ v_e \f$ [m/s]; must be non-negative.
     * @param initial_mass Initial (wet) mass \f$ m_0 \f$ [kg]; must satisfy \f$ m_0 \ge m_f \f$.
     * @param final_mass Final (dry) mass \f$ m_f \f$ [kg]; must be positive.
     * @return Ideal velocity change \f$ \Delta v \f$ [m/s].
     * @throws std::invalid_argument If @p final_mass is not positive or exceeds @p initial_mass.
     */
    double rocket_delta_v(double exhaust_velocity, double initial_mass, double final_mass);

    /**
     * @brief Mass ratio required for a velocity change, \f$ m_0 / m_f = e^{\Delta v / v_e} \f$.
     * @param delta_v Required velocity change \f$ \Delta v \f$ [m/s]; must be non-negative.
     * @param exhaust_velocity Effective exhaust velocity \f$ v_e \f$ [m/s]; must be positive.
     * @return Dimensionless mass ratio \f$ m_0 / m_f \f$.
     * @throws std::invalid_argument If @p delta_v is negative or @p exhaust_velocity is not positive.
     */
    double rocket_mass_ratio(double delta_v, double exhaust_velocity);

    /**
     * @brief Propellant mass for a velocity change, \f$ m_p = m_f (e^{\Delta v / v_e} - 1) \f$.
     *
     * Evaluated with @c expm1, so it stays accurate for small \f$ \Delta v / v_e \f$.
     *
     * @param delta_v Required velocity change \f$ \Delta v \f$ [m/s]; must be non-negative.
     * @param exhaust_velocity Effective exhaust velocity \f$ v_e \f$ [m/s]; must be positive.
     * @param final_mass Final (dry) mass \f$ m_f \f$ [kg]; must be positive.
     * @return Propellant mass \f$ m_p \f$ [kg].
     * @throws std::invalid_argument If any precondition above is violated.
     */
    double rocket_propellant_mass(double delta_v, double exhaust_velocity, double final_mass);

    /**
     * @brief Effective exhaust velocity from specific impulse, \f$ v_e = I_{sp}\, g_0 \f$.
     * @param specific_impulse Specific impulse \f$ I_{sp} \f$ [s]; must be non-negative.
     * @param standard_gravity Standard gravity \f$ g_0 \f$ [m/s²]; must be positive.
     * @return Effective exhaust velocity [m/s].
     * @throws std::invalid_argument If @p specific_impulse is negative or @p standard_gravity is not positive.
     */
    double exhaust_velocity_from_isp(double specific_impulse, double standard_gravity);

    /**
     * @brief Specific impulse from effective exhaust velocity, \f$ I_{sp} = v_e / g_0 \f$.
     * @param exhaust_velocity Effective exhaust velocity \f$ v_e \f$ [m/s]; must be non-negative.
     * @param standard_gravity Standard gravity \f$ g_0 \f$ [m/s²]; must be positive.
     * @return Specific impulse [s].
     * @throws std::invalid_argument If @p exhaust_velocity is negative or @p standard_gravity is not positive.
     */
    double specific_impulse(double exhaust_velocity, double standard_gravity);

    /**
     * @brief Rocket thrust, \f$ F = \dot m v_e + (p_e - p_a) A_e \f$.
     * @param mass_flow_rate Propellant mass flow rate \f$ \dot m \f$ [kg/s]; must be non-negative.
     * @param exhaust_velocity Exhaust velocity \f$ v_e \f$ [m/s].
     * @param exit_pressure Nozzle exit pressure \f$ p_e \f$ [Pa]; must be non-negative.
     * @param ambient_pressure Ambient pressure \f$ p_a \f$ [Pa]; must be non-negative.
     * @param exit_area Nozzle exit area \f$ A_e \f$ [m²]; must be non-negative.
     * @return Thrust \f$ F \f$ [N].
     * @throws std::invalid_argument If a flow rate, pressure or area is negative.
     */
    double rocket_thrust(double mass_flow_rate, double exhaust_velocity, double exit_pressure, double ambient_pressure,
                         double exit_area);

    /**
     * @brief Burn time at constant mass flow, \f$ t = m_p / \dot m \f$.
     * @param propellant_mass Propellant mass \f$ m_p \f$ [kg]; must be non-negative.
     * @param mass_flow_rate Mass flow rate \f$ \dot m \f$ [kg/s]; must be positive.
     * @return Burn time [s].
     * @throws std::invalid_argument If @p propellant_mass is negative or @p mass_flow_rate is not positive.
     */
    double burn_time(double propellant_mass, double mass_flow_rate);

} // namespace physicsforge::mechanics

#endif // PHYSICSFORGE_MECHANICS_ROCKETRY_H
