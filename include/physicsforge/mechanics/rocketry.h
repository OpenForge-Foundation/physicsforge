// include/physicsforge/mechanics/rocketry.h
//
// Rocket propulsion: the Tsiolkovsky rocket equation and related quantities.
//
// The ideal rocket equation assumes a constant effective exhaust velocity and
// neglects gravity and drag losses:
//
//     delta_v = v_e ln(m0 / mf)
//
// where m0 is the initial (wet) mass and mf the final (dry) mass.

#ifndef PHYSICSFORGE_MECHANICS_ROCKETRY_H
#define PHYSICSFORGE_MECHANICS_ROCKETRY_H

namespace physicsforge::mechanics {

    // Tsiolkovsky rocket equation: delta_v = v_e ln(m0 / mf).
    // Requires v_e >= 0 and m0 >= mf > 0.
    double rocket_delta_v(double exhaust_velocity, double initial_mass, double final_mass);

    // Mass ratio m0 / mf required to achieve delta_v: R = exp(delta_v / v_e).
    // Requires delta_v >= 0, v_e > 0.
    double rocket_mass_ratio(double delta_v, double exhaust_velocity);

    // Propellant mass needed for delta_v given the final (dry) mass:
    // m_p = mf (exp(delta_v / v_e) - 1). Uses expm1 for accuracy at small delta_v.
    double rocket_propellant_mass(double delta_v, double exhaust_velocity, double final_mass);

    // Effective exhaust velocity from specific impulse: v_e = I_sp g0. Requires I_sp >= 0, g0 > 0.
    double exhaust_velocity_from_isp(double specific_impulse, double standard_gravity);

    // Specific impulse from effective exhaust velocity: I_sp = v_e / g0. Requires v_e >= 0, g0 > 0.
    double specific_impulse(double exhaust_velocity, double standard_gravity);

    // Rocket thrust: F = mdot v_e + (p_e - p_a) A_e. Requires mdot >= 0, p_e, p_a >= 0, A_e >= 0.
    double rocket_thrust(double mass_flow_rate, double exhaust_velocity, double exit_pressure, double ambient_pressure,
                         double exit_area);

    // Burn time at constant mass flow: t = m_p / mdot. Requires m_p >= 0, mdot > 0.
    double burn_time(double propellant_mass, double mass_flow_rate);

} // namespace physicsforge::mechanics

#endif // PHYSICSFORGE_MECHANICS_ROCKETRY_H
