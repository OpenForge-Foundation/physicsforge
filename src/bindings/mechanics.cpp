// src/bindings/mechanics.cpp - exposes physicsforge::mechanics
//
// Scalar laws are wrapped with py::vectorize, so each accepts any mix of Python
// scalars and NumPy arrays and broadcasts in C++ (no Python-level loop).
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include "physicsforge/mechanics/dynamics.h"
#include "physicsforge/mechanics/gravitation.h"
#include "physicsforge/mechanics/kinematics.h"
#include "physicsforge/mechanics/rocketry.h"
#include "physicsforge/mechanics/rotation.h"

namespace py = pybind11;

void bind_mechanics(py::module_ &m) {
    namespace pm = physicsforge::mechanics;

    // Kinematics
    m.def("displacement", py::vectorize(&pm::displacement), py::arg("initial_velocity"), py::arg("acceleration"),
          py::arg("time"));
    m.def("final_velocity", py::vectorize(&pm::final_velocity), py::arg("initial_velocity"), py::arg("acceleration"),
          py::arg("time"));
    m.def("final_velocity_from_displacement", py::vectorize(&pm::final_velocity_from_displacement),
          py::arg("initial_velocity"), py::arg("acceleration"), py::arg("displacement"));
    m.def("average_velocity", py::vectorize(&pm::average_velocity), py::arg("displacement"), py::arg("time"));
    m.def("average_acceleration", py::vectorize(&pm::average_acceleration), py::arg("velocity_change"),
          py::arg("time"));
    m.def("projectile_range", py::vectorize(&pm::projectile_range), py::arg("initial_speed"), py::arg("launch_angle"),
          py::arg("gravity"));
    m.def("projectile_time_of_flight", py::vectorize(&pm::projectile_time_of_flight), py::arg("initial_speed"),
          py::arg("launch_angle"), py::arg("gravity"));
    m.def("projectile_max_height", py::vectorize(&pm::projectile_max_height), py::arg("initial_speed"),
          py::arg("launch_angle"), py::arg("gravity"));

    // Dynamics
    m.def("net_force", py::vectorize(&pm::net_force), py::arg("mass"), py::arg("acceleration"));
    m.def("acceleration_from_force", py::vectorize(&pm::acceleration_from_force), py::arg("force"), py::arg("mass"));
    m.def("weight", py::vectorize(&pm::weight), py::arg("mass"), py::arg("gravity"));
    m.def("momentum", py::vectorize(&pm::momentum), py::arg("mass"), py::arg("velocity"));
    m.def("impulse", py::vectorize(&pm::impulse), py::arg("force"), py::arg("duration"));
    m.def("kinetic_energy", py::vectorize(&pm::kinetic_energy), py::arg("mass"), py::arg("velocity"));
    m.def("potential_energy", py::vectorize(&pm::potential_energy), py::arg("mass"), py::arg("height"),
          py::arg("gravity"));
    m.def("work", py::vectorize(&pm::work), py::arg("force"), py::arg("distance"), py::arg("angle"));
    m.def("power", py::vectorize(&pm::power), py::arg("work"), py::arg("time"));
    m.def("spring_force", py::vectorize(&pm::spring_force), py::arg("spring_constant"), py::arg("extension"));
    m.def("spring_potential_energy", py::vectorize(&pm::spring_potential_energy), py::arg("spring_constant"),
          py::arg("extension"));
    m.def("friction_force", py::vectorize(&pm::friction_force), py::arg("coefficient"), py::arg("normal_force"));

    // Gravitation
    m.def("gravitational_force", py::vectorize(&pm::gravitational_force), py::arg("mass1"), py::arg("mass2"),
          py::arg("distance"));
    m.def("gravitational_force_vector", &pm::gravitational_force_vector, py::arg("mass1"), py::arg("mass2"),
          py::arg("position1"), py::arg("position2"));
    m.def("gravitational_field_strength", py::vectorize(&pm::gravitational_field_strength), py::arg("mass"),
          py::arg("distance"));
    m.def("gravitational_potential_energy", py::vectorize(&pm::gravitational_potential_energy), py::arg("mass1"),
          py::arg("mass2"), py::arg("distance"));
    m.def("escape_velocity", py::vectorize(&pm::escape_velocity), py::arg("central_mass"), py::arg("distance"));
    m.def("orbital_velocity", py::vectorize(&pm::orbital_velocity), py::arg("central_mass"), py::arg("distance"));
    m.def("orbital_period", py::vectorize(&pm::orbital_period), py::arg("semi_major_axis"), py::arg("central_mass"));

    // Rotation
    m.def("torque", py::vectorize(&pm::torque), py::arg("lever_arm"), py::arg("force"), py::arg("angle"));
    m.def("torque_vector", &pm::torque_vector, py::arg("lever_arm"), py::arg("force"));
    m.def("angular_acceleration", py::vectorize(&pm::angular_acceleration), py::arg("torque"),
          py::arg("moment_of_inertia"));
    m.def("angular_momentum", py::vectorize(&pm::angular_momentum), py::arg("moment_of_inertia"),
          py::arg("angular_velocity"));
    m.def("angular_momentum_vector", &pm::angular_momentum_vector, py::arg("position"), py::arg("linear_momentum"));
    m.def("rotational_kinetic_energy", py::vectorize(&pm::rotational_kinetic_energy), py::arg("moment_of_inertia"),
          py::arg("angular_velocity"));
    m.def("tangential_velocity", py::vectorize(&pm::tangential_velocity), py::arg("angular_velocity"),
          py::arg("radius"));
    m.def("centripetal_acceleration", py::vectorize(&pm::centripetal_acceleration), py::arg("speed"),
          py::arg("radius"));
    m.def("centripetal_force", py::vectorize(&pm::centripetal_force), py::arg("mass"), py::arg("speed"),
          py::arg("radius"));
    m.def("moment_of_inertia_point_mass", py::vectorize(&pm::moment_of_inertia_point_mass), py::arg("mass"),
          py::arg("radius"));
    m.def("moment_of_inertia_solid_sphere", py::vectorize(&pm::moment_of_inertia_solid_sphere), py::arg("mass"),
          py::arg("radius"));
    m.def("moment_of_inertia_hollow_sphere", py::vectorize(&pm::moment_of_inertia_hollow_sphere), py::arg("mass"),
          py::arg("radius"));
    m.def("moment_of_inertia_solid_cylinder", py::vectorize(&pm::moment_of_inertia_solid_cylinder), py::arg("mass"),
          py::arg("radius"));
    m.def("moment_of_inertia_hoop", py::vectorize(&pm::moment_of_inertia_hoop), py::arg("mass"), py::arg("radius"));
    m.def("moment_of_inertia_rod_center", py::vectorize(&pm::moment_of_inertia_rod_center), py::arg("mass"),
          py::arg("length"));
    m.def("moment_of_inertia_rod_end", py::vectorize(&pm::moment_of_inertia_rod_end), py::arg("mass"),
          py::arg("length"));
    m.def("parallel_axis", py::vectorize(&pm::parallel_axis), py::arg("moment_of_inertia_cm"), py::arg("mass"),
          py::arg("distance"));

    // Rocketry
    m.def("rocket_delta_v", py::vectorize(&pm::rocket_delta_v), py::arg("exhaust_velocity"), py::arg("initial_mass"),
          py::arg("final_mass"));
    m.def("rocket_mass_ratio", py::vectorize(&pm::rocket_mass_ratio), py::arg("delta_v"), py::arg("exhaust_velocity"));
    m.def("rocket_propellant_mass", py::vectorize(&pm::rocket_propellant_mass), py::arg("delta_v"),
          py::arg("exhaust_velocity"), py::arg("final_mass"));
    m.def("exhaust_velocity_from_isp", py::vectorize(&pm::exhaust_velocity_from_isp), py::arg("specific_impulse"),
          py::arg("standard_gravity"));
    m.def("specific_impulse", py::vectorize(&pm::specific_impulse), py::arg("exhaust_velocity"),
          py::arg("standard_gravity"));
    m.def("rocket_thrust", py::vectorize(&pm::rocket_thrust), py::arg("mass_flow_rate"), py::arg("exhaust_velocity"),
          py::arg("exit_pressure"), py::arg("ambient_pressure"), py::arg("exit_area"));
    m.def("burn_time", py::vectorize(&pm::burn_time), py::arg("propellant_mass"), py::arg("mass_flow_rate"));
}
