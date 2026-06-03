@doc raw"""

    GLOBAL_MEASUREMENTS = (
        "density",
        "double_occ",
        "spin-z",
        "energy_per_site"
    )

List of all the global measurements that are made.

"""
const GLOBAL_MEASUREMENTS = (
    "density",
    "double_occ",
    "spin-z",
    "energy_per_site"
)


@doc raw"""

    LOCAL_MEASUREMENTS = Base.ImmutableDict(
        "density"                  => "ORBITAL_ID",
        "site-dependent_density"   => "ORBITAL_ID",
        "double_occ"               => "ORBITAL_ID",
        "spin-z"                   => "ORBITAL_ID",
        "site-dependent_density"   => "ORBITAL_ID"
        "hopping_energy"           => "HOPPING_ID",
        "hubbard_energy"           => "HUBBARD_ID",
        "ext_hub_energy"           => "EXT_HUB_ID",
        "phonon_kin_energy"        => "PHONON_ID",
        "phonon_pot_energy"        => "PHONON_ID",
        "X"                        => "PHONON_ID",
        "holstein_energy"          => "HOLSTEIN_ID",
        "ssh_energy"               => "SSH_ID",
    )

List of all the local measurements than can be made, with a mapping to the
corresponding type of ID each measurement is reported in terms of.

"""
const LOCAL_MEASUREMENTS = Base.ImmutableDict(
    "density"                  => "ORBITAL_ID",
    "site-dependent_density"   => "ORBITAL_ID",
    "double_occ"               => "ORBITAL_ID",
    "spin-z"                   => "ORBITAL_ID",
    "site-dependent_spin-z"    => "ORBITAL_ID",
    "hopping_energy"           => "HOPPING_ID",
    "hubbard_energy"           => "HUBBARD_ID",
    "ext_hub_energy"           => "EXT_HUB_ID",
    "phonon_kin_energy"        => "PHONON_ID",
    "phonon_pot_energy"        => "PHONON_ID",
    "X"                        => "PHONON_ID",
    "holstein_energy"          => "HOLSTEIN_ID",
    "ssh_energy"               => "SSH_ID"
)



@doc raw"""

    CORRELATION_FUNCTIONS = Base.ImmutableDict(
        "density"          => "ORBITAL_ID",
        "spin-x"           => "ORBITAL_ID",
        "spin-z"           => "ORBITAL_ID",
        "pair"             => "BOND_ID",
    )

List of all the correlation functions that can be measured, along with the corresponding
type of ID the correlation measurement is reported in terms of.
Correlation functions are well defined in both position and momentum space.

"""
const CORRELATION_FUNCTIONS = Base.ImmutableDict(
    "density"          => "ORBITAL_ID",
    "spin-x"           => "ORBITAL_ID",
    "spin-z"           => "ORBITAL_ID",
    "pair"             => "BOND_ID"
)