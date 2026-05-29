@doc raw"""
    DETERMINANTAL_PARAMETERS = Base.ImmutableDict(
        "μ"             => "charge", 
        "density"       => "charge", 
        "spin-x"        => "spin", 
        "spin-z"        => "spin", 
        "s-wave"        => "pair", 
        "d-wave"        => "pair"
    )

List of all determinantal parameters that can be added to the wavefunction, along with the corresponding
ordering type.

"""
const DETERMINANTAL_PARAMETERS = Base.ImmutableDict(
        "charge" => ["μ", "density"],
        "spin"   => ["spin-x", "spin-z"],
        "pair"   => ["s-wave", "d-wave"]
)


@doc raw"""

    JASTROW_PARAMETERS = Base.ImmutableDict(
        "electron-electron" => ["density-density", "spin-spin"],
        "phonon-phonon"   => ["density-density", "displacement-displacement"],
        "electron-phonon"   => ["density-density", "density-displacement"]
    )

List of all Jastrow parameters that can be added, along with the corresponding
ordering type.

"""
const JASTROW_PARAMETERS = Base.ImmutableDict(
    "electron-electron"     => ["density-density", "spin-spin"],
    "phonon-phonon"         => ["density-density", "displacement-displacement"],
    "electron-phonon"       => ["density-density", "density-displacement"]
)