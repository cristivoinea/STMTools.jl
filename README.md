# STMTools.jl

`STMTools.jl` is a Julia package for modelling scanning tunnelling microscopy
(STM) experiments in quantum Hall systems.


## Installation

Install `STMTools.jl` directly from its Git repository into the currently
active Julia environment:

```julia
import Pkg
Pkg.add(url="https://github.com/cristivoinea/STMTools.jl.git")
```

Once installed, the package can be loaded from any Julia session using that
environment:

```julia
using STMTools
```


## Basic use

The following example computes the STM extraction spectrum for a Laughlin state with 5 particles and 3 additional quasiparticles. First, one needs to define the setup, including the tunneling direction (bias), impurity charge, and gate/impurity/tip distance from the sample.

```julia
using STMTools, Plots

ne = 5
nbr_qp = 3
charge_qp = 1//3
nm = 3*ne - 2 -nbr_qp

bias = -1
impurity_charge = -1.
d_i = 0.
d_g = 7.
nbr_g = 1.
d_t  = 0.2
rpa = true
field = 0.

width = 0.01
enrg_res = 0.002
phi = pi/2
thetas = 0:0.02:pi
```

The interaction pseudopotentials and the one-body potentials from the impurity, tip, and linear field need to be retrieved first:
```julia
interaction_pspot = get_pseudopotentials(nm, nbr_g, d_g; rpa=true)
impurity_pot = get_impurity_potentials(nm, nbr_g, d_g, q_i, d_i; rpa=true)
tip_pot = get_tip_potentials(nm, thetas, nbr_g, d_g, q_i, d_i, r_t, d_t; rpa=true)
```

The LDOS can now be computed and visualised:
```julia
enrg_range, dist_range, ldos = ldos_anisotropic(
    ne, nm, bias, thetas, phi,
    interaction_pspot, impurity_pot; tip_pot = tip_pot, field=field,
    nbr_qp=nbr_qp, charge_qp=charge_qp,
    width=width, enrg_res=enrg_res)


p = heatmap(dist_range, enrg_range, ldos, title="LDOS using ED",
            yguide=L"E\,\,(E_C)", xguide=L"x\,\,(\ell_B)", ylims=(-1.5, 0.))
```
<img src="misc/ldos_quick_start.png"
     width="500">

If one wants to trade resolution over speed, then the LDOS can be also computed using the kernel polynomial method instead of full exact diagonalization.
```julia
enrg_range, dist_range, ldos = ldos_anisotropic_kpm(
    ne, nm, bias, thetas, phi,
    interaction_pspot, impurity_pot; tip_pot = tip_pot, field=field,
    nbr_qp=nbr_qp, charge_qp=charge_qp,
    width=width, enrg_res=enrg_res,
    num_moments=2048)


p = heatmap(dist_range, enrg_range, ldos, title="LDOS using KPM",
            yguide=L"E\,\,(E_C)", xguide=L"x\,\,(\ell_B)", ylims=(-1.5, 0.))
```
<img src="misc/ldos_kpm_quick_start.png"
     width="500">