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

The following example evaluates a single-gate screened electron-electron
interaction at momentum `q = 1`:

```julia
using STMTools

q = 1.0
eps_xy = 6.6
eps_z = 3.25
beta = sqrt(eps_xy / eps_z)
eps = sqrt(eps_xy * eps_z)
gate_distance = 7.0

Vq = v_interaction_singlegate_twobody(
    q, beta, eps, gate_distance,
)

println(round(Vq; digits=3))
# Result: 2.231
```