# Development setup

In Nuvolos
- create `/files/normann`
- create `/files/normann/venvs`

**Run the following commands:**
```
git clone https://github.com/NormannR/Tasmanian.jl
git checkout v8
git clone https://github.com/NormannR/Dynare.jl
git checkout sparsegrids
```

**In Julia REPL**
```julia
] activate /files/normann/venvs/sparsegrids
] develop /files/normann/Tasmanian.jl
] develop /files/normann/Dynare.jl
using Dynare
```
N.B.: For the commands above to work, you may need to amend the `Tasmanian.jl` and `Dynare.jl` lines in the `files/normann/venvs/sparsegrids/Manifest.toml` file, which might be set to suit my desktop configuration.

**Running the examples**
```julia
; cd /files/normann/Dynare.jl/test/models/global
context = @dynare "irbc" "-DN=2" "-DMETHODS=\"NewtonRaphson()\"";
```