


julia --project=../docs/ -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd()*"/../")); Pkg.instantiate()'
julia --color=yes --project=../docs/ make.jl


