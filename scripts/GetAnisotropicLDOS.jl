using STMTools
using ArgParse

using FuzzifiED
using DelimitedFiles
using LinearAlgebra, SparseArrays


settings = ArgParseSettings()
@add_arg_table! settings begin
    "-n"
        help = "number of electrons in the system"
        arg_type = Int
        required = true
    "-s"
        help = "twice the number of flux quanta in the system"
        required = true
        arg_type = Int
    "--bias"
        help = "type of tunneling: +1/-1 for electron removal/addition"
        arg_type = Int
        required = true
    "--impurity-dist"
        help = "distance from sample to impurity"
        arg_type = Float64
        default = nothing
    "--impurity-charge"
        help = "charge of the impurity"
        arg_type = Float64
        default = 1.
    "--gate-dist"
        help = "distance from sample to gate"
        arg_type = Float64
        required = nothing 
    "--nbr-gates"
        help = "number of gates in the setup. options: 0/1/2 (equally spaced)"
        arg_type = Int
        default = 0
    "--tip-dist"
        help = "distance from sample to tip"
        arg_type = Float64
        default = nothing
    "--rpa"
        help = "turn on RPA screening effect"
        action = :store_true
    "--field"
        help = "additional linear field"
        arg_type = Float64
        default = 0.
    "--width"
        help = "experimental resolution (width of the Gaussian applied)"
        arg_type = Float64
        default = 0.01
    "--enrg-res"
        help = "LDOS sampling along energy axis"
        arg_type = Float64
        default = 0.001
    "--theta-res"
        help = "LDOS sampling along theta axis"
        arg_type = Float64
        default = 0.02
    "--phi-cut"
        help = "LDOS direction along phi axis"
        arg_type = Float64
        default = π/2
    "--input"
        help = "directory containing the pseudopotential files"
        default = pwd()
    "--output"
        help = "directory for the output file"
        default = pwd()
    "--threads"
        help = "threads used by FuzzifiED"
        arg_type = Int
        default = Threads.nthreads()
end

# parse arguments
args = parse_args(settings)
bias = args["bias"]
ne = args["n"]
nm = args["s"] + 1
impurity_charge = args["impurity-charge"]
d_i = args["impurity-dist"]
d_g = args["gate-dist"]
nbr_g = args["nbr-gates"]
d_t  = args["tip-dist"]
rpa = args["rpa"]
field = args["field"]
width = args["width"]
enrg_res = args["enrg-res"]
theta_res = args["theta-res"]
phi = args["phi-cut"]

input_dir = abspath(args["input"])
output_dir = abspath(args["output"])

FuzzifiED.NumThreads=args["threads"]
FuzzifiED.ElementType = Float64

tip_flag = (d_t !== nothing)
rpa_tag = rpa ? "_rpa" : ""
if  d_i === nothing
    throw(ArgumentError("Impurity needs to be specified in the current code."))
end
if nbr_g > 0 && d_g === nothing
    throw(ArgumentError("Distance to gate needs to be specified in this setup."))
end
screening_tag = (nbr_g == 0) ? "bare" : "screened$(nbr_g)"
if nbr_g > 1 
    throw(ArgumentError("Function not implemented for the desired gate setup."))
end


thetas = 0:theta_res:π
if ne != nm # initial number of anyons assuming we work with Laughlin 1/3
    nbr_qp = (3*ne - 3) - (nm - 1)
    charge_qp = 1//3
else # exception for nu=1
    nbr_qp = 0
    charge_qp = 1//1
end


interaction_pspot = vec(readdlm(joinpath(input_dir, "$(screening_tag)_coulomb_dg_$(fmt(d_g))$(rpa_tag)_disk.dat")))
impurity_pot = vec(readdlm(joinpath(input_dir, "$(screening_tag)_impurity_qi_$(fmt(impurity_charge))_di_$(fmt(d_i))_dg_$(fmt(d_g))$(rpa_tag)_disk.dat")))
println("Retrieved pseudopotentials and impurity one-body potentials.")

if tip_flag
    tip_data = readdlm(joinpath(
            input_dir, "$(screening_tag)_tip_2s_$(nm-1)_qi_$(fmt(impurity_charge))_di_$(fmt(d_i))_dg_$(fmt(d_g))_dt_$(fmt(d_t))$(rpa_tag)_disk.dat"), 
            ',', Float64, '\n'; header=true)[1]
    tip_data[:,1] ≈ thetas || error("Tip potentials are calculated at different tip locations than requested.")
    tip_pot = tip_data[:, 2:end]
    println("Retrieved tip one-body potentials.")
end


enrg_range, dist_range, ldos = ldos_anisotropic(
    ne, nm, bias, thetas, phi,
    interaction_pspot, impurity_pot; tip_pot = (tip_flag ? tip_pot : nothing), field=field,
    nbr_qp=nbr_qp, charge_qp=charge_qp,
    width=width, enrg_res=enrg_res)

open(joinpath(output_dir, "ldos_n_$(ne)_2s_$(nm-1)_bias_$(bias)_qi_$(fmt(impurity_charge))_di_$(fmt(d_i))_dg_$(fmt(d_g))_dt_$(fmt(d_t))$(rpa_tag)_field_$(fmt(field)).txt"); write=true) do f
         write(f, "#d/l_B #E/E_c #LDOS\n")
         writedlm(f, hcat(repeat(dist_range, inner=[length(enrg_range)]), vec(repeat(enrg_range, outer=[length(dist_range)])), vec(ldos) ), ',')
end
