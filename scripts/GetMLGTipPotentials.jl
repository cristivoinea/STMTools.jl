using STMTools
using ArgParse, DelimitedFiles
using QuadGK: quadgk

args = ArgParseSettings()
@add_arg_table! args begin
    "-s"
        help = "twice the number of flux quanta in the system"
        required = true
        arg_type = Int
    "--tip-dist"
        help = "distance from sample to tip (in units of the magnetic length)"
        arg_type = Float64
        default = 0.2
    "--tip-radius"
        help = "radius of the tip (in units of the magnetic length)"
        arg_type = Float64
        default = 0.073
    "--impurity-dist"
        help = "distance from sample to impurity (in units of the magnetic length)"
        arg_type = Float64
        default = 0.
    "--gate-dist"
        help = "distance from sample to gate (in units of the magnetic length)"
        arg_type = Float64
        default = 7.
    "--nbr-gates"
        help = "number of gates in the setup. options: 0/1/2 (equally spaced)"
        arg_type = Int
        default = 1
    #"--tip-theta"
    #    help = "range of polar angles for the tip in spherical geometry; accepted inputs: theta (for single value) or start:step(1):stop"
    #    arg_type = String
    #    default = "0"
    #"--tip-phi"
    #    help = "range of azimuth angles for the tip in spherical geometry; accepted inputs: phi (for single value) or start:step(1):stop"
    #    arg_type = String
    #    default = "0"
    "--theta-res"
        help = "step along the theta direction for the tip position"
        arg_type = Float64
        default = 0.1
    "--phi"
        help = "value of the phi angle of the tip (fixed)"
        arg_type = Real
        default = π/2
    "--rpa"
        help = "turn on RPA screening effect"
        action = :store_true
    "--impurity-charge"
        help = "charge of the impurity"
        arg_type = Float64
        default = 1.
    "--b-field"
        help = "strength of the magnetic field"
        arg_type = Float64
        default = 13.9
    "--eps-xy"
        help = "in-plane (parallel) dielectric constant"
        arg_type = Float64
        default = 6.6
    "--eps-z"
        help = "out-of-plane (perp.) dielectric constant"
        arg_type = Float64
        default = 3.25
    "--alpha"
        help = "fine structure constant of graphene"
        arg_type = Float64
        default = 1.85
    "--nu-vec"
        help = "filling of the N=1, N=2, LLs. for example, at nu = 3+eps, nuN = [1+eps]"
        nargs = '+'
        arg_type = Float64
        default = Float64[]
    "--output"
        help = "directory for the output file"
        default = pwd()
end



# parse arguments
parsed_args = parse_args(args)
nm = parsed_args["s"] + 1
d_i = parsed_args["impurity-dist"]
q_i = parsed_args["impurity-charge"]
d_g = parsed_args["gate-dist"]
nbr_g = parsed_args["nbr-gates"]
rpa = parsed_args["rpa"]

d_t = parsed_args["tip-dist"]
r_t = parsed_args["tip-radius"]

#thetas = parse_float_range(parsed_args["tip-theta"])
#phis = parse_float_range(parsed_args["tip-phi"])
thetas = 0:parsed_args["theta-res"]:π
phi = parsed_args["phi"]

eps_z = parsed_args["eps-z"]
eps_xy = parsed_args["eps-xy"]
beta = sqrt(eps_xy/eps_z)
eps = sqrt(eps_z*eps_xy)
alpha = parsed_args["alpha"]/eps 

nuN = parsed_args["nu-vec"]

output_dir = abspath(parsed_args["output"])

if nbr_g == 0
    screening_tag = "bare"
elseif nbr_g == 1
    screening_tag = "screened1"
else
    throw(ArgumentError("Number of gates in the setup cannot take values other than 0, 1."))
end


matrix_el = get_tip_potentials(nm, thetas, beta, eps, nbr_g, d_g, q_i, d_i, r_t, d_t; rpa=rpa, nuN=nuN, alpha=alpha)

open(joinpath(output_dir, "$(screening_tag)_tip_2s_$(nm-1)_qi_$(fmt(q_i))_di_$(fmt(d_i))_dg_$(fmt(d_g))_dt_$(fmt(d_t))$(rpa ? "_rpa" : "")_disk.dat"); write=true) do f
         write(f, "#theta #V_m\n")
         writedlm(f, vcat(thetas', matrix_el)', ',')
end
    

