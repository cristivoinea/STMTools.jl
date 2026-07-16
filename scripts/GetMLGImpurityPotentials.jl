using STMTools
using ArgParse, DelimitedFiles

args = ArgParseSettings()
@add_arg_table! args begin
    "--impurity-dist"
        help = "distance from sample to impurity"
        arg_type = Float64
        default = 0.
    "--impurity-charge"
        help = "charge of the impurity"
        arg_type = Float64
        default = 1.
    "--gate-dist"
        help = "distance from sample to gate"
        arg_type = Float64
        default = 7.
    "--nbr-gates"
        help = "number of gates in the setup. options: 0/1/2 (equally spaced)"
        arg_type = Int
        default = 1
    "--rpa"
        help = "turn on RPA screening effect"
        action = :store_true
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
    "--m-max"
        help = "maximum number of orbitals in the system"
        arg_type = Int
        default = 40
end

# parse arguments
parsed_args = parse_args(args)
d_i = parsed_args["impurity-dist"]
q_i = parsed_args["impurity-charge"]
d_g = parsed_args["gate-dist"]
nbr_g = parsed_args["nbr-gates"]
rpa = parsed_args["rpa"]

eps_z = parsed_args["eps-z"]
eps_xy = parsed_args["eps-xy"]
beta = sqrt(eps_xy/eps_z)
eps = sqrt(eps_z*eps_xy)
alpha = parsed_args["alpha"]/eps 

nuN = parsed_args["nu-vec"]
m_max = parsed_args["m-max"]

if nbr_g == 0
    screening_tag = "bare"
elseif nbr_g == 1
    screening_tag = "screened1"
else
    throw(ArgumentError("Number of gates in the setup cannot take values other than 0, 1."))
end


matrix_el = get_impurity_potentials(m_max, beta, eps, nbr_g, d_g, q_i, d_i; rpa=rpa, nuN=nuN, alpha=alpha)

writedlm("./$(screening_tag)_impurity_qi_$(fmt(q_i))_di_$(fmt(d_i))_dg_$(fmt(d_g))$(rpa ? "_rpa" : "")_disk.dat", matrix_el)