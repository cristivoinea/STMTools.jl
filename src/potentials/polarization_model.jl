export PolarizationTable, 
       MLGScreeningModel, 
       load_mlg_screening_model,
       pol_mlg,
       v_rpa_mlg


struct PolarizationTable
    q::Vector{Float64}
    pol::Vector{Float64}
    tail_slope::Union{Nothing, Float64}
end


struct MLGScreeningModel
    neutrality_table::PolarizationTable
    neutrality_itp::Any
    filled_levels_table::Vector{PolarizationTable}
    filled_levels_itp::Vector{Any}
    nuN::Vector{Float64}
    alpha::Float64
end


function _read_scalar(path::String)
    raw = strip(read(path, String))
    return isempty(raw) ? nothing : parse(Float64, raw)
end


"""Load a polarization table stored as `<prefix>_q.csv`, `<prefix>_pi.csv`,
and an optional `<prefix>_tail_slope.txt`."""
function load_polarization_table(prefix::String)
    q = vec(readdlm(prefix * "_q.csv", ',', Float64))
    pol = vec(readdlm(prefix * "_pi.csv", ',', Float64))
    if length(q) != length(pol)
        error("Polarization table length mismatch for prefix $prefix")
    end
    tail_path = prefix * "_tail_slope.txt"
    tail_slope = isfile(tail_path) ? _read_scalar(tail_path) : nothing
    return PolarizationTable(q, pol, tail_slope)
end


function load_mlg_screening_model(
    data_dir::String;
    nuN::Vector{Float64} = Float64[],
    alpha::Float64,
)
    neutrality_table = load_polarization_table(joinpath(data_dir, "MLG_pi_nu_0"))
    neutrality_itp = QuadraticSpline(neutrality_table.pol, neutrality_table.q; extrapolation = ExtrapolationType.Linear)
    filled_levels_table = PolarizationTable[]
    filled_levels_itp = Any[]
    for n in 1:length(nuN)
        push!(filled_levels_table, load_polarization_table(joinpath(data_dir, "MLG_piN_$(n)")))
        push!(filled_levels_itp, QuadraticSpline(filled_levels_table[end].pol, filled_levels_table[end].q))
    end
    return MLGScreeningModel(neutrality_table, neutrality_itp,
                            filled_levels_table, filled_levels_itp,
                            copy(nuN), alpha)
end



"""Evaluate the dimensionless graphene polarization model used in the RPA kernel."""
function pol_mlg(model::MLGScreeningModel, q::Float64)
    qabs = abs(q)
    if qabs > model.neutrality_table.q[end]
        if model.neutrality_table.tail_slope == nothing
            pol_q = model.neutrality_itp(qabs)
        end
        q_last = model.neutrality_table.q[end]
        pol_last = model.neutrality_table.pol[end]
        slope = model.neutrality_table.tail_slope
        pol_q =  qabs * slope + (pol_last - q_last * slope) * q_last^3 / qabs^3
    else
        pol_q = model.neutrality_itp(qabs)
    end
    pol_q *= 4.0
    for (weight, filled_itp) in zip(model.nuN, model.filled_levels_itp)
        pol_q += weight * filled_itp(qabs)
    end
    return pol_q * model.alpha / sqrt(2.0)
end


"""Apply the notebook RPA dressing `V -> V / (1 - V Π)`."""
function v_rpa_mlg(v::Function, pol_fn::Function, v_ee::Function, q::Float64)
    return v(q) / (1 - v_ee(q) * pol_fn(q))
end


