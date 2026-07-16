export eigensystem_corrected, chemical_potential, ldos_anisotropic

using FuzzifiED


function eigensystem_corrected(
    ne::Integer,
    nm::Integer,
    bs::Basis,
    tms_hmt::Terms,
    interaction::AbstractVector{<:Real},
    nbr_qp::Integer,
    charge_qp::Rational;
    full::Bool=false,
)
    hmt_mat = OpMat{Float64}(Operator(bs, tms_hmt); disp_std=false)
    if bs.dim < 100 || full
        enrg, st = eigen(Hermitian(Matrix(hmt_mat)))
    else
        enrg, st = GetEigensystem(hmt_mat, 10; disp_std=false)
    end

    enrg_corrected = enrg .- BackgroundCorrection(ne, nm, interaction, nbr_qp, charge_qp)
    p = sortperm(enrg_corrected)
    enrg_corrected = enrg_corrected[p]
    st = st[:,p];

    return enrg_corrected, st
end


function chemical_potential(
    ne::Integer,
    nm::Integer,
    tms_hmt_mu::Terms,
    interaction::AbstractVector{<:Real},
    nbr_qp::Integer,
    charge_qp::Rational,
)
    if ne == nm
        return 0.59 # fitted value for IQH (in absence of full spin/valley simulation)
    end

    qnd = [GetNeQNDiag(nm)]
    cfs_extract = Confs(1 * nm, [ne - 1], qnd)
    bs_extract = Basis(cfs_extract)
    println("Electron extraction HS dim = $(bs_extract.dim)")

    e0_extract = eigensystem_corrected(
        ne-1, nm, bs_extract, tms_hmt_mu, interaction, nbr_qp - denominator(charge_qp), charge_qp)[1][1]

    cfs_add = Confs(1 * nm, [ne + 1], qnd)
    bs_add = Basis(cfs_add)
    println("Electron addition HS dim = $(bs_add.dim)")

    e0_add = eigensystem_corrected(
        ne+1, nm, bs_add, tms_hmt_mu, interaction, nbr_qp + denominator(charge_qp), charge_qp)[1][1]
        
    return (e0_add - e0_extract)/2
end


function ldos_anisotropic(
    ne::Integer,
    nm::Integer,
    bias::Integer,
    thetas::AbstractVector{<:Real},
    phi::Real,
    interaction_pspot::AbstractVector{<:Real},
    impurity_pot::AbstractVector{<:Real};
    tip_pot::Union{Nothing,AbstractMatrix{<:Real}}=nothing,
    field::Real=0.0,
    nbr_qp::Integer=0,
    charge_qp::Rational=1,
    width::Real=0.01,
    enrg_res::Real=0.001,
)
    qnd = [GetNeQNDiag(1 * nm)]
    cfs = Confs(1 * nm, [ne], qnd)
    bs = Basis(cfs)
    bias in (-1,1) || throw(ArgumentError("Bias only takes input in (-1,1) and signifies the direction of the tunnelling electron."))
    cfs_tunneling = Confs(1 * nm, [ne + bias], qnd)
    bs_tunneling = Basis(cfs_tunneling)

    ms = (nm-1):-2:-(nm-1)
    wf_matrix = zeros(Complex, length(ms), length(thetas))
    fermion_op = []
    for k in eachindex(ms)
        wf_matrix[k,:] = computeYQQm.(nm, nm - k,thetas, phi)
        push!(fermion_op, Operator(bs, bs_tunneling, [Term(1, [(bias == 1) ? 1 : 0 , nm - k + 1])]; red_q=0))
    end

    tms_impurity_isotropic = get_onebody_terms(nm, impurity_pot)
    tms_linear_field = get_onebody_terms(nm, field*linear_field_pot(nm); theta=pi/2, phi=0)

    tms_hmt_isotropic = SimplifyTerms(
        GetDenIntTerms(nm, 1, 0.5 .* interaction_pspot) + tms_impurity_isotropic)
    mu = bias*chemical_potential(ne, nm, tms_hmt_isotropic, interaction_pspot, nbr_qp, charge_qp)

    tms_hmt_notip = SimplifyTerms(tms_hmt_isotropic + tms_linear_field);
    
    enrg_tunneling_corrected = zeros(Real, bs_tunneling.dim, length(thetas))
    overlaps = zeros(Complex, bs_tunneling.dim, length(thetas))
    if tip_pot === nothing
        tms_hmt = tms_hmt_notip

        enrg_initial, st_initial = eigensystem_corrected(
                ne, nm, bs, tms_hmt, interaction_pspot, nbr_qp, charge_qp)

        enrg_tunneling, st_tunneling = eigensystem_corrected(
                ne, nm, bs_tunneling, tms_hmt, interaction_pspot, nbr_qp + bias*denominator(charge_qp), charge_qp, full=true)

        enrg_tunneling_corrected .= enrg_tunneling .- enrg_initial[1] .- mu

        GS_tunneling = zeros(bs_tunneling.dim, length(ms))
        for j in eachindex(ms)    
            GS_tunneling[:,j] = fermion_op[j] * st_initial[:,1]
        end
        overlaps = (st_tunneling' * GS_tunneling) * wf_matrix
    else
        for i in eachindex(thetas)
            tms_tip = get_onebody_terms(nm, tip_pot[i, :]; theta=thetas[i], phi=phi)
            tms_hmt = SimplifyTerms(tms_hmt_notip + tms_tip)

            enrg_initial, st_initial = eigensystem_corrected(
                ne, nm, bs, tms_hmt, interaction_pspot, nbr_qp, charge_qp)

            enrg_tunneling, st_tunneling = eigensystem_corrected(
                ne, nm, bs_tunneling, tms_hmt, interaction_pspot, nbr_qp + bias*denominator(charge_qp), charge_qp, full=true)

            GS_tunneling = zeros(bs_tunneling.dim, length(ms))
            for j in eachindex(ms)    
                GS_tunneling[:,j] = fermion_op[j] * st_initial[:,1]
            end

            overlaps[:, i] = (st_tunneling' * GS_tunneling) * wf_matrix[:, i]
            enrg_tunneling_corrected[:, i] = enrg_tunneling .- enrg_initial[1] .- mu
        end
    end

    enrg_min = minimum(enrg_tunneling_corrected)
    enrg_max = maximum(enrg_tunneling_corrected)
    enrg_spread = (enrg_max - enrg_min)/10

    enrg_range = -(enrg_max + 1enrg_spread):enrg_res:-(enrg_min - 1enrg_spread)
    dist_range = sqrt((nm-1)/2)*thetas

    LDOS = 1/(width*sqrt(2*π)) .* dropdims(sum(
            exp.(-1/(2width^2) .*(reshape(enrg_range, 1, 1, :) .+ enrg_tunneling_corrected).^2 
                ) .* abs.(overlaps).^2, dims=1), dims=1)'

    return enrg_range, dist_range, LDOS
end