export ldos_anisotropic_kpm

using FuzzifiED
using LinearAlgebra
using SparseArrays


function _kpm_scaling(bounds, safety::Real)
    emin, emax = float.(bounds)
    isfinite(emin) && isfinite(emax) ||
        throw(ArgumentError("KPM spectral bounds must be finite."))
    emax > emin || throw(ArgumentError("KPM requires emax > emin."))
    0 < safety < 1 ||
        throw(ArgumentError("safety must lie in the open interval (0, 1)."))

    return ((emax - emin)/(2*(1 - safety)), (emax + emin)/2)
end


function _kpm_gershgorin_bounds(hmt::SparseMatrixCSC)
    dim = size(hmt, 1)
    centers = zeros(Float64, dim)
    radii = zeros(Float64, dim)
    rows, cols, vals = findnz(hmt)
    @inbounds for k in eachindex(vals)
        row = rows[k]
        col = cols[k]
        if row == col
            centers[row] = real(vals[k])
        else
            radii[row] += abs(vals[k])
        end
    end

    emin = minimum(centers .- radii)
    emax = maximum(centers .+ radii)
    if emax == emin
        padding = max(abs(emin), 1.0)*sqrt(eps(Float64))
        return (emin - padding, emax + padding)
    end
    return (emin, emax)
end


function _kpm_moments(
    hmt::AbstractMatrix,
    state::AbstractVector,
    num_moments::Integer,
    scaling,
)
    num_moments >= 2 ||
        throw(ArgumentError("num_moments must be at least 2."))
    length(state) == size(hmt, 1) ||
        throw(DimensionMismatch("State and Hamiltonian dimensions do not match."))

    a, b = scaling
    work_type = promote_type(eltype(hmt), eltype(state), typeof(a), typeof(b))
    state_work = Vector{work_type}(state)

    # The state is deliberately not normalised: its norm is the integrated
    # local tunnelling spectral weight, as in ldos_anisotropic.
    alpha = copy(state_work)
    beta = similar(alpha)
    mul!(beta, hmt, alpha)
    @. beta = (beta - b*alpha)/a

    moments = Vector{typeof(dot(alpha, alpha))}(undef, num_moments)
    moments[1] = dot(alpha, alpha)
    moments[2] = dot(alpha, beta)

    scratch = similar(alpha)
    out_index = 3
    while out_index <= num_moments
        mul!(scratch, hmt, beta)
        @. scratch = 2*(scratch - b*beta)/a - alpha
        alpha, beta, scratch = beta, scratch, alpha

        moments[out_index] = 2*dot(alpha, alpha) - moments[1]
        if out_index + 1 <= num_moments
            moments[out_index + 1] = 2*dot(alpha, beta) - moments[2]
        end
        out_index += 2
    end

    return moments
end


function _kpm_lorentz_kernel(num_moments::Integer, lambda::Real)
    lambda > 0 || throw(ArgumentError("Lorentz lambda must be positive."))
    kernel = Vector{Float64}(undef, num_moments)

    if lambda < 100
        denom = sinh(lambda)
        @inbounds for n in 0:(num_moments-1)
            kernel[n+1] = sinh(lambda*(1 - n/num_moments))/denom
        end
    else
        # Stable form of sinh(lambda*(1-t))/sinh(lambda).
        denom = 1 - exp(-2lambda)
        @inbounds for n in 0:(num_moments-1)
            t = n/num_moments
            kernel[n+1] = exp(-lambda*t)*
                (1 - exp(-2lambda*(1 - t)))/denom
        end
    end
    kernel[1] = 1.0
    return kernel
end


function _kpm_jackson_kernel(num_moments::Integer)
    kernel = Vector{Float64}(undef, num_moments)
    denom = num_moments + 1
    cot_term = cot(pi/denom)
    @inbounds for n in 0:(num_moments-1)
        angle = pi*n/denom
        kernel[n+1] = ((num_moments - n + 1)*cos(angle) +
                       sin(angle)*cot_term)/denom
    end
    return kernel
end


function _kpm_kernel(
    kernel,
    num_moments::Integer,
    width::Real,
    scaling_a::Real,
    kernel_parameter,
)
    if kernel isa AbstractVector
        length(kernel) == num_moments ||
            throw(DimensionMismatch("Kernel length must equal num_moments."))
        return collect(float.(kernel))
    elseif kernel === :lorentz
        # a*lambda/N is the approximate linewidth in physical energy.
        lambda = kernel_parameter === nothing ?
            num_moments*width/scaling_a : float(kernel_parameter)
        return _kpm_lorentz_kernel(num_moments, lambda)
    elseif kernel === :jackson
        return _kpm_jackson_kernel(num_moments)
    elseif kernel === :none
        return ones(Float64, num_moments)
    else
        throw(ArgumentError(
            "kernel must be :lorentz, :jackson, :none, or a vector."))
    end
end


function _kpm_reconstruct(
    moments::AbstractVector,
    enrg_range::AbstractVector,
    scaling,
    kernel,
)
    a, b = scaling
    coeffs = moments .* kernel
    coeffs[2:end] .*= 2
    density = zeros(Float64, length(enrg_range))

    @inbounds for i in eachindex(enrg_range)
        x = (enrg_range[i] - b)/a
        abs(x) < 1 || continue

        t_prev = one(x)
        value = coeffs[1]
        if length(coeffs) > 1
            t_curr = x
            value += coeffs[2]*t_curr
            for n in 3:length(coeffs)
                t_next = 2x*t_curr - t_prev
                value += coeffs[n]*t_next
                t_prev, t_curr = t_curr, t_next
            end
        end

        rho = real(value)/(pi*sqrt(1 - x^2)*a)
        density[i] = max(rho, 0.0)
    end
    return density
end


function _kpm_hamiltonian(
    bs_tunneling::Basis,
    tms_hmt::Terms,
    energy_shift::Real,
)
    raw_hmt = OpMat{Float64}(
        Operator(bs_tunneling, tms_hmt); disp_std=false)
    sparse_hmt = SparseMatrixCSC_fix(raw_hmt)
    return -sparse_hmt +
        spdiagm(0 => fill(float(energy_shift), bs_tunneling.dim))
end


"""
    ldos_anisotropic_kpm(ne, nm, bias, thetas, phi,
                         interaction_pspot, impurity_pot; kwargs...)

KPM replacement for `ldos_anisotropic`. The positional arguments and the
return value are identical:

```julia
enrg_range, dist_range, ldos = ldos_anisotropic_kpm(...)
```

`ldos` has shape `(length(enrg_range), length(thetas))`. The KPM recurrence
keeps each tunnelling state unnormalised, so integrated spectral weight has
the same meaning as in the exact implementation.

The existing `width` keyword is interpreted as the target physical linewidth
of the default Lorentz kernel. Set `kernel_parameter` explicitly to control
the dimensionless Lorentz parameter instead. `bounds`, when supplied, must be
the `(emin, emax)` bounds of the final plotted energy axis; otherwise safe
Gershgorin bounds are computed from each sparse Hamiltonian.
"""
function ldos_anisotropic_kpm(
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
    charge_qp::Rational=1//1,
    width::Real=0.01,
    enrg_res::Real=0.001,
    num_moments::Integer=512,
    bounds=nothing,
    safety::Real=0.2,
    kernel=:lorentz,
    kernel_parameter=nothing,
)
    bias in (-1, 1) || throw(ArgumentError(
        "Bias only takes input in (-1,1) and signifies the direction of the tunnelling electron."))
    width > 0 || throw(ArgumentError("width must be positive."))
    enrg_res > 0 || throw(ArgumentError("enrg_res must be positive."))
    num_moments >= 2 ||
        throw(ArgumentError("num_moments must be at least 2."))
    isempty(thetas) && throw(ArgumentError("thetas must not be empty."))
    if tip_pot !== nothing
        size(tip_pot, 1) == nm ||
            throw(DimensionMismatch("tip_pot must have nm rows."))
        size(tip_pot, 2) == length(thetas) ||
            throw(DimensionMismatch("tip_pot must have one column per theta."))
    end

    num_positions = length(thetas)
    qnd = [GetNeQNDiag(nm)]
    bs = Basis(Confs(nm, [ne], qnd))
    bs_tunneling = Basis(Confs(nm, [ne + bias], qnd))

    wf_matrix = zeros(ComplexF64, nm, num_positions)
    fermion_op = Vector{Operator}(undef, nm)
    for k in 1:nm
        wf_matrix[k, :] = computeYQQm.(nm, nm-k, thetas, phi)
        fermion_op[k] = Operator(
            bs,
            bs_tunneling,
            [Term(1, [(bias == 1) ? 1 : 0, nm-k+1])];
            red_q=0,
        )
    end

    tms_impurity_isotropic = get_onebody_terms(nm, impurity_pot)
    tms_linear_field = get_onebody_terms(
        nm, field*linear_field_pot(nm); theta=pi/2, phi=0)
    tms_hmt_isotropic = SimplifyTerms(
        GetDenIntTerms(nm, 1, 0.5.*interaction_pspot) +
        tms_impurity_isotropic)
    mu = bias*chemical_potential(
        ne, nm, tms_hmt_isotropic, interaction_pspot,
        nbr_qp, charge_qp)
    tms_hmt_notip = SimplifyTerms(tms_hmt_isotropic + tms_linear_field)

    final_background = BackgroundCorrection(
        ne + bias,
        nm,
        interaction_pspot,
        nbr_qp + bias*denominator(charge_qp),
        charge_qp,
    )

    moments = Vector{Vector{ComplexF64}}(undef, num_positions)
    scalings = Vector{Tuple{Float64,Float64}}(undef, num_positions)
    column_bounds = Vector{Tuple{Float64,Float64}}(undef, num_positions)

    function register_column!(i, hmt, initial_state)
        probe_state = zeros(ComplexF64, bs_tunneling.dim)
        for j in eachindex(fermion_op)
            probe_state .+= wf_matrix[j, i].*
                ComplexF64.(fermion_op[j]*initial_state)
        end
        local_bounds = bounds === nothing ?
            _kpm_gershgorin_bounds(hmt) : Tuple(float.(bounds))
        scaling = _kpm_scaling(local_bounds, safety)
        moments[i] = _kpm_moments(hmt, probe_state, num_moments, scaling)
        scalings[i] = scaling
        column_bounds[i] = local_bounds
        return nothing
    end

    if tip_pot === nothing
        enrg_initial, st_initial = eigensystem_corrected(
            ne, nm, bs, tms_hmt_notip, interaction_pspot,
            nbr_qp, charge_qp)
        energy_shift = final_background + enrg_initial[1] + mu
        hmt = _kpm_hamiltonian(
            bs_tunneling, tms_hmt_notip, energy_shift)

        for i in eachindex(thetas)
            register_column!(i, hmt, st_initial[:, 1])
        end
    else
        for i in eachindex(thetas)
            tms_tip = get_onebody_terms(
                nm, tip_pot[:, i]; theta=thetas[i], phi=phi)
            tms_hmt = SimplifyTerms(tms_hmt_notip + tms_tip)
            enrg_initial, st_initial = eigensystem_corrected(
                ne, nm, bs, tms_hmt, interaction_pspot,
                nbr_qp, charge_qp)
            energy_shift = final_background + enrg_initial[1] + mu
            hmt = _kpm_hamiltonian(
                bs_tunneling, tms_hmt, energy_shift)
            register_column!(i, hmt, st_initial[:, 1])
        end
    end

    enrg_min = minimum(first.(column_bounds))
    enrg_max = maximum(last.(column_bounds))
    enrg_spread = (enrg_max - enrg_min)/10
    enrg_range = (enrg_min - enrg_spread):enrg_res:(enrg_max + enrg_spread)
    dist_range = sqrt((nm-1)/2).*thetas

    ldos = zeros(Float64, length(enrg_range), num_positions)
    for i in eachindex(thetas)
        kernel_values = _kpm_kernel(
            kernel, num_moments, width, scalings[i][1], kernel_parameter)
        ldos[:, i] = _kpm_reconstruct(
            moments[i], enrg_range, scalings[i], kernel_values)
    end

    return enrg_range, dist_range, ldos
end
