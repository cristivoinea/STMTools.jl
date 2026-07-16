export get_pseudopotentials, get_impurity_potentials, get_tip_potentials


function get_pseudopotentials(nm, nbr_g, d_g; eps_z=3.25, eps_xy=6.6, rpa=false, nuN=[], alpha=1.85)
    beta = sqrt(eps_xy/eps_z)
    eps = sqrt(eps_z*eps_xy)
    alpha /= eps 

    v0 = if nbr_g == 0
        q -> v_interaction_bare_twobody(q)
    elseif nbr_g == 1
        q -> v_interaction_singlegate_twobody(q, beta, eps, d_g)
    elseif nbr_g == 2
        q -> v_interaction_doublegate_twobody(q, beta, eps, d_g)
    else
        throw(ArgumentError("Number of gates in the setup cannot take values other than 0, 1, or 2."))
    end

    matrix_el = zeros(nm)
    if rpa
        model = load_mlg_screening_model(joinpath(@__DIR__, "../../data/mlg/"); nuN = nuN, alpha = alpha)
        pol_fn(q) = pol_mlg(model, q)
    end
    v = if rpa
        q -> v_rpa_mlg(v0, pol_fn, v0, q)
    else
        q -> v0(q)
    end

    Threads.@threads for m in 0:nm-1
        matrix_el[m+1] = disk_matrixel(v, m; type="twobody")
    end
    return matrix_el * eps
end


function get_impurity_potentials(nm, nbr_g, d_g, q_i, d_i; eps_z=3.25, eps_xy=6.6, rpa=false, nuN=[], alpha=1.85)
    beta = sqrt(eps_xy/eps_z)
    eps = sqrt(eps_z*eps_xy)
    alpha /= eps 

    v0 = if nbr_g == 0
        q -> v_impurity_bare_onebody(q, beta, eps, q_i, d_i)
    elseif nbr_g == 1
        q -> v_impurity_singlegate_onebody(q, beta, eps, q_i, d_i, d_g)
    end
    v_ee = if nbr_g == 0
        q -> v_interaction_bare_twobody(q)
    elseif nbr_g == 1
        q -> v_interaction_singlegate_twobody(q, beta, eps, d_g)#
    end
    nbr_g in (0,1) || throw(ArgumentError("Number of gates in the setup cannot take values other than 0, 1."))

    matrix_el = zeros(nm)
    if rpa
        model = load_mlg_screening_model(joinpath(@__DIR__, "../../data/mlg/"); nuN = nuN, alpha = alpha)
        pol_fn(q) = pol_mlg(model, q)
    end
    v = if rpa
        q -> v_rpa_mlg(v0, pol_fn, v_ee, q)
    else
        q -> v0(q)
    end

    Threads.@threads for m in 0:nm-1
        matrix_el[m+1] = disk_matrixel(v, m; type="onebody")
    end
    return matrix_el * eps
end


function get_tip_potentials(nm, thetas, nbr_g, d_g, q_i, d_i, r_t, d_t; eps_z=3.25, eps_xy=6.6, rpa=false, nuN=[], alpha=1.85)
    beta = sqrt(eps_xy/eps_z)
    eps = sqrt(eps_z*eps_xy)
    alpha /= eps 

    v0_imp = if nbr_g == 0
        q -> v_impurity_bare_onebody(q, beta, eps, q_i, d_i)
    elseif nbr_g == 1
        q -> v_impurity_singlegate_onebody(q, beta, eps, q_i, d_i, d_g)
    end
    v_ee = if nbr_g == 0
        q -> v_interaction_bare_twobody(q)
    elseif nbr_g == 1
        q -> v_interaction_singlegate_twobody(q, beta, eps, d_g)#
    end
    nbr_g in (0,1) || throw(ArgumentError("Number of gates in the setup cannot take values other than 0, 1."))

    
   if rpa
        model = load_mlg_screening_model(joinpath(@__DIR__, "../../data/mlg/"); nuN = nuN, alpha = alpha)
        pol_fn(q) = pol_mlg(model, q)
    end
    v_imp = if rpa
        q -> v_rpa_mlg(v0_imp, pol_fn, v_ee, q)
    else
        q -> v0(q)
    end


    matrix_el = zeros(nm, length(thetas))
    for i in eachindex(thetas)
        #for phi in phis
        theta = thetas[i]
        x = sqrt((nm-1)/2)*theta
        charge_integrand(q) = charge_tip_singlegate_onebody(q, x, v_imp, r_t, d_t)
        impurity_charge_t = quadgk(charge_integrand, 0, Inf, atol = 1e-10, rtol = 1e-10)[1]

        v0_tip = if nbr_g == 0
            q -> v_tip_bare_onebody(q, beta, eps, impurity_charge_t, d_t)
        elseif nbr_g == 1
            q -> v_tip_singlegate_onebody(q, beta, eps, impurity_charge_t, d_t, d_g)
        end

        v_tip = if rpa
            q -> v_rpa_mlg(v0_tip, pol_fn, v_ee, q)
        else
            q -> v0_tip(q)
        end

            
        Threads.@threads for m in 0:nm-1
            matrix_el[m+1, i] = disk_matrixel(v_tip, m; type="onebody")
        end
        #end
    end
    return matrix_el * eps
end