export disk_matrixel,
       charge_tip_singlegate_onebody,
       v_tip_bare_onebody,
       v_tip_singlegate_onebody,
       v_impurity_bare_onebody,
       v_impurity_singlegate_onebody,
       v_interaction_bare_twobody,
       v_interaction_singlegate_twobody, 
       v_interaction_doublegate_twobody


function disk_matrixel(
    vq::Function,
    m::Int;
    type:: String = "twobody",
    atol::Float64 = 1e-10,
    rtol::Float64 = 1e-10,
)
    integrand(q) = begin
            x = q * q * ((type == "onebody") ? 0.5 : 1)
            q * laguerrel(m, x) * exp(-x) * vq(q) / (2π)
    end
        return quadgk(integrand, 0, Inf, atol = atol, rtol = rtol)[1]
end


function v_impurity_bare_onebody(q::Float64, beta::Float64, eps::Float64, Z::Float64, d_i::Float64)
    return (4π*Z/q) * exp(-beta*q*d_i) / (1 + eps)
end


function v_impurity_singlegate_onebody(q::Float64, beta::Float64, eps::Float64, Z::Float64, d_i::Float64, d_g::Float64)
    if d_i == 0
        return (4π*Z/q) / (1 + eps*coth(beta*q*d_g))
    else
        return (4π*Z/q) * exp(-beta*q*d_i) * (expm1(-2*beta*q*(d_g-d_i))/expm1(-2*beta*q*d_g)) / (1 + eps*coth(beta*q*d_g))
    end
end


function charge_tip_singlegate_onebody(q::Float64, r::Float64, vq_imp::Function, R_t::Float64, d_t::Float64)
    return -R_t*vq_imp(q) * exp(-q*d_t) * q * besselj0(q*r) / (2*π)
end


function v_tip_bare_onebody(q::Float64, beta::Float64, eps::Float64, Z_t::Float64, d_t::Float64)
    return (4π*Z_t/q) * exp(-q*d_t) / (1 + eps)
end


function v_tip_singlegate_onebody(q::Float64, beta::Float64, eps::Float64, Z_t::Float64, d_t::Float64, d_g::Float64)
    return (4π*Z_t/q) * exp(-q*d_t) / (1 + eps*coth(beta*q*d_g))
end

    
function v_interaction_bare_twobody(q::Float64)
    return 2π/q
end


function v_interaction_singlegate_twobody(q::Float64, beta::Float64, eps::Float64, d_g::Float64) 
    return (4π/q) / (1 + eps*coth(beta*q*d_g))
end

    
function v_interaction_doublegate_twobody(q::Float64, beta::Float64, eps::Float64, d_g::Float64)
    return (2/q) * tanh(beta*q*d_g)
end