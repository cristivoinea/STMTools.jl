export parse_float_range
export BackgroundCorrection
export computeYQQm
export linear_field_pot
export get_onebody_terms
export fmt

using Printf

"""
    _canonicalize_csc!(colptr, rowval, nzval)

Sort the row indices and corresponding nonzero values within each CSC column.
FuzzifiED's `OpMat` storage is not guaranteed to be in the canonical row
order required by `SparseMatrixCSC` operations.

This is the local equivalent of the fix proposed in
https://github.com/FuzzifiED/FuzzifiED.jl/pull/7.
"""
function _canonicalize_csc!(
    colptr::Vector{Int},
    rowval::Vector{Int},
    nzval,
)
    for col in 1:(length(colptr)-1)
        indices = colptr[col]:(colptr[col+1]-1)
        length(indices) <= 1 && continue

        permutation = sortperm(view(rowval, indices))
        rowval[indices] = rowval[indices][permutation]
        nzval[indices] = nzval[indices][permutation]
    end
    return nothing
end


"""
    SparseMatrixCSC_fix(mat::OpMat)

Convert a FuzzifiED `OpMat` to a canonical `SparseMatrixCSC` and then apply
the symmetry convention recorded by `mat.sym_q`.

This local implementation follows FuzzifiED.jl pull request #7 and can be
removed once STMTools requires a FuzzifiED release containing that fix.
"""
function SparseMatrixCSC_fix(mat::OpMat)
    colptr = copy(mat.colptr)
    rowval = copy(mat.rowid)
    nzval = copy(mat.elval)
    _canonicalize_csc!(colptr, rowval, nzval)

    sparse_mat = SparseMatrixCSC(
        mat.dimf, mat.dimd, colptr, rowval, nzval)
    if mat.sym_q == 0
        return sparse_mat
    elseif mat.sym_q == 1
        return sparse_mat + adjoint(sparse_mat) -
            spdiagm(0 => diag(sparse_mat))
    elseif mat.sym_q == 2
        return sparse_mat + transpose(sparse_mat) -
            spdiagm(0 => diag(sparse_mat))
    else
        throw(ArgumentError("Unknown OpMat symmetry flag: $(mat.sym_q)."))
    end
end


function parse_float_range(s::AbstractString)
    parts = split(s, ":")
    nums = parse.(Float64, parts)
    if length(nums) == 1
        return nums
    elseif length(nums) == 2
        start, stop = nums
        return range(start, stop; step=1.0)
    elseif length(nums) == 3
        start, step, stop = nums
        return start:step:stop
    else
        error("Argument must be a single number, or a range like start:stop / start:step:stop, got $s")
    end
end


function computeYQQm(nm :: Int64, m :: Int64, θ :: Union{Int64, Float64, Irrational}, ϕ :: Union{Int64, Float64, Irrational})
    return sqrt(nm*binomial(nm-1, m)/(4*π))*(cos(θ/2)^m)*(sin(θ/2)^(nm-1-m))*exp(1im*(m-(nm-1)/2)*ϕ) *(-1)^(nm-1-m)
end


function rotate_potential(nm, θ, ϕ, mat_el)
    rotated_mat_el = zeros(length(mat_el), length(mat_el))
    rotated_mat_el[diagind(rotated_mat_el)] = mat_el

    phase = Matrix{Complex}(I,nm,nm)
    phase[diagind(phase)] .*= [exp(-1im*j*ϕ) for j in 1:nm]
    if θ == 0 
        return rotated_mat_el
    else
        D = wignerd((nm-1)/2, θ)    
        if ϕ != 0
            D = phase * D
        end
        return D' * rotated_mat_el * D
    end
end


function linear_field_pot(nm)
    pot = zeros(nm)
    for i in 1:nm
        pot[i] = (2*i - 2 - (nm-1))/(nm + 1)
    end
    return pot
end


function get_onebody_terms(nm, lzcons_pot; theta=0, phi=0)
    if theta == 0
        return Terms([Term(lzcons_pot[i], [1, nm - i + 1, 0, nm - i + 1]) for i in 1:nm])
    else
        aniso_pot = rotate_potential(nm, theta, phi, lzcons_pot)
        return Terms([Term(aniso_pot[i,j], [1, nm - i + 1, 0, nm - j + 1]) for i in 1:nm for j in 1:nm])
    end
end

function BackgroundCorrection(ne, nm, pspot, n_q=0, e_q=0)
    C = pspot[1:nm]' * (2*nm .- 1 .- 2 .*(0:(nm-1)))
    return (ne^2 - (n_q*e_q)^2)* C / (2 * nm^2)
end


function fmt(x)
    return @sprintf("%.3f", x)
end
