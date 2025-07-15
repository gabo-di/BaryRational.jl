# so x can be real while f and w are complex
using Infiltrator
struct MovingAAAapprox{T,W} <: BRInterp
    m::Int # m <= M amount of support points
    X::AbstractVector{T} # support points, fixed size M
    F::AbstractVector{W} # function values at support points, size M
    W::AbstractVector{W} # baricentric weights, size m
    A::AbstractMatrix{W} # M × m Loewner matrix 
    J::AbstractVector{Int} # permutation vector, size M
end

function Base.getproperty(obj::MovingAAAapprox, sym::Symbol)
    if sym === :f
        return view(obj.F, obj.J[1:obj.m])
    elseif sym === :x
        return view(obj.X, obj.J[1:obj.m])
    elseif sym === :w
        return view(obj.W, 1:obj.m)
    else 
        return getfield(obj, sym)
    end
end

(a::MovingAAAapprox)(zz) = bary(zz, a)
(a::MovingAAAapprox)(zz::T) where {T <: AbstractVector} = bary.(zz, a)

function movingAAAapprox(Z, F; tol=1e-13, mmax=150,
             verbose=false, clean=false, do_sort=true, cleanup_tol=1e-13)

    U = eltype(Z)
    S = eltype(F)
    # filter out any NaN's or Inf's in the input
    keep = isfinite.(F)
    F = F[keep]
    Z = Z[keep]

    # Remove repeated elements of Z and corresponding elements of F
    ii = unique(i -> Z[i], eachindex(Z))
    Z = Z[ii]
    F = F[ii]

    M = length(Z)                    # number of sample points
    F, Z = promote(F, Z)
    T = promote_type(S, U)

    # call usual aaa
    r = aaa(Z, F; tol=tol, mmax=mmax, verbose=verbose, clean=clean, do_sort=do_sort, cleanup_tol=cleanup_tol )

    m = length(r.x) # number of points used on the interpolation  

    A  = Matrix{T}(undef, M, m)
    C  = Matrix{T}(undef, M, m)
    J = zeros(Int,M)
    j_notsupport = m


    for ii in eachindex(Z) 
        z = Z[ii]
        idx_0, idx_1 = _nearby(z, r.x)
        if isnothing(idx_1)
            # the point is on the support
            J[idx_0] = ii
            # next column of Cauchy matrix
            C[:, idx_0] = T(1) ./ (Z .- z)

            # Loewner matrix
            A[:, idx_0] = (F .- F[ii]) .* C[:, idx_0]
        else
            # the point is not on the support
            j_notsupport += 1
            J[j_notsupport] = ii
        end
    end

    w = compute_weights(m, J[m+1:end], A)


    return MovingAAAapprox(m, 
        Z,
        F,
        w,
        A,
        J
        )
end


function update_aaa!(a::MovingAAAapprox, Z, F)
    # assume that z has the same order than a.X
    m = a.m
    _j = a.J[1:m]
    for ii in eachindex(Z)
        a.X[ii] = Z[ii]
        a.F[ii] = F[ii]
        k = findfirst(==(ii), _j)
        if !isnothing(k)
            # this is a support point so actualize stuff
            a.A[:,k] = (F .- F[ii]) ./ (Z .- Z[ii])
        end
    end
    a.W .= compute_weights(m, a.J[m+1:end], a.A)
    return nothing
end

# TODO
# 1. see fit of gaussian, polynomial (data only in x, on x y, which point is more probable)
#       pol fit is nice with few points, gaussian not, it uses the one that is further from mean(f)
# 1. see if I can add and delete elements of vector on struct
#       I can
#
# 1. do not move the points order on movingAAAapprox, use array permutations and views
    # i. we need J vector on the struct 
        #  done but maybe change A so the W is the same as in aaa
# 1. maybe check bary and derivative functions so they use each index and are compatible with movingAAAapprox 
#       maybe not necesary, seems to have good loops?
# 1. re fit of movingAAAapprox
