# so x can be real while f and w are complex
struct MovingAAAapprox{T,W} <: BRInterp
    m::Int # m <= M amount of support points
    x::AbstractVector{T} # support points, fixed size m
    f::AbstractVector{W} # function values at support points, size m
    w::AbstractVector{W} # baricentric weights, size m
    A::AbstractMatrix{W} # M × m Loewner matrix 
    _j::AbstractVector{Int} # permutation vector, size M
end


(a::MovingAAAapprox)(zz) = bary(zz, a)
(a::MovingAAAapprox)(zz::T) where {T <: AbstractVector} = bary.(zz, a)

function movingAAAapprox(Z, F; tol=1e-13, mmax=150,
             verbose=false, clean=0, do_sort=true, cleanup_tol=1e-13)

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
    J = zeros(Int,M)
    j_notsupport = m


    for ii in eachindex(Z) 
        z = Z[ii]
        idx_0 = findfirst(==(ii), r._j)
        if !isnothing(idx_0)
            # the point is on the support
            J[idx_0] = ii

            # Loewner matrix
            A[:, idx_0] = (F .- F[ii]) ./ (Z .- z)
        else
            # the point is not on the support
            j_notsupport += 1
            J[j_notsupport] = ii
        end
    end

    w = compute_weights(m, J[m+1:end], A)


    return MovingAAAapprox(m, 
        Z[J[1:m]],
        F[J[1:m]],
        w,
        A,
        J
        )
end


# note that is not standard, we some times modify a and return nothing and 
# some times modify a and return the real response
function update_movingaaa!(a::MovingAAAapprox, Z, F; tol=1e-13, mmax=150,
             verbose=false, clean=0, do_sort=true, cleanup_tol=1e-13)
    # assume that z has the same order than a.X
    m = a.m
    _j = a._j[1:m]
    for ii in eachindex(Z)
        k = findfirst(==(ii), _j)
        if !isnothing(k)
            a.x[k] = Z[ii]
            a.f[k] = F[ii]
            # this is a support point so actualize stuff
            a.A[:,k] = (F .- F[ii]) ./ (Z .- Z[ii])
        end
    end
    a.w .= compute_weights(m, a._j[m+1:end], a.A)
    
    # compute error
    abstol = tol * norm(F, Inf)
    err = norm(a(Z[a._j[m+1:end]]) - F[a._j[m+1:end]])
    if err > abstol
        verbose && println("err = ",err," bigger than tolerance, calling aaa algorihtm")
        r = movingAAAapprox(Z, F; tol=tol, mmax=mmax,
             verbose=verbose, clean=clean, do_sort=do_sort, cleanup_tol=cleanup_tol)
        return r
    else
        verbose && println("err = ",err," less than tolerance, not calling aaa algorihtm")
    end

    return nothing
end
