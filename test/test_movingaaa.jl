function test_moving_aaa_spiral(T=Float64, verbose=false)
    clean = false
    # verbose = true
    zz = exp.(range(-0.5, complex(0.5, 15*pi), length=20));

    yy = tan.(pi*zz/2);
    dyy = (pi/2) ./ (cos.(pi*zz/2) .^ 2);
    f = aaa(zz, yy, clean=clean, verbose=verbose)
    fm = movingAAAapprox(zz, yy, clean=clean, verbose=verbose)
    p1 = norm(abs.(f(zz) - yy)) < 1e-8
    p2 = norm(abs.(fm(zz) - yy)) < 1e-8
    # this has error
    # p1_d = norm(abs.(deriv.(f, zz) - dyy)) < 1e-8
    # p2_d = norm(abs.(deriv.(fm, zz) - dyy)) < 1e-8

    v0 = exp(im*0.3);
    dt = 0.5;
    zz_t = zz .+ v0 .* dt ;
    yy_t = tan.(pi*zz_t*(1+dt)/2);
    dyy_t = (pi*(1+dt)/2) ./ (cos.(pi*zz_t/2 * (1+dt)) .^ 2);

    f_t = aaa(zz_t, yy_t, clean=clean, verbose=verbose)
    r = update_movingaaa!(fm, zz_t, yy_t; verbose=verbose, clean=clean)
    if !isnothing(r)
        fm = r
    end

    p3 = norm(abs.(f_t(zz_t) - yy_t)) < 1e-8
    p4 = norm(abs.(fm(zz_t) - yy_t)) < 1e-8
    # this has error
    # p3_d = norm(abs.(deriv.(f_t, zz_t) - dyy_t)) < 1e-8
    # p4_d = norm(abs.(deriv.(fm_t, zz_t) - dyy_t)) < 1e-8
    


    return all((p1, p2, p3, p4))
end

function test_moving_aaa_refit_aaa(T=Float64, verbose=false)
    clean = false
    # verbose = true
    zz = exp.(range(-0.5, complex(0.5, 15*pi), length=20));

    yy = (pi*zz/2);
    dyy = (pi/2);
    f = aaa(zz, yy, clean=clean, verbose=verbose)
    fm = movingAAAapprox(zz, yy, clean=clean, verbose=verbose)
    p1 = norm(abs.(f(zz) - yy)) < 1e-8
    p2 = norm(abs.(fm(zz) - yy)) < 1e-8
    # this has error
    # p1_d = norm(abs.(deriv.(f, zz) - dyy)) < 1e-8
    # p2_d = norm(abs.(deriv.(fm, zz) - dyy)) < 1e-8

    v0 = exp(im*0.3);
    dt = 0.5;
    zz_t = zz .+ v0 .* dt ;
    # this is a completely different function so we need to  re fit
    yy_t = cos.(pi*zz_t*(1+dt)/2);
    dyy_t = (pi*(1+dt)/2) * (-sin.(pi*zz_t/2 * (1+dt)));

    f_t = aaa(zz_t, yy_t, clean=clean, verbose=verbose)
    r = update_movingaaa!(fm, zz_t, yy_t; verbose=verbose, clean=clean)
    if !isnothing(r)
        fm = r
    end

    p3 = norm(abs.(f_t(zz_t) - yy_t)) < 1e-8
    p4 = norm(abs.(fm(zz_t) - yy_t)) < 1e-8
    # this has error
    # p3_d = norm(abs.(deriv.(f_t, zz_t) - dyy_t)) < 1e-8
    # p4_d = norm(abs.(deriv.(fm_t, zz_t) - dyy_t)) < 1e-8
    


    return all((p1, p2, p3, p4))
end
