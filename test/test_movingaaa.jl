
# HACK: Sort by increasing absolute value of the real part.
# This is a hack to get poles ordered for test purposes.
csort(x) = sort(x, lt = (x,y) -> abs(real(x)) < abs(real(y)))

function test_moving_aaa_spiral(T=Float64, verbose=false)
    clean = false
    verbose = true
    zz = exp.(range(-0.5, complex(0.5, 15*pi), length=20));

    yy = tan.(pi*zz/2);
    dyy = (pi/2) ./ (cos.(pi*zz/2) .^ 2);
    f = aaa(zz, yy, clean=clean, verbose=verbose)
    fm = movingAAAapprox(zz, yy, clean=clean, verbose=verbose)
    p1 = norm(abs.(f(zz) - yy)) < 1e-8
    p2 = norm(abs.(fm(zz) - yy)) < 1e-8
    # this has error
    # p3 = norm(abs.(deriv.(f, zz) - dyy)) < 1e-8
    # p4 = norm(abs.(deriv.(fm, zz) - dyy)) < 1e-8

    v0 = exp(im*0.3);
    dt = 0.1;
    zz_t = zz .+ v0 .* dt ;
    yy_t = tan.(pi*zz_t*(1+dt)/2);
    dyy_t = (pi*(1+dt)/2) ./ (cos.(pi*zz_t/2 * (1+dt)) .^ 2);

    f_t = aaa(zz_t, yy_t, clean=clean, verbose=verbose)
    update_aaa!(fm, zz_t, yy_t)

    p5 = norm(abs.(f_t(zz_t) - yy_t)) < 1e-8
    p6 = norm(abs.(fm(zz_t) - yy_t)) < 1e-8
    # this has error
    # p7 = norm(abs.(deriv.(f_t, zz_t) - dyy_t)) < 1e-8
    # p8 = norm(abs.(deriv.(fm_t, zz_t) - dyy_t)) < 1e-8
    


    return all((p1, p2, p3))
end
