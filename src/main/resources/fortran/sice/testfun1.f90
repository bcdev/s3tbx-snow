program testfun1

    EXTERNAL fun1, fun2
    common as, bs, cs, csza, cvza, r400, r865, r1020
    common xa(168), ya(168), NSOLO, thv

    open(898, file = 'ice_index.dat')
    do 67 j = 1, 168
        read(898, *) xa(j), an, ya(j)
    67       continue

    NSOLO = 1
    thv = 0.97

    x1 = 0.4425
    x2 = 0.70875
    x3 = 1.020

    r400 = 0.95
    r865 = 0.7
    r1020 = 0.45

    csza = 0.5
    cvza = 0.8

    y1 = 0.8
    y2 = 0.6
    y3 = 0.4

    d1 = (y3 - y1) * (x2 - x1) - (y2 - y1) * (x3 - x1)
    d2 = (x3 * x3 - x1 * x1) * (x2 - x1) - (x2 * x2 - x1 * x1) * (x3 - x1)

    !             second order polynomial coefficients  for planar albedo:
    as = d1 / d2
    bs = (y3 - y2 - as * (x3 * x3 - x2 * x2)) / (x3 - x2)
    cs = y3 - as * x3 * x3 - bs * x3

    !             limits of integration
    at = 0.3
    bt = 2.4
    aat = 0.7

    x = 0.4
    result = fun1(x)

    write(*,*) "result: ", result

    as = -1.6304738716662046
    bs = 1.6587701807575812
    cs = 0.3699402668538485
    r400 = 0.78
    csza = 0.7395136208713016
    x = 1.54
    result = fun1(x)

    write(*,*) "result: ", result

end

real function fun1(x)

    real dif(168)
    integer kss(2)
    common as, bs, cs, csza, cvza, r400, r865, r1020
    common xa(168), ya(168), NSOLO, thv

    ks = 1

    pi = acos(-1.)

    alam3 = 0.865
    alam4 = 1.02
    akap3 = 2.4e-7
    akap4 = 2.25e-6

    alpha3 = 4. * pi * akap3 / alam3
    alpha4 = 4. * pi * akap4 / alam4

    eps = sqrt(alpha3 / alpha4)
    ax1 = 1. / (1. - eps)
    ax2 = 1. / (1. - 1. / eps)

    !r0 = r865**ax1 * r1020**ax2
    r0 = 0.9092759022694742

    um1 = 3. * (1. + 2. * csza) / 7.
    um2 = 3. * (1. + 2. * cvza) / 7.

    xx = um1 * um2 / r0

    !al = (alog(r1020 / r0))**2. / xx / xx / alpha4
    al = 15517.04843100975;

    write(*,*) "as, bs, cs: ", as, bs, cs
    write(*,*) "csza, cvza: ", csza, cvza
    write(*,*) "rrr: ", r400, r865, r1020
    write(*,*) "xx: ", xx
    write(*,*) "thv: ", thv
    write(*,*) "al: ", al

    if (x.lt.0.4)astra = 2.e-11

    do 68 k = 1, 168

        dif(k) = abs(x - xa(k))
        write(*,*) "dif: ", k, x, xa(k), dif(k)

        if (dif(k).le.1.e-6) ks = k
        if (dif(k).le.1.e-6) go to 69
    68                continue

    kss(1) = minloc(dif, 1)
    l = kss(1)

    delta = x - xa(l)
    if (delta.le.0) LS1 = L - 1
    if (delta.le.0) LS2 = L
    if (delta.ge.0) LS1 = L
    if (delta.ge.0) LS2 = L + 1
    x0 = xa(ls1)
    x1 = xa(ls2)
    y0 = ya(ls1)
    y1 = ya(ls2)
    y = y0 + (x - x0) * (y1 - y0) / (x1 - x0)
    go to 18
    69   continue
    y = ya(ks)
    18         continue

    ASTRA = y
    write(*,*) "ks: ", ks, y, ya(ks)
    write(*,*) "astra: ", ASTRA

    dega = al * 4. * pi * astra / x
    pow = sqrt(dega)
    if (pow.gt.1.e-6)rsd = exp(-pow)
    if (pow.le.1.e-6) rsd = 1.

    if (NSOLO.eq.0) f1 = rsd**um1
    if (NSOLO.eq.1) f1 = rsd

    write(*,*) "f1: ", f1

    if (r400.le.thv.and.x.le.1.02) f1 = as * x * x + bs * x + cs

    p0 = 32.38
    p1 = -160140.33
    p2 = 7959.53

    t1 = 85.34 * 1.e-3
    t2 = 401.79 * 1.e-3
    funcs = p0 + p1 * exp(-x / t1) + p2 * exp(-x / t2)
    if (x.lt.0.4) funcs = p0 + p1 * exp(-0.4 / t1) + p2 * exp(-0.4 / t2)
    write(*,*) "funcs: ", funcs

    fun1 = f1 * funcs
    write(*,*) "fun1: ", fun1

    return
END


