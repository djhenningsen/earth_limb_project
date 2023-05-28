def hstdlondt(lat):
    a=2.00778275e-05
    b=3.45996369e-07
    c=5.06082712e-02

    return a*lat**2+b*lat+c #degrees longitude per second
