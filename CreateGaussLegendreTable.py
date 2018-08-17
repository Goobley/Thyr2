import numpy.polynomial.legendre as leg

GauLegDegree = 20
x, w = leg.leggauss(GauLegDegree)

with open('GaussLegendreTable.ixx', 'w') as f:
    f.write('namespace GaussLegendre\n{\n')
    f.write('    constexpr f64 x[20] = {')
    for i, xVal in enumerate(x):
        f.write('%.17f' % xVal) 
        if i < len(x) - 1:
            f.write(',')
    f.write('};\n')

    f.write('    constexpr f64 w[20] = {')
    for i, wVal in enumerate(w):
        f.write('%.17f' % wVal) 
        if i < len(w) - 1:
            f.write(',')
    f.write('};\n')
    f.write('}')
