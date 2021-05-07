from math import pi, gcd
from ket import plugins
from ket.lib import qft

def period():
    reg1 = quant(4)
    h(reg1)
    reg2 = plugins.pown(7, reg1, 15)
    qft(reg1)
    return measure(reg1).get()

r = period()
results = [r]
for _ in range(4):
    results.append(period())
    r = gcd(r, results[-1])

print(results)
r = 2**4/r

print('measurements =', results)
print('r =', r)
p = gcd(int(7**(r/2))+1, 15)
q = gcd(int(7**(r/2))-1, 15)
print(15, '=', p , "x", q)