from basic import *


qm = QuineMcCluskey(['A','B','C','D'])
ones = [2,6,8,9,10,11,14,15]
dontcares = []
primes = qm.simplify(ones, dontcares)[1]
print(qm.get_function(primes))
