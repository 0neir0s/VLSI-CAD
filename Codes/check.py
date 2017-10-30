from basic import *


qm = QuineMcCluskey(['A','B','C','D'])
ones = [4,8,10,11,12,15]
dontcares = [9,14]
primes = qm.simplify(ones, dontcares)[1]
print(qm.get_function(primes))
