from basic import QuineMcCluskey
from QM import QM

# qm = QuineMcCluskey()
# ones = [2,3,7,9,11,13]
# dontcares = [1,10,15]
# print(qm.simplify(ones, dontcares))

q = QM(['A','B','C','D'])
ones = [2,3,7,9,11,13]
dontcares = [1,10,15]
print q.get_function(q.solve(ones,dontcares)[1])
