from basic import QuineMcCluskey

qm = QuineMcCluskey()
ones = [4,8,10,11,12,15]
dontcares = [9,14]
print(qm.simplify(ones, dontcares))