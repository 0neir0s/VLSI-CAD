#!/usr/bin/env python

"""An implementation of the Quine McCluskey algorithm.

This implementation of the Quine McCluskey algorithm has no inherent limits
(other than the calculation time) on the size of the inputs.

Also, in the limited tests of the author of this module, this implementation is
considerably faster than other public Python implementations for non-trivial
inputs.

Another unique feature of this implementation is the possibility to use the XOR
and XNOR operators, in addition to the normal AND operator, to minimise the
terms. This slows down the algorithm, but in some cases it can be a big win in
terms of complexity of the output.
"""

from __future__ import print_function
import math

class QuineMcCluskey:
	"""The Quine McCluskey class.

	The QuineMcCluskey class minimises boolean functions using the Quine
	McCluskey algorithm.

	If the class was instantiiated with the use_xor set to True, then the
	resulting boolean function may contain XOR and XNOR operators.
	"""

	def __init__(self,variables):
		"""The class constructor"""
		self.variables = variables	 
		self.numvars = len(variables)
		self.n_bits = 0         # number of bits (i.e. self.n_bits == len(ones[i]) for every i).

	def updateVariables(self,variables):
		"""The class constructor"""
		self.variables = variables	 
		self.numvars = len(variables)
		self.n_bits = 0         # number of bits (i.e. self.n_bits == len(ones[i]) for every i).

	def __num2str(self, i):
		"""
		Convert an integer to its bit-representation in a string.

		Args:
		    i (int): the number to convert.

		Returns:
		    The binary string representation of the parameter i.
		"""
		x = ['1' if i & (1 << k) else '0' for k in range(self.n_bits - 1, -1, -1)]
		return "".join(x)



	def simplify(self, ones, dc = []):
		"""Simplify a list of terms.

		Args:
		    ones (list of int): list of integers that describe when the output
		    function is '1', e.g. [1, 2, 6, 8, 15].

		Kwargs:
		    dc (list of int): list of numbers for which we don't care if they
		    have one or zero in the output.

		Returns:
		    see: simplify_los.

		Example:
		    ones = [2, 6, 10, 14]
		    dc = []

		    This will produce the ouput: ['--10']
		    This means x = b1 & ~b0, (bit1 AND NOT bit0)

		Example:
		    ones = [1, 2, 5, 6, 9, 10, 13, 14]
		    dc = []

		    This will produce the ouput: ['--^^'].
		    In other words, x = b1 ^ b0, (bit1 XOR bit0).
		"""
		terms = ones + dc
		if len(terms) == 0:
		    return None

		# Calculate the number of bits to use
		# Needed internally by __num2str()
		self.n_bits = int(math.ceil(math.log(max(terms) + 1, 2)))

		# Generate the sets of ones and dontcares
		ones = set(self.__num2str(i) for i in ones)
		dc = set(self.__num2str(i) for i in dc)

		return self.simplify_los(ones, dc)



	def simplify_los(self, ones, dc = []):
		"""The simplification algorithm for a list of string-encoded inputs.

		Args:
			ones (list of str): list of strings that describe when the output
			function is '1', e.g. ['0001', '0010', '0110', '1000', '1111'].

		Kwargs:
			dc: (list of str)set of strings that define the don't care
			combinations.

		Returns:
			Returns a set of strings which represent the reduced minterms.
			'-' don't care: this bit can be either zero or one.
			'1' the bit must be one.
			'0' the bit must be zero.

		Example:
			ones = ['0010', '0110', '1010', '1110']
			dc = []

			This will produce the ouput: ['--10'].
			In other words, x = b1 & ~b0, (bit1 AND NOT bit0).
		"""
		self.profile_cmp = 0    # number of comparisons (for profiling)

		terms = ones | dc
		if len(terms) == 0:
			return None

		# Calculate the number of bits to use
		self.n_bits = max(len(i) for i in terms)
		if self.n_bits != min(len(i) for i in terms):
			return None

		# First step of Quine-McCluskey method : prime implicants
		prime_implicants = self.__get_prime_implicants(terms)
		# Convert prime implicants to a proper form (value+mask)
		prime_implicants = self.__get_prime_tuples(prime_implicants)
		return prime_implicants
	
	def __get_prime_tuples(self,primes):
		prime_tuples = list()
		for prime in primes:
			value = int(prime.replace('-','0'),2)
			mask = int(prime.replace('1','0').replace('-','1'),2)
			prime_tuples.append((value,mask))
		prime_tuples = set(prime_tuples)
		return prime_tuples

	def possible_covers(self, primes, ones):
		"""
		Use the prime implicants to find the essential prime implicants of the
		function, as well as other prime implicants that are necessary to cover
		the function. This method uses the Petrick's method, which is a technique
		for determining all minimum sum-of-products solutions from a prime implicant
		chart.

		primes: the prime implicants that we want to minimize.
		ones: a list of indices for the minterms for which we want the function to
		evaluate to 1.
		"""
		if type(ones[0]) is int:
			ones = set(self.__num2str(i) for i in ones)
		chart = []
		for one in ones:
			column = []
			for i in xrange(len(primes)):
				if (int(one,2) & (~primes[i][1])) == primes[i][0]:
					column.append(i)
			chart.append(column)

		covers = []
		if len(chart) > 0:
			covers = [set([i]) for i in chart[0]]
		for i in xrange(1,len(chart)):
			new_covers = []
			for cover in covers:
				for prime_index in chart[i]:
					x = set(cover)
					x.add(prime_index)
					append = True
					for j in xrange(len(new_covers)-1,-1,-1):
						if x <= new_covers[j]:
							del new_covers[j]
						elif x > new_covers[j]:
							append = False
					if append:
						new_covers.append(x)
			covers = new_covers
		return covers

	def petricks_method(self, primes, ones):
		"""
		Use the prime implicants to find the essential prime implicants of the
		function, as well as other prime implicants that are necessary to cover
		the function. This method uses the Petrick's method, which is a technique
		for determining all minimum sum-of-products solutions from a prime implicant
		chart.

		primes: the prime implicants that we want to minimize.
		ones: a list of indices for the minterms for which we want the function to
		evaluate to 1.
		"""

		chart = []
		for one in ones:
			column = []
			for i in xrange(len(primes)):
				if (int(one,2) & (~primes[i][1])) == primes[i][0]:
					column.append(i)
			chart.append(column)

		covers = []
		if len(chart) > 0:
			covers = [set([i]) for i in chart[0]]
		for i in xrange(1,len(chart)):
			new_covers = []
			for cover in covers:
				for prime_index in chart[i]:
					x = set(cover)
					x.add(prime_index)
					append = True
					for j in xrange(len(new_covers)-1,-1,-1):
						if x <= new_covers[j]:
							del new_covers[j]
						elif x > new_covers[j]:
							append = False
					if append:
						new_covers.append(x)
			covers = new_covers

		min_complexity = 99999999
		for cover in covers:
			primes_in_cover = [primes[prime_index] for prime_index in cover]
			complexity = self.calculate_complexity(primes_in_cover)
			if complexity < min_complexity:
				min_complexity = complexity
				result = primes_in_cover
		return min_complexity,result

	
	def calculate_complexity(self, minterms):
		"""
		Calculate the complexity of the given function. The complexity is calculated
		based on the following rules:
		A NOT gate adds 1 to the complexity.
		A n-input AND or OR gate adds n to the complexity.

		minterms: a list of minterms that form the function

		returns: an integer that is the complexity of the function

		>>> qm = QM(['A','B','C'])

		>>> qm.calculate_complexity([(1,6)])
		0
		>>> qm.calculate_complexity([(0,6)])
		1
		>>> qm.calculate_complexity([(3,4)])
		2
		>>> qm.calculate_complexity([(7,0)])
		3
		>>> qm.calculate_complexity([(1,6),(2,5),(4,3)])
		3
		>>> qm.calculate_complexity([(0,6),(2,5),(4,3)])
		4
		>>> qm.calculate_complexity([(0,6),(0,5),(4,3)])
		5
		>>> qm.calculate_complexity([(0,6),(0,5),(0,3)])
		6
		>>> qm.calculate_complexity([(3,4),(7,0),(5,2)])
		10
		>>> qm.calculate_complexity([(1,4),(7,0),(5,2)])
		11
		>>> qm.calculate_complexity([(2,4),(7,0),(5,2)])
		11
		>>> qm.calculate_complexity([(0,4),(7,0),(5,2)])
		12
		>>> qm.calculate_complexity([(0,4),(0,0),(5,2)])
		15
		>>> qm.calculate_complexity([(0,4),(0,0),(0,2)])
		17
		"""

		complexity = len(minterms)
		if complexity == 1:
			complexity = 0
		mask = (1<<self.numvars)-1
		for minterm in minterms:
			masked = ~minterm[1] & mask
			term_complexity = bitcount(masked)
			if term_complexity == 1:
				term_complexity = 0
			complexity += term_complexity
			complexity += bitcount(~minterm[0] & masked)
		return complexity

	def __get_prime_implicants(self, terms):
		"""Simplify the set 'terms'.

		Args:
		    terms (set of str): set of strings representing the minterms of
		    ones and dontcares.

		Returns:
		    A list of prime implicants. These are the minterms that cannot be
		    reduced with step 1 of the Quine McCluskey method.

		This is the very first step in the Quine McCluskey algorithm. This
		generates all prime implicants, whether they are redundant or not.
		"""

		# Sort and remove duplicates.
		n_groups = self.n_bits + 1
		marked = set()

		# Group terms into the list groups.
		# groups is a list of length n_groups.
		# Each element of groups is a set of terms with the same number
		# of ones.  In other words, each term contained in the set
		# groups[i] contains exactly i ones.
		groups = [set() for i in range(n_groups)]
		for t in terms:
			n_bits = t.count('1')
			groups[n_bits].add(t)

		done = False
		while not done:
		    # Group terms into groups.
		    # groups is a list of length n_groups.
		    # Each element of groups is a set of terms with the same
		    # number of ones.  In other words, each term contained in the
		    # set groups[i] contains exactly i ones.
			groups = dict()
			for t in terms:
				n_ones = t.count('1')
				key = n_ones
				if key not in groups:
					groups[key] = set()
				groups[key].add(t)

			terms = set()           # The set of new created terms
			used = set()            # The set of used terms

			# Find prime implicants
			for key in groups:
				key_next = key+1
				if key_next in groups:
					group_next = groups[key_next]
					for t1 in groups[key]:
						# Optimisation:
						# The Quine-McCluskey algorithm compares t1 with
						# each element of the next group. (Normal approach)
						# But in reality it is faster to construct all
						# possible permutations of t1 by adding a '1' in
						# opportune positions and check if this new term is
						# contained in the set groups[key_next].
						for i, c1 in enumerate(t1):
							if c1 == '0':
								self.profile_cmp += 1
								t2 = t1[:i] + '1' + t1[i+1:]
								if t2 in group_next:
									t12 = t1[:i] + '-' + t1[i+1:]
									used.add(t1)
									used.add(t2)
									terms.add(t12)

			# Add the unused terms to the list of marked terms
			for g in list(groups.values()):
				marked |= g - used

			if len(used) == 0:
				done = True

		# Prepare the list of prime implicants
		pi = marked
		for g in list(groups.values()):
			pi |= g
		return pi

	def permutations(self, value = ''):
		"""Iterator to generate all possible values out of a string.

		Args:
		    value (str): A string containing any of the above characters.

		Returns:
		    The output strings contain only '0' and '1'.

		Example:
		    from qm import QuineMcCluskey
		    qm = QuineMcCluskey()
		    for i in qm.permutations('1--^^'):
		        print(i)

		The operation performed by this generator function can be seen as the
		inverse of binary minimisation methonds such as Karnaugh maps, Quine
		McCluskey or Espresso.  It takes as input a minterm and generates all
		possible maxterms from it.  Inputs and outputs are strings.

		Possible input characters:
		    '0': the bit at this position will always be zero.
		    '1': the bit at this position will always be one.
		    '-': don't care: this bit can be zero or one.

		Algorithm description:
		    This lovely piece of spaghetti code generates all possibe
		    permutations of a given string describing logic operations.
		    This could be achieved by recursively running through all
		    possibilities, but a more linear approach has been preferred.
		    The basic idea of this algorithm is to consider all bit
		    positions from 0 upwards (direction = +1) until the last bit
		    position. When the last bit position has been reached, then the
		    generated string is yielded.  At this point the algorithm works
		    its way backward (direction = -1) until it finds an operator
		    like '-', '^' or '~'.  The bit at this position is then flipped
		    (generally from '0' to '1') and the direction flag again
		    inverted. This way the bit position pointer (i) runs forth and
		    back several times until all possible permutations have been
		    generated.
		    When the position pointer reaches position -1, all possible
		    combinations have been visited.
		"""
		n_bits = len(value)
		res = ['0' for i in range(n_bits)]
		i = 0
		direction = +1
		while i >= 0:
			# binary constant
			if value[i] == '0' or value[i] == '1':
				res[i] = value[i]
			# dontcare operator
			elif value[i] == '-':
				if direction == +1:
					res[i] = '0'
				elif res[i] == '0':
					res[i] = '1'
					direction = +1
			# unknown input
			else:
				res[i] = '#'

			i = i + direction
			if i == n_bits:
				direction = -1
				i = n_bits - 1
				yield "".join(res)
	def get_function(self, minterms):
		"""
		Return in human readable form a sum of products function.

		minterms: a list of minterms that form the function

		returns: a string that represents the function using operators AND, OR and
		NOT.
		"""

		if isinstance(minterms,str):
			return minterms

		def parentheses(glue, array):
			if len(array) > 1:
				return ''.join([glue.join(array)])
			else:
				return glue.join(array)

		or_terms = []
		for minterm in minterms:
			and_terms = []
			for j in xrange(len(self.variables)):
				if minterm[0] & 1<<j:
					and_terms.append(self.variables[j])
				elif not minterm[1] & 1<<j:
					and_terms.append('%s\'' % self.variables[j])
			or_terms.append(parentheses('.', and_terms))
		return parentheses('+', or_terms)

def bitcount(i):
	""" Count set bits of the input. """
	res = 0
	while i > 0:
		res += i&1
		i>>=1
	return res

qm = QuineMcCluskey(['A','B','C'])
ones1 = [2,3,4,5,6,7]
ones2 = [0,1,2,3,4,5]

# Convert prime implicants to a proper form (value+mask)
prime_implicants1 = list(qm.simplify(ones1,[]))
covers1 = qm.possible_covers(prime_implicants1,ones1)
prime_implicants2 = list(qm.simplify(ones2,[]))
covers2 = qm.possible_covers(prime_implicants2,ones2)	

min_complexity = 99999999
for cover1 in covers1:
	for cover2 in covers2:
			implicant1 = [prime_implicants1[i] for i in list(cover1)]
			implicant2 = [prime_implicants2[i] for i in list(cover2)]
			cover = list(set().union(implicant1,implicant2))
			complexity = qm.calculate_complexity(cover)
			print(cover)
			if complexity < min_complexity:
				min_complexity = complexity
				result = cover
print(min_complexity,result)
'''
		# Second step of Quine McCluskey method : prime implicant chart and Petrick's Method.
		final_terms = qm.petricks_method(list(prime_implicants),ones)
		#print(final_terms)
'''