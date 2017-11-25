#!/usr/bin/env python

"Library code"
from __future__ import print_function
import math


class QuineMcCluskey:
	def __init__(self,variables):
		"""The class constructor"""
		self.variables = variables	 
		self.numvars = len(variables)
		self.n_bits = 0         

	def updateVariables(self,variables):
		"""The class constructor"""
		self.variables = variables	 
		self.numvars = len(variables)
		self.n_bits = 0         

	def __num2str(self, i):
		"""Convert an integer to its bit-representation in a string"""
		x = ['1' if i & (1 << k) else '0' for k in range(self.n_bits - 1, -1, -1)]
		return "".join(x)



	def minimize(self, ones, dc = []):
		"""Simplify a list of term_mins"""
		term_mins = ones + dc
		if len(term_mins) == 0:
		    return None

		# Calculate the number of bits to use
		self.n_bits = int(math.ceil(math.log(max(term_mins) + 1, 2)))
		if(term_mins == [0]):
			self.n_bits = 1
		# Generate the sets of ones and dontcares
		ones = set(self.__num2str(i) for i in ones)
		dc = set(self.__num2str(i) for i in dc)

		return self.minimize_los(ones, dc)



	def minimize_los(self, ones, dc = []):
		"""The main algorithm takes place here"""
		self.global_prof = 0   
		term_mins = ones | dc
		if len(term_mins) == 0:
			return None

		# Calculate the number of bits to use
		self.n_bits = max(len(i) for i in term_mins)
		if self.n_bits != min(len(i) for i in term_mins):
			return None
		# First step of Quine-McCluskey method : prime implicants
		prime_implicants = self.__find_PI(term_mins)
		# Convert prime implicants to a proper form (value+mask)
		prime_implicants = self.__get_prime_tuples(prime_implicants)
		# Second step of Quine McCluskey method : prime implicant chart and Petrick's Method.
		final_term_mins = self.__petricks_method(list(prime_implicants),ones)
		return list(prime_implicants),final_term_mins
	
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
		Use the prime implicants to find the exhastive covers
		the function.

		primes: the prime implicants that we want to minimize.
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

		#this loop gets a list of the covers that can be used to represent the function
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

	def __petricks_method(self, primes, ones):
		"""Find the essential prime implicants reusing the cover function written above"""

		chart = []
		print(primes)
		for one in ones:
			column = []
			for i in xrange(len(primes)):
				if (int(one,2) & (~primes[i][1])) == primes[i][0]:
					column.append(i)
			chart.append(column)

		covers = []
		if len(chart) > 0:
			covers = [set([i]) for i in chart[0]]
		#this loop iterates over the covers and applies the petrick's method to them
		print(chart)
		for i in xrange(1,len(chart)):
			new_covers = []
			for cover in covers:
				for prime_index in chart[i]:
					print(chart[i])
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
			print('Covers:',covers)
			print('--------------------------------------------------------')
		# Now we have to find the most efficient cover among the ones we've found
		# The idea is to find the min cover adhering to the standard cost rules
		min_cost = 99999999
		for cover in covers:
			primes_in_cover = [primes[prime_index] for prime_index in cover]
			cost = self.calculate_cost(primes_in_cover)
			if cost < min_cost:
				min_cost = cost
				result = primes_in_cover
		return result

	
	def calculate_cost(self, minterm_mins):
		"""
		Calculates cost of the given cover of the function.
		A NOT gate adds 1 to the cost.
		AND or OR gate adds 2 to the cost.
		"""

		cost = len(minterm_mins)
		if cost == 1:
			cost = 0
		mask = (1<<self.numvars)-1
		for minterm in minterm_mins:
			masked = ~minterm[1] & mask
			term_cost = no_of_bits(masked)
			if term_cost == 1:
				term_cost = 0
			cost += term_cost
			cost += no_of_bits(~minterm[0] & masked)
		return cost

	def __find_PI(self, term_mins):
		"""
		Returns prime implicants
		"""
		n_groups = self.n_bits + 1
		marked = set()
		# Group term_mins into the list groups.
		groups = [set() for i in range(n_groups)]
		for t in term_mins:
			n_bits = t.count('1')
			groups[n_bits].add(t)

		done = False
		while not done:
		    # Group term_mins into groups.
			groups = dict()
			for t in term_mins:
				n_ones = t.count('1')
				key = n_ones
				if key not in groups:
					groups[key] = set()
				groups[key].add(t)
			term_mins = set()           # The set of new created term_mins
			used = set()            # The set of used term_mins

			# Find prime implicants
			for key in groups:
				key_next = key+1
				if key_next in groups:
					group_next = groups[key_next]
					for t1 in groups[key]:
						for i, c1 in enumerate(t1):
							if c1 == '0':
								self.global_prof += 1
								t2 = t1[:i] + '1' + t1[i+1:]
								if t2 in group_next:
									t12 = t1[:i] + '-' + t1[i+1:]
									used.add(t1)
									used.add(t2)
									term_mins.add(t12)

			# Add the unused term_mins to the list of marked term_mins
			for g in list(groups.values()):
				marked |= g - used

			if len(used) == 0:
				done = True

		# Prepare the list of prime implicants
		pi = marked
		for g in list(groups.values()):
			pi |= g
		return pi

	def readable_format(self, minterm_mins):
		"""function in readable form given the essential prime implicants as the inputs"""

		if isinstance(minterm_mins,str):
			return minterm_mins

		def parentheses(glue, array):
			if len(array) > 1:
				return ''.join([glue.join(array)])
			else:
				return glue.join(array)

		or_term_mins = []
		for minterm in minterm_mins:
			and_term_mins = []
			for j in xrange(len(self.variables)):
				if minterm[0] & 1<<j:
					and_term_mins.append(self.variables[j])
				elif not minterm[1] & 1<<j:
					and_term_mins.append('%s\'' % self.variables[j])
			or_term_mins.append(parentheses('.', and_term_mins))
		return parentheses(' + ', or_term_mins)

def no_of_bits(i):
	""" Count set bits of the input. """
	res = 0
	while i > 0:
		res += i&1
		i>>=1
	return res