class monomial:
	def __init__(self, coeff, degree):
		self.degree = degree
		self.coeff = coeff

	def __str__(self):
		if self.coeff == 1:
			return "x^{}".format(self.degree)
		elif self.coeff == -1:
			return "-x^{}".format(self.degree)
		else:
			return "{}x^{}".format(self.coeff, self.degree)

	def __neg__(a):
		return monomial(-a.coeff, a.degree)

	def __add__(a, b):
		# monomial + monomial case
		if isinstance(b, monomial): 
			# if the monomials are of the same degree, add coefficients
			if a.degree == b.degree:
				return monomial(a.coeff + b.coeff, a.degree)
			# otherwise, just return a polynomial made of the monomials
			else:
				return polynomial([a, b])
		# monomial + polynomial case
		elif isinstance(b, polynomial):
			# if there is already a term of the same degree, add their coefficients
			if a.degree in [a.degree for a in b.terms]:
				newterm = a + b[a.degree] # recursively call monomial + monomial case
				return polynomial(list(filter(lambda x: x.degree != a.degree, 
					b.terms)) + [newterm])
			# otherwise, just add the monomial as another term
			else:
				return polynomial(b.terms + [a])

	def __sub__(a, b):
		return a + (-b)

	def __mul__(a, b):
		return monomial(a.coeff * b.coeff, a.degree + b.degree)

	def eval(self, x):
		return self.coeff * x**self.degree


class polynomial:
	def __init__(self, terms):
		self.terms = terms # list of monomials
		self.degree = max(a.degree for a in terms)
		self.clean()

	def clean(self):
		""" 
		A cleanup helper function.
		Ensures that the terms are ordered by degree, and that
		terms with coefficient 0 are eliminated.
		TODO: consolidate terms of same degree
		"""
		self.terms = list(filter(lambda x: x.coeff != 0, self.terms))
		self.terms.sort(key = lambda x: x.degree, reverse = True)

	def __str__(self):
		return ' + '.join(str(a) for a in self.terms)

	def __getitem__(self, degree):
		if degree in [i.degree for i in self.terms]:
			return next(x for x in self.terms if x.degree == degree)
		else:
			#raise IndexError("No term of that degree")
			return 0
			# still don't know whether it should throw an error or give
			# the technically mathematically correct answer

	def __neg__(a):
		return polynomial([-term for term in a.terms])

	def __add__(a, b):
		# polynomial + monomial case
		if isinstance(b, monomial):
			# recursively call monomial + polynomial case
			return b + a
		# polynomial + polynomial case
		elif isinstance(b, polynomial):
			out = polynomial([i for i in a.terms]) # necessary to create a deep copy of a
			# recursively call monomial + polynomial case
			for term in b.terms:
				out = term + out
			return out

	def __sub__(a, b):
		return a + (-b)

	def __mul__(a, b):
		return polynomial([i*j for i in a.terms for j in b.terms])

	def eval(self, x):
		return sum(a.eval(x) for a in self.terms)


a1 = monomial(2, 1) # 2x
a2 = monomial(-1, 2) # -x^2
a3 = monomial(1, 2) # x^2
a5 = monomial(1, 5) # x^5

f = a1 + a2 # -x^2 + 2x
g = a3 + a5 # x^5 + x^2

h = f - g
print(h.degree)
print(str(h))
print("h(1) =", h.eval(1))