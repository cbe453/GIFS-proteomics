# Small class for containing peptide information and a couple general methods.
class peptide:
	def __init__(self, content, name):
		self.query_name = name
		self.title = content[0][1]
		self.rt_in_seconds = content[1][1]
		self.index = content[2][1]
		self.charge = content[3][1]
		self.mass_min = content[4][1]
		self.mass_max = content[5][1]
		self.int_min = content[6][1]
		self.int_max = content[7][1]
		self.num_vals = content[8][1]
		self.num_used = content[9][1]
		self.ions = content[10][1].strip().split(",")
		self.counts = {126: 0, 127: 0, 128: 0, 129: 0, 130:  0,
						131: 0, 229: 0}
	
	# Return the list of ions associated with the peptide
	def getIons(self):
		if self.ions != None:
			return self.ions
		else:
			print("Ion list for " + self.query_name + " is empty!")
			return None
	
	# Prints all information associated with the peptide in a labelled line by line format.
	def toString(self):
		return ("Query: " + self.query_name + "\nTitle: " + self.title + "\nRtInSeconds: " + self.rt_in_seconds + 
			  "\nIndex: " + self.index + "\nCharge: " + self.charge + "\nMassMin: " + 
			  self.mass_min + "\nMassMax: " + self.mass_max + "\nIntMin: " + self.int_min +
			  "\nIntMax: " + self.int_max + "\nNumVals: " + self.num_vals + "\nNumUsed: " +
			  self.num_used + "\nIons: " + str(self.ions) + "\nCounts: " + str(self.counts))
	
	# Prints all information associated with a peptide in a tabulated format. Can easily
	# be converted to commas or whatever		  
	def tabFormat(self):
		return (self.query_name + "\t" + self.title + "\t" + self.rt_in_seconds + 
			  "\t" + self.index + "\t" + self.charge + "\t" + 
			  self.mass_min + "\t" + self.mass_max + "\t" + self.int_min +
			  "\t" + self.int_max + "\t" + self.num_vals + "\t" +
			  self.num_used + "\t" + str(self.ions) + "\t" + str(self.counts))

