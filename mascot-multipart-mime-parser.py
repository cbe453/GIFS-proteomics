#!/usr/bin/env python3

__version__ = '0.1'


__disclaimer__ = 'TBD'

import email

import xml.etree.ElementTree as et

import sys, re

from collections import defaultdict

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
	
	def getIons(self):
		return self.ions
	
	def toString(self):
		return ("Query: " + self.query_name + "\nTitle: " + self.title + "\nRtInSeconds: " + self.rt_in_seconds + 
			  "\nIndex: " + self.index + "\nCharge: " + self.charge + "\nMassMin: " + 
			  self.mass_min + "\nMassMax: " + self.mass_max + "\nIntMin: " + self.int_min +
			  "\nIntMax: " + self.int_max + "\nNumVals: " + self.num_vals + "\nNumUsed: " +
			  self.num_used + "\nIons: " + str(self.ions))
			  
	def tabFormat(self):
		return (self.query_name + "\t" + self.title + "\t" + self.rt_in_seconds + 
			  "\t" + self.index + "\t" + self.charge + "\t" + 
			  self.mass_min + "\t" + self.mass_max + "\t" + self.int_min +
			  "\t" + self.int_max + "\t" + self.num_vals + "\t" +
			  self.num_used + "\t" + str(self.ions))

def findRecur(root):
    global indent
    if (root.tag.title() == '{Http://Www.Unimod.Org/Xmlns/Schema/Unimod_2}Misc_Notes' and root.text != None):
    	if (root.text != None and bool(re.search('Sixplex', root.text))):
    		global masses
    		masses = root.text.split(':')[1].split(' ')
    		return
    for elem in root.getchildren():
        findRecur(elem)

def parse_xml(part):
	content = et.fromstring(part.get_payload())
	findRecur(content)
	return content

def parse_key_value_pairs(part):
	pairs = list()
	for line in part.get_payload().split('\n'):
		things = line.split('=')
		key = things[0].strip()
		value = ''.join(things[1:])
		pairs.append((key, value))
	return pairs

def parse_enzyme(part):
	return part.get_payload()

def choose_handler(name):
	handler = None
	key = None
	for k in mime_parts.keys():
		if name.startswith(k):
			handler = mime_parts[k]
			key = k

	if key is None or key != 'query' and key != name:
		msg = 'MIME part name "{}" does not correspond to a valid key'.format(name)
		raise KeyError(msg)
       
	return key, handler

def part_iterator(infile):

	mascot = email.message_from_file(infile)

	for part in mascot.walk():
		if part.get_content_type() == 'multipart/mixed':
			continue

		assert(part.get_content_type() == 'application/x-mascot')
		
		name = part.get_param('name')
		try: 
			key, handler = choose_handler(name)
		except KeyError as err:
			print(err, file=sys.stderr)
		
		content = handler(part)
		yield (key, name, content)

def main(infile):
	parts = part_iterator(infile)
	peptides = defaultdict()
	matches = []
	
	for i, (kind, name, content) in enumerate(parts, 1):
		if kind == 'query':
			trunc_name = re.sub('uery', '', name)
			peptides[trunc_name] = peptide(content, trunc_name)
		elif kind == 'peptides':
			for key, item in content:
				if (re.match('q[0-9]+_p[0-9]+$', key) and item != '-1' ):
					matches.append(item)
					print(key + "\t" +item)
	
	for i in range(len(masses)):
		masses[i] = float(masses[i])
	
	tag_dict = defaultdict(list)
	multi_tag_count = 0
	
	for key in peptides.keys():
		tag_flag = False
		delta_flag = False
		tag_hits = set([])
		
		for ion in peptides[key].getIons():
			ion_float = float(ion.split(':')[0])
			ion_count = int(float(ion.split(':')[1]))
			ion_int = int(ion_float)
			if (ion_int == 229 and ion_count > 1000):
				delta_flag = True
			for mass in masses:
				if ((mass - 0.001) < ion_float < (mass + 0.001)):
					tag_flag = True
					tag_hits.add(int(mass))
		
		if tag_flag and delta_flag:
			if len(tag_hits) > 1:
				multi_tag_count += 1
			for tag in tag_hits:
				tag_dict[tag].append(peptides[key])
	
	for key in tag_dict.keys():
		for pep in tag_dict[key]:
			print(str(key) + "\t" + pep.tabFormat())			 	
	
	print(masses)
	print(multi_tag_count)
	
	#for item in matches:
		#print(item)
    
mime_parts = {'parameters': parse_key_value_pairs,
               'masses' : parse_key_value_pairs,
               'quantitation': parse_xml,
               'unimod': parse_xml,
               'enzyme': parse_enzyme,
               'header': parse_key_value_pairs,
               'summary': parse_key_value_pairs,
               'peptides': parse_key_value_pairs,
               'proteins': parse_key_value_pairs,
               'query': parse_key_value_pairs,
               'index': parse_key_value_pairs}

if __name__ == '__main__':
    with open(sys.argv[1], 'r') as input:
        sys.exit(main(input))


