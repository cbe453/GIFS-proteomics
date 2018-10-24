#!/usr/bin/env python3

__version__ = '0.1'


__disclaimer__ = 'TBD'

import email

import xml.etree.ElementTree as et

import sys

class peptide:
	def __init__(self, content):
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
	
	def toString(self):
		print("Title: " + self.title + "\nRtInSeconds: " + self.rt_in_seconds + 
			  "\nIndex: " + self.index + "\nCharge: " + self.charge + "\nMassMin: " + 
			  self.mass_min + "\nMassMax: " + self.mass_max + "\nIntMin: " + self.int_min +
			  "\nIntMax: " + self.int_max + "\nNumVals: " + self.num_vals + "\nNumUsed: " +
			  self.num_used + "\nIons: " + str(self.ions))
			  
	def tabFormat(self):
		return (self.title + "\t" + self.rt_in_seconds + 
			  "\t" + self.index + "\t" + self.charge + "\t" + 
			  self.mass_min + "\t" + self.mass_max + "\t" + self.int_min +
			  "\t" + self.int_max + "\t" + self.num_vals + "\t" +
			  self.num_used + "\t" + str(self.ions))


def parse_xml(part):
	return et.fromstring(part.get_payload())

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
	
	for i, (kind, name, content) in enumerate(parts, 1):
		print(i, kind, name)
		if kind == 'query':
			#print(content)
			new_peptide = peptide(content)
			print(new_peptide.tabFormat())
				#sys.stdout.write(new_peptide.tabFormat())
	return 0

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


