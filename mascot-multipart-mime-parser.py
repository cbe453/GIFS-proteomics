#!/usr/bin/env python3

__version__ = '0.1'


__disclaimer__ = 'TBD'

import email

import xml.etree.ElementTree as et

import sys, re, os

from collections import defaultdict
from peptide import peptide


# Probably unnecessarily complicated recursive method of extracting the secondary masses 
#from the XML portion of the .dat file. It works
def findRecur(root):
	global found_masses
	if found_masses:
		return
	if (root.tag.title() == '{Http://Www.Unimod.Org/Xmlns/Schema/Unimod_2}Misc_Notes' and root.text != None):
		if (root.text != None and bool(re.search('Sixplex', root.text))):
			global masses
			masses = root.text.split(':')[1].split(' ')
			found_masses = True
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

def main():

	peptides = defaultdict()
	matches = defaultdict(dict)
	global found_masses 
	found_masses = False
	outfile = open(sys.argv[1], 'w')
	outfile.write("Protein\tQuery-File\tSequence\t126\t127\t128\t129\t130\t131\t229\n")
	#debug = open('debug.out', 'w')
	
	for file in sys.argv[2:]:
		read_file = open(file, 'r')
		parts = part_iterator(read_file)
		file_basename = os.path.basename(file)
	
		# Iterate through input file and pull out relevant data. (mainly peptide information 
		# query information as well as the mass tags outlined in the XML portion of the file)
		# Isoforms that are matched in queries are collapsed down to a single entry in 
		# the dictionary.
		for i, (kind, name, content) in enumerate(parts, 1):
			if kind == 'query':
				trunc_name = re.sub('uery', '', name)
				new_key = trunc_name + "-" + file_basename
				#debug.write(content + "\n")
				if len(content) < 11:
					print(name + " from " + file_basename + " is missing content! Skipping this record...")
					continue
				else:
					peptides[new_key] = peptide(content, trunc_name)
			elif kind == 'peptides':
				for key, item in content:
					# Two assumptions are made here. We are assuming that the first hit i.e. 'p1'
					# reported by mascot is the top protein match. This prevents some 
					# ambiguous cases. We also remove cases where a single peptide query
					# matches both Cucumber and Watermelon as these cases are ambiguous. 
					if (re.match('q[0-9]+_p1$', key) and item != '-1'):
						if(re.search('Cla', item) and re.search('Csa', item)):
							continue
						else:
							sequence = item.split(';')[0].split(',')[4]
							for protein in item.split(';')[1].split(','):
								new_item = (key.split('_')[0] + "-" + str(file_basename))
								temp_protein = re.sub('\.[0-9]+', '', protein)
								trunc_protein = re.sub('\"', '', temp_protein)
								if sequence in matches[trunc_protein.split(':')[0]]:
									matches[trunc_protein.split(':')[0]][sequence].add(new_item)
								else:
									matches[trunc_protein.split(':')[0]][sequence] = set([new_item])
		read_file.close()
		
	# Convert string bases masses to floats
	for i in range(len(masses)):
		masses[i] = float(masses[i])
	
	# Iterate through the dictionary of dictionaries of lists  structured as 
	# proteins -> sequences -> query_lists. Decoys are removed, proteins with fewer than
	# 3 supporting queries are ignored. Prints out proteins, sequences, file names and 
	# tag counts for proteins with more than 3 supporting queries along with the 229 isobaric
	# mass tags and secondary masses
	for protein in matches.keys():
		if (re.match('DECOY', protein)):
			continue
		else:
			if (sum(len(list) for list in matches[protein].values())) >= 3:
				for sequence in matches[protein].keys():
					for query in matches[protein][sequence]:
						cur_peptide = peptides[query]
						primary_flag = False
						
						for ion in cur_peptide.getIons():
							ion_int = int(float(ion.split(':')[0]))
							ion_count = int(float(ion.split(':')[1]))
							ion_float = float(ion.split(':')[0])

							if (ion_int == 229):
								primary_flag = True
								cur_peptide.counts[229] = ion_count
								continue
							for mass in masses:
								if ((mass - 0.001) < ion_float < (mass + 0.001)):
									cur_peptide.counts[int(mass)] = ion_count
						
						if primary_flag:
							#print(protein + "\t" + query + "\t" + sequence
								  #+ "\t" + str(cur_peptide.counts))
							outfile.write(protein + "\t" + query + "\t" + sequence + "\t"
								  + '\t'.join(str(x) for x in cur_peptide.counts.values()) + "\n")
	outfile.close()
    
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
    sys.exit(main())


