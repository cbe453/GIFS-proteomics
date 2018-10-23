#!/usr/bin/env python3

__version__ = '0.1'


__disclaimer__ = 'TBD'

import email

import xml.etree.ElementTree as et

import sys

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
           print(content[10][1])
        
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


