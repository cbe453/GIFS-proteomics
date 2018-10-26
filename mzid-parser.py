#!/usr/bin/env python3

import sys
from argparse import ArgumentParser as _ArgParser
from argparse import ArgumentError
from collections import namedtuple
from collections import Counter

from pyteomics import mzid
from pyteomics import mzml
from pyteomics import mzxml
from pyteomics import pepxml
from pyteomics import tandem

xmlns = {'mzIdentML': 'http://psidev.info/psi/pi/mzIdentML/1.1'}



Format = namedtuple('Format', ['name', 'reader'])

supported_formats = {'mzIdentML': mzid,
                    'mzML': mzml,
                    'mzXML': mzxml,
                    'pepXML': pepxml,
                    'X!Tandem': tandem}

class ArgParser(_ArgParser):

   def __init__(self):
       super(ArgParser, self).__init__(description='Extract hits into CSV file')

       self.add_argument('-i', '--input-file', dest='infile', type=str, required=True,
                         help='Search results in XML format')
       self.add_argument('-f', '--format', dest='format',
                         required='True',
                         choices=list(supported_formats.keys()),
                         help='Input file format')

def main(ns):

   unprocessed = supported_formats[ns.format].DataFrame(ns.infile)
   #filtered = mzid.filter("170712_UCD_K_batch1_TMT_F10_500ng_180min_ReZipTip_01_p.mzid", fdr=0.1)
   print(list(unprocessed.columns.values))
   seq = list(unprocessed['PeptideSequence'])
   mod = list(unprocessed['Modification'])
   assert(len(seq) == len(mod))

   #for k in range(len(seq)):
       #print(k, seq[k])
       #try:
           #for m in mod[k]:
               #print('\t', m)
       #except:
           #pass
           
   #for row in list(unprocessed.values):
      #print(row)
   #for k in list(unprocessed['accession']):
       #print(k)

   return 0

if __name__ == '__main__':
   parser = ArgParser()
   try:
       ns = parser.parse_args()
   except ArgumentError as err:
       print(err)
       sys.exit(1)

   sys.exit(main(ns))
