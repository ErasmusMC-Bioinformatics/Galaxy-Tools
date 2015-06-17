import logging
log = logging.getLogger(__name__)

from galaxy import util
import galaxy
import galaxy.model
import galaxy.datatypes
import galaxy.datatypes.data

from galaxy.datatypes.metadata import MetadataElement

from galaxy.datatypes.sequence import Sequence
from galaxy.datatypes.tabular import Tabular
from galaxy.datatypes.xml import GenericXml

from galaxy.datatypes.data import Data


import re

class DotBracket ( Sequence ):
    edam_format = "format_1457"
    file_ext = "dbn"
    
    sequence_regexp = re.compile( "^[ACGTURYKMSWBDHVN]+$", re.I)
    structure_regexp = re.compile( "^[\(\)\.\[\]]*" )
    
    def set_meta( self, dataset, **kwd ):
        """
        Set the number of sequences and the number of data lines
        in dataset.
        """
        if self.max_optional_metadata_filesize >= 0 and dataset.get_size() > self.max_optional_metadata_filesize:
            dataset.metadata.data_lines = None
            dataset.metadata.sequences = None
            dataset.metadata.seconday_structures = None
            return
        
        data_lines = 0
        sequences = 0
        
        for line in file( dataset.file_name ):
            line = line.strip()
            data_lines += 1
            
            if line and line.startswith( '>' ):
                sequences += 1
        
        dataset.metadata.data_lines = data_lines
        dataset.metadata.sequences = sequences
    
    def sniff(self, filename):
        """
        The Dot-Bracket format is as follows:
        
        >sequenceName1
        CCCaaaGGG
        (((...)))
        >sequenceName2
        GGGuuuCCC
        (((...)))
        
        Because it remains unclear whether the Dot-Bracket format may
        contain multiple sequences per file, sniffing is only applied
        on the first 3 lines.
        """
    
    i = 0
    pairs = False
    
    with open( filename ) as handle:
        for line in handle:
            line = line.strip()
            
            state = i % 3
            
            #header line
            if state == 0:
                if(line[0] != '>'):
                    return False
            
            #sequence line
            elif state == 1:
                if not sequence_regexp.match(line):
                    return False
                else:
                    sequence_size = len(line)
            
            #dot-bracket structure line
            elif state == 2:
                if (sequence_size != len(line)) or (not structure_regexp.match(line)) or (line.count('(') != line.count(')')) or (line.count('[') != line.count(']')):
                    return False
                else:
                    return True
            
            i += 1
        
        # Number of lines is less than 3
        return False

class ConnectivityTable( Tabular ):
    edam_format = "format_3309"
    file_ext = "ct"
    
    header_regexp = re.compile( "^[0-9]+" + "(?:\t|[ ]+)" + ".*?" + "(?:ENERGY|energy|dG)")
    structure_regexp = re.compile( "^[0-9]+" + "(?:\t|[ ]+)" +  "[ACGTURYKMSWBDHVN]+" + "(?:\t|[ ]+)" + "[^\t]+" + "(?:\t|[ ]+)" + "[^\t]+" + "(?:\t|[ ]+)" + "[^\t]+" + "(?:\t|[ ]+)" + "[^\t]+")
    
    def __init__(self, **kwd):
        Tabular.__init__( self, **kwd )
        
        self.columns = 6
        self.column_names = ['base_index', 'base', 'neighbor_left', 'neighbor_right', 'partner', 'natural_numbering']
        self.column_types = ['int', 'str', 'int', 'int', 'int', 'int']

    def set_meta( self, dataset, **kwd ):
        data_lines = 0
        
        for line in file( dataset.file_name ):
            data_lines += 1
        
        dataset.metadata.data_lines = data_lines
    
    def sniff(self, filename):
        """
        The ConnectivityTable (CT) basic format is defined as follows:
        
5	energy = -12.3	sequence name
1	G	0	2	0	1
2	A	1	3	0	2
3	A	2	4	0	3
4	A	3	5	0	4
5	C	4	6	1	5
        
        The links given at the edam ontology page do not indicate what
        type of separator is used (space or tab) while different
        implementations exist. The implementation by RNAStructure is as
        follows:

   10    ENERGY = -34.8  seqname
    1 G       0    2    9    1
    2 G       1    3    8    2
    3 G       2    4    7    3
    4 a       3    5    0    4
    5 a       4    6    0    5
    6 a       5    7    0    6
    7 C       6    8    3    7
    8 C       7    9    2    8
    9 C       8   10    1    9
   10 a       9    0    0   10

        """
        filename = filename.file_name
        
        i = 0
        with open( filename ) as handle:
            for line in handle:
                line = line.strip()
                
                if(i == 0):
                    if not self.header_regexp.match(line):
                        return False
                else:
                    if not self.structure_regexp.match(line.upper()):
                        return False
                i += 1
        return True


class RNAML( GenericXml ):
    edam_format = "format_3311"
    file_ext = "rnaml"
