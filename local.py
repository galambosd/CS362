from Bio import SeqIO
from Bio.Seq import Seq
import numpy
from sys import argv

fasta1_name = argv[1]
fasta2_name = argv[2]
scoringFile_name = argv[3]

# put fasta/scoring parsing in another file?
#parse the fasta files
with open(fasta1_name, 'rU'), open(fasta2_name, 'rU') as handle1, handle2:
    fasta1_record= SeqIO.parse(handle1, 'fasta')
    fasta2_record= SeqIO.parse(handle2, 'fasta')
    fasta1Seq = fasta1_record.seq
    fasta2Seq = fasta2_record.seq
    fasta1 = str(fasta1Seq.upper())
    fasta2 = str(fasta2Seq.upper())

# parse the scoring file, ignoring the gap extend
with scoringFile_name as handle:
    scoringFile = open(scoringFile_name, 'r')
    lines = scoringFile.readlines()
    cols = lines[2].split('\t')
    score_params = []
    for col in cols:
        score_params.append(col.rstrip('\n'))
    # return score_params here

# parse the scoring file -- ignore the gap extend
# pass along an array of scoring parameters

# can python have objects as elements of arrays?

# alignment table class
class alignmentTable(object):
    """docstring for alignmentTable."""
    def __init__(self, size_m, size_n):
        super(alignmentTable, self).__init__()
        self.size_m = size_m
        self.size_n = size_n
        array = [][]

    class cell(alignmentTable):
        """docstring for cell."""
        # input i,j as if it's a matrix -- putting something in (i,j) will put in array[j,i]
        def __init__(self, value, prev):
            super(cell, self).__init__()
            self.value = value
            self.prev = prev

    # global init method
    def globalInit(self):
        pass
    # local init method
    def localInit(self):
        i,j = 0
        while i < size_n:
            array[0,i]= cell(0, None)
        while j< size_m:
            array[j,0] = cell(0, None)

    # fill out the table -- maybe have separate helper methods that actually return largest values for local/global?
        # track largest value as we go
        # update backtrace pointer (could come from helper method too)
    def fillOutHelper(self, i, j, scoring):
        pass
    def fillOut(self, i, j, scoring):
        pass
    def localBacktrace(self, arg):
        pass
    def globalBacktrace(self, arg):
        # return the indicies of the
        pass
    # local backtrace method (largest --> closest 0)
    # global backtrace method (sink --> source)

# make a class for the table with an init method (hard-coding init might be easier)
# make another class inside this that's each cell -- has a value, and backtrace pointers
# set up a table (a 2-D array w/ each cell is the cell object)
# initialize it
# fill out the values
# keep updating the max value
# get max value -- store 'hangover' part of alignment
# traceback from max value until you get to 0
# get the 'hangover' parts of v and w from the beginning and end
# print the alignment and the max value

if __name__=="__main__":
    print 'Test:':
    print fasta1
    print fasta2
    print score_params
    i, j = 0
    while i < size_n:
        print array[0,i].value ,
    while j< size_m:
        print array[j,0].value ,
