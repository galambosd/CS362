from Bio import SeqIO
from Bio.Seq import Seq
import numpy
from sys import argv

fasta1_name = argv[1]
fasta2_name = argv[2]
scoringFile_name = argv[3]

# put fasta/scoring parsing in another file?
#parse the fasta files
with open(fasta1_name, 'rU') as handle1, open(fasta2_name, 'rU') as handle2:
    for record in SeqIO.parse(handle1, 'fasta'):
        fasta1_record = record
    for record in SeqIO.parse(handle2, 'fasta'):
        fasta2_record = record
    fasta1Seq = fasta1_record.seq
    fasta2Seq = fasta2_record.seq
    fasta1 = str(fasta1Seq.upper())
    fasta2 = str(fasta2Seq.upper())

# parse the scoring file, ignoring the gap extend penalty 
    scoringFile=open(scoringFile_name, 'r')
    lines = scoringFile.readlines()
    cols = lines[1].split('\t')
    score_params = []
    for col in cols:
        score_params.append(int(col.rstrip('\n')))
    # return score_params here

# parse the scoring file -- ignore the gap extend
# pass along an array of scoring parameters


# alignment table class
class LocalAlignmentTable(object):
    """docstring for LocalAlignmentTable."""
    def __init__(self, size_m, size_n):
        super(LocalAlignmentTable, self).__init__()
        self.size_m = size_m
        self.size_n = size_n
        self.array = [[0 for x in range(0,size_m+1)] for y in range(0,size_n+1)]
        self.max_loc = None

    class Cell:
        """docstring for cell."""
        # input i,j as if it's a matrix -- putting something in (i,j) will put in array[j,i]
        def __init__(self,i,j,value, prev):
            self.value = value
            self.prev = prev
            self.prev_type = None
            self.i=i
            self.j=j
            self.optimal_score = None

    # global init method
    # local init method
    def localInit(self):
        i= 0
        j= 0
        # go to i<= size_n+1 -- want to initialize values for empty string + all of v
        while i < self.size_n+1:
            self.array[i][0]= self.Cell(i,0,0, None)
            i+=1
        while j< self.size_m+1:
            self.array[0][j] = self.Cell(0,j,0, None)
            j+=1

    # fill out the table -- maybe have separate helper methods that actually return largest values for local/global?
        # track largest value as we go
        # update backtrace pointer (could come from helper method too)

    def match(self, char1, char2, scoring):
        if char1==char2:
            return scoring[0]
        else:
            return scoring[1]

    def fillOutHelper(self, i, j, scoring):
        # a list of possible values to choose from to make indexing possible
        options = [self.get_val(i-1, j-1) + self.match(fasta1[i-1], fasta2[j-1],scoring),
            self.get_val(i, j-1) + scoring[2], self.get_val(i-1,j)+scoring[2], 0]
        max_index = numpy.argmax(options)
        max_val = options[max_index]
        return max_val, max_index

    def fillOut(self, scoring):
        # location of the max value across the whole local alignment table
        self.max_loc= (0,0)
        # table is filled out row-by-row
        for j in numpy.arange(1, self.size_m+1):
            for i in numpy.arange(1, self.size_n+1):
                # options for previous cells that match up with the values in fillOutHelper(). Down, left and up in the matrix
                options = [(i-1,j-1), (i,j-1), (i-1,j)]
                descriptive_location = ['d', 'l', 'u']
                value, prev_loc= self.fillOutHelper(i,j, scoring)
                # if the max value returned was a 0, no pointer to the previous cell
                if value ==0:
                    prev = None, None
                else:
                    # which cell is the previous cell in relation to the current one?
                    prev = options[prev_loc]
                    # actually point the prev variable to the right cell
                    prev = self.array[prev[0]][prev[1]], descriptive_location[prev_loc]
                    # TEST
                self.set_cell(self.Cell(i,j,value, prev),i,j)
                if value > self.get_val(self.max_loc[0], self.max_loc[1]):
                    self.max_loc = (i,j)
        self.optimal_score = self.get_val(self.max_loc[0], self.max_loc[1])
        return self.optimal_score

    def localBacktrace(self, fasta1, fasta2):
        align1=''
        align2 = ''
        # identify the max value cell. Don't want the value, but want to access the object in that cell directly
        start_cell = self.array[self.max_loc[0]][self.max_loc[1]]
        # print self.array
        #print start_cell
        # add the parts to the alignment strings that aren't part of the
        align1 = align1 + fasta1[start_cell.i:]
        align2 = align2 + fasta2[start_cell.j:]
        prev_cell = (start_cell.prev)[0]
        prev_align = (start_cell.prev)[1]
        while prev_cell.value > 0:
            print align1
            print align2
            # use the value of prev_align and various if-tests to determine what to insert
            if prev_align ==  'd':
                align1 = fasta1[prev_cell.i]+align1
                align2 = fasta2[prev_cell.j]+align2
            elif prev_align == 'l':
                align2 = fasta2[prev_cell.j]+align2
                align1= '-'+align1
            else:
                align1 = fasta1[prev_cell.i]+align1
                align2= '-'+align2

            prev_align = (prev_cell.prev)[1]
            self.last_i=prev_cell.i
            self.last_j=prev_cell.j
            print self.last_i, self.last_j
            prev_cell = (prev_cell.prev)[0]
        align1= fasta1[:self.last_i] + align1
        align2= fasta2[:self.last_j] + align2
        return align1, align2

    def get_val(self, i, j):
        return self.array[i][j].value

    def set_cell(self, value, i, j):
        self.array[i][j] = value

    def format_alignment(self, align1, align2):
        length_dif = abs(self.last_i-self.last_j)
        if self.last_i>self.last_j:
            align2=' '*length_dif + align2
        else:
            align1=' '*length_dif + align1
        return '*'*10+'\n'+align1+'\n'+align2+'\n'+'*'*10


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
    print 'Test:'
    print fasta1, len(fasta1)
    print fasta2, len(fasta2)
    print score_params
    divider =  '*'*20
    print divider
    i = 0
    j=0
    table = LocalAlignmentTable(len(fasta2), len(fasta1))
    table.localInit()
    print table.fillOut(score_params)
    print table.max_loc
    print divider
    printout1, printout2 = table.localBacktrace(fasta1,fasta2)
    print table.format_alignment(printout1, printout2)
