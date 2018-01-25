from Bio import SeqIO
from Bio.Seq import Seq
import numpy as np
from sys import argv

class LocalAlignmentTable(object):
    """docstring for LocalAlignmentTable."""
    def __init__(self, fasta1, fasta2, scoring):
        # initialize a 2-d array for the alignment, initially setting all the cells to 0
        # get the right size
        self.size_m = len(fasta2)
        self.size_n = len(fasta1)
        if (self.size_m == 0) or (self.size_n==0):
            raise IOError('One of the FASTA files that you entered is empty.')
        self.array = [[0 for x in np.arange(0,self.size_m+1)] for y in np.arange(0,self.size_n+1)]
        # input the sequences to align
        self.fasta1= fasta1
        self.fasta2=fasta2
        # get the scoring matrix
        self.scoring = scoring
        self.array[0][0]=self.LocalCell(0, None)
        i=1
        j=1
        try:
            self.gap_extend = self.scoring[3]
        except IndexError as e:
            raise e('The scoring file you entered is empty or improperly formatted.')
        # set the first column, first row to the right values
        while i < self.size_n+1:
            self.set_cell(self.LocalCell(0, None),i,0)
            i+=1
        while j< self.size_m+1:
            self.set_cell(self.LocalCell(0, None),0,j)
            j+=1

    class LocalCell:
        def __init__(self, value, prev):
            # create variables for the whole
            self.value = value
            self.prev = prev
    def match(self, char1, char2):
        if char1==char2:
            return self.scoring[0]
        else:
            return self.scoring[1]

    def calcAdjacent(self, i,j):
        E = self.get_cell(i,j-1).value + self.gap_extend
        F = self.get_cell(i-1,j).value + self.gap_extend
        G = self.get_cell(i-1,j-1).value + self.match(self.fasta1[i-1], self.fasta2[j-1])
        # return the list of the values of E, F and G
        return [E,F,G,0]

    def fillOut(self):
        # fill out the table, while setting 'previous cell' pointer to the right value
        # fill it out row by row
        # for each cell, compute E, F and G
        options = [(0, -1),(-1,0),(-1,-1)]
        max_loc = (0,0)
        max_val = 0
        for i in range(1, self.size_n+1):
            for j in range(1, self.size_m+1):
                # decide which one is the max
                # store value to S
                # adjust the prev pointer accordingly
                EFG = self.calcAdjacent(i,j)
                val = max(EFG)
                if val == 0:
                    prev = None
                else:
                # here, prev is tuple of i and j that shows relative location of
                # prev cell to current cell
                    prev = options[np.argmax(EFG)]
                self.set_cell(self.LocalCell(val, prev), i, j)
                if val > max_val:
                    max_loc = (i,j)
                    max_val = val
        # return the optimal score
        return max_loc

    def get_cell(self, i, j):
        return self.array[i][j]
    def set_cell(self, cell, i, j):
        self.array[i][j] = cell
    def format_alignment(self, align1, align2,i,j, per_line):
        length_dif = abs(i-j)
        if i>j:
            align2=' '*length_dif + align2
        else:
            align1=' '*length_dif + align1
        output_string = ''
        # i=0
        # left = len(align1)
        # while left >=per_line:
        #     output_string+=align1[i:i+per_line]+'\n'+align2[i:i+per_line]+'\n'
        #     output_string+='\n'
        #     i+=per_line
        #     left-=per_line
        # output_string+=align1[-left:]+'\n'+align2[-left:]
        return align1, align2

    def backTraceLocal(self, max_loc):
        # prepare two empty strings:
        align1 = ''
        align2 = ''
        # go to the max cell
        i = max_loc[0]
        j = max_loc[1]
        current_cell = self.get_cell(i,j)
        align1 = align1 + self.fasta1[i:]
        align2 = align2 + self.fasta2[j:]
        # while the current cell's value isn't 0:
        while current_cell.prev != None:
            # if the last cell is diagonal, add to both strings:
            if current_cell.prev ==  (-1,-1):
                align1 = self.fasta1[i-1]+align1
                align2 = self.fasta2[j-1]+align2
                j-=1
                i-=1
            # if it's to the left, put a gap in v
            elif current_cell.prev == (0,-1):
                align2 = self.fasta2[j-1]+align2
                align1= '-'+align1
                j-=1
            # if it isn't either of the above, put a gap in w
            else:
                align1 = self.fasta1[i-1]+align1
                align2= '-'+align2
                i-=1
            current_cell = self.get_cell(i,j)
        align1 = self.fasta1[0:i] + align1
        align2 = self.fasta2[0:j] + align2
        return align1, align2, i, j

def fileOpen(fasta1_name, fasta2_name, scoringFile_name):
    with open(fasta1_name, 'rU') as handle1, open(fasta2_name, 'rU') as handle2:
        for record in SeqIO.parse(handle1, 'fasta'):
            fasta1_record = record
        for record in SeqIO.parse(handle2, 'fasta'):
            fasta2_record = record
        fasta1Seq = fasta1_record.seq
        fasta2Seq = fasta2_record.seq
        fasta1 = str(fasta1Seq.upper())
        fasta2 = str(fasta2Seq.upper())
        handle1.close()
        handle2.close()

    scoringFile=open(scoringFile_name, 'r')
    lines = scoringFile.readlines()
    cols = lines[1].split('\t')
    scoringFile.close()
    score_params = []
    for col in cols:
        score_params.append(int(col.rstrip('\n')))

    return fasta1, fasta2, score_params

def main():
    fasta1_name = argv[1]
    fasta2_name = argv[2]
    scoringFile_name = argv[3]
    fasta1, fasta2, score_params = fileOpen(fasta1_name, fasta2_name, scoringFile_name)
    print 'Local alignment of \'{0}\' and \'{1}\' with scoring parameters \'{2}\':\n'.format(fasta1_name, fasta2_name, scoringFile_name)
    table = LocalAlignmentTable(fasta1, fasta2, score_params)
    max_loc = table.fillOut()
    max_score = table.get_cell(max_loc[0], max_loc[1]).value
    print 'Optimal score:{0}.'.format(max_score)
    printout1, printout2, i, j = table.backTraceLocal(max_loc)
    align1, align2 = table.format_alignment(printout1, printout2, i, j, 150)
    print align1
    print align2
    print 'END OF ALIGNMENT.'

main()
