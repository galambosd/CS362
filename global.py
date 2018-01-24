from Bio import SeqIO
from Bio.Seq import Seq
import numpy as np
from sys import argv

class GlobalAlignmentTable(object):
    """docstring for GlobalAlignmentTable."""
    def __init__(self, fasta1, fasta2, scoring):
        # initialize a 2-d array for the alignment, initially setting all the cells to 0
        # get the right size
        self.size_m = len(fasta2)
        self.size_n = len(fasta1)
        if (self.size_m == 0) or (self.size_n==0):
            raise IOException('One of the FASTA files that you entered is empty.')
        self.array = [[0 for x in range(0,size_m+1)] for y in range(0,size_n+1)]
        # input the sequences to align
        self.fasta1= fasta1
        self.fasta2=fasta2
        # get the scoring matrix
        self.scoring = scoring
        self.array[0][0]=self.AffineCell(0,-3,-3,-1*np.infty)
        i=1
        j=1
        try:
            gap_open = self.scoring[2]
            gap_extend = self.scoring[3]
        except IndexError as e:
            raise e('The scoring file you entered is empty or improperly formatted.')
        # set the first column, first row to the right values
        while i < self.size_n+1:
            self.array[i][0] = self.AffineCell(gap_open+i*gap_extend,-1*np.infty,gap_open+i*gap_extend,-1*np.infty, None)
            i+=1
        while j< self.size_m+1:
            self.array[0][j] = self.AffineCell(gap_open+i*gap_extend,gap_open+i*gap_extend,-1*np.infty,-1*np.infty, None)
            j+=1

    class AffineCell:
        def __init__(self, S, E, F, G, prev)
            # create variables for the whole
            self.S=S
            self.E=E
            self.F=F
            self.G=G
            self.prev = prev
    def match(self, char1, char2):
        if char1==char2:
            return self.scoring[0]
        else:
            return self.scoring[1]

    def calcEFG(self, i,j):
        E = max(self.get_cell(i,j-1).E + gap_extend, self.get_cell(i, j-1).S+gap_open+gap_extend)
        F = max(self.get_cell(i-1,j).F + gap_extend, self.get_cell(i-1, j).S+gap_open+gap_extend)
        G = self.get_cell(i-1,j-1).S + self.match(self.fasta1[i-1], self.fasta2[j-1])
        # return the list of the values of E, F and G
        return [E,F,G]

    def fillOut(self):
        # fill out the table, while setting 'previous cell' pointer to the right value
        # fill it out row by row
        # for each cell, compute E, F and G
        options = [(i, j-1),(i-1,j),(i-1,j-1)]
        for i in range(1, size_n):
            for j in range(1, size_m):
                # decide which one is the max
                # store value to S
                # adjust the prev pointer accordingly
                EFG = calcEFG(i,j)
                S_val = max(EFG)
                prev_loc = options[np.argmax(EFG)]
                prev_cell = get_cell(prev_loc[0], prev_loc[1])
                set_cell(AffineCell(S_val, E, F, G, get_cell(prev_cell)), i, j)
        # return the optimal score
        return get_cell(self.size_m, self.size_n).S

    def get_cell(self, i, j):
        return self.array[i][j]
    def set_cell(self, cell, i, j):
        self.array[i][j] = cell
    def format_alignment(self, align1, align2):
        pass
    def backTraceGlobal(self):
        # prepare two empty strings:
        align1 = ''
        align2 = ''
        # go to cell (m,n)
        current_cell = get_cell(self.size_m, self.size_n)
        # while the current cell isn't the sink:
        while current_cell.prev != get_cell(0,0):
            
        # keep following the traceback arrows until the source is reached
        # return a correctly formatted alignment, using the helper method format_alignment()
        pass



def main(arg):
    fasta1_name = argv[1]
    fasta2_name = argv[2]
    scoringFile_name = argv[3]

    with open(fasta1_name, 'rU') as handle1, open(fasta2_name, 'rU') as handle2:
        for record in SeqIO.parse(handle1, 'fasta'):
            fasta1_record = record
        for record in SeqIO.parse(handle2, 'fasta'):
            fasta2_record = record
        fasta1Seq = fasta1_record.seq
        fasta2Seq = fasta2_record.seq
        fasta1 = str(fasta1Seq.upper())
        fasta2 = str(fasta2Seq.upper())

    scoringFile=open(scoringFile_name, 'r')
    lines = scoringFile.readlines()
    cols = lines[1].split('\t')
    score_params = []
    for col in cols:
        score_params.append(int(col.rstrip('\n')))


    def run_tests():
        pass

if __name__ == '__main__':
    main()
    run_tests()
