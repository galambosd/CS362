global.py
DEPENDENCIES: BioPython, numpy
OUTPUT: Optimal global alignment with affine gap penalties and maximum score.
Output is formatted to print successive, matched portions of
the alignment that fit comfortably on the Terminal screen.
BUGS: No known major bugs.

local_2.py
DEPENDENCIES: BioPython, numpy
OUTPUT: Optimal local alignment and maximum score.
Output is formatted to print the whole aligned first string
followed by the whole aligned second string (ie, the alignments match
if the Terminal screen is wide enough).
BUGS: No known major bugs.
