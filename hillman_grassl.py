from sage.combinat.tableau import Tableau
from sage.rings.integer import Integer

def Hillman_Grassl(part):
    '''
    Converts from a plane partition to a Hillman Grassl Tableau

    Returns the Hillman Grassl Tableau

    Note: See Hillman-Grassl_Correspondence.txt for more information on the Algorithm

    Args:
        part: a weak reverse plane partition

    Returns:
        A Hillman Grassl Tableau
    '''
    ##1. Creating the empty Tableau and setup variables

    #A copy of part stored as a list of lists
    ppart = [list(r) for r in part]

    #A copy of part with all zero entries stored as a list of lists
    HG = [[Integer(0)] * len(r) for r in ppart]

    #Holder variables for later
    row = len(ppart)
    col = 0
    prow = row
    pcol = 0
    pcolstore = 0

    ##5. Loop steps 2-4

    #Continues to run while part still has non-zero entires
    while any(x != 0 for r in ppart for x in r):

        ##2. Finding the southwest entry

        #stores the lowest row with a non-zero entry and the entry
        negi, j = min((-i, j) for i in range(row)
                      for j in range(len(ppart[i]))
                      if ppart[i][j] != 0)

        #fixes the negative and gives friendlier names
        #prow: present row, pcol: present column, pcolstore: storing
        #the present column for step 4
        prow, pcol = -negi, j
        pcolstore = pcol

        ##3.Develop the path

        #gets the length of the present row
        col = len(ppart[prow])

        #Loops until our present entry is on the right side of the partition
        while pcol < col:
            #checks to see if the northern element is the same as the present one
            #if so it moves up and decrements the old element
            #also resets the length of the new row
            if prow > 0 and ppart[prow][pcol] == ppart[prow-1][pcol]:
                ppart[prow][pcol] -= 1
                prow -= 1
                col = len(ppart[prow])

            #if the northern element was not the same it moves east
            #and still decrements the old element
            else:
                ppart[prow][pcol] -= 1
                pcol += 1

        ##4 Increment HG
        HG[prow][pcolstore] += 1
    #Convert HG to a Tableau
    HGT = Tableau(HG)

    return HGT

def Inverse_HG(HG):
    '''
    Convert from a Hillman Grassl Tableau to a weak reverse plane
    partition.

    Args:
        HG: The Hillman-Grassl Tableau.  Must be a Tableau, list 
            of lists, or any other iterator of iterators.

    Returns:
        A weak reverse plane partition, stored as a Tableau.
    '''
    # Store the HG and partition as mutable lists
    HGL = [list(r) for r in HG]
    ppart = [[Integer(0)] * len(r) for r in HGL]

    # Continue adding paths while there are non-zero entries
    while any(x != 0 for x in r for r in HGL):
        # Find the north-eastern most entry in HGL
        # Take advantage of how python orders tuples to do this
        negj, i = min((-j, i) for i in range(len(HGL))
                      for j in range(len(HGL[i]))
                      if HGL[i][j] != 0)
        
        # row and col are the coords of the entry in HGL
        row, col = i, -negj
        HGL[row][col] -= 1
        
        # prow and pcol are the current entry in the path
        # start with prow == row, end with pcol == col
        prow = row
        pcol = len(HGL[row]) - 1

        # Loop over each entry in the path
        while pcol >= col:
            # checks to see if the southern element is the same as the present one
            # if so it moves down and increments the old element
            # otherwise, move west
            if prow < len(ppart) - 1 and \
               pcol < len(ppart[prow+1]) and \
               ppart[prow][pcol] == ppart[prow + 1][pcol]:
                ppart[prow][pcol] += 1
                prow += 1
            else:
                ppart[prow][pcol] += 1
                pcol -= 1
    # done, return a tableau
    return Tableau(ppart)
