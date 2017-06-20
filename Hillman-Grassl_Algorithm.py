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
    HG = [[0]* len(ppart) for r in ppart]

    #Holder variables for later
    row = len(ppart)
    col = 0
    prow = row
    pcol = 0
    pcolstore = 0

    ##5. Loop steps 2-4

    #Continues to run while part still has non-zero entires
    any(x != 0 for x in ppart.entries):

        ##2. Finding the southwest entry

        max(i for i in range(rows)
            if ppart[i][0])

        #stores the lowest row with a non-zero entry
        negi, j = min((-i, j) for i in range(rows)
                      for j in range(len(x[i]))
                      if x[i][j] != 0)

        #fixes the negative and gives friendlier names
        #prow: present row, pcol: present column, pcolstore: storing the present column for step 4
        prow, pcol = -negi, j
        pcolstore = pcol

        ##3.Develop the path

        #gets the length of the present row
        col = len(ppart[prow])

        #Loops until our present entry is on the right side of the partition
        while pcol <= col:
            #checks to see if the northern element is the same as the present one
            #if so it moves up and decrements the old element
            #also resets the length of the new row
            if ppart[prow][pcol] == ppart[prow-1][pcol]:
                ppart[prow][pcol] -= 1
                prow-=1
                col = len(ppart[prow])

            #if the northern element was not the same it moves east
            #and still decrements the old element
            else:
                ppart[prow][pcol] -= 1
                pcol+=1

        ##4 Increment HG
        HG[prow][pcolstore] += 1
    #Convert HG to a Tableau
    return HG

def Inverse_HG(HG):
    '''
    Converts from a Hillman Grassl Tableau to a plane partition

    Returns the plane partition

    Args:
        HG: a Hillman Grassl

    Returns:
        a plane partition
    '''
    #1. Creating the empty partition
    part = partition([0]*x for x in part)
    #setup variables
    row, col = HG.shape
    prow = 0
    pcol = 0
    #6. Loop steps 2-5
    any(x != 0 for x in ppart.entries):
        #2. Finding the northeast entry
        for x in row:
            for y in col:
                if HG[row][col] != 0:
                    prow = row
                    pcol = col
                    break
            if HG[row][col] != 0:
                break
        #3.Develop the path
        while pcol > 0:
            if HG[prow][pcol] == HG[prow+1][pcol]:
                HG[prow][pcol] -= 1
                prow+=1
            else:
                HG[prow][pcol] -= 1
                pcol-=1
        #5 Increment part
        part[prow][pcol] += 1
    return part
