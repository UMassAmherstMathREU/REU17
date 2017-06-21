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

        #stores the lowest row with a non-zero entry and the entry
        negi, j = min((-i, j) for i in range(row)
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
    HGT = Tableau(HG)

    return HGT

def Inverse_HG(HG):
    '''
    Converts from a Hillman Grassl Tableau to a plane partition

    Returns the plane partition

    Args:
        HG: a Hillman Grassl

    Returns:
        a plane partition
    '''
    ##1. Creating the empty Tableau and setup variables

    #A copy of HG stored as a list of lists
    HGL = [list(r) for r in part]

    #A copy of HG with all zero entries stored as a list of lists
    ppart = [[0]* len(HGL) for r in HGL]

    #Holder variables for later
    row = len(HGL)
    col = 0
    prow = 0
    pcol = 0
    pcolstore = 0

    ##5. Loop steps 2-4

    #Continues to run while HGL still has non-zero entires
    any(x != 0 for x in HGL.entries):

        ##2. Finding the northeast entry

        max(i for i in range(row)
            if HGL[i][0])

        #stores the higest row with a non-zero entry and the entry
        #prow: present row, pcol: present column, pcolstore: storing the present column for step 4
        prow, pcol = min((i, j) for i in range(row)
                      for j in range(len(x[i]))
                      if x[i][j] != 0)
        pcolstore = pcol

        ##3.Develop the path

        #assigns col to the bounding column and pcol to the eastern edge of the HGL
        col = pcol
        pcol = len(HGL[prow])

        #Loops until our present entry reaches the bounding column
        while pcol >= col:
            #checks to see if the southern element is the same as the present one
            #if so it moves down and increments the old element
            #also resets the length of the new row
            if ppart[prow][pcol] == ppart[prow+1][pcol]:
                ppart[prow][pcol] += 1
                prow+=1
                col = len(HGL[prow])

            #if the southern element was not the same it moves west
            #and still increments the old element
            else:
                ppart[prow][pcol] -= 1
                pcol+=1

        ##4 Decrement HGL
        HGL[prow][pcolstore] -= 1
    #Convert ppart to a Tableau
    part = Tableau(ppart)

    return part
