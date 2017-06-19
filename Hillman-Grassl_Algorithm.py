def Hillman_Grassl(part):
    '''
    Converts from a plane partition to a Hillman Grassl Tableau

    Returns the Hillman Grassl Tableau

    Args:
        part: a weak reverse plane partition

    Returns:
        A Hillman Grassl Tableau
    '''
    #1. Creating the empty Tableau
    HG = Tableau([0]*x for x in part)
    #Setup variables
    ppart = part
    row, col = ppart.shape
    prow = row
    pcol = 0
    pcolstore = 0
    #6. Loop steps 2-5
    any(x != 0 for x in ppart.entries):
        #2. Finding the southwest entry
        for x in col:
            if ppart[prow][x] != 0:
                pcol = x
                pcolstore = x
                break
            if x == col and ppart[prow][x] == 0:
                prow-=1
        #3.Develop the path
        while pcol <= col:
            if ppart[prow][pcol] == ppart[prow-1][pcol]:
                ppart[prow][pcol] -= 1
                prow-=1
            else:
                ppart[prow][pcol] -= 1
                pcol+=1
        #5 Increment HG
        HG[prow][pcolstore] += 1
    #Fix this later
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
