def Hillman_Grassl(part):
    '''
    Converts from a plane partition of shape lambda to an positive integer valued
    function over lambda. Information on plane partitions and their underlying
    shapes can be found with explain(plane_partition) or on Wikipedia.

    Args:
        part: a weak reverse plane partition with shape lambda

    Returns:
        a function from lambda to the positive integers
    '''
    #1. Creating the empty Tableau
    HG = Tableau([0]*x for x in part)
    #Setup variables
    ppart = part
    row, col = ppart.shape
    prow = row
    path = {}
    pcol = 0
    pcolstore
    #2. Finding the southwest entry
    for x in col:
        if ppart[prow][x] != 0:
            path[(prow,x)] = ppart[prow][x]
            pcol = x
            pcolstore = x
            break
        if x == col and ppart[prow][x] == 0:
            prow-=1
    #3.Develop the path
    while pcol <= col:
        if ppart[prow-1][pcol] == ppart[prow-1][pcol]:
            path[(prow-1,pcol)] = ppart[prow-1][pcol]
            prow-=1
        else:
            path[(prow,pcol+1)] = ppart[prow][pcol+1]
            pcol+=1
    #4 Decrement all elements in the path of PP by 1
    for i in path(x,y):
        ppart[x][y] -= 1
    #5 Increment HG
    HG[prow][pcolstore] += 1
    #6. Loop steps 2-5
    #Fix this later
