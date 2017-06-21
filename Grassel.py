def Hillman_Grassl(part):
    '''
    Converts from a plane partition of shape lambda to an positive integer valued
    function over lambda. Information on plane partitions and their underlying
    shapes can be found with explain(plane_partition) or on Wikipedia.

    Args:
        part: a weak reverse plane partition with shape lambda (sage class ?)

    Returns:
        a function from lambda to the positive integers
    '''
    #1. Creating the empty Tableau
    rows = part.shape
    res = [[0]*row for row in rows]
    #Setup variables
    path = {}
    pcol = 0
    pcolstore

    left_col = len(part.shape)
    #2. Finding the southwest entry
    for j in xrange(left_col-1,0,-1):
        if part[j][0] != 0:
            sw_entry = (j,0)
            path[(j,0)] = part[j][0]
            break

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
