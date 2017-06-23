def Test_HG_RSK():
    length = randint(3,4)
    width = randint(3,4)
    height = randint(3,4)
    P = PlanePartitions((length,width,height))
    P = P.random_element()
    RP = PlanePartition([list(reversed(r)) for r in reversed(P)])
    print pp_to_mat(P)
    print
    print Hillman_Grassl(RP).pp()
