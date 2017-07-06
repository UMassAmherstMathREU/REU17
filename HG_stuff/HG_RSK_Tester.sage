from collections import Counter

def Test_HG_RSK():
    x = 0
    while x < 10000:
        length = randint(1,5)
        width = randint(1,5)
        height = randint(1,5)
        P = PlanePartitions((length,width,height))
        P = P.random_element()
        RP = PlanePartition([list(reversed(r)) for r in reversed(P)])
        #print P
        #print
        #print pp_to_mat(P)
        #print
        #print HG.pp()
        HG = Hillman_Grassl(RP)
        tot1 = sum(HG.entries())
        tot2 = sum(1 for i,j,k in P.cells() if i == j)
        if tot1 != tot2:
            print "It failed: %r, %r" % tot1, tot2
            return

        x += 1
        if x % 500 == 0:
            print x
    print "No failure!"

def Test_Trans():
    x = 0
    while x < 1000:

        #test case 1

        length = randint(1,10)
        width = randint(1,10)
        height = randint(1,10)
        '''
        #test1
        P = PlanePartitions((length,width,height))
        P = P.random_element()
        PT = P.transpose()
        RP = PlanePartition([list(reversed(r)) for r in reversed(P)])
        RPT = PlanePartition([list(reversed(r)) for r in reversed(PT)])
        HG = Hillman_Grassl(RP)
        HGT = HG.conjugate()

        if (HG != HGT and P == PT) or (HG == HGT and P != PT):
            print "It failed:"
            return
        '''
        #test2
        P = PlanePartitions((length,width,height))
        P = P.random_element()
        APT = P.transpose()
        ARPT = PlanePartition([list(reversed(r)) for r in reversed(APT)])
        AHG = Hillman_Grassl(ARPT)
        BRP = PlanePartition([list(reversed(r)) for r in reversed(P)])
        BHG = Hillman_Grassl(BRP)
        BHGT = BHG.conjugate()

        if AHG != BHGT:
            print "It failed:"
            return

        x += 1
        if x % 100 == 0:
            print x
    print "No failure!"
