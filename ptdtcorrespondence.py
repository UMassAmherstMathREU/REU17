from sage.structure.sage_object import SageObject
from hillman_grassl import Hillman_Grassl, Inverse_HG
def reverse_hook(shape, i, j):
    conj = shape.conjugate()
    vertical = shape

class PTDTCorrespondence(SageObject):
    def __init__(self, shape):
        shape = Partition(shape)
        
        self._dt2pt = dict()
        self._pt2dt = dict()
        self._shape = shape

        if not shape:
            # empty partition, do nothing
            return
        
        max_length = shape.hook_length(0, 0)

        # Calculate the shape of the tail, possibly useless
        outer1 = Partition(shape[i] + max_length - i
                          for i in range(len(shape)))
        conj = shape.conjugate()
        outer2 = Partition(conj[i] + max_length - i
                          for i in range(len(conj))).conjugate()
        outer = Partition(outer1[i] if i < len(shape) else outer2[i]
                          for i in range(len(outer2)))
        tail = SkewPartition([outer, shape])
        self._tail = tail

        cells_by_hook = [[] for _ in range(max_length)]
        for i, j in tail.cells():
            vert = conj[j] if j < len(conj) else 0
            horiz = shape[i] if i < len(shape) else 0
            hook = i - horiz + j - vert + 1
            cells_by_hook[hook - 1].append((i, j))

        # match the inside of the shape (pt_l) to the tail (dt_l)
        for i, j in shape.cells():
            hook = shape.hook_length(i, j)
            match = cells_by_hook[hook - 1][-1]
            del cells_by_hook[hook - 1][-1]
            self._dt2pt[match] = i, j, True
            self._pt2dt[i, j, True] = match

        # match the remaining 
        for hook in range(max_length):
            num_matched = 0
            for i, j in cells_by_hook[hook]:
                p, q = hook - num_matched, num_matched
                self._dt2pt[i, j] = p, q, False
                self._pt2dt[p, q, False] = i, j
                num_matched += 1
            assert num_matched == hook + 1

        # Check everything worked as expected
        assert all((i, j) in self._dt2pt for i, j in tail.cells())
        assert all((i - j, j, False) in self._pt2dt
                   for i in range(0, max_length)
                   for j in range(i + 1))
        assert all((i, j, True) in self._pt2dt
                   for i, j in shape.cells())
        assert all(self._pt2dt[self._dt2pt[k]] == k
                   for k in self._dt2pt)
        assert all(self._dt2pt[self._pt2dt[k]] == k
                   for k in self._pt2dt)

    def dt2pt(self, dt_l):
        # Check our input is good
        dt_l = SkewTableau(dt_l)
        if dt_l.shape().inner() != self._shape:
            raise ValueError("DT has incorrect inner shape")

        # If the shape is empty, DT_0 = DT_l
        if not self._shape:
            return Tableau([]), dt_l
        
        # fill dt_l with zeros and flip it
        max_len = max(len(r) for r in dt_l)
        new_rows = [[0] * (max_len - len(r)) +
                    [x for x in reversed(r) if x is not None]
                    for r in reversed(dt_l)]

        # convert to hillman-grassl
        dt_l_hg = Hillman_Grassl(new_rows)

        # create new hillman-grassl's
        width = len(dt_l_hg[0])
        height = len(dt_l_hg)
        pt_l_hg = [[0] * l for l in self._shape]
        dt_0_hg = [[0] * width for i in range(height)]

        # fill in pt and dt_0.  This is the messy part
        for r in range(len(dt_l_hg)):
            for c in range(len(dt_l_hg[r])):
                value = dt_l_hg[r][c]
                if value == 0:
                    continue
                # r and c count from far corner, but our board drawings
                # count from the near corner
                i = height - r - 1
                j = width - c - 1
                if i >= len(self._shape) and j >= self._shape[0]:
                    # we're in the big infinite box
                    dt_0_hg[r][c] = value
                elif (i, j) in self._dt2pt:
                    # we're in the tail
                    x, y, use_pt = self._dt2pt[i, j]
                    if use_pt:
                        pt_l_hg[x][y] = value
                    else:
                        dt_0_hg[height - x - 1][width - y - 1] = value
                elif i >= len(self.shape):
                    # left strip
                    offset = self._shape.conjugate()[j]
                    dt_0_hg[r + offset][c] = value
                else:
                    # top strip
                    offset = self._shape[i]
                    dt_0_hg[r][c + offset] = value

        # convert back using inverse hg
        dt_0_wrpp = Inverse_HG(dt_0_hg)
        pt_l = Inverse_HG(pt_l_hg)

        # flip dt_0 around so it is a normal plane partition
        dt_rows = [[x for x in reversed(r) if x != 0]
                   for r in reversed(dt_0_wrpp)]
        dt_0 = Tableau(r for r in dt_rows if r)
        return pt_l, dt_0
