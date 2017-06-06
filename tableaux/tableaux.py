# Storage for Young Tableaux.
# Standardness or semi-standardness not implied by data structure, but
# can be checked by supplied methods

import util.order as order

class YoungTable:
    def __init__(self):
        self.rows = []
        self.row_pos = 0 # allows for iteration across table
        self.flipped = False # determines order of iteration

    # construct iterator to iterate over rows (or columns if flipped) of table
    # each time next is called, this iterator will return an iterator over the
    # current row (or column)
    def __iter__(self):
        if self.flipped:

        else:
            for row in self.rows:
                yield row

    def standard(self):
        # check if each row is non-decreasing and each column is increasing
        # check rows
        for row in self.rows:
            if not order.ordered(row, increasing):
                return False
        for i in enumerate(rows):
            
