class LabelledBoxConfiguration(ClonableElement):
    def __init__(self, leg1, leg2, leg3, boxes = None, labels = None, check = True):
        self._parent = LabelledBoxConfigurations(leg1, leg2, leg3)
        self._boxes = dict(boxes) if boxes is not None else dict()
        self._labels = list(labels) if labels is not None else list()
        self._is_immutable = True
        if check:
            self.check()

    def __copy__():
        t = type(self)
        result = type(t).__new__(t)
        result._parent = self._parent
        result._boxes = self._boxes
        result._labels = self._labels
        return result

    def _box_type(i, j, k):
        return sum(1 for p in [(j, k) in self._parent.leg(0).cells(),
                               (i, k) in self._parent.leg(1).cells(),
                               (i, j) in self._parent.leg(2).cells()]
                   if p)

    def _is_valid_box(i, j, k):
        num_negative = sum(1 for l in [i, j, k] if l < 0)
        box_type = _box_type(i, j, k)
        return box_type == 2 or box_type == 3 or box_type == 1 and num_negative == 1

    def _is_filled_if_valid(i, j, k, gen):
        if not _is_valid_box(i, j, k):
            return True
        if (i, j, k) not in self._boxes:
            return False
        label = self._boxes[(i, j, k)]
        if label == -1:
            return True
        assert (label >= 0 and label < len(self._labels),
                "%s is not the index of any label" % label)
        value = self._labels[value]
        return value == gen

    def check():
        """Verify that this object defines a valid box configuration.

        The rules for a valid box configuration are defined in section
        2.5 of [PT2008], starting on page 12.

        """
        for (i, j, k) in self._boxes:
            assert all(l in ZZ for l in (i, j, k)), "Box coordinates must be integers"
            assert _is_valid_box(i, j, k), "%s is not a valid box" % ((i, j, k),)
            box_type = _box_type(i, j, k)
            # TODO: fix up these rules.  Simplify how we check neighbors to make it shorter.
            if box_type == 1:
                assert (self._boxes[(i, j, k)] == -1,
                        "%s is type I and should have no label" % ((i, j, k),))
                assert (_is_filled_if_valid(i + 1, j, k),
                        "%s is filled but %s is not" % ((i, j, k), (i + 1, j, k)))
                assert (_is_filled_if_valid(i, j + 1, k)
                        "%s is filled but %s is not" % ((i, j, k), (i, j + 1, k)))
                assert ((i, j, k + 1) in self._boxes
                        or (j, k + 1) not in self._parent.leg(0).cells(),
                        "%s is filled but %s is not" % ((i, j, k), (i, j, k + 1)))
                assert (self._boxes[(i + 1, j, k)] == -1
                        or self._labels[self._boxes[(i + 1, j, l)]] == 0,
                        "%s is filled but %s has a different label" % ((i, j, k), (i + 1, j, k)))
            elif j < 0:
                assert (self._boxes[(i, j, k)] == -1,
                        "%s is type I and should have no label" % (i, j, k))
                assert ((i, k) in self._parent.leg(1).cells(),
                        "%s is not a valid box" % (i, j, k))
                assert ((i, j + 1, k) in self._boxes
                        or j == -1 and (0, k) not in self._parent.leg(0).cells()
                        or j == -1 and (i, 0) not in self._parent.leg(2).cells(),
                        "%s is filled but %s is not" % ((i, j, k), (i, j + 1, k)))
                assert ((i + 1, j, k) in self._boxes
                        or (i + 1, k) not in self._parent.leg(1).cells(),
                        "%s is filled but %s is not" % ((i, j, k), (i + 1, j, k)))
                assert ((i, j, k + 1) in self._boxes
                        or (i, k + 1) not in self._parent.leg(0).cells(),
                        "%s is filled but %s is not" % ((i, j, k), (i, j, k + 1)))
                assert (self._boxes[(i, j + 1, k)] == -1
                        or self._labels[self._boxes[(i, j + 1, l)]] == 1,
                        "%s is filled but %s has a different label" % ((i, j, k), (i, j + 1, k)))
            elif k < 0:
                assert (self._boxes[(i, j, k)] == -1,
                        "%s is type I and should have no label" % (i, j, k))
                assert ((i, j) in self._parent.leg(2).cells(),
                        "%s is not a valid box" % (i, j, k))
                assert ((i, j, k + 1) in self._boxes
                        or k == -1 and (j, 0) not in self._parent.leg(0).cells()
                        or k == -1 and (i, 0) not in self._parent.leg(1).cells(),
                        "%s is filled but %s is not" % ((i, j, k), (i, j, k + 1)))
                assert ((i + 1, j, k) in self._boxes
                        or (i + 1, j) not in self._parent.leg(2).cells(),
                        "%s is filled but %s is not" % ((i, j, k), (i + 1, j, k)))
                assert ((i, j + 1, k) in self._boxes
                        or (i, j + 1) not in self._parent.leg(2).cells(),
                        "%s is filled but %s is not" % ((i, j, k), (i, j + 1, k)))
                assert (self._boxes[(i, j, k + 1)] == -1
                        or self._labels[self._boxes[(i, j, k + 1)]] == 2,
                        "%s is filled but %s has a different label" % ((i, j, k), (i, j, k + 1)))
            else:
                # All co-ordinates are positive, so we're either II or III
                # Count the number of inclusions
                box_type = sum(1 for p in [(j, k) in self._parent.leg(0).cells(),
                                           (i, k) in self._parent.leg(1).cells(),
                                           (i, j) in self._parent.leg(2).cells()]
                               if p)
                assert box_type == 2 or box_type == 3, "%s is not a valid box" % (i, j, k)

                if box_type == 2:
                    assert (self._boxes[(i, j, k)] == -1,
                            "%s is type II and should have no label" % (i, j, k))
                    assert ((i 
                else:
                    pass

class LabelledBoxConfigurations(Parent):
    Element = LabelledBoxConfiguration

    def __init__(self, leg1, leg2, leg3):
        self.legs = tuple(Partition(leg) for leg in [leg1, leg2, leg3])
        Parent.__init__(self, category=Sets())
