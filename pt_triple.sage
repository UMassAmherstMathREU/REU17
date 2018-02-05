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
        return gen is not None and value == gen

    def _has_matching_label_or_filled_if_valid(i, j, k, lab):
        if not _is_valid_box(i, j, k):
            return True
        if (i, j, k) not in self._boxes:
            return False
        label = self._boxes[(i, j, k)]
        if label == -1:
            return True
        return label == lab

    def check():
        """Verify that this object defines a valid box configuration.

        The rules for a valid box configuration are defined in section
        2.5 of [PT2008], starting on page 12.

        """
        for (i, j, k) in self._boxes:
            assert all(l in ZZ for l in (i, j, k)), "Box coordinates must be integers"
            assert _is_valid_box(i, j, k), "%s is not a valid box" % ((i, j, k),)
            box_type = _box_type(i, j, k)
            
            if i < 0:
                gen = 0
            elif j < 0:
                gen = 1
            elif k < 0:
                gen = 2
            else:
                gen = None
            
            if box_type == 1 or box_type == 2:
                assert (self._boxes[(i, j, k)] == -1,
                        "%s is type I or II and should have no label" % ((i, j, k),))
                assert (_is_filled_if_valid(i + 1, j, k, gen),
                        "%s is filled but %s is not" % ((i, j, k), (i + 1, j, k)))
                assert (_is_filled_if_valid(i, j + 1, k, gen)
                        "%s is filled but %s is not" % ((i, j, k), (i, j + 1, k)))
                assert (_is_filled_if_valid(i, j, k + 1, gen)
                        "%s is filled but %s is not" % ((i, j, k), (i, j, k + 1))) 
            else:
                label = self._boxes[(i, j, k)]
                assert (label >= 0 and label < len(self._lables),
                        "%s is not the index of any label" % label)
                assert (_has_matching_label_or_filled_if_valid(i + 1, j, k, label),
                        "%s and %s do not have matching labels" % ((i, j, k), (i + 1, j, k)))
                assert (_has_matching_label_or_filled_if_valid(i, j + 1, k, label),
                        "%s and %s do not have matching labels" % ((i, j, k), (i, j + 1, k)))
                assert (_has_matching_label_or_filled_if_valid(i, j, k + 1, label),
                        "%s and %s do not have matching labels" % ((i, j, k), (i, j, k + 1)))

class LabelledBoxConfigurations(Parent):
    Element = LabelledBoxConfiguration

    def __init__(self, leg1, leg2, leg3):
        self.legs = tuple(Partition(leg) for leg in [leg1, leg2, leg3])
        Parent.__init__(self, category=Sets())
