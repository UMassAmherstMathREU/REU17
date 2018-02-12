from sage.structure.list_clone import ClonableElement
from sage.structure.unique_representation import CachedRepresentation
from collections import deque

class LabelledBoxConfiguration(ClonableElement):
    def __init__(self, leg1, leg2, leg3, boxes = None, labels = None, check = True):
        self._parent = LabelledBoxConfigurations(leg1, leg2, leg3)
        self._boxes = dict(boxes) if boxes is not None else dict()
        self._labels = list(labels) if labels is not None else list()
        self._is_immutable = True
        if check:
            self.check()
        ClonableElement.__init__(self, self._parent)

    def __copy__(self):
        t = type(self)
        result = t.__new__(t)
        result._parent = self._parent
        result._boxes = dict(self._boxes)
        result._labels = list(self._labels)
        return result

    def _is_filled_if_valid(self, i, j, k, gen):
        if not self._parent.is_valid_box(i, j, k):
            return True
        if (i, j, k) not in self._boxes:
            return False
        label = self._boxes[(i, j, k)]
        if label == -1:
            return True
        assert label >= 0 and label < len(self._labels), \
            "%s is not the index of any label" % label
        value = self._labels[label]
        return gen is not None and value == gen

    def _has_matching_label_or_filled_if_valid(self, i, j, k, lab):
        if not self._parent.is_valid_box(i, j, k):
            return True
        value = self._labels[lab]
        # If child is type II and independent from us, return true
        if (value == 0 and (j, k) not in self._parent.leg(0).cells() or
            value == 1 and (i, k) not in self._parent.leg(1).cells() or
            value == 2 and (i, j) not in self._parent.leg(2).cells()):
            return True
        # Otherwise, there had better be a box there
        if (i, j, k) not in self._boxes:
            return False
        # The child must be filled or matching
        label = self._boxes[(i, j, k)]
        return label == -1 or label == lab

    def check(self):
        """Verify that this object defines a valid box configuration.

        The rules for a valid box configuration are defined in section
        2.5 of [PT2008], starting on page 12.

        This method does not check that there are no unused labels, or
        that the same label is not used twice for no good reason.  If
        either of these issues are present, then
        unrestricted_components() will not function properly.
        However, if the configuration is otherwise valid except for
        one of these issues, then there is another configuration
        without these issues that is the same up to renaming the
        unrestricted components, and possibly removing unused labels,
        to which this configuration will compare equal.
        """
        for label in self._labels:
            assert label in ZZ, "%s is not an integer and cannot be a label" % label
            assert label >= -1 and label <= 2, "%s is not a valid label" % label

        for (i, j, k) in self._boxes:
            assert all(l in ZZ for l in (i, j, k)), "Box coordinates must be integers"
            assert self._parent.is_valid_box(i, j, k), "%s is not a valid box" % ((i, j, k),)
            box_type = self._parent.box_type(i, j, k)

            if i < 0:
                gen = 0
            elif j < 0:
                gen = 1
            elif k < 0:
                gen = 2
            else:
                gen = None

            if box_type == 1 or box_type == 2:
                assert self._boxes[(i, j, k)] == -1, \
                    "%s is type I or II and should have no label" % ((i, j, k),)
                assert self._is_filled_if_valid(i + 1, j, k, gen), \
                    "%s is filled but %s is not" % ((i, j, k), (i + 1, j, k))
                assert self._is_filled_if_valid(i, j + 1, k, gen), \
                    "%s is filled but %s is not" % ((i, j, k), (i, j + 1, k))
                assert self._is_filled_if_valid(i, j, k + 1, gen), \
                    "%s is filled but %s is not" % ((i, j, k), (i, j, k + 1))
            else:
                label = self._boxes[(i, j, k)]
                assert label >= -1 and label < len(self._labels), \
                    "%s is not the index of any label" % label
                assert self._has_matching_label_or_filled_if_valid(i + 1, j, k, label), \
                    "%s is filled but %s is not, or has the wrong label" % ((i, j, k), (i + 1, j, k))
                assert self._has_matching_label_or_filled_if_valid(i, j + 1, k, label), \
                    "%s is filled but %s is not, or has the wrong label" % ((i, j, k), (i, j + 1, k))
                assert self._has_matching_label_or_filled_if_valid(i, j, k + 1, label), \
                    "%s is filled but %s is not, or has the wrong label" % ((i, j, k), (i, j, k + 1))

    def length(self):
        """Return the length of the configuration.

        The length of the configuration, as defined in [PT2008].  This
        corresponds to the length of a submodule.

        """
        return sum(2 if self._parent.box_type(i, j, k) == 3 and self._boxes[(i, j, k)] == -1 else 1
                   for (i, j, k) in self._boxes)

    def unrestricted_components(self):
        """Return the number of unrestricted components.

        An unrestricted component is a connected component of labelled
        type III boxes, which are not forced to have any specific
        value for their label.  If there are $n$ unrestricted
        components, then the configuration corresponds to a family of
        submodules parameterized by $(\mathbb P^1)^n$.

        """
        return sum(1 for i in self._labels if i == -1)

    def fill_box(self, i, j, k):
        self._require_mutable()
        self._boxes[(i, j, k)] = -1

    def set_box_reference(self, i, j, k, label = None):
        self._require_mutable()
        if label is None:
            self._labels.append(-1)
            label = len(self._labels) - 1
        self._boxes[(i, j, k)] = label
        return label

    def set_box_value(self, i, j, k, value):
        self._require_mutable()
        self._labels.append(value)
        self._boxes[(i, j, k)] = len(self._labels) - 1

    def propogate_labels(self, i, j, k):
        self._require_mutable()
        q = deque([(i + 1, j, k), (i, j + 1, k), (i, j, k + 1)])
        label = self._boxes[(i, j, k)]
        while q:
            i, j, k = q.popleft()
            if (i, j, k) not in self._boxes:
                continue
            if (i, j, k) in self._boxes and self._boxes[(i, j, k)] != -1 \
               and self._boxes[(i, j, k)] != label:
                self._boxes[(i, j, k)] = label
                q.append((i + 1, j, k))
                q.append((i, j + 1, k))
                q.append((i, j, k + 1))

    def _update_constraints(self, i, j, k, constraints, axis):
        new_box = [i, j, k]
        new_box[axis] += 1
        new_box = tuple(new_box)
        i, j, k = new_box

        if new_box not in self._boxes:
            t = self._parent.box_type(i, j, k)
            if t == 3:
                constraints[:] = [True, True, True]
            if t == 2:
                # Determine the box that the type II cell is not in
                notin = [(j, k) in self._parent.leg(0).cells(),
                         (i, k) in self._parent.leg(1).cells(),
                         (i, j) in self._parent.leg(2).cells()].index(False)
                constraints[notin] = True
        elif self._boxes[new_box] != -1:
            value = self._labels[self._boxes[(i, j, k)]]
            if value != -1:
                constraints[value] = True

    def children(self):
        """Find all possible configurations which can be obtained by adding a box

        More precisely, return a list of labelled box configurations
        which can be obtained from this one by one of the following
        operations:

        - Adding a type I or type II box
        - Adding a labelled type III box
        - Transforming a labelled type III box into an unlabelled type III box

        """
        result = []
        # Add a type I
        # Add a type II
        for (i, j, k) in self._parent.type_II_boxes():
            if (i, j, k) not in self._boxes and \
               self._is_filled_if_valid(i + 1, j, k, None) and \
               self._is_filled_if_valid(i, j + 1, k, None) and \
               self._is_filled_if_valid(i, j, k + 1, None):
                with self.clone() as child:
                    child.fill_box(i, j, k)
                    # NEED TO RELAX LABELS IF POSSIBLE
                    # This could get tricky with chain reaction stuff happening
                    # Maybe just re-compute all type III labels on creating a new child?
                result.append(child)
        # Add a labelled type III
        for (i, j, k) in self._parent.type_III_boxes():
            if (i, j, k) not in self._boxes:
                # We need to determine the possible labels at a given position
                constraints = [False, False, False]
                for axis in range(3):
                    self._update_constraints(i, j, k, constraints, axis)
                num_constraints = sum(1 for c in constraints if c)
                if num_constraints == 0:
                    # Add a free label (unrestricted component)
                    with self.clone() as child:
                        new_value = child.set_box_reference(i, j, k)
                        to_replace = [child._boxes[box]
                                      for box in [(i + 1, j, k), (i, j + 1, k), (i, j, k + 1)]
                                      if box in child._boxes and child._boxes[box] != -1]
                        for box in child._boxes:
                            if child._boxes[box] in to_replace:
                                child._boxes[box] = new_value
                    result.append(child)
                elif num_constraints == 1:
                    # Add a labelled box with specified label
                    c = constraints.index(True)
                    with self.clone() as child:
                        child.set_box_value(i, j, k, c)
                        child.propogate_labels(i, j, k)
                    result.append(child)

        # Remove a label from a type III
        for (i, j, k) in self._parent.type_III_boxes():
            if (i, j, k) in self._boxes and self._boxes[(i, j, k)] != -1 and \
               self._is_filled_if_valid(i + 1, j, k, None) and \
               self._is_filled_if_valid(i, j + 1, k, None) and \
               self._is_filled_if_valid(i, j, k + 1, None):
                with self.clone() as child:
                    child.fill_box(i, j, k)
                result.append(child)
        # TODO
        return result

    def __eq__(self, other):
        """Check if two box configurations are equivalent

        Two box configurations are equivalent if they have the same
        set of labelled and unlabelled boxes.  Note that this
        information is enough to uniquely determine the "most generic"
        equivalent box configuration, up to relabelling of the
        unrestricted components.

        """
        # If all of our boxes are contained in theirs, and they have the same length,
        # then they're the same
        return (isinstance(other, type(self)) and
                self._parent == other._parent and
                len(self._boxes) == len(other._boxes) and
                all((i, j, k) in other._boxes and
                    (self._boxes[(i, j, k)] == -1) == (other._boxes[(i, j, k)] == -1)
                    for (i, j, k) in self._boxes))

    def __hash__(self):
        unlabelled_boxes = frozenset(b for b in self._boxes if self._boxes[b] == -1)
        labelled_boxes = frozenset(b for b in self._boxes if self._boxes[b] != -1)
        return hash((labelled_boxes, unlabelled_boxes))

    def _repr_(self):
        return "Labelled box configuration of length %s with outgoing partitions %s" \
            % (self.length(), self._parent.legs())

    def _ascii_art_(self):
        """Print cross-sections of the configuration.

        """
        # This is the messiest code in the world.  There's definately a better way to do this.
        # Working on this is giving me a headache.  I think it works reasonably well.
        legs = self._parent.legs()
        top = max(legs[0][0] if legs[0] else 0,
                  legs[1][0] if legs[1] else 0) - 1
        front = max(len(legs[1]), len(legs[2])) - 1
        right = max(len(legs[0]), legs[2][0] if legs[2] else 0) - 1
        back = min([i for (i, j, k) in self._boxes] + [-1])
        left = min([j for (i, j, k) in self._boxes] + [-1])
        bottom = min([k for (i, j, k) in self._boxes] + [-1])
        lastCol = False
        result = ""
        for layer in range(top, bottom - 1, -1):
            result += "Layer %s:\n" % layer
            for row in range(back, front + 2):
                # Top line
                for col in range(left, right + 1):
                    if (not self._parent.is_valid_box(row, col, layer)
                        and not self._parent.is_valid_box(row - 1, col, layer)
                        and not self._parent.is_valid_box(row, col - 1, layer)
                        and not self._parent.is_valid_box(row - 1, col - 1, layer)):
                        result += "      "
                    elif (row, col, layer) not in self._boxes and (row - 1, col, layer) not in self._boxes:
                        result += "+     "
                    else:
                        result += "+-----"
                if (self._parent.is_valid_box(row, right, layer) or
                    self._parent.is_valid_box(row - 1, right, layer)):
                    result += "+"
                result += "\n"
                # Next Line
                for col in range(left, right + 1):
                    if not self._parent.is_valid_box(row, col, layer):
                        if (row, col - 1, layer) in self._boxes:
                            result += "|     "
                        else:
                            result += "      "
                    else:
                        if (row, col, layer) not in self._boxes and (row, col - 1, layer) not in self._boxes:
                            result += " "
                        else:
                            result += "|"
                        result += { 1: "  I  ", 2: " II  ", 3: " III "}[self._parent.box_type(row, col, layer)]
                if (row, right, layer) in self._boxes:
                    result += "|"
                result += "\n"
                # Last line
                for col in range(left, right + 1):
                    if (row, col, layer) not in self._boxes:
                        if (row, col - 1, layer) not in self._boxes:
                            result += "      "
                        else:
                            result += "|     "
                    elif self._boxes[(row, col, layer)] == -1:
                        result += "|     "
                    else:
                        # Labelled box.  What to do?
                        label = self._boxes[(row, col, layer)]
                        if self._labels[label] == -1:
                            result += "| L_%-2d" % label
                        else:
                            result += "| *%d* " % (self._labels[label] + 1)
                if (row, right, layer) in self._boxes:
                    result += "|"
                result += "\n"
        return ascii_art(result)

class LabelledBoxConfigurations(Parent, CachedRepresentation):
    Element = LabelledBoxConfiguration

    @staticmethod
    def __classcall__(cls, leg1, leg2, leg3):
        leg1 = Partition(leg1)
        leg2 = Partition(leg2)
        leg3 = Partition(leg3)
        return super(LabelledBoxConfigurations, cls).__classcall__(cls, leg1, leg2, leg3)

    def __init__(self, leg1, leg2, leg3):
        self._legs = (leg1, leg2, leg3)
        Parent.__init__(self, category=Sets())

    def legs(self):
        return self._legs

    def leg(self, i):
        return self._legs[i]

    def box_type(self, i, j, k):
        return sum(1 for p in [(j, k) in self._legs[0].cells(),
                               (i, k) in self._legs[1].cells(),
                               (i, j) in self._legs[2].cells()]
                   if p)

    def is_valid_box(self, i, j, k):
        num_negative = sum(1 for l in [i, j, k] if l < 0)
        box_type = self.box_type(i, j, k)
        return box_type == 2 or box_type == 3 or box_type == 1 and num_negative == 1

    def _type_n_boxes(self, n):
        max_x = max(len(self._legs[1]), len(self._legs[2]))
        max_y = max(self._legs[0][0] if self._legs[0] else 0,
                    len(self._legs[2]))
        max_z = max(self._legs[0][0] if self._legs[0] else 0,
                    self._legs[1][0] if self._legs[1] else 0)
        return [(i, j, k)
                for i in range(max_x)
                for j in range(max_y)
                for k in range(max_z)
                if self.box_type(i, j, k) == n]

    def type_II_boxes(self):
        return self._type_n_boxes(2)

    def type_III_boxes(self):
        return self._type_n_boxes(3)

    def of_size_n(self, n):
        S = {LabelledBoxConfiguration(self._legs[0], self._legs[1], self._legs[2])}
        for i in range(n):
            S = {y for x in S for y in x.children()}
        return S

