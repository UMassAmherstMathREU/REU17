from sage.structure.list_clone import ClonableElement
from sage.structure.unique_representation import CachedRepresentation
from collections import deque

def cy_pt_vertex_series(mu1, mu2, mu3, prec=5):
    """Compute the Calibi-Yau vertex.

    This computes the series \(W^P_{\vec \mu} |_{s_1+s_2+s_3=0}\),
    using the conjectured equation (5.3) in [PT2008].

    The arguments mu1,mu2,mu3 are iterables defining the outgoing
    partitions.  The argument prec determines the number of terms to
    compute.

    """

    R.<q> = LaurentSeriesRing(ZZ)
    cfgs = LabelledBoxConfigurations(mu1, mu2, mu3)
    v = cfgs.renormalized_volume()
    # If I understand correctly:
    # We need the euler characteristic of (P^1)^n, which is 2^n
    return sum(2^p.unrestricted_components() * (-q)^(p.length() + v)
               for p in cfgs.up_to_size_n(prec - 1)) + O(q^(prec + v))
    
    

### Helper functions used internally

def higher_boxes(w):
    return [tuple(w[j] + 1 if i == j else w[j] for j in range(3))
            for i in range(3)]

def lower_boxes(w):
    return [tuple(w[j] - 1 if i == j else w[j] for j in range(3))
            for i in range(3)]

def adjacent_boxes(w):
    return higher_boxes(w) + lower_boxes(w)

def remove_nth(w, n):
    w = list(w)
    del w[n]
    return tuple(w)

class LabelledBoxConfiguration(ClonableElement):

    """Class to represent a labelled box configuration.

    Internally, (0,1,2) represent labels fixed in each direction, and
    3 represents a "free label", and -1 represents represents a
    supported box with no label.

    """
    
    def __init__(self, leg1, leg2, leg3, boxes = None, labels = None, check = True):
        """Create a new labelled box configuration with specified data.

        """
        self._parent = LabelledBoxConfigurations(leg1, leg2, leg3)
        self._boxes = dict(boxes) if boxes is not None else dict()
        self._is_immutable = True
        if check:
            self.check()
        ClonableElement.__init__(self, self._parent)

    def __copy__(self):
        """Create a copy of this box configuration.

        """
        t = type(self)
        result = t.__new__(t)
        result._parent = self._parent
        result._boxes = dict(self._boxes)
        return result

    def _is_filled_if_valid(self, w, gen = -1):
        """Check if the box at w contains the subspace induced by gen.

        More precisely, check that \(x_1^i x_2^j x_3^k \cdot
        \mathbf{gen}_0\) where \(w = (i,j,k)\) is in the corresponding
        module.  If (i,j,k) is an invalid location, this product is
        \(0\) and we always return true.

        If gen is -1, then check if (i,j,k) is supported by an
        unlabelled box, or is not a valid location for a box.

        If gen is 0, 1, or 2, then check that (i,j,k) is either
        unlabelled, labelled by gen, is a type II box not in the
        cylinder corresponding to gen, or is not a valid location for a box.

        If gen is 3, then check that (i,j,k) is unlabelled, contains
        an unrestricted label, or is not a valid location for a box.

        """
        if not self._parent.is_valid_box(*w):
            return True
        if w in self._boxes and self._boxes[w] in (-1, gen):
            return True
        if gen in (0,1,2) and remove_nth(w, gen) not in self._parent.leg(gen).cells():
            return True
        return False

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
        for w in self._boxes:
            assert all(l in ZZ for l in w), "Box coordinates must be integers"
            assert self._parent.is_valid_box(*w), "%s is not a valid box" % (w,)
            box_type = self._parent.box_type(*w)

            if box_type == 1 or box_type == 2:
                assert self._boxes[w] == -1, \
                    "%s is type I or II and should have no label" % (w,)

            if box_type == 1:
                gen = [i < 0 for i in w].index(True)
            else:
                gen = self._boxes[w]
            
            for v in higher_boxes(w):
                assert self._is_filled_if_valid(v, gen), \
                    "%s is filled but %s is not" % (w, v)

    def length(self):
        """Return the length of the configuration.

        The length of the configuration, as defined in [PT2008].  This
        corresponds to the length of a submodule.

        """
        return sum(2 if self._parent.box_type(*w) == 3 and self._boxes[w] == -1 else 1
                   for w in self._boxes)

    def unrestricted_components(self):
        """Return the number of unrestricted components.

        An unrestricted component is a connected component of labelled
        type III boxes, which are not forced to have any specific
        value for their label.  If there are $n$ unrestricted
        components, then the configuration corresponds to a family of
        submodules parameterized by $(\mathbb P^1)^n$.

        """
        unvisited = {w for w in self._boxes
                     if self._boxes[w] == 3}
        to_visit = deque()
        components = 0
        while unvisited:
            components += 1
            w = iter(unvisited).next()
            unvisited.remove(w)
            to_visit.append(w)
            while to_visit:
                w = to_visit.pop()
                for v in adjacent_boxes(w):
                    if v in unvisited:
                        unvisited.remove(v)
                        to_visit.append(v)
        return components

    def fill_box(self, w, value = -1):
        """Update the value in a box.

        By default, use an unlabelled box.

        """
        self._require_mutable()
        self._boxes[w] = value

    def force(self, w, lab):
        """Force the value of a type III box and propogate changes.

        If the box at w is an unrestricted label, fill it in with the
        value of lab.  If it is empty, or unlabelled, do nothing.

        """
        self._require_mutable()
        if w in self._boxes:
            if self._boxes[w] == 3:
                self._boxes[w] = lab
                for v in adjacent_boxes(w):
                    self.force(v, lab)

    def relax(self, w):
        """Relax the value of a type III box to an unrestricted label.

        Check the entire connected component of this box for any
        external constraints from type I and II boxes, and then relax
        if possible.

        """
        self._require_mutable()
        component = set()
        to_check = deque([w])
        while to_check:
            v = to_check.pop()
            if v not in self._boxes:
                if self._parent.box_type(*v) == 2:
                    return
            elif self._parent.box_type(*v) == 1:
                return
            elif self._boxes[v] == self._boxes[w]:
                component.add(v)
                for u in adjacent_boxes(v):
                    if u not in component:
                        to_check.append(u)
        # If we made it this far, it's safe to convert the whole
        # component
        for v in component:
            self._boxes[v] = 3

    def children(self):
        """Find all possible configurations which can be obtained by adding a
        box

        More precisely, return a list of labelled box configurations
        which can be obtained from this one by one of the following
        operations:

        - Adding a type I or type II box
        - Adding a labelled type III box
        - Transforming a labelled type III box into an unlabelled type III box

        """
        result = []

        # Add a type I
        for i in range(3):
            for u in self._parent.leg(i).cells():
                # Try to insert into the column corresponding to u in
                # the ith leg
                w = list(u)
                w.insert(i, -1)
                while tuple(w) in self._boxes:
                    w[i] -= 1
                w = tuple(w)
                # We found a candidate spot for a new box
                # Check if neighbors are unrestricted labels, or i-labels
                if all(self._is_filled_if_valid(v, 3) or
                       self._is_filled_if_valid(v, i)
                       for v in higher_boxes(w)):
                    with self.clone() as child:
                        child.fill_box(w)
                        # We may need to restrict an unrestricted label
                        v = tuple(w[j] + 1 if j == i else w[j] for j in range(3))
                        child.force(v, i)
                    result.append(child)
                    
        # Add a type II
        for w in self._parent.type_II_boxes():
            if w not in self._boxes and all(self._is_filled_if_valid(v)
                                            for v in higher_boxes(w)):
                with self.clone() as child:
                    child.fill_box(w)
                    for v in lower_boxes(w):
                        child.relax(v)
                result.append(child)

        # Add labelled type III or remove a label
        for w in self._parent.type_III_boxes():
            if w not in self._boxes:
                # Try to add a labelled box
                choices = [i for i in range(3)
                           if all(self._is_filled_if_valid(v, i) for v in higher_boxes(w))]
                assert len(choices) != 2, "there should always be 0,1, or 3 choices"
                
                if len(choices) == 1:
                    # Theres only one possible fixed label.
                    # Upper blocks may be free, so we need to force them to be fixed.
                    with self.clone() as child:
                        child.fill_box(w, choices[0])
                        for v in higher_boxes(w):
                            child.force(v, choices[0])
                    result.append(child)
                elif len(choices) == 3:
                    # Use a free label
                    with self.clone() as child:
                        child.fill_box(w, 3)
                    result.append(child)
                    
            elif self._boxes[w] != -1:
                # Try to remove a label
                if all(self._is_filled_if_valid(v) for v in higher_boxes(w)):
                    with self.clone() as child:
                        child.fill_box(w)
                        for v in lower_boxes(w):
                            child.relax(v)
                    result.append(child)
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
                self._boxes == other._boxes and
                all((i, j, k) in other._boxes and
                    (self._boxes[(i, j, k)] == -1) == (other._boxes[(i, j, k)] == -1)
                    for (i, j, k) in self._boxes))

    def __hash__(self):
        return hash((self._parent, frozenset(self._boxes.items())))

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
                    else:
                        value = self._boxes[(row, col, layer)]
                        if value == -1:
                            result += "|     "
                        elif value == 3:
                            result += "| *?* "
                        else:
                            result += "| *%d* " % (value + 1)
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

    def renormalized_volume(self):
        return -len(self.type_II_boxes()) - 2 * len(self.type_III_boxes())

    def of_size_n(self, n):
        S = {LabelledBoxConfiguration(self._legs[0], self._legs[1], self._legs[2])}
        for i in range(n):
            S = {y for x in S for y in x.children()}
        return S

    def up_to_size_n(self, n):
        S = {LabelledBoxConfiguration(self._legs[0], self._legs[1], self._legs[2])}
        result = S
        for i in range(n):
            S = {y for x in S for y in x.children()}
            result = result.union(S)
        return result
        

