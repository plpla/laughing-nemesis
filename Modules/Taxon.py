#! /usr/bin/python

"""
Class used to construct the taxonomic tree. Use number ID.
Taxon represent a node in the tree.
"""

import FileUtility
import sys
try:
    import cPickle as pickle
except:
    import pickle

"""
class TaxonName():
    def __init__(self):
        self.Name = ""
        self.Rank = ""

    def getName(self):
        return self.Name

    def getRank(self):
        return self.Rank

    def setAll(self, name, rank):
        self.Name = name
        self.Rank = rank

    def setName(self, name):
        self.Name = name

    def setRank(self, rank):
        self.Rank = rank
"""


class Taxon():
    def __init__(self, taxonid):
        self._ID = taxonid
        self._parent = -1
        self._kids = []

    @property
    def id(self):
        return self._ID

    @property
    def parent(self):
        return self._parent

    @property
    def kids(self):
        return self._kids

    def has_parent(self):
        if self.parent != -1:
            return True
        return False

    def is_kid(self, taxon):
        if taxon in self.kids:
                return True
        return False

    def add_kid(self, taxon):
        if not taxon in self.kids:
            self._kids.append(taxon)

    def set_parent(self, taxon):
        if self.parent == -1 or self.parent == taxon:
            self._parent = taxon
        else:
            raise ValueError("Can't set a second parent")


"""
The taxonomic tree
"""
class Taxonomic_tree():
    def __init__(self):
        self._nodes = {}

    @property
    def nodes(self):
        return self._nodes

    def node_exist(self, node):
        detect = False
        try:
            self.nodes[node]
            detect = True
        except KeyError:
            detect = False
        return detect

    def create_node(self, node):
        if self.node_exist(node.id):
            raise ValueError("add_node:node already exist. Node number:"+str(node.id))
        self._nodes[node.id] = node

    def get_node(self, node_id):
        try:
            return self.nodes[node_id]
        except:
            raise ValueError("get_node: required node does not exist "+str(node_id))

    """
    def dumpTree(self, file_name):
        sys.stderr.write("Dumping tree to disk\n")
        f = open(file_name, 'wb')
        pickle.dump(self.getNodes(), f, protocol=2)
    """

    # def loadTree(self, file_name):
    #    """
    #    Load a tree from a previously pickle dumped file
    #    :param file_name: the pickle dump file to load
    #    :return: Nothing
    #    """
    #   sys.stderr.write("Loading tree\n")
    #    f = open(file_name, 'rb')
    #    self.Nodes = pickle.load(f)


    def find_LCA(self, node1, node2):
        """
        Find the last common ancester of 2 nodes in a given tree
        :param node1: Node #1
        :param node2: Node #2
        :return: The id of the node corresponding to the LCA. If not found raise ValueError since root is the ancestor
        of every node
        """
        if self.node_exist(node1) and self.node_exist(node2):
            pass
        else:
            message = "findLCA: "
            if not self.node_exist(node1):
                message = "Node %s does not exist " % node1
            if not self.node_exist(node2):
                message += "Node %s does not exist " % node2
            raise ValueError(message)

        family_node1 = [node1]
        family_node2 = [node2]

        lca = -1
        while self.get_node(node1).has_parent(): #TODO: use root instead of "hasParents()"
            parent = self.get_node(node1).parent
            family_node1.append(parent)
            node1 = parent
        while self.get_node(node2).has_parent():
            parent = self.get_node(node2).parent
            family_node2.append(parent)
            node2 = parent
        for i in family_node1:
            if i in family_node2:
                lca = i
                break
        if lca == -1:
            raise ValueError("Have not found the LCA for nodes %s and %s" % (node1, node2))
        return lca

    """
    def addTaxonName(self, file):
        sys.stderr.write("Adding taxon name to tree\n")
        FileUtility.isValid(file)
        for line in open(file):
            id = int(line.split('\t')[0])
            name = line.split('\t')[1]
            rank = line.split('\t')[2].rstrip('\n')
            if self.nodeExist(id):
                self.getNodes()[id].getTaxonName().setAll(name, rank)
            else:
                sys.stderr.write("Warning: Taxon %s is missing from tree\n" % id)
    """


    def check_tree(self):
        """
        Confirm the validity of a tree. This can be long and could be used only as a test. Raise ValueError if not valid
        :return: True if tree is valid. Else False
        """
        num_nodes = len(self.nodes)
        sys.stderr.write(str(num_nodes)+ " to be checked\n")
        checked = 0
        root = []
        for i in self.nodes:
            if self.nodes[i].parent >= 1:
                if self.node_exist(self.nodes[i].parent):
                    pass
                else:
                    raise ValueError("tree is not valid because parent %s does not exist" % self.nodes[i].parent)
            else:
                sys.stderr.write("Found a root\n")
                root.append(self.nodes[i].id)
            #else:
            #    raise ValueError("Node "+str(self.nodes[i].id)+" has no parents")

            for kid in self.nodes[i].kids:
                if self.node_exist(kid):
                    pass
                else:
                    raise ValueError("Logic error: tree is not valid because kid %s does not exist:" % kid)
            checked += 1
            if checked % 50000 == 0:
                sys.stderr.write("%s nodes checked\n" % checked)
        if len(root) != 1:
            sys.stderr.write("Num of root: "+str(len(root)))
            print("Roots id:")
            for i in root:
                print(i)
            raise ValueError("There is more than one root or no root to tree \n")
        sys.stderr.write("Taxonomic tree is considered valid\n")
        return True


    def read_tree_of_life(self, file_in):
        FileUtility.isValid(file_in)
        to_read = FileUtility.countLines(file_in)
        sys.stderr.write(str(to_read)+" lines to read to construct tree\n")
        num_line = 0
        for line in open(file_in):
            parent = int(line.split()[0])
            kid = int(line.split()[1])
            if kid == parent:
                sys.stderr.write("Warning: I can't create a link from %s to %s (line %s)\n" % (parent, kid, num_line))
                continue
            if self.node_exist(parent):
                if self.node_exist(kid):
                    self.get_node(parent).add_kid(kid)
                    self.get_node(kid).set_parent(parent)
                else:
                    new_node = Taxon(kid)
                    self.create_node(new_node)
                    self.get_node(parent).add_kid(kid)
                    self.get_node(kid).set_parent(parent)
            else:
                if self.node_exist(kid):
                    new_node = Taxon(parent)
                    self.create_node(new_node)
                    self.get_node(kid).set_parent(parent)
                    self.get_node(parent).add_kid(kid)
                else:
                    new_kid = Taxon(kid)
                    new_parent = Taxon(parent)
                    new_kid.set_parent(parent)
                    new_parent.add_kid(kid)
                    self.create_node(new_kid)
                    self.create_node(new_parent)
            num_line += 1
            if num_line % 100000 == 0:
                sys.stderr.write(str(num_line)+" lines read\n")
        sys.stderr.write("Tree base constructed\n")

#test main	
if __name__=="__main__":
    file_in=sys.argv[1]
    sys.stderr.write("Loading taxonomic tree\n")
    tree = Taxonomic_tree()
    tree.read_tree_of_life(file_in)
    sys.stderr.write("Taxonomic tree is loaded\n")
    sys.stderr.write("Checking if tree is valid\n")
    tree.check_tree()
    sys.stderr.write("Tree is ready")
