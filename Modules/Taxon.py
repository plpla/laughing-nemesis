#! /usr/bin/python

"""
Class used to construct the taxonomic tree. Use number ID.
Taxon represent a node in the tree.
"""

import threading
import Modules.Error
import FileUtility
import sys
try:
    import cPickle as pickle
except:
    import pickle


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

class Taxon():
    def __init__(self, taxonid):
        self.ID = taxonid
        self.Parent = ""
        self.Kids = []
        self.TaxonName = TaxonName()

    def hasParent(self):
        a = False
        if self.Parent != "":
            a = True
        return a

    def hasKid(self, taxonId):
        a = False
        if taxonId in self.Kids:
                a = 1
        return a

    def getKids(self):
        return self.Kids

    def getTaxonName(self):
        return self.TaxonName

    def getParents(self):
        return self.Parent

    def getId(self):
        return self.ID

    def addKid(self, taxonId):
        if not self.hasKid(taxonId):
            self.Kids.append(taxonId)
        else:
            print(taxonId+" is already a kid")

    def setParent(self, parent):
        if parent != "":
            self.Parent = parent
"""
The taxonomic tree
"""
class TaxonomicTree():
    def __init__(self):
        self.Nodes = {}

    def getNodes(self):
        return self.Nodes

    def nodeExist(self, nodeId):
        detect = False
        try:
            self.Nodes[nodeId]
            detect = True
        except:
            detect = False
        return detect

    def addNode(self, node):
        if self.nodeExist(node.getId()):
            raise ValueError("add_node:node already exist. Node number:"+str(node.getId()))
        self.Nodes[node.getId()] = node

    def getNode(self, node_id):
        if self.nodeExist(node_id):
            return self.Nodes[node_id]
        else:
            raise ValueError("get_node: required node does not exist %s\n" % node_id)

    def dumpTree(self, file_name):
        sys.stderr.write("Dumping tree to disk\n")
        f = open(file_name, 'wb')
        pickle.dump(self.getNodes(), f, protocol=2)

    def loadTree(self, file_name):
        """
        Load a tree from a previously pickle dumped file
        :param file_name: the pickle dump file to load
        :return: Nothing
        """
        sys.stderr.write("Loading tree\n")
        f = open(file_name, 'rb')
        self.Nodes = pickle.load(f)


    def findLCA(self, node1, node2):
        """
        Find the last common ancester of 2 nodes in a given tree
        :param node1: Node #1
        :param node2: Node #2
        :return: The id of the node corresponding to the LCA. If not found raise ValueError since root is the ancestor
        of every node
        """
        if self.nodeExist(node1) and self.nodeExist(node2):
            pass
        else:
            message = "findLCA: "
            if not self.nodeExist(node1):
                message = "Node %s does not exist " % node1
            if not self.nodeExist(node2):
                message += "Node %s does not exist " % node2
            raise ValueError(message)
        parentsOfNode1 = []
        parentsOfNode1.append(node1)
        parentsOfNode2 = []
        parentsOfNode2.append(node2)
        lca = -1
        while self.getNode(node1).hasParent(): #TODO: use root instead of "hasParents()"
            parentToAdd = self.getNode(node1).getParents()
            parentsOfNode1.append(parentToAdd)
            node1 = parentToAdd
        while self.getNode(node2).hasParent():
            parentToAdd = self.getNode(node2).getParents()
            parentsOfNode2.append(parentToAdd)
            node2 = parentToAdd
        for i in parentsOfNode1:
            if i in parentsOfNode2:
                lca = i
                break
        if lca == -1:
            raise ValueError("Have not found the LCA for nodes %s and %s" % (node1, node2))
        return lca

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



    def checkTree(self):
        """
        Confirm the validity of a tree. This can be long and could be used only as a test. Raise ValueError if not valid
        :return: True if tree is valid. Else False
        """
        numOfNode = len(self.getNodes())
        sys.stderr.write(str(numOfNode)+ " to be checked\n")
        checked = 0
        root = []
        for i in self.Nodes:
            #print(self.Nodes[i].getId()==1);
            #if not(self.Nodes[i].getId()):
            #	print("my parents")
            #	print(self.Nodes[i].getParents());
            #	print("my kids")
            #	print(self.Nodes[i].getKids());
            if self.Nodes[i].getParents():
                if self.nodeExist(self.Nodes[i].getParents()):
                    pass
                #print("One parent checked")
                else:
                    raise ValueError("tree is not valid because parent %s does not exist" % self.Nodes[i].getParents())
            else:
                sys.stderr.write("Found a root\n")
                root.append(self.Nodes[i])
            for kid in self.Nodes[i].getKids():
                if self.nodeExist(kid):
                    pass;
                #print("One kid checked");
                else:
                    raise ValueError("Logic error: tree is not valid because kid %s does not exist:" % kid)
            checked += 1
            if checked % 50000 == 0:
                sys.stderr.write("%s nodes checked\n" % checked)
        if len(root) > 1:
            raise ValueError("There is more than one root to tree \n")
        sys.stderr.write("Taxonomic tree is considered valid\n")
        return True


    def readTreeOfLife(self, file):
        FileUtility.isValid(file)
        toBeRead = FileUtility.countLines(file)
        sys.stderr.write(str(toBeRead)+" lines to read to construct tree\n")
        numOfLine = 0
        for line in open(file):
            parent=int(line.split()[0])
            kid=int(line.split()[1])
            if kid == parent:
                sys.stderr.write("Warning: I can't create a link from %s to %s (line %s)\n" % (parent, kid, numOfLine))
                continue
            if self.nodeExist(parent):
                if self.nodeExist(kid):
                    self.getNode(parent).addKid(kid)
                    self.getNode(kid).setParent(parent)
                else:
                    newNode = Taxon(kid)
                    self.addNode(newNode)
                    self.getNode(parent).addKid(kid)
                    self.getNode(kid).setParent(parent)
            else:
                if self.nodeExist(kid):
                    newNode = Taxon(parent)
                    self.addNode(newNode)
                    self.getNode(kid).setParent(parent)
                    self.getNode(parent).addKid(kid)
                else:
                    newKid = Taxon(kid)
                    newParent = Taxon(parent)
                    newKid.setParent(parent)
                    newParent.addKid(kid)
                    self.addNode(newKid)
                    self.addNode(newParent)
            numOfLine += 1
            if numOfLine%100000 == 0:
                sys.stderr.write(str(numOfLine)+" lines read\n")
        sys.stderr.write("Tree base constructed\n")

#test main	
if __name__=="__main__":
    files=sys.argv[1]
    sys.stderr.write("Loading taxonomic tree\n")
    tree=TaxonomicTree()
    tree.readTreeOfLife(files)
    sys.stderr.write("Taxonomic tree is loaded\n")
    sys.stderr.write("Adding taxon name to tree\n")
    tree.addTaxonName(sys.argv[2])
    sys.stderr.write("Taxon name added to tree\n")
    sys.stderr.write("Checking if tree is valid\n")
    tree.checkTree()
    sys.stderr("Writing tree to file")
    tree.dumpTree("Data/Tree.bin")
    tree=TaxonomicTree()
    sys.stderr.write("Checking tree file")
    tree.loadTree("Data/Tree.bin")
    sys.stderr.write("Tree is ready")
    tree.checkTree()
