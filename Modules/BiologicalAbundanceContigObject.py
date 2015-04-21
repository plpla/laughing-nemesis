#!/usr/bin/python

import logging
from ContigObject import *
from ContigIdentificationObject import *
from collections import defaultdict


#possible_taxonomic_level=["species", "phylum", "order", "species group", "kingdom", "class", "genus", "family",
                          #"subspecies"]

possible_taxonomic_level=["phylum", "class", "order", "family", "genus", "species"]

def level_up(branch, tree, name=""):
    """

    :param branch: Branch is  dict[taxon_id]= proportion
    :param tree: The taxonomic tree
    :return: The branch after a level up
    """
    #logging.debug("Initial branch\n")
    #logging.debug(str(branch) + "\n")
    result = defaultdict(float)
    for i in branch:
        if i == "":
            raise ValueError("THAT IS NOT NORMAL!!!")
        if i > 1:
            parent = tree.get_node(i).parent
            #print("Parent of i %s is %s" % (i, parent))
        else:
            parent = i
        result[parent] += branch[i]
    #if name == "contig-328000361":
    #logging.debug("Results from level up:\n")
    #logging.debug(result+"\n")
    return result

"""
class Taxonomy_history():
    def __init__(self):
        self.history = {}
        for level in possible_taxonomic_level:
            self.history[level] = [-1, -1]    #will contain [id, score]
"""

class BiologicalAbundanceContigObject(ContigObject):
    "Class representing a contig and the information about it from RayMeta"
    def __init__(self, name="", lengthInKmer=0, coloredKmers=0, coverage_depth=0):
        ContigObject.__init__(self, name, lengthInKmer, coloredKmers, coverage_depth)
        self.contigIdentifications = []
        self.LCA_id = ""
        self.LCA_name = ""
        self.LCA_score = 0.0
        self.LCA_rank = ""
        self.history = {}
        self.init_history()

    def is_history_set(self):
        v = False
        for level in self.history:
            if self.history[level][0] != -1 and self.history[level][1] != -1:
                v = True
                break #need only one to prove is has been set.
        return v


    def init_history(self):
        for level in possible_taxonomic_level:
            self.history[level] = [-1, -1, ""]  #id, score, name

    def addNewContigIdentification(self, Name, Length, Matches):
        newContigIdentification = ContigIdentificationObject(Name, Length, Matches)
        self.contigIdentifications.append(newContigIdentification)

    def calculatePLvalues(self):
        for entry in self.contigIdentifications:
            entry.calculatePLvalue(ContigObject.getLengthInKmer(self), ContigObject.getColoredKmers(self))

    def selectBestIdentifications(self, numberOfBestId):
        #New step: removePL-value of 0;
        nonNulContigIdentification = []
        for entry in self.contigIdentifications:
            entry.calculatePLvalue(ContigObject.getLengthInKmer(self), ContigObject.getColoredKmers(self))
            if entry.getPLvalue() > 0:
                nonNulContigIdentification.append(entry)
        self.contigIdentifications = nonNulContigIdentification
        #print("Select after select non-=null");
        #print(len(self.contigIdentifications));
        if len(self.contigIdentifications) <= 0:
            return
        else:
            self.contigIdentifications.sort(key=lambda x: x.PLvalue, reverse=True)
            #print("Select after sort");
            #print(len(self.contigIdentifications));
            if len(self.contigIdentifications)>numberOfBestId:
                self.contigIdentifications = self.contigIdentifications[0:numberOfBestId]

    def clean_contig_identifications(self):
        new_contig_id = []
        for i in self.contigIdentifications:
            if i.taxon_name != "Unknown" and i.taxon_name != "":
                new_contig_id.append(i)
        self.contigIdentifications = new_contig_id



    def find_lca(self, tree, converter, verbosity):
        numberOfId = len(self.contigIdentifications)
        #3 cases: 0 id, 1 id and 2 id or more
        if numberOfId == 0:
            if verbosity:
                sys.stderr.write("Case where there is only 0 match:\n")
            self.LCA_id = "No color"
            self.LCA_name = "Unknown"
        if numberOfId == 1:
            if verbosity:
                sys.stderr.write("Case where there is only 1 match:\n")
            if self.contigIdentifications[0].getSequenceName().split('|')[0] == "gi":
                if verbosity:
                    sys.stderr.write(self.contigIdentifications[0].getSequenceName()+"\n")
                id = int(self.contigIdentifications[0].getSequenceName().split('|')[1])
                if converter.genome_is_valid(id):
                    id = converter.convert_to_taxon(id)
                    self.LCA_id = id
                    name = converter.get_taxon_name(id)
                    self.LCA_name = name
                    self.contigIdentifications[0].taxon_name = name
                    self.contigIdentifications[0].taxon_rank = converter.get_taxon_rank(id)
                    self.contigIdentifications[0].taxon_id = id
                else:
                    self.LCA_id = "Invalid id ("+str(id)+")"
                    self.LCA_name = "Unknown"
                    if verbosity:
                        sys.stderr.write("NOT VALID\n")
            else:
                self.LCA_id = "Unknown"
                self.LCA_name = "Unknown"
        if numberOfId > 1:
            if verbosity:
                sys.stderr.write("Case where there is %s match\n" % numberOfId)
            index = 0
            lca = 0
            #1- We set lca to a value that is know to exist in the tree
            while lca == 0 and index < numberOfId:
                if verbosity:
                    sys.stderr.write("Searching LCA starting point with genome: ")
                    sys.stderr.write(self.contigIdentifications[index].getSequenceName()+"\n")
                try:
                    lca = int(self.contigIdentifications[index].getSequenceName().split('|')[1])
                except (IndexError, ValueError) as e:
                    if verbosity:
                        sys.stderr.write("An error caused by the genome name was handled: %s" % e)
                    lca = 0
                    index += 1
                    continue
                validity = converter.genome_is_valid(lca)
                if verbosity:
                    sys.stderr.write("id %s is valid: %s\n" % (lca, validity))
                if not validity:
                    lca = 0
                    index += 1
            if not converter.genome_is_valid(lca):
                self.LCA_id = "Unknown"
                self.LCA_name = "Unknown"
                self.contigIdentifications[0].taxon_name = "Unknown"
                self.contigIdentifications[0].taxon_rank = "Unknown"
                self.contigIdentifications[0].taxon_id = "Unknown"
                # If you are here, you have reached the last possible identification and
                #  have not found one that is valid.
            else:
            #LCA has now a valid value. We can iterate on each entry to find the true lca
                id1 = converter.convert_to_taxon(lca)
                name = converter.get_taxon_name(id1)
                self.contigIdentifications[index].taxon_name = name
                self.contigIdentifications[index].taxon_rank = converter.get_taxon_rank(id1)
                self.contigIdentifications[index].taxon_id = id1

                for identification in self.contigIdentifications:
                    if identification.getSequenceName().split('|')[0] == "gi":
                        if verbosity:
                            sys.stderr.write("Found a GI. Fetching name\n")
                        id2 = int(identification.getSequenceName().split('|')[1])
                        if verbosity:
                            sys.stderr.write("Converting name to taxon id. Checking if valid\n")
                        if converter.genome_is_valid(id2):
                            id2 = converter.convert_to_taxon(id2)
                            name = converter.get_taxon_name(id2)
                            identification.taxon_name = name
                            identification.taxon_rank = converter.get_taxon_rank(id2)
                            identification.taxon_id = id2
                            if verbosity:
                                sys.stderr.write("Searching the LCA in the tree for %s and %s\n" % (id1, id2))
                            lca = tree.find_LCA(id1, id2)
                            if verbosity:
                                sys.stderr.write("Resulting LCA is %s\n" %lca)
                            id1 = lca
                        else:
                            if verbosity:
                                sys.stderr.write("Sequence %s cant be converted\n" % id2)
                            identification.taxon_name = "Unknown"
                            identification.taxon_rank = "Unknown"
                            identification.taxon_id = "Unknown"
                            continue
                    else:
                        if verbosity:
                            sys.stderr.write("%s , %s there is no GI \n" % (self.name, identification.getSequenceName()))
                        identification.taxon_name = "Unknown"
                        identification.taxon_rank = "Unknown"
                        identification.taxon_id = "Unknown"
                        continue #We should do something about it...
            self.LCA_id = lca
            if tree.node_exist(self.LCA_id):
                node = tree.get_node(self.LCA_id)
                if verbosity:
                    tax_name = converter.get_taxon_name(node.id)
                    sys.stderr.write("Setting %s LCA_name to %s by taxon %s\n" % (self.name, tax_name, node.id))
                self.LCA_name = converter.get_taxon_name(node.id)
            else:
                self.LCA_name = "Unknown"


    def compute_history(self, converter, tree):
        if self.LCA_id == "No color":
            #In this case, there is nothing we can do so we set everything to 0.
            for level in self.history:
                self.history[level] = [0, 0, ""]
            return

        self.clean_contig_identifications()
        branch = defaultdict(float)
        matches_sum = 0
        for identification in self.contigIdentifications:
            #print("Match are %s and id are: %s" %(identification.sequenceName, identification.taxon_id))
            if identification.taxon_id != "Unknown":  # Just a precaution
                matches_sum += identification.matchesInContig
        for identification in self.contigIdentifications:
            if identification.taxon_id != "Unknown":
                branch[identification.taxon_id] += float(identification.matchesInContig)/matches_sum
        #We are now ready!
        at_root = False
        while not at_root:
            at_root = True
            for b in branch:
                if b > 1:
                    at_root = False
                try:
                    level = converter.get_taxon_rank(b)
                except ValueError:
                    continue
                if level in possible_taxonomic_level:
                        if branch[b] > self.history[level][1]:
                            self.history[level] = [b, branch[b], ""]
            branch = level_up(branch, tree)
        return


    def set_lca_by_level(self, level):
        """
        requires computed history. Does not set LCA_name
        """
        if self.is_history_set():
            self.LCA_id = self.history[level][0]
            self.LCA_score = self.history[level][1]





    def get_taxonomic_lca(self, threshold=0.9):
        """"
        requires possible_taxonomic_level to be ordered. Does not set LCA_name
        """
        if not self.is_history_set():
            self.LCA_id = "unknown"
            self.LCA_score = 0
            self.LCA_name = "unknown"
            return
        detect = 0
        for level in possible_taxonomic_level:
            if self.history[level][1] >= threshold:
                detect=1
                self.LCA_id = self.history[level][0]
                self.LCA_score = self.history[level][1]
        if detect == 0:
            i = len(possible_taxonomic_level)-1
            level = possible_taxonomic_level[i]
            score = self.history[level][1]
            while i >= 0:
                i = i-1
                level = possible_taxonomic_level[i]
                #sys.stderr.write("%s is current level %s\n" %(level, self.name))
                if score == self.history[level][1]:
                    i += 1
                    level = possible_taxonomic_level[i]
                    self.LCA_id = self.history[level][0]
                    self.LCA_score = self.history[level][1]
                    break
                score = self.history[level][1]







    def refine_lca(self, tree, converter, threshold=0.9):
        #1 for each leaf, check proportion if there is one or many having thresh.
        # dict[taxon] = node, match
        #2 Go up one level. Add the leafs and check.
        #3 repeat!
        sys.stderr.write("Working with %s\n" %self.name)
        if self.LCA_id == "No color":
            self.LCA_score = 0.0
            return
        #print("Contig id len %s" %len(self.contigIdentifications))
        branch = defaultdict(float)
        matches_sum = 0
        for identification in self.contigIdentifications:
            #print("Match are %s and id are: %s" %(identification.sequenceName, identification.taxon_id))
            if identification.taxon_id != "Unknown":  # Just a precaution
                matches_sum += identification.matchesInContig
        #print("%s colored kmers over %s matches" %(matches_sum, len(self.contigIdentifications)))
        for identification in self.contigIdentifications:
            if identification.taxon_id != "Unknown":
                branch[identification.taxon_id] += float(identification.matchesInContig)/matches_sum
                #if self.name == "contig-328000361":
                #print("First Taxon id: %s and score %s" %(identification.taxon_id, float(identification.matchesInContig)/matches_sum))
                if branch[identification.taxon_id] >= threshold:
                    #if self.name == "contig-328000361":
                    #print("Final Taxon id: %s and score %s and treshold: %s" %(identification.taxon_id, branch[identification.taxon_id], threshold))
                    self.LCA_id = identification.taxon_id
                    #The easiest case. One leaf is the best id.
                    self.LCA_score = branch[identification.taxon_id]
                    return
        non_root = True
        num_of_level = 0
        while non_root:
            if num_of_level >= 90:
                sys.stderr.write(branch)
                raise ValueError("SHOULD NOT DO IT!")
            #if self.name == "contig-328000361":
            #print("At level %s" %num_of_level)
            num_of_level += 1
            non_root = False
            #print("Leveling up %s" % self.name)
            #if self.name == "contig-328000361":
            #    branch = level_up(branch, tree, self.name)
            #else:
            branch = level_up(branch, tree)
            for i in branch:
                if i > 1:
                    non_root = True
                if branch[i] >= threshold:
                    self.LCA_id = i
                    #if self.name == "contig-328000361":
                    #print("Second Taxon id: %s and score %s" % (identification.taxon_id, branch[i]))
                    self.LCA_score = branch[i]
                    #print("Found LCA for %s" % self.name)
                    return
        #if self.name == "contig-328000361":
        #print("Reached the end")

    def removeNulIdentification(self):
        contigIdModified = []
        for entry in self.contigIdentifications:
            if entry.getPLvalue() > 0:
                contigIdModified.append(entry)
        return contigIdModified


    def update_LCA_name_and_rank(self, converter):
        try:
            self.LCA_rank = converter.get_taxon_rank(self.LCA_id)
        except ValueError:
            #print(self.name)
            self.LCA_rank = "Unknown rank"
        try:
            self.LCA_name = converter.get_taxon_name(self.LCA_id)
        except ValueError:
            self.LCA_name = "No name"
        for level in self.history:
            try:
                self.history[level][2] = converter.get_taxon_name(self.history[level][0])
            except ValueError:
                if self.history[level][0] > 1:
                    raise ValueError("BROKEN!")
                else:
                    continue

    def showTSV(self):
        #self.calculatePLvalues();
        if len(self.contigIdentifications) != 0:
            #print("Longueur contigIdentifications:"+str(len(self.contigIdentifications)));
            for entry in self.contigIdentifications:
                line = ContigObject.showTSV(self)
                line = line + "\t" + entry.getSequenceName() + "\t" + str(entry.getSequenceLengthInKmer()) + "\t" + str(
                    entry.getMatchesInContig()) + "\t" + str(entry.getPLvalue())
                print(line)

