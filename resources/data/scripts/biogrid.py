import networkx as nx


class BioGRID():

    """Implements Network_parser for BioGrid
       Is a template class, should be innitialized for PPI: BioGRID_PPI or GI : BioGRID_GI"""

    def __init__(self, file_path, tax_id, required_experimental_evidence, experimental_system_type):
        """init Network_parser subclass

        :file_path: path to biogrid file FOR A SPECIFIC ORGANISM 

        """
        self._file_path = file_path
        self._tax_id = tax_id
        self.__required_experimental_evidence = set([ ev.lower() for ev in  required_experimental_evidence])
        self.__required_experimental_system = experimental_system_type

    def get_network(self):
        """get an undirected network from the biogrid file for set experimental codes
        :returns: undirected graph 

        """

        print ("allowed evidence:\t", self.__required_experimental_evidence)
        print ("required system:\t", self.__required_experimental_system)
        G = nx.Graph()
        evidence_accepted = set()
        with open(self._file_path, 'r') as i_stream:
            for line in i_stream:
                line = line.strip().split('\t')
                experimental_system_type = line[12]
                if experimental_system_type.lower() == self.__required_experimental_system:
                    experimental_evidence = line[11]
                    if experimental_evidence.lower() in self.__required_experimental_evidence \
                            or len(self.__required_experimental_evidence) == 0:
                        evidence_accepted.add(experimental_evidence)
                        tax_id_1 = line[15]
                        tax_id_2 = line[16]
                        if tax_id_1 == self._tax_id and tax_id_2 == self._tax_id: 
                            gene_id_1 = line[1]
                            gene_id_2 = line[2]
                            if gene_id_1 != gene_id_2:
                                G.add_edge(gene_id_1, gene_id_2)

        print(evidence_accepted)
        #print ("check for ubiquitin: ", '7316' in G)
        print ("density:\t\t", nx.density(G))
        print ("nodes:\t\t", G.number_of_nodes())
        print ("edges:\t\t", G.number_of_edges())
        #components = list(nx.connected_component_subgraphs(G))
        #print ("nr of components:\t\t", len(components))
        #print ("diameter:\t\t", nx.diameter(max(components, key=len)))

        if len(evidence_accepted) != len(self.__required_experimental_evidence) and len(self.__required_experimental_evidence) != 0:
            raise Warning("did not use all possible evidence codes")

        return G


class BioGRID_PPI(BioGRID):

    """Class to parse experimental PPI"""

    def __init__(self, file_path, tax_id):
        """Initializes the parent BioGRID class to parse experimental PPI"""
        required_experimental_evidence = {'Two-hybrid', 'Affinity Capture-Luminescence',
                                          'Affinity Capture-MS', 'Affinity Capture-RNA', 'Affinity Capture-Western'}
        experimental_system_type = "physical"
        BioGRID.__init__(
            self, file_path, tax_id, required_experimental_evidence, experimental_system_type)


class BioGRID_GI(BioGRID):

    """Class to parse experimental GI"""

    def __init__(self, file_path, tax_id):
        """Initializes the parent BioGRID class to parse experimental GI"""
        required_experimental_evidence = {}
        experimental_system_type = "genetic"
        BioGRID.__init__(
            self, file_path, tax_id, required_experimental_evidence, experimental_system_type)


def main():
    """run the BioGRID_PPI parser, dev purposes"""

    human_biogrid_file_path = "./dummy_data/bioGRID/BIOGRID-ORGANISM-Homo_sapiens-3.5.166.tab2.txt"
    tax_id = '9606'

    biogrid_db = BioGRID_PPI(human_biogrid_file_path, tax_id)
    G_PPI = biogrid_db.get_network()

    biogrid_gi_db = BioGRID_GI(human_biogrid_file_path, tax_id)
    G_GI = biogrid_gi_db.get_network()


if __name__ == "__main__":
    main()
