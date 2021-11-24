import pandas as pd
import os
import networkx as nx

class CoexDB(object):

    """Parse coexpress db data dump"""

    def __init__(self, coexdb_rootdir):
        """TODO: to be defined1.

        :coexdb_filepath: (str). String pointing to the root -dir of the Coexdb
        data dump.

        """
        self._coexdb_rootdir = coexdb_rootdir
        self.genes = self._genes()

    def _genes(self):
        """Get all entrez ids

        :returns: list(str)

        """
        return os.listdir(self._coexdb_rootdir)

    def coexpressed(self, gene, threshold, genes_subset=None, legacy=False):
        """ For a given gene, get top x% coexpressed genes.

        :gene: (str). Entrez ID.
        :threshold: (float). Determines top number of genes to keep.
        :returns: list(str). List of top coexpressed gene IDs.

        legacy = version <=6
        """
        
        if legacy:
            df = pd.read_csv('{}{}{}'.format(self._coexdb_rootdir, '/', gene),
                         header=None, names=['Entrez ID','legacy column', 'value'], sep='\t',
                         dtype={'Entrez ID': 'str'})

            df = df.sort_values(by='value', ascending=False)

        else:
            df = pd.read_csv('{}{}{}'.format(self._coexdb_rootdir, '/', gene),
                         header=None, names=['Entrez ID','value'], sep='\t',
                         dtype={'Entrez ID': 'str'})

            df = df.sort_values(by='value', ascending=True)


        if len(df) != len(self.genes):
            raise Warning('{} lines in file. Should match nr of genes: {}.'
                          .format(len(df), len(self.genes)))

        if genes_subset is not None:
           df = df.loc[df['Entrez ID'].isin(genes_subset)]


        nr_rows_to_keep = int(round(len(df)*threshold))
        return list(df.iloc[0:nr_rows_to_keep]['Entrez ID'])

    def get_network(self, threshold, genes_subset=None, legacy=False):

        G = nx.Graph()
        if genes_subset is None:
            genes_intersection = self.genes
        else:
            genes_intersection = list(set(genes_subset).intersection(self.genes))
            print(len(self.genes))
            print(genes_subset)
            print(genes_intersection)
        for i, gene in enumerate(genes_intersection):
            print(i+1, '/', len(genes_intersection))
            for coex_gene in self.coexpressed(gene, threshold, genes_intersection, legacy):
                if gene != coex_gene:
                    G.add_edge(gene, coex_gene)
        return G


def main():
    """Get coexpression network at given threshold (networkx undirected graph).
    :returns: TODO

    """

    coexdb = CoexDB('dummy_data/coexdb/Sce')
    # coexdb = CoexDB('../../data_collection_and_preprocessing/collect_networks/raw_data/Hsa-r.v18-12.G22897-S22897.combat_pca_subagging.mrgeo.d/')
    subset = ['850532', '850594', '850651', '850708','850761', '850813']
    G = coexdb.get_network(0.21, genes_subset=subset)
    print(G.number_of_nodes())
    print(G.number_of_edges())
    nx.write_edgelist(G, 'test2.edgelist')
    pass

if __name__ == "__main__":
    main()
