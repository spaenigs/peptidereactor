__author__ = "Siba Mohsen"
__email__ = "mohsen@students.uni-marbug.de"

from Bio import PDB
import numpy as np
import pandas as pd
from scipy.spatial import Delaunay
from .utility import Utils


class ProteinStructure:
    """reads the protein structure file and extracts 5 functions from the Delaunay graph."""

    def __init__(self, protein, max_angstroms=15):
        self.protein = protein
        self.max_angstroms = max_angstroms

    def __read(self):
        """
        This function reads a .pdb files and extracts the residues and the coordinates of the C-alpha atoms.
        Then, the Delaunay triangulation is employed on the coordinates. A graph of tetrahedrons is made.
        From the graph, the edges of the tetrahedrons are extracted. The maximum length to be consider can be manually
        set by the user or taken as default (15 Angstroms).
        The distance between every two residues/points is calculated using distance_3d().
        :return: data frame containing two columns: column 'residues' contains the tuples of residues and column
        'distance' contains the distance between the two residues of a tuple.
        """
        parser = PDB.PDBParser()
        io = PDB.PDBIO()
        struct = parser.get_structure(self.protein, str(self.protein))
        residues = []
        coordinates = []
        coord_to_list = []
        for model in struct:
            for chain in model:
                for residue in chain:
                    try:
                        residues.append(residue.get_resname())
                        atom = residue['CA']
                        coordinates.append(atom.get_coord())
                    except:
                        pass

        # Extracting all the residues in the .pdb file will return a list of three letters strings. Some of them are not
        # included in the 20 amino acid group, so I just wrote this line to remove the uninteresting residues for
        # our study.
        residues = [resid for resid in residues if resid not in ["HOH", "PO4", "HEM"]]
        for point in coordinates:
            coord_to_list.append(point.tolist())

        aa_and_coordinates = list(zip(residues, coord_to_list))

        # In favor of applying the Delaunay triangulation on the coordinates, these have to be stored in an array not
        # a list.
        # That's why we transform our list of arrays to an array of 3D arrays.
        coo_delaunay = np.asarray(coordinates)

        # To apply the Delaunay function, simply use the Delaunay function from the scipy.spatial library.
        try:
            tri = Delaunay(coo_delaunay)
        except Exception as e:
            return pd.DataFrame()
            # import pydevd_pycharm
            # pydevd_pycharm.settrace('localhost', port=8889, stdoutToServer=True, stderrToServer=True)

        self.tri = tri
        # Now, we have a graph with vertices and edges and every 4 vertices/residues/(x,y,z)-coordinates form a
        # tetrahedron.

        # The first step to make toward the extraction of vertices and edges is to extract the indices of the
        # vertices of each four points/vertices bounded through edges to form the tetrahedron.
        indices = tri.simplices
        # The output is an array of arrays, each containing 4 indices.
        # To extract the vertices coordinates, find the coordinates by their indices.
        vertices = coo_delaunay[indices]
        # The output is an array, containing arrays of 4 arrays indicating the coordinates of each vertex/residue.
        edges = set()
        for (i0, i1, i2, i3) in indices:
            edges.add(Utils.sorted_tuple(i0, i1))
            edges.add(Utils.sorted_tuple(i0, i2))
            edges.add(Utils.sorted_tuple(i0, i3))
            edges.add(Utils.sorted_tuple(i1, i2))
            edges.add(Utils.sorted_tuple(i1, i3))
            edges.add(Utils.sorted_tuple(i2, i3))

        indices_to_coordinates = []
        for edge in edges:
            if Utils.distance_3d(coo_delaunay[edge[1]], coo_delaunay[edge[0]]) <= self.max_angstroms:
                indices_to_coordinates.append((coo_delaunay[edge[0]].tolist(), coo_delaunay[edge[1]].tolist()))

        len_edge = []
        for edge in edges:
            if Utils.distance_3d(coo_delaunay[edge[1]], coo_delaunay[edge[0]]) <= self.max_angstroms:
                len_edge.append(Utils.distance_3d(coo_delaunay[edge[1]], coo_delaunay[edge[0]]))

        coordinates_to_residue = []
        for t in indices_to_coordinates:
            for tuple1 in aa_and_coordinates:
                for tuple2 in aa_and_coordinates:
                    if t[0] == tuple1[1] and t[1] == tuple2[1]:
                        coordinates_to_residue.append((tuple1[0], tuple2[0]))

        # To create a data frame of each residue tuple/edge and its distance.
        # The data frame will make it easier for us to generate the 210 feature vector later.
        data = list(zip(coordinates_to_residue, len_edge))
        matrix = pd.DataFrame(data, columns=['residues', 'distance'])

        return matrix

    def average_distance(self):
        """
        Function 1.
        Given the data frame containing all the edges and distances, they must be first grouped by their residue's name.
        Then, the mean of each distance is calculated.
        A new data frame is generated, its values are then stored in a list of 210 elements,
        following the generate_aa_matrix() sorted list.
        :return: list of mean distances of each edge (210 features instead of 400).
        """
        data_frame = self.__read()
        if data_frame.empty:
            return []
        output_distance = []
        df = Utils.reverse_imp_tuples(data_frame)
        df = df.groupby('residues', as_index=False)['distance'].mean()
        for residue_edge in Utils.generate_aa_matrix():
            if residue_edge in set(df['residues']):
                output_distance.append(df.loc[df['residues'] == residue_edge, 'distance'].iloc[0])
            else:
                output_distance.append(0)
        return output_distance

    def total_distance(self):
        """
        Function 2.
        Given the data frame containing all the edges and distances, they must be first grouped by their residue's name.
        Then, the total distance's sum of each kind of edge is calculated.
        A new data frame is generated, its values are then stored in a list of 210 elements,
        following the generate_aa_matrix() sorted list.
        :return: list of distance's sum of each edge (210 features instead of 400).
        """
        data_frame = self.__read()
        if data_frame.empty:
            return []
        output_distance = []
        df = Utils.reverse_imp_tuples(data_frame)
        df = df.groupby('residues', as_index=False)['distance'].sum()
        for residue_edge in Utils.generate_aa_matrix():
            if residue_edge in set(df['residues']):
                output_distance.append(df.loc[df['residues'] == residue_edge, 'distance'].iloc[0])
            else:
                output_distance.append(0)
        return output_distance

    def number_instances(self):
        """
        Function 3.
        Given the data frame containing all the edges and distances, they must be first grouped by their residue's name.
        Then, the occurrence of each edge/residue's tuple is counted.
        A new data frame is generated, its values are then stored in a list of 210 elements,
        following the generate_aa_matrix() sorted list.
        :return: list of the occurrence of each edge (210 features instead of 400).
        """
        data_frame = self.__read()
        if data_frame.empty:
            return []
        output_distance = []
        df = Utils.reverse_imp_tuples(data_frame)
        df = df.groupby('residues', as_index=False)['distance'].count()
        for residue_edge in Utils.generate_aa_matrix():
            if residue_edge in set(df['residues']):
                output_distance.append(df.loc[df['residues'] == residue_edge, 'distance'].iloc[0])
            else:
                output_distance.append(0)
        return output_distance

    def frequency_instances(self):
        """
        Function 3.
        Given the data frame containing all the edges and distances, they must be first grouped by their residue's name.
        Then, the occurrence of each edge/residue's tuple is counted. After that, the occurrence is divided by the
        number of all residue edges in the data frame, so we get the frequency of each edge.
        A new data frame is generated, its values are then stored in a list of 210 elements,
        following the generate_aa_matrix() sorted list.
        :return: list of frequency of the occurrence of each edge (210 features instead of 400).
        """
        data_frame = self.__read()
        if data_frame.empty:
            return []
        output_distance = []
        df = Utils.reverse_imp_tuples(data_frame)
        all_len = len(df)
        df = df.groupby('residues', as_index=False)['distance'].count()
        df['distance'] = df.distance / all_len
        for residue_edge in Utils.generate_aa_matrix():
            if residue_edge in set(df['residues']):
                output_distance.append(df.loc[df['residues'] == residue_edge, 'distance'].iloc[0])
            else:
                output_distance.append(0)
        return output_distance

    def cartesian_product(self):
        """
        Function 5.
        Given the data frame containing all the edges and distances, the average of the distance of a certain edge
        as well as the number of instances of each kind of edge are calculated using the functions average_distance()
        and number_instances().
        The outputs of each of these functions are 210-dimensional lists.
        Then, the cartesian product of the two 210-dimensional lists is calculated by adding their contents element-wise
        to a new list like described in the publication 'Encoding protein structures with functions of graph' under the
        paragraph 'Machine Learning'.
        :return: list of the average distances followed by the number of instances (210 + 210 = 420 features).
        """
        average_distance = self.average_distance()
        number_of_instances = self.number_instances()
        cartesian_prod = [element for tup in zip(average_distance, number_of_instances) for element in tup]
        return cartesian_prod
