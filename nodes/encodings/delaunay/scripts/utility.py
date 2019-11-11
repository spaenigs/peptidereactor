__author__ = "Siba Mohsen"
__email__ = "mohsen@students.uni-marbug.de"

from Bio.Data import IUPACData


class Utils:

    """groups the variables and static methods used in the package AminoAcidComposition in order to be used freely"""

    # The tuple of 20 amino acids should be immutable
    aa_20 = tuple(IUPACData.protein_letters_1to3.keys())

    # The tuple of 20 twins amino acids
    aa_twins = ('AA', 'CC', 'DD', 'EE', 'FF', 'GG', 'HH', 'II', 'KK', 'LL',
                'MM', 'NN', 'PP', 'QQ', 'RR', 'SS', 'TT', 'VV', 'WW', 'YY')

    # The set of amino acid residues encoded in three letters code.
    # Since the .pdb dataset includes residues encoded in three letters code, we will be working
    # with them in this paper.
    aa_residues = tuple(a.upper() for a in IUPACData.protein_letters_1to3.values())

    @staticmethod
    def frequency(seq, aa_group):
        """
        frequency(): calculates the frequency of a specific amino acid group in a peptide sequence (e.g. hydrophobic
        amino acid group) and returns the distance between the closest amino acids in order to classify them later
        in 6 distance classes.
        Example: Matsuda et al.(2005), § 'Distance frequency' will return after calling this function: [3,2,3,2,3].
        INPUT:
        @:param seq: String. The sequence of amino acid.
        @:param aa_group: String. The group for which the frequency distance should be calculated, e.g. basic,
        hydrophobic.
        OUTPUT:
        List<int>. List of the distances between two successive specific amino acids.
        """
        freq_index = [pos for pos, char in enumerate(seq) if char in aa_group]
        distance = [(freq_index[d + 1] - freq_index[d]) for d in range(0, len(freq_index) - 1)]
        return distance

    @staticmethod
    def reverse(tuples):
        """
        Given a tuple, this function reverses its arguments.
        e.g. Input = tuple('a','b','c') -----reverse----> Output = tuple('c','b','a').
        :param tuples: a tuple of any size
        :return: the input tuple reversed
        NOTE: Applying -1 on the usual slice notation of lists/strings/tuples in python list[<start>:<stop>:<step>]
        allows us to have the reverse of this lists/string/tuple. In this case, step = -1, so the iterable will start
        slicing from the last element till the first one.
        """
        new_tup = tuples[::-1]
        return new_tup

    @staticmethod
    def sorted_tuple(vertex_index_1, vertex_index_2):
        """
        This function allows us to extract the same edge of the tetrahedron only once, even if it comes in reverse
        order.
        This will limit the number of edges to exactly six edges per tetrahedron rather than counting the same edge
        twice: e.g. (ALA, SER) and (SER, ALA) in a tetrahedron containing the vertices [ALA SER CYS TYR].
        :param vertex_index_1: index of the vertex 1.
        :param vertex_index_2: index of the vertex 2.
        :return: one sorted tuple for each edge in a tetrahedron.
        """
        return (vertex_index_1, vertex_index_2) if vertex_index_1 < vertex_index_2 else \
            (vertex_index_2, vertex_index_1)

    @staticmethod
    def distance_3d(point1, point2):
        """
        To calculate the distance between two 3D points using the Pythagorean theorem.
        :param point1: coordinates of the first point as [x1,y1,z1] or (x1,y1,z1).
        :param point2: coordinates of the second point as [x2,y2,z2] or (x2,y2,z2).
        :return: distance between two points in 3D graph.
        """
        return (((point2[0] - point1[0]) ** 2) + ((point2[1] - point1[1]) ** 2) + (
                    (point2[2] - point1[2]) ** 2)) ** (1 / 2)

    @staticmethod
    def generate_aa_matrix():
        """
        Extracting the edges after applying the Delaunay triangulation on the coordinates of C-alpha (CA) atoms of a
        .pdb file consists of calculating the vertices(residues) of the Delaunay graph.
        Each of these vertices/residues represents a point with coordinates (x,y,z) in the graph and is displayed as a
        three letters code of amino acids.
        The edge (A_i, A_j) is the connection between two residues, A_i(x1,y1,z1) and A_j(x2,y2,z2).
        Since we have 20 amino acids, we will have 20x20 kinds of connection between edges.
        This is exactly a 20x20 matrix, or a so called adjacency matrix.
        In this paper, the authors focused on the upper part of the matrix, as well as on the diagonal (A_i, A_i).
        The part of the matrix below the diagonal is ignored, because it is just the reverse display of the upper
        tuples. Therefore their distances remain the same as for the upper one.
        At the end, instead of having 20x20 = 400 elements in the matrix, we will only have 210 features.
        210 = 400 - (380/2) + 20.
        :return: a list of the 210 ordered elements of the matrix.
        The output looks like this:
        [('ALA', 'ALA'), ('ALA', 'CYS'), ('ALA', 'ASP'), ('ALA', 'GLU'), ('ALA', 'PHE'), ... , ('TYR', 'TYR')]
        """
        matrix_elements = []
        for i in Utils.aa_residues:
            for j in Utils.aa_residues:
                if Utils.reverse((i, j)) not in matrix_elements:
                    matrix_elements.append((i, j))
        return matrix_elements

    @staticmethod
    def reverse_imp_tuples(data_frame):
        """
        This function reverses the important tuples in order to be considered with the 210 residue'S tuples instead
        of the 400 tuples.
        :param data_frame: data frame with a column 'residues'.
        :return: returns the same data frame with reversed tuples.
        """
        for index, row in data_frame.iterrows():
            tuple_value = row['residues']
            if tuple_value not in Utils.generate_aa_matrix():
                data_frame.at[index, 'residues'] = tuple_value[::-1]  # or reverse(tuple_value)
        return data_frame

