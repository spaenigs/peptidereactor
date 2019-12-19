__author__ = "Siba Mohsen"
__email__ = "mohsen@students.uni-marbug.de"

from .utility import Utils


class DistanceFrequency:
    """Calculates the distance frequency-feature vector of a peptide sequence."""

    def __init__(self, seq, dn=13, dc=10):
        self.seq = seq
        self.dn = dn
        self.dc = dc

    # Define the three groups of amino acids to work with:
    # basic, hydrophobic and others amino acids
    # I chose to write them as variables, not dictionaries since they are fixed by the authors of Matsuda et. al.
    basic_amino_acid = 'RKH'
    hydrophobic_amino_acids = 'IVLFMAGWP'
    other_amino_acids = 'DNEQYST'

    def __fragmentation(self):

        """
        fragmentation(): breaks an antimicrobial peptide into three fragments, namely the N-, M- and C- Terminal.
        Note: The values dn and dc are set per default to 13 and 10 respectively. These values were deduced from:
        Abdullah M. Khamis, Magbubah Essack, Xin Gao and Vladimir B. Bajic, Distinct profiling of antimicrobial
        peptide families, Bioinformatics, 31(6), 2015, 849–856.
        From table 1 (Khamis et al.), column 1 shows different types of antimicrobial peptides,column 2 shows the dn
        values (N-Terminal = 4*dn) of each type and column 3 shows the dc values: the numbers of amino acids located
        at the end of the peptide.
        In column 2 dn can take 4 different values: 10,12,14,16. If we calculate the average of these values,
        dn = 13 will result.
        In column 3 dc can only take 2 values: 8 and 10. Here, dc = 10 was selected, because 10 is most common for all
        the types. Since the length of antimicrobial peptides differs from 5 to up to 100 amino acids(Khamis et al.),
        it is difficult to determine one uniform value for dn and dc.
        In this method, dn and dc can be specified by the user. If they are not, the method defaults to dn = 13
        and dc = 10.
        INPUT:
        @:param dn: int. The N-terminal of the sequence of amino acids to be fragmented
        @:param dc: int. The C-terminal of the sequence of amino acids to be fragmented
        OUTPUT:
        dictionary: keys = [N,M,C] and values = fragments of the sequence seq as strings.
        """

        length_sequence = len(self.seq)
        if length_sequence >= (4*self.dn + 20 + self.dc):
            return {'N': self.seq[:(4*self.dn)],
                    'M': self.seq[(4*self.dn): -self.dc],
                    'C': self.seq[-self.dc:length_sequence]}
        elif length_sequence <= (4*self.dn + self.dc):
            return {'N': self.seq[:((length_sequence - self.dc)//2)],
                    'M': self.seq[((length_sequence - self.dc)//2): (length_sequence - self.dc)],
                    'C': self.seq[-self.dc:length_sequence]}
        else:
            return {'N': self.seq[:(4*self.dn)],
                    'M': self.seq[(4*self.dn): -self.dc],
                    'C': self.seq[-self.dc:length_sequence]}

    def __fragmentation_n_terminal(self):

        """
        fragmentation_n_terminal(): breaks the N-Terminal of an antimicrobial peptide into four fragments: n1, n2, n3
        and n4.
        Note: As mentioned in Matsuda et al.(2005) the N-Terminal must be fragmented in 4 portions dn.
        In fact, this fragmentation can only take place if the length of the N-Terminal = 4*dn.
        In this case, the function will return the 4 parts of the N-terminal, packed in a dictionary and each
        fragment represents a different sub sequence of the N-terminal. Otherwise, the function will return the
        same sub sequences for n1, n2, n3 and n4, as assigned in Matsuda et al.(2005) under the paragraph
        "Feature Vector".
        In this method, dn can be specified. If it is not, the method defaults to dn = 13.
        INPUT:
        @:param dn: int. The N-terminal of the sequence of amino acids to be fragmented
        OUTPUT:
        dictionary: keys = [n1,n2,n3,n4] and values = fragments of the N-Terminal as strings.
        """
        fragment_n = self.__fragmentation()['N']
        if len(fragment_n) == 4*self.dn:
            return {'n1': fragment_n[:self.dn],
                    'n2': fragment_n[self.dn:2*self.dn],
                    'n3': fragment_n[2*self.dn:3*self.dn],
                    'n4': fragment_n[3*self.dn:4*self.dn]}
        else:
            return {'n1': fragment_n[:len(fragment_n)],
                    'n2': fragment_n[:len(fragment_n)],
                    'n3': fragment_n[:len(fragment_n)],
                    'n4': fragment_n[:len(fragment_n)]}

    @staticmethod
    def __seq_to_list(seq):
        """
        seq_to_list(): makes from each amino acid sequence a list of 20 elements which refer to the 20 amino acids
        introduced above as aa_20. It counts the occurrence of an amino acid then place it in the relevant list's
        position.
        Note: Since the fragmentation functions above return fragments of the sequence as strings and these strings are
        supposed to represent a list of 20 amino acid to be included in the feature vector, this function will map the
        occurrence of each amino acid in a list of 20 elements.
        INPUT:
        @:param seq: String. The sequence of amino acid to be transformed to a tuple.
        OUTPUT:
        result: List<int>. Ordered list counting the occurrence of each amino acid of the input sequence.
        """
        return [seq.count(aa) if aa in seq else 0 for aa in Utils.aa_20]

    @staticmethod
    def __detect_twin_aa(seq):
        """
        detect_twin_aa(): detects the occurrence of two same, consecutive amino acids and store their occurrence in a
        vector of 20 amino acids.
        Note: As mentioned in the paragraph "Feature Vector" in Matsuda et al.(2005), the vector should also include
        'the composition of 20 twin amino acids (e.g., RR, KK) in the middle part.'
        INPUT:
        @:param seq: String. The sequence of amino acids.
        OUTPUT:
        result: list<int>. Ordered list counting the occurrence of each twin successive amino acids from the input
        sequence.
        """
        return [seq.count(aa) if aa in seq else 0 for aa in Utils.aa_twins]

    @staticmethod
    def __distance_frequency(seq, aa_group):
        """
        distance_frequency(): counts the occurrence of a distance frequency then classify it in one of the 6 frequency
        classes introduced in Matsuda et al.(2005).
        The distance classes are represented in the following tuple:(H=1, 1<H<=6, 6<H<=11, 11<H<=16, 16<H<=21, H>21),
        where H is the distance.
        Example: For the same example of the paper, this function will first calculate the frequency() [3,2,3,2,3] then
        return their classification in a list of 6 elements representing the distance classes above [0,5,0,0,0,0].
        INPUT:
        @:param seq: String. The sequence of amino acid.
        @:param aa_group: String. The group for which the frequency distance should be calculated, e.g. basic,
        hydrophobic.
        OUTPUT:
        List<int>. List of the occurrence of each distance frequency class in the sequence.
        """
        h1, h2, h3, h4, h5, h6 = 0, 0, 0, 0, 0, 0
        freq = Utils.frequency(seq, aa_group)
        for x in freq:
            if x == 1:
                h1 += 1
            elif x in [2, 3, 4, 5, 6]:
                h2 += 1
            elif x in [7, 8, 9, 10, 11]:
                h3 += 1
            elif x in [12, 13, 14, 15, 16]:
                h4 += 1
            elif x in [17, 18, 19, 20, 21]:
                h5 += 1
            elif x > 21:
                h6 += 1
        result = [h1, h2, h3, h4, h5, h6]
        return result

    def feature_vector(self):
        """
        feature_vector(): delivers the feature vector to be used in Machine Learning models later on. Given a peptide
        sequence, it cuts the amino acid sequence in several fragments using the functions defined above then stores all
        the returned lists in a list. Every element of this list takes its position exactly like described in
        Matsuda et al.(2005). At the end, the list of lists will be flatten and the function will return a list of 184
        integers, ordered as a tuple called 'fragments' within the function.
        OUTPUT:
        List<List<int>>. Contains 184 components grouped in lists.
        """

        n1 = self.__fragmentation_n_terminal()['n1']
        n2 = self.__fragmentation_n_terminal()['n2']
        n3 = self.__fragmentation_n_terminal()['n3']
        n4 = self.__fragmentation_n_terminal()['n4']
        n_terminal = self.__fragmentation()['N']
        m_terminal = self.__fragmentation()['M']
        twins_m_terminal = self.__detect_twin_aa(m_terminal)
        c_terminal = self.__fragmentation()['C']
        dist_freq_n = self.__distance_frequency(n_terminal, self.basic_amino_acid)
        dist_freq_m = self.__distance_frequency(m_terminal, self.basic_amino_acid)
        dist_freq_m_hydro = self.__distance_frequency(n_terminal, self.hydrophobic_amino_acids)
        dist_freq_m_others = self.__distance_frequency(n_terminal, self.other_amino_acids)
        entire_seq = self.seq
        fragments = (n1, n2, n3, n4,    # Distribution of amino acids in the four parts of the N-Terminal
                     m_terminal,        # Distribution of amino acids in the Middle-Terminal
                                        # twins_m_terminal, # Distribution of twin amino acids in the Middle-Terminal
                     c_terminal,        # Distribution of amino acids in the C-Terminal
                                        # Distribution of amino acids frequency in 6 distance classes:
                                        # dist_freq_n, dist_freq_m, dist_freq_m_hydro, dist_freq_m_others,
                     entire_seq)        # Distribution of amino acids along the entire peptide length

        vector = []
        for fragment in fragments:
            vector.append(self.__seq_to_list(fragment))
        vector.insert(5, twins_m_terminal)
        vector.insert(7, dist_freq_n)
        vector.insert(8, dist_freq_m)
        vector.insert(9, dist_freq_m_hydro)
        vector.insert(10, dist_freq_m_others)
        # to transform the list of lists to one list containing 184 successive components
        vector = [item for sublist in vector for item in sublist]
        return vector
