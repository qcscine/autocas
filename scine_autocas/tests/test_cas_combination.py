# pylint: disable=C0114, C0115, C0116
__copyright__ = """ This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details.
"""

import unittest

from scine_autocas.autocas_utils.cas_combination import combine_active_spaces


class TestCombineActiveSpaces(unittest.TestCase):
    def test_cas_combination(self):
        """
        Idea of the test: Provide artificial orbital map, cas indices, and occupations.
        Check the combined CAS.
        """
        occ_1 = [2, 2, 0, 0]
        occ_2 = [2, 0]
        indices_1 = [5, 6, 7, 8]
        indices_2 = [6, 7]
        orbital_groups = [
            [[2, 3, 4], [2, 4, 5]],  # This group is never selected in any CAS
            [[7, 8, 9], [7, 8, 10]],  # Selected in both CAS
            [[0], [0]],  # never selected
            [[1], [1]],  # never selected
            [[10], [9]],  # never selected
            [[5], [6]],  # selected in both CAS
            [[6], [3]]  # selected only in the first CAS
        ]
        combined_occupations, combined_indices = combine_active_spaces([occ_1, occ_2], [indices_1, indices_2],
                                                                       orbital_groups)
        reference_occupations = [2, 2, 0, 0, 0]
        reference_indices_1 = [5, 6, 7, 8, 9]
        reference_indices_2 = [6, 3, 7, 8, 10]
        self.assertEqual(reference_occupations, combined_occupations[0])
        self.assertEqual(reference_occupations, combined_occupations[1])
        self.assertEqual(reference_indices_1, combined_indices[0])
        self.assertEqual(reference_indices_2, combined_indices[1])


if __name__ == "__main__":
    unittest.main()
