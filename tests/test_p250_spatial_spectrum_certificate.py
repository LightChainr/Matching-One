from __future__ import annotations

import copy
import json
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))

from p250_spatial_spectrum_certificate import (
    GaussianQuotient, NN, MATCHING, alias_audit, build_certificate,
    charge_coefficients, parent_label, phase_multiply, projection_control, totient,
)
from verify_p250_spatial_spectrum import verify_certificate


class SpatialSpectrumTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data = build_certificate()

    def test_checked_certificate_reproduces_exactly(self):
        path = ROOT / 'results' / 'p250-spatial-spectrum' / 'certificate.json'
        self.assertEqual(self.data, json.loads(path.read_text(encoding='utf-8')))

    def test_independent_union_find_verifies_bfs_certificate(self):
        result = verify_certificate(self.data)
        self.assertTrue(result['verified'])
        self.assertEqual(result['exact_spatial_rank_lower_bound'], 100)
        self.assertFalse(result['production_data_read'])

    def test_empty_and_full_topology_controls(self):
        geometry = GaussianQuotient(3, 2)
        self.assertEqual(geometry.component_ranks([], NN), ([], {}))
        for steps in (NN, MATCHING):
            components, ranks = geometry.component_ranks(geometry.representatives, steps)
            self.assertEqual(components, [{'size': 13, 'ambient_rank': 2}])
            self.assertEqual(set(ranks.values()), {2})

    def test_contractible_site_control(self):
        geometry = GaussianQuotient(3, 2)
        components, ranks = geometry.component_ranks([geometry.key(0, 0)], NN)
        self.assertEqual(components, [{'size': 1, 'ambient_rank': 0}])

    def test_actual_witnesses_have_two_rank_one_carriers(self):
        for row, occupied in zip(self.data['witnesses'], (31, 29)):
            self.assertEqual(row['occupied_count'], occupied)
            self.assertEqual(row['black_NN_components'], [{'size': occupied, 'ambient_rank': 1}])
            self.assertEqual(row['white_matching_components'], [{'size': 505-occupied, 'ambient_rank': 1}])

    def test_each_charge_has_gauge_robust_zero_and_nonzero_sites(self):
        for row in self.data['witnesses']:
            for charge in row['charges'].values():
                self.assertEqual(charge['zero_parents'], 72)
                self.assertEqual(charge['nonzero_parents'], 29)
                self.assertEqual(charge['zero_coefficients'], [0, 0, 0, 0])
                self.assertEqual(charge['nonzero_coefficients'], [2, 0, 0, 0])

    def test_constant_fibers_have_no_nontrivial_deck_charge(self):
        for r in (1, 2, 3, 4):
            self.assertEqual(charge_coefficients([1]*5, r), [0]*4)
            self.assertEqual(charge_coefficients([-1]*5, r), [0]*4)

    def test_deck_charge_conjugacy_exact(self):
        def conjugate(value):
            out = [0]*4
            for power, coefficient in enumerate(value):
                term = phase_multiply([coefficient, 0, 0, 0], -power)
                out = [x+y for x, y in zip(out, term)]
            return out
        for values in ([1, -1, 0, 1, 0], [1, -1, -1, -1, -1]):
            for r in (1, 2):
                self.assertEqual(conjugate(charge_coefficients(values, r)),
                                 charge_coefficients(values, 5-r))

    def test_gauge_action_is_invertible_and_order_five(self):
        for values in ([0, 0, 0, 0], [2, 0, 0, 0], [1, -2, 3, -4]):
            for t in range(5):
                self.assertEqual(phase_multiply(phase_multiply(values, t), -t), list(values))
            self.assertEqual(phase_multiply(values, 5), list(values))

    def test_parent_periods_and_rotation(self):
        self.assertEqual(parent_label(10, 1), 0)
        self.assertEqual(parent_label(-1, 10), 0)
        self.assertEqual(len({parent_label(j, 0) for j in range(101)}), 101)
        for x in range(-6, 7):
            for y in range(-6, 7):
                self.assertEqual(parent_label(-y, x), (-10*parent_label(x, y)) % 101)

    def test_radius_six_has_eight_repeated_parent_classes(self):
        self.assertEqual([(alias_audit(r)['displacement_labels'],
                           alias_audit(r)['distinct_parent_vertices']) for r in (4, 5, 6)],
                         [(41, 41), (61, 61), (85, 77)])
        self.assertEqual(parent_label(5, 0), parent_label(-5, -1))
        self.assertEqual(len(alias_audit(6)['repeated_classes']), 8)

    def test_cyclotomic_degree_arithmetic(self):
        self.assertEqual((totient(5), totient(101), totient(505)), (4, 100, 400))
        self.assertEqual(totient(505)//totient(5), 100)

    def test_projection_leakage_without_microscopic_noncommutation(self):
        result = projection_control()
        self.assertTrue(all(result['checks'].values()))
        self.assertEqual(result['commutator'][0][1], '-4/9')

    def test_corrupt_occupied_configuration_rejected(self):
        data = copy.deepcopy(self.data)
        data['witnesses'][0]['occupied_lifts'][0] = [99, 99]
        with self.assertRaisesRegex(ValueError, 'staircase'):
            verify_certificate(data)

    def test_corrupt_topology_rejected(self):
        data = copy.deepcopy(self.data)
        data['witnesses'][1]['white_matching_components'][0]['ambient_rank'] = 2
        with self.assertRaisesRegex(ValueError, 'topology'):
            verify_certificate(data)

    def test_corrupt_charge_rejected(self):
        data = copy.deepcopy(self.data)
        data['witnesses'][0]['charges']['1']['nonzero_coefficients'] = [0]*4
        with self.assertRaisesRegex(ValueError, 'nonzero witness'):
            verify_certificate(data)

    def test_corrupt_rank_or_ensemble_rejected(self):
        for key, value in [('complete_spatial_rank_lower_bound', 101), ('ensemble', 'all p')]:
            data = copy.deepcopy(self.data)
            data['conclusion'][key] = value
            with self.assertRaises(ValueError):
                verify_certificate(data)

    def test_corrupt_alias_rejected(self):
        data = copy.deepcopy(self.data)
        data['endpoint_aliases'][2]['distinct_parent_vertices'] = 85
        with self.assertRaisesRegex(ValueError, 'alias'):
            verify_certificate(data)

    def test_invalid_inputs_fail_explicitly(self):
        for call in (lambda: GaussianQuotient(0, 0), lambda: alias_audit(-1),
                     lambda: charge_coefficients([1]*4, 1),
                     lambda: charge_coefficients([1]*5, 0), lambda: totient(0)):
            with self.assertRaises(ValueError):
                call()


if __name__ == '__main__':
    unittest.main()
