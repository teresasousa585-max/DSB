import unittest
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

from n32_beam_params import (
    MAX_DISTANCE_MM,
    MIN_DISTANCE_MM,
    MAX_STEER_DEG,
    build_packet,
    calculate_beam_params,
)


class N32BeamParamsTest(unittest.TestCase):
    def assert_valid_result(self, result):
        self.assertEqual(len(result.amplitudes), 32)
        self.assertEqual(len(result.phases), 32)
        self.assertTrue(all(0 <= amp <= 255 for amp in result.amplitudes))
        self.assertTrue(all(0 <= phase < 2500 for phase in result.phases))
        self.assertLessEqual(result.combined_angle_deg, MAX_STEER_DEG + 1e-6)

    def test_broadside_1m(self):
        result = calculate_beam_params(distance_mm=1000, az_deg=0, el_deg=0)
        self.assert_valid_result(result)
        self.assertEqual(result.phases[0], 0)

    def test_steering_and_distance_limits(self):
        for distance, az, el in [
            (MIN_DISTANCE_MM, 15, 0),
            (MIN_DISTANCE_MM, -15, 0),
            (MAX_DISTANCE_MM, 15, 15),
            (10, 90, 90),
        ]:
            result = calculate_beam_params(distance_mm=distance, az_deg=az, el_deg=el)
            self.assert_valid_result(result)
            self.assertGreaterEqual(result.distance_mm, MIN_DISTANCE_MM)
            self.assertLessEqual(result.distance_mm, MAX_DISTANCE_MM)

    def test_packet_format_and_checksum(self):
        result = calculate_beam_params(distance_mm=1000, az_deg=15, el_deg=0)
        packet = build_packet(result.amplitudes, result.phases)
        self.assertEqual(len(packet), 102)
        self.assertEqual(packet[0:3], b"\xAA\xBB\x01")
        self.assertEqual(packet[-2:], b"\x0D\x0A")
        self.assertEqual(packet[-3], sum(packet[2:-3]) & 0xFF)


if __name__ == "__main__":
    unittest.main()
