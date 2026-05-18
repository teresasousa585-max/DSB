import csv
import math
import os
import struct
from collections import namedtuple


CHANNEL_COUNT = 32
SOUND_SPEED_M_S = 343.0
CARRIER_HZ = 40000.0
FPGA_CLK_HZ = 100_000_000
PWM_PERIOD_TICKS = 2500
MIN_DISTANCE_MM = 100
MAX_DISTANCE_MM = 4999
MAX_STEER_DEG = 15.0

MODE_FITTED = "Fitted RLS formula"
MODE_FULL = "Full amplitude"
MODE_HAMMING = "Radial Hamming"
AMPLITUDE_MODES = (MODE_FITTED, MODE_FULL, MODE_HAMMING)

BeamResult = namedtuple(
    "BeamResult",
    [
        "amplitudes",
        "phases",
        "raw_amplitudes",
        "distance_mm",
        "actual_az_deg",
        "actual_el_deg",
        "combined_angle_deg",
        "target_m",
    ],
)


DEFAULT_COORDS_MM = (
    (0.00, 0.00),
    (19.50, 0.00),
    (12.16, 15.25),
    (-4.34, 19.01),
    (-17.57, 8.46),
    (-17.57, -8.46),
    (-4.34, -19.01),
    (12.16, -15.25),
    (37.00, 0.00),
    (28.34, 23.78),
    (6.42, 36.44),
    (-18.50, 32.04),
    (-34.77, 12.65),
    (-34.77, -12.65),
    (-18.50, -32.04),
    (6.42, -36.44),
    (28.34, -23.78),
    (58.00, 0.00),
    (52.99, 23.59),
    (38.81, 43.10),
    (17.92, 55.16),
    (-6.06, 57.68),
    (-29.00, 50.23),
    (-46.92, 34.09),
    (-56.73, 12.06),
    (-56.73, -12.06),
    (-46.92, -34.09),
    (-29.00, -50.23),
    (-6.06, -57.68),
    (17.92, -55.16),
    (38.81, -43.10),
    (52.99, -23.59),
)


# Refit from the N32 simulation model over 100..4999 mm and 0..15 degrees.
# Each row is [a0, a1*r_m, a2*theta_n, a3*r_m^2, a4*theta_n^2, a5*r_m*theta_n].
AMPLITUDE_COEFF_MODELS = (
    (0.961311464201, -0.37294021301, -0.0472811300656, 0.0329300289827, 0.0515598103784, 0.168986631635),
    (-1.02536467443, -0.543031154246, 6.43112151055, 0.0539052143166, -5.71611393016, 0.874272452205),
    (1.37582109018, 3.73172277494, -16.0706989262, -0.386406277334, 14.1904781851, -2.8525104199),
    (-0.74444565781, -2.58864698391, 9.36033912493, 0.269574727217, -8.31808457918, 1.72526861009),
)


def clamp(value, low, high):
    return max(low, min(high, value))


def default_layout_csv():
    here = os.path.dirname(os.path.abspath(__file__))
    return os.path.abspath(os.path.join(here, "..", "analysis", "n32_array_coordinates.csv"))


def load_layout(path=None):
    csv_path = path or default_layout_csv()
    if not os.path.exists(csv_path):
        return list(DEFAULT_COORDS_MM)

    coords = []
    with open(csv_path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            coords.append((float(row["X_mm"]), float(row["Y_mm"])))

    if len(coords) != CHANNEL_COUNT:
        raise ValueError("N32 layout must contain exactly 32 elements")
    return coords


def limit_direction(az_deg, el_deg, max_angle_deg=MAX_STEER_DEG):
    sx = math.sin(math.radians(float(az_deg)))
    sy = math.sin(math.radians(float(el_deg)))
    max_s = math.sin(math.radians(max_angle_deg))
    sxy = math.hypot(sx, sy)

    if sxy > max_s and sxy > 0.0:
        scale = max_s / sxy
        sx *= scale
        sy *= scale
        sxy = max_s

    sz = math.sqrt(max(0.0, 1.0 - sx * sx - sy * sy))
    actual_az = math.degrees(math.asin(clamp(sx, -1.0, 1.0)))
    actual_el = math.degrees(math.asin(clamp(sy, -1.0, 1.0)))
    combined_angle = math.degrees(math.asin(clamp(sxy, 0.0, 1.0)))
    return sx, sy, sz, actual_az, actual_el, combined_angle


def fitted_amplitude(distance_mm, combined_angle_deg, rho_mm):
    r_m = clamp(float(distance_mm), MIN_DISTANCE_MM, MAX_DISTANCE_MM) / 1000.0
    theta_n = clamp(abs(float(combined_angle_deg)) / MAX_STEER_DEG, 0.0, 1.0)
    rho_n = clamp(float(rho_mm) / 58.004, 0.0, 1.0)
    features = (1.0, r_m, theta_n, r_m * r_m, theta_n * theta_n, r_m * theta_n)

    coeff = []
    for model in AMPLITUDE_COEFF_MODELS:
        coeff.append(sum(a * x for a, x in zip(model, features)))

    amp = coeff[0] + coeff[1] * rho_n + coeff[2] * rho_n * rho_n + coeff[3] * rho_n * rho_n * rho_n
    return clamp(amp, 0.0, 1.0)


def radial_hamming_amplitude(rho_mm):
    rho_n = clamp(float(rho_mm) / 58.004, 0.0, 1.0)
    return clamp(0.54 + 0.46 * math.cos(math.pi * rho_n), 0.0, 1.0)


def calculate_beam_params(distance_mm=1000, az_deg=0.0, el_deg=0.0, amp_mode=MODE_FITTED, coords_mm=None):
    coords = list(coords_mm or load_layout())
    if len(coords) != CHANNEL_COUNT:
        raise ValueError("N32 layout must contain exactly 32 elements")

    distance_mm = int(round(clamp(distance_mm, MIN_DISTANCE_MM, MAX_DISTANCE_MM)))
    distance_m = distance_mm / 1000.0
    sx, sy, sz, actual_az, actual_el, combined_angle = limit_direction(az_deg, el_deg)
    target = (distance_m * sx, distance_m * sy, distance_m * sz)

    phases = []
    raw_amplitudes = []
    amplitudes = []

    for x_mm, y_mm in coords:
        x_m = x_mm / 1000.0
        y_m = y_mm / 1000.0
        path_m = math.sqrt((target[0] - x_m) ** 2 + (target[1] - y_m) ** 2 + target[2] ** 2)
        phase_cycles = (path_m - distance_m) * CARRIER_HZ / SOUND_SPEED_M_S
        phase_ticks = int(round((phase_cycles % 1.0) * PWM_PERIOD_TICKS)) % PWM_PERIOD_TICKS
        phases.append(phase_ticks)

        rho_mm = math.hypot(x_mm, y_mm)
        if amp_mode == MODE_FULL:
            amp = 1.0
        elif amp_mode == MODE_HAMMING:
            amp = radial_hamming_amplitude(rho_mm)
        else:
            amp = fitted_amplitude(distance_mm, combined_angle, rho_mm)
        raw_amplitudes.append(amp)
        amplitudes.append(int(round(clamp(amp, 0.0, 1.0) * 255.0)))

    return BeamResult(
        amplitudes=amplitudes,
        phases=phases,
        raw_amplitudes=raw_amplitudes,
        distance_mm=distance_mm,
        actual_az_deg=actual_az,
        actual_el_deg=actual_el,
        combined_angle_deg=combined_angle,
        target_m=target,
    )


def build_packet(amplitudes, phases):
    if len(amplitudes) != CHANNEL_COUNT or len(phases) != CHANNEL_COUNT:
        raise ValueError("Packet requires 32 amplitudes and 32 phases")

    packet = bytearray([0xAA, 0xBB, 0x01])
    for amp in amplitudes:
        packet.append(int(clamp(amp, 0, 255)) & 0xFF)

    for phase in phases:
        ph = int(phase)
        if ph < 0 or ph > 0x0FFF:
            raise ValueError("Phase must fit in the parser's 12-bit field")
        packet.extend(struct.pack("<H", ph))

    packet.append(sum(packet[2:]) & 0xFF)
    packet.extend([0x0D, 0x0A])
    if len(packet) != 102:
        raise AssertionError("N32 UART packet length must be 102 bytes")
    return bytes(packet)


def format_packet_hex(packet):
    return " ".join("{:02X}".format(b) for b in packet)
