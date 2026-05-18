# UART Beam Control Protocol

UART format: `115200` baud, 8 data bits, no parity, 1 stop bit.

This protocol is the current N=32 host-to-FPGA frame used by
`upper/Phased_Array_Control.py` and `verilog/uart_protocol_parser.sv`.

## Frame Layout

```text
AA BB 01 amp[0] ... amp[31] phase0_l phase0_h ... phase31_l phase31_h checksum 0D 0A
```

Total length is 102 bytes:

```text
2 byte  header       AA BB
1 byte  command      01
32 byte amplitude    uint8, channel E00..E31
64 byte phase        uint16 little-endian, channel E00..E31
1 byte  checksum     8-bit sum from command through phase payload
2 byte  tail         0D 0A
```

## Channel Order

The host uses `analysis/n32_array_coordinates.csv` order directly:

```text
E00 -> transducer_io[0]
E01 -> transducer_io[1]
...
E31 -> transducer_io[31]
```

There is no 5x5 snake mapping in the N=32 path.

## Field Encoding

`amplitude[i]` is an unsigned 8-bit value:

```text
0   = off
255 = full channel weight
```

`phase[i]` is a little-endian unsigned 16-bit value, but the RTL parser keeps
only the low 12 bits:

```text
phase[i] = {phase_h[3:0], phase_l}
valid host range for the current 40 kHz PWM carrier: 0..2499 ticks
```

With the current 100 MHz PWM clock:

```text
2500 ticks = 25 us = one 40 kHz carrier period
```

## Checksum

The checksum is the low 8 bits of the byte sum starting at the command byte:

```text
checksum = sum(packet[2 : 2 + 1 + 32 + 64]) & 0xFF
```

It does not include `AA BB`, the checksum byte itself, or `0D 0A`.

## Host Parameter Limits

The N=32 host UI clamps or limits runtime parameters before packet generation:

```text
focus distance: 100..4999 mm
combined steering angle: <= 15 degrees
carrier frequency: 40 kHz
sound speed: 343 m/s
```

The FPGA parser updates the 32-channel shadow values only after the checksum
and `0D 0A` tail are valid.
