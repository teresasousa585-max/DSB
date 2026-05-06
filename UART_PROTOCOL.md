# UART Beam Control Protocol

UART format: `921600` baud, 8 data bits, no parity, 1 stop bit.

Frame layout:

```text
AA 55 cmd seq len_l len_h payload... crc_l crc_h
```

CRC is `CRC-16/CCITT-FALSE` over `cmd, seq, len_l, len_h, payload...`:

```text
poly = 0x1021
init = 0xFFFF
bit order = MSB first
crc bytes on wire = low byte, then high byte
```

Commands:

```text
0x10 WRITE_PARAMS
  len = 75
  payload = 25 channels * {phase_l, phase_h, amplitude}
  phase uses low 10 bits: phase = {phase_h[1:0], phase_l}
  amplitude is 0..255
  The frame is first written into a shadow buffer, then committed on the next
  40 kHz carrier boundary.

0x11 START
  len = 0
  Enables ultrasound output.

0x12 STOP
  len = 0
  Disables ultrasound output immediately and forces PWM low.

0x13 SOFT_RESET
  len = 0
  Disables ultrasound output, resets the carrier/envelope/PWM chain, then
  reloads the active phase/amplitude table on a carrier boundary.
```

`cmd_status[7:4]` is the last latched parser error:

```text
0 = no error / last accepted frame OK
1 = unsupported command
2 = invalid length
3 = CRC mismatch
4 = UART frame error
```
