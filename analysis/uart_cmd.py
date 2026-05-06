#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
FPGA超声阵列DSB系统的UART命令帧生成器

UART协议参数:
- 波特率: 921600, 8N1
- 帧格式: AA 55 cmd seq len_l len_h payload... crc_l crc_h
- CRC: CRC-16/CCITT-FALSE
  - poly = 0x1021
  - init = 0xFFFF
  - MSB first
  - 输出字节顺序: low byte first, then high byte
"""

import struct


# =============================================================================
# CRC-16/CCITT-FALSE 实现
# =============================================================================
def crc16_ccitt_false(data: bytes) -> int:
    """
    计算 CRC-16/CCITT-FALSE.

    参数:
        data: 待计算CRC的字节序列

    返回:
        16位CRC值 (0x0000 ~ 0xFFFF)
    """
    crc = 0xFFFF
    for byte in data:
        crc ^= byte << 8
        for _ in range(8):
            if crc & 0x8000:
                crc = (crc << 1) ^ 0x1021
            else:
                crc <<= 1
        crc &= 0xFFFF
    return crc


# =============================================================================
# 帧构建函数
# =============================================================================
def _build_frame(cmd: int, seq: int, payload: bytes = b"") -> bytes:
    """
    构建UART命令帧的通用函数.

    帧格式: AA 55 cmd seq len_l len_h payload... crc_l crc_h

    参数:
        cmd: 命令字节 (0x10, 0x11, 0x12, 0x13)
        seq: 序列号 (0~255)
        payload: 载荷数据

    返回:
        完整的帧字节序列
    """
    length = len(payload)
    len_l = length & 0xFF
    len_h = (length >> 8) & 0xFF

    # 头部: 同步字 + 命令 + 序列号 + 长度
    header = bytes([0xAA, 0x55, cmd, seq, len_l, len_h])

    # 计算CRC: 对 header[2:] + payload 进行计算
    # 即 cmd + seq + len_l + len_h + payload
    crc_data = header[2:] + payload
    crc = crc16_ccitt_false(crc_data)

    crc_l = crc & 0xFF
    crc_h = (crc >> 8) & 0xFF

    frame = header + payload + bytes([crc_l, crc_h])
    return frame


def build_write_params_frame(seq: int, phases: list, amplitudes: list) -> bytes:
    """
    构建 WRITE_PARAMS (0x10) 命令帧.

    参数:
        seq: 序列号
        phases: 25路相位值列表 (0~1023, 低10位有效)
        amplitudes: 25路幅度值列表 (0~255)

    返回:
        完整的帧字节序列 (长度 = 6 + 75 + 2 = 83 字节)

    异常:
        ValueError: 通道数不为25或数值超出范围
    """
    if len(phases) != 25 or len(amplitudes) != 25:
        raise ValueError("phases 和 amplitudes 必须各包含25个元素")

    payload = bytearray()
    for i in range(25):
        phase = int(phases[i]) & 0x3FF  # 低10位有效
        amplitude = int(amplitudes[i]) & 0xFF

        phase_l = phase & 0xFF
        phase_h = (phase >> 8) & 0x03

        payload.append(phase_l)
        payload.append(phase_h)
        payload.append(amplitude)

    return _build_frame(0x10, seq, bytes(payload))


def build_start_frame(seq: int = 0) -> bytes:
    """构建 START (0x11) 命令帧."""
    return _build_frame(0x11, seq)


def build_stop_frame(seq: int = 0) -> bytes:
    """构建 STOP (0x12) 命令帧."""
    return _build_frame(0x12, seq)


def build_soft_reset_frame(seq: int = 0) -> bytes:
    """构建 SOFT_RESET (0x13) 命令帧."""
    return _build_frame(0x13, seq)


# =============================================================================
# Hex dump 打印
# =============================================================================
def hex_dump(data: bytes, title: str = ""):
    """
    以 xxd 风格打印字节序列.

    格式示例:
        0000: AA 55 10 00 4B 00 ...  ..U..K.

    参数:
        data: 待打印的字节序列
        title: 可选的标题
    """
    if title:
        print(f"\n=== {title} ===")

    for i in range(0, len(data), 16):
        chunk = data[i:i + 16]

        # 偏移量
        offset_str = f"{i:04X}:"

        # Hex 部分
        hex_parts = []
        for j in range(16):
            if j < len(chunk):
                hex_parts.append(f"{chunk[j]:02X}")
            else:
                hex_parts.append("  ")
        # 分成两组8字节，中间加空格
        hex_str = " ".join(hex_parts[:8]) + "  " + " ".join(hex_parts[8:])

        # ASCII 部分
        ascii_str = ""
        for b in chunk:
            if 32 <= b <= 126:
                ascii_str += chr(b)
            else:
                ascii_str += "."

        print(f"{offset_str} {hex_str}  {ascii_str}")


# =============================================================================
# pyserial 发送 stub
# =============================================================================
def send_frame(port: str, baud: int, frame: bytes):
    """
    通过串口发送帧 (可选功能，不强制依赖 pyserial).

    参数:
        port: 串口名称 (如 "COM3" 或 "/dev/ttyUSB0")
        baud: 波特率
        frame: 待发送的帧字节序列
    """
    try:
        import serial
        with serial.Serial(port, baud, timeout=1) as ser:
            ser.write(frame)
            print(f"Sent {len(frame)} bytes to {port}")
    except ImportError:
        print("pyserial not installed. Frame not sent.")


# =============================================================================
# 示例运行
# =============================================================================
def run():
    """打印4种命令的hex dump示例，并输出CRC值供交叉验证."""

    # -------------------------------------------------------------------------
    # 1. WRITE_PARAMS: 25路全部 phase=0, amplitude=255
    # -------------------------------------------------------------------------
    phases = [0] * 25
    amplitudes = [255] * 25

    write_frame = build_write_params_frame(seq=0, phases=phases, amplitudes=amplitudes)
    hex_dump(write_frame, "WRITE_PARAMS (0x10) - 25ch phase=0, amp=255")

    # 提取并打印CRC值 (最后两个字节)
    crc_in_frame = write_frame[-2] | (write_frame[-1] << 8)
    print(f"CRC-16/CCITT-FALSE: 0x{crc_in_frame:04X} ({crc_in_frame})")
    print(f"  -> 可在 https://crccalc.com/ 选择 'CRC-16/CCITT-FALSE' 验证")
    print(f"  -> CRC数据域 (cmd+seq+len+payload): {write_frame[2:-2].hex()}")

    # -------------------------------------------------------------------------
    # 2. START
    # -------------------------------------------------------------------------
    start_frame = build_start_frame(seq=0)
    hex_dump(start_frame, "START (0x11)")
    crc_start = start_frame[-2] | (start_frame[-1] << 8)
    print(f"CRC-16/CCITT-FALSE: 0x{crc_start:04X} ({crc_start})")

    # -------------------------------------------------------------------------
    # 3. STOP
    # -------------------------------------------------------------------------
    stop_frame = build_stop_frame(seq=0)
    hex_dump(stop_frame, "STOP (0x12)")
    crc_stop = stop_frame[-2] | (stop_frame[-1] << 8)
    print(f"CRC-16/CCITT-FALSE: 0x{crc_stop:04X} ({crc_stop})")

    # -------------------------------------------------------------------------
    # 4. SOFT_RESET
    # -------------------------------------------------------------------------
    reset_frame = build_soft_reset_frame(seq=0)
    hex_dump(reset_frame, "SOFT_RESET (0x13)")
    crc_reset = reset_frame[-2] | (reset_frame[-1] << 8)
    print(f"CRC-16/CCITT-FALSE: 0x{crc_reset:04X} ({crc_reset})")

    print("\n" + "=" * 60)
    print("所有帧构建完成。")
    print("提示: 使用 send_frame('COM3', 921600, frame) 可通过串口发送。")


if __name__ == "__main__":
    run()
