import argparse
import sys

from analysis.lut_generator import run as run_lut
from analysis.array_sim import run as run_array
from analysis.dsb_sim import run as run_dsb
from analysis.uart_cmd import run as run_uart


def main():
    parser = argparse.ArgumentParser(description="FPGA超声阵列数据分析与仿真工具")
    subparsers = parser.add_subparsers(dest="command", required=True)

    subparsers.add_parser("lut", help="LUT计算与对拍")
    subparsers.add_parser("array", help="波阵面与方向图仿真")
    subparsers.add_parser("dsb", help="DSB信号与频谱仿真")
    subparsers.add_parser("uart", help="UART命令生成")
    subparsers.add_parser("all", help="运行全部")

    args = parser.parse_args()

    if args.command == "lut":
        run_lut()
    elif args.command == "array":
        run_array()
    elif args.command == "dsb":
        run_dsb()
    elif args.command == "uart":
        run_uart()
    elif args.command == "all":
        for name, fn in [("lut", run_lut), ("array", run_array), ("dsb", run_dsb), ("uart", run_uart)]:
            print(f"\n{'='*70}")
            print(f"Running: {name}")
            print('='*70)
            fn()
    else:
        parser.print_help()
        sys.exit(1)


if __name__ == "__main__":
    main()
