import sys
import os
import serial
import serial.tools.list_ports
import datetime
import struct
import math
from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QVBoxLayout, QHBoxLayout,
    QWidget, QPushButton, QComboBox, QTextEdit, QLabel,
    QMessageBox, QSlider, QGroupBox, QGridLayout
)
from PyQt6.QtCore import QTimer, Qt
from PyQt6.QtGui import QFont, QIcon

# ==========================================
# 资源路径处理函数
# ==========================================
def resource_path(relative_path):
    if hasattr(sys, '_MEIPASS'):
        return os.path.join(sys._MEIPASS, relative_path)
    return os.path.join(os.path.abspath("."), relative_path)

# ==========================================
# 🌌 V18 极客深空·相控阵专属样式表
# ==========================================
STYLESHEET = """
QMainWindow {
    background-color: #1E1E1E; 
    border-image: url('bg.jpg') 0 0 0 0 stretch stretch;
}
QWidget {
    font-family: "Microsoft YaHei UI", "Consolas";
    font-size: 14px;
    color: #E0E0E0;
}
QWidget#LeftPanel {
    background-color: rgba(20, 20, 25, 0.7); 
    border-right: 1px solid rgba(64, 196, 255, 0.2); 
}
QTextEdit, QComboBox {
    background-color: rgba(0, 0, 0, 0.4); 
    border: 1px solid rgba(255, 255, 255, 0.15);
    border-radius: 4px; 
    padding: 6px;
    color: #00E5FF; 
}
QPushButton {
    background-color: rgba(255, 255, 255, 0.08); 
    border: 1px solid rgba(255, 255, 255, 0.1);
    color: #E0E0E0;
    padding: 8px;
    border-radius: 6px;
}
QPushButton:hover {
    background-color: rgba(64, 196, 255, 0.25); 
    border: 1px solid #40C4FF;
    color: #40C4FF;
}
QPushButton#btn_action {
    background-color: rgba(0, 188, 212, 0.5); 
    border: 1px solid rgba(0, 229, 255, 0.6);
    font-weight: bold;
}
QPushButton#btn_danger {
    background-color: rgba(255, 87, 34, 0.5); 
    border: 1px solid rgba(255, 87, 34, 0.6);
}
QGroupBox {
    border: 1px solid rgba(64, 196, 255, 0.3);
    border-radius: 8px;
    margin-top: 15px;
    padding-top: 20px;
    color: #40C4FF; 
    font-weight: bold;
    background-color: rgba(0, 0, 0, 0.2); 
}
QGroupBox::title { subcontrol-origin: margin; left: 10px; padding: 0 5px; }
QSlider::groove:horizontal {
    border: 1px solid #999999;
    height: 8px;
    background: rgba(0, 0, 0, 0.5);
    border-radius: 4px;
}
QSlider::handle:horizontal {
    background: #00E5FF;
    border: 1px solid #00E5FF;
    width: 18px;
    margin: -5px 0; 
    border-radius: 9px;
}
/* 阵列探头 Label 样式 */
QLabel.transducer {
    background-color: rgba(0, 229, 255, 0.05);
    border: 1px solid rgba(0, 229, 255, 0.3);
    border-radius: 35px; /* 变成圆形 */
    color: #FFF;
    font-family: "Consolas";
    font-size: 12px;
}
"""

class SerialLogic:
    def __init__(self):
        self.ser = None

    def open(self, port, baud):
        if self.ser and self.ser.is_open:
            return False, "已打开"
        try:
            self.ser = serial.Serial(port=port.split(" ")[0], baudrate=int(baud), timeout=0)
            return True, "成功"
        except Exception as e:
            return False, str(e)

    def close(self):
        if self.ser:
            self.ser.close()

    def send(self, data):
        if not self.ser or not self.ser.is_open:
            return False, "未连接"
        try:
            self.ser.write(data)
            return True, "成功"
        except Exception as e:
            return False, str(e)


class PhasedArrayWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("相控阵波束控制终端 - Phased Array Controller")
        self.resize(1100, 800)
        
        icon_path = resource_path("icon.jpg")
        if os.path.exists(icon_path):
            self.setWindowIcon(QIcon(icon_path))

        self.logic = SerialLogic()
        
        # 物理常数与 FPGA 参数
        self.d = 0.0172           # 阵元间距 17.2mm
        self.c = 340.0            # 声速 340m/s
        self.fpga_clk = 100_000_000 # 100MHz 
        self.period = 2500        # 2500 计数值 = 25us (40kHz)
        
        # 【新增】物理网格(0~24) 到 FPGA硬件通道(0~24) 的蛇形映射表
        # 对应硬件 PCB 的 S 型走线 CH1 ~ CH25
        self.grid_to_hw = [
             4,  3,  2,  1,  0,
             5,  6,  7,  8,  9,
            14, 13, 12, 11, 10,
            15, 16, 17, 18, 19,
            24, 23, 22, 21, 20
        ]
        
        # 权重矩阵预置
        self.windows = {
            "高斯窗 (Gaussian)": [
                5, 21, 35, 21, 5,
                21, 94, 155, 94, 21,
                35, 155, 255, 155, 35,
                21, 94, 155, 94, 21,
                5, 21, 35, 21, 5
            ],
            "汉明窗 (Hamming)": [
                2, 11, 20, 11, 2,
                11, 74, 138, 74, 11,
                20, 138, 255, 138, 20,
                11, 74, 138, 74, 11,
                2, 11, 20, 11, 2
            ],
            "矩形窗 (Rectangular, 漏音大)": [255] * 25
        }
        
        # 保存发送给硬件的 25 个真实参数槽位
        self.current_amps = [255] * 25
        self.current_phases = [0] * 25
        self.transducer_labels = []

        self.init_ui()
        self.refresh_ports()
        self.calculate_array() # 初始化时计算一次默认界面

    def init_ui(self):
        main_widget = QWidget()
        self.setCentralWidget(main_widget)
        main_layout = QHBoxLayout(main_widget)
        main_layout.setContentsMargins(0, 0, 0, 0)
        main_layout.setSpacing(0)

        # ===============================================
        # 1. 左侧面板 (串口连接与日志)
        # ===============================================
        left_panel = QWidget()
        left_panel.setObjectName("LeftPanel")
        left_panel.setFixedWidth(300)
        left_layout = QVBoxLayout(left_panel)
        
        lbl_logo = QLabel("📡 硬件连接 (LINK)")
        lbl_logo.setStyleSheet("color: #00E5FF; font-weight: bold; font-size: 18px;")
        left_layout.addWidget(lbl_logo)
        
        self.cmb_port = QComboBox()
        self.cmb_baud = QComboBox()
        self.cmb_baud.addItems(["115200", "921600", "9600"])
        
        left_layout.addWidget(QLabel("端口号 (Port):"))
        left_layout.addWidget(self.cmb_port)
        btn_refresh = QPushButton("🔄 刷新端口")
        btn_refresh.clicked.connect(self.refresh_ports)
        left_layout.addWidget(btn_refresh)
        
        left_layout.addWidget(QLabel("波特率 (Baudrate):"))
        left_layout.addWidget(self.cmb_baud)
        
        self.btn_open = QPushButton("🔗 连接 FPGA (CONNECT)")
        self.btn_open.setObjectName("btn_action")
        self.btn_open.setMinimumHeight(45)
        self.btn_open.clicked.connect(self.toggle_serial)
        left_layout.addWidget(self.btn_open)
        
        left_layout.addSpacing(20)
        left_layout.addWidget(QLabel("通信日志 (LOG):"))
        self.txt_log = QTextEdit()
        self.txt_log.setReadOnly(True)
        left_layout.addWidget(self.txt_log)
        
        btn_clear = QPushButton("🗑 清空日志")
        btn_clear.clicked.connect(self.txt_log.clear)
        left_layout.addWidget(btn_clear)

        # ===============================================
        # 2. 右侧面板 (阵列解算与控制)
        # ===============================================
        right_panel = QWidget()
        right_layout = QVBoxLayout(right_panel)
        right_layout.setContentsMargins(20, 20, 20, 20)
        
        # --- 顶部控制区 ---
        grp_ctrl = QGroupBox("🎛️ 波束偏转与赋形控制 (Beam Steering & Apodization)")
        ctrl_l = QVBoxLayout(grp_ctrl)
        
        # 窗函数选择
        h_win = QHBoxLayout()
        h_win.addWidget(QLabel("空间加窗模式 (降低漏音):"))
        self.cmb_win = QComboBox()
        self.cmb_win.addItems(self.windows.keys())
        self.cmb_win.currentIndexChanged.connect(self.calculate_array)
        h_win.addWidget(self.cmb_win)
        h_win.addStretch()
        ctrl_l.addLayout(h_win)
        
        # 水平角度
        h_az = QHBoxLayout()
        self.lbl_az = QLabel("水平偏转角 (Azimuth θ): 0°")
        self.lbl_az.setFixedWidth(200)
        self.sld_az = QSlider(Qt.Orientation.Horizontal)
        self.sld_az.setRange(-60, 60)
        self.sld_az.setValue(0)
        self.sld_az.valueChanged.connect(self.calculate_array)
        h_az.addWidget(self.lbl_az)
        h_az.addWidget(self.sld_az)
        ctrl_l.addLayout(h_az)

        # 垂直角度
        h_el = QHBoxLayout()
        self.lbl_el = QLabel("垂直偏转角 (Elevation φ): 0°")
        self.lbl_el.setFixedWidth(200)
        self.sld_el = QSlider(Qt.Orientation.Horizontal)
        self.sld_el.setRange(-60, 60)
        self.sld_el.setValue(0)
        self.sld_el.valueChanged.connect(self.calculate_array)
        h_el.addWidget(self.lbl_el)
        h_el.addWidget(self.sld_el)
        ctrl_l.addLayout(h_el)
        
        right_layout.addWidget(grp_ctrl)

        # --- 中间阵列显示区 ---
        grp_array = QGroupBox("🚀 5x5 物理阵列实时映射图 (硬件视角)")
        array_l = QGridLayout(grp_array)
        array_l.setSpacing(10)
        
        # 创建 5x5 的 Label 阵列 (UI 严格按物理空间排列)
        for row in range(5):
            for col in range(5):
                lbl = QLabel("CH--\nA: 0\nP: 0") # 加入了 CH 标号
                lbl.setProperty("class", "transducer")
                lbl.setAlignment(Qt.AlignmentFlag.AlignCenter)
                lbl.setFixedSize(70, 70)
                array_l.addWidget(lbl, row, col, alignment=Qt.AlignmentFlag.AlignCenter)
                self.transducer_labels.append(lbl)
                
        right_layout.addWidget(grp_array, stretch=1)
        
        # --- 底部发送按钮 ---
        self.btn_send = QPushButton("⚡ 下发参数至 FPGA (SYNC)")
        self.btn_send.setObjectName("btn_action")
        self.btn_send.setMinimumHeight(55)
        self.btn_send.setFont(QFont("Microsoft YaHei UI", 14, QFont.Weight.Bold))
        self.btn_send.clicked.connect(self.send_to_fpga)
        right_layout.addWidget(self.btn_send)

        main_layout.addWidget(left_panel)
        main_layout.addWidget(right_panel)

    def refresh_ports(self):
        self.cmb_port.clear()
        for p in serial.tools.list_ports.comports():
            self.cmb_port.addItem(f"{p.device} ({p.description})")

    def toggle_serial(self):
        if self.btn_open.text() == "🔗 连接 FPGA (CONNECT)":
            ok, msg = self.logic.open(self.cmb_port.currentText(), self.cmb_baud.currentText())
            if ok:
                self.btn_open.setText("🚫 断开连接 (DISCONNECT)")
                self.btn_open.setObjectName("btn_danger")
                self.btn_open.setStyleSheet("background-color: rgba(255, 87, 34, 0.6);")
                self.log(f"成功连接至 {self.cmb_port.currentText()}", "green")
            else:
                QMessageBox.critical(self, "错误", f"打开串口失败: {msg}")
        else:
            self.logic.close()
            self.btn_open.setText("🔗 连接 FPGA (CONNECT)")
            self.btn_open.setObjectName("btn_action")
            self.btn_open.setStyleSheet("")
            self.log("连接已断开", "red")

    def calculate_array(self):
        az_deg = self.sld_az.value()
        el_deg = self.sld_el.value()
        self.lbl_az.setText(f"水平偏转角 (Azimuth θ): {az_deg}°")
        self.lbl_el.setText(f"垂直偏转角 (Elevation φ): {el_deg}°")
        
        theta = math.radians(az_deg)
        phi = math.radians(el_deg)
        
        win_key = self.cmb_win.currentText()
        base_amps = self.windows[win_key] # 理论物理网格幅度 (0~24)
        
        raw_ticks_grid = []
        # 物理坐标系：行(y)从上到下，列(x)从左到右 (左上角为 -2, 2)
        y_coords = [2, 1, 0, -1, -2]
        x_coords = [-2, -1, 0, 1, 2]
        
        # 步骤 1：遍历 5x5 物理空间网格，算出理想的延迟时间
        for r in range(5):
            for c in range(5):
                y_idx = y_coords[r]
                x_idx = x_coords[c]
                dt = (x_idx * self.d * math.sin(theta) + y_idx * self.d * math.sin(phi)) / self.c
                ticks = dt * self.fpga_clk
                raw_ticks_grid.append(ticks)
                
        min_tick = min(raw_ticks_grid)
        
        # 步骤 2：通过 S 型映射表，将物理参数填入 FPGA 硬件通道槽位
        for r in range(5):
            for c in range(5):
                grid_idx = r * 5 + c               # UI 和物理网格索引 (0~24)
                hw_idx = self.grid_to_hw[grid_idx] # FPGA 硬件通道索引 (0~24)
                
                # 计算相位并装入正确的硬件发包槽位
                phase_val = int(round(raw_ticks_grid[grid_idx] - min_tick)) % self.period
                self.current_phases[hw_idx] = phase_val
                
                # 将幅度装入正确的硬件发包槽位
                amp_val = base_amps[grid_idx]
                self.current_amps[hw_idx] = amp_val
                
                # 更新 UI 阵列气泡
                lbl = self.transducer_labels[grid_idx]
                # CH 编号通常从 1 开始，所以显示 hw_idx + 1
                lbl.setText(f"CH{hw_idx+1}\nA: {amp_val}\nP: {phase_val}")
                
                # 高亮最早发波的波阵面前沿 (P=0)
                if phase_val == 0:
                    lbl.setStyleSheet("background-color: rgba(255, 87, 34, 0.6); border: 2px solid #FF5722; border-radius: 35px;")
                else:
                    lbl.setStyleSheet("background-color: rgba(0, 229, 255, 0.1); border: 1px solid rgba(0, 229, 255, 0.3); border-radius: 35px;")

    def send_to_fpga(self):
        # 组装 79 字节通信协议帧
        packet = bytearray([0xAA, 0xBB, 0x01]) # 帧头 + 指令码
        
        # 25字节 幅度 (uint8)，注意这里发的是映射过后的硬件槽位数组
        for amp in self.current_amps:
            packet.append(amp & 0xFF)
            
        # 50字节 相位 (uint16 小端)，同理，发硬件槽位数组
        for ph in self.current_phases:
            packet.extend(struct.pack('<H', ph))
            
        # 校验和 (单字节, 从指令码开始累加)
        checksum = sum(packet[2:]) % 256
        packet.append(checksum)
        
        # 帧尾: 0x0D 0x0A (\r\n)
        packet.extend([0x0D, 0x0A])
        
        # 硬件发送
        ok, msg = self.logic.send(packet)
        if ok:
            # 格式化打印十六进制
            hex_str = " ".join([f"{b:02X}" for b in packet])
            self.log(f"成功下发参数 [{len(packet)} Bytes]<br>{hex_str}", "#00E5FF")
        else:
            self.log(f"发送失败: {msg}", "red")

    def log(self, text, color="#FFF"):
        ts = datetime.datetime.now().strftime('%H:%M:%S.%f')[:-3]
        self.txt_log.append(f"<span style='color:{color};'>[{ts}] {text}</span>")


if __name__ == "__main__":
    app = QApplication(sys.argv)
    
    # 动态加载样式表中的背景图片
    bg_path = resource_path("bg.jpg").replace("\\", "/") 
    new_stylesheet = STYLESHEET.replace("bg.jpg", bg_path)
    app.setStyleSheet(new_stylesheet)
    
    win = PhasedArrayWindow()
    win.show()
    sys.exit(app.exec())