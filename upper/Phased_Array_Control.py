import datetime
import os
import sys

try:
    import serial
    import serial.tools.list_ports
except ImportError:
    serial = None

from PyQt6.QtCore import QTimer, Qt
from PyQt6.QtGui import QFont, QIcon
from PyQt6.QtWidgets import (
    QApplication,
    QCheckBox,
    QComboBox,
    QGridLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QMainWindow,
    QMessageBox,
    QPushButton,
    QSlider,
    QSpinBox,
    QTextEdit,
    QVBoxLayout,
    QWidget,
)

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
if SCRIPT_DIR not in sys.path:
    sys.path.insert(0, SCRIPT_DIR)

from n32_beam_params import (  # noqa: E402
    AMPLITUDE_MODES,
    MAX_DISTANCE_MM,
    MAX_STEER_DEG,
    MIN_DISTANCE_MM,
    MODE_FITTED,
    build_packet,
    calculate_beam_params,
    format_packet_hex,
    load_layout,
)


def resource_path(relative_path):
    if hasattr(sys, "_MEIPASS"):
        return os.path.join(sys._MEIPASS, relative_path)
    return os.path.join(os.path.abspath("."), relative_path)


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
    background-color: rgba(20, 20, 25, 0.72);
    border-right: 1px solid rgba(64, 196, 255, 0.24);
}
QTextEdit, QComboBox, QSpinBox {
    background-color: rgba(0, 0, 0, 0.42);
    border: 1px solid rgba(255, 255, 255, 0.16);
    border-radius: 4px;
    padding: 6px;
    color: #00E5FF;
}
QPushButton {
    background-color: rgba(255, 255, 255, 0.08);
    border: 1px solid rgba(255, 255, 255, 0.12);
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
    background-color: rgba(0, 188, 212, 0.50);
    border: 1px solid rgba(0, 229, 255, 0.60);
    font-weight: bold;
}
QPushButton#btn_danger {
    background-color: rgba(255, 87, 34, 0.55);
    border: 1px solid rgba(255, 87, 34, 0.70);
}
QGroupBox {
    border: 1px solid rgba(64, 196, 255, 0.30);
    border-radius: 8px;
    margin-top: 15px;
    padding-top: 20px;
    color: #40C4FF;
    font-weight: bold;
    background-color: rgba(0, 0, 0, 0.22);
}
QGroupBox::title {
    subcontrol-origin: margin;
    left: 10px;
    padding: 0 5px;
}
QSlider::groove:horizontal {
    border: 1px solid #777777;
    height: 8px;
    background: rgba(0, 0, 0, 0.50);
    border-radius: 4px;
}
QSlider::handle:horizontal {
    background: #00E5FF;
    border: 1px solid #00E5FF;
    width: 18px;
    margin: -5px 0;
    border-radius: 9px;
}
QLabel.transducer {
    background-color: rgba(0, 229, 255, 0.08);
    border: 1px solid rgba(0, 229, 255, 0.30);
    border-radius: 6px;
    color: #FFFFFF;
    font-family: "Consolas";
    font-size: 11px;
}
"""


class SerialLogic:
    def __init__(self):
        self.ser = None

    def open(self, port, baud):
        if serial is None:
            return False, "pyserial is not installed"
        if not port:
            return False, "no serial port selected"
        if self.ser and self.ser.is_open:
            return False, "already open"
        try:
            self.ser = serial.Serial(port=port.split(" ")[0], baudrate=int(baud), timeout=0)
            return True, "ok"
        except Exception as exc:
            return False, str(exc)

    def close(self):
        if self.ser:
            self.ser.close()

    def send(self, data):
        if not self.ser or not self.ser.is_open:
            return False, "not connected"
        try:
            self.ser.write(data)
            return True, "ok"
        except Exception as exc:
            return False, str(exc)


class PhasedArrayWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("N32 Phased Array Controller")
        self.resize(1180, 840)

        icon_path = resource_path("icon.jpg")
        if os.path.exists(icon_path):
            self.setWindowIcon(QIcon(icon_path))

        self.logic = SerialLogic()
        self.coords_mm = load_layout()
        self.current_amps = [255] * 32
        self.current_phases = [0] * 32
        self.transducer_labels = []
        self.auto_timer = QTimer(self)
        self.auto_timer.setSingleShot(True)
        self.auto_timer.timeout.connect(self.send_to_fpga)

        self.init_ui()
        self.refresh_ports()
        self.calculate_array()

    def init_ui(self):
        main_widget = QWidget()
        self.setCentralWidget(main_widget)
        main_layout = QHBoxLayout(main_widget)
        main_layout.setContentsMargins(0, 0, 0, 0)
        main_layout.setSpacing(0)

        left_panel = QWidget()
        left_panel.setObjectName("LeftPanel")
        left_panel.setFixedWidth(315)
        left_layout = QVBoxLayout(left_panel)

        lbl_logo = QLabel("FPGA LINK")
        lbl_logo.setStyleSheet("color: #00E5FF; font-weight: bold; font-size: 18px;")
        left_layout.addWidget(lbl_logo)

        self.cmb_port = QComboBox()
        self.cmb_baud = QComboBox()
        self.cmb_baud.addItems(["115200", "921600", "9600"])

        left_layout.addWidget(QLabel("Port:"))
        left_layout.addWidget(self.cmb_port)
        btn_refresh = QPushButton("Refresh ports")
        btn_refresh.clicked.connect(self.refresh_ports)
        left_layout.addWidget(btn_refresh)

        left_layout.addWidget(QLabel("Baudrate:"))
        left_layout.addWidget(self.cmb_baud)

        self.btn_open = QPushButton("Connect FPGA")
        self.btn_open.setObjectName("btn_action")
        self.btn_open.setMinimumHeight(44)
        self.btn_open.clicked.connect(self.toggle_serial)
        left_layout.addWidget(self.btn_open)

        self.chk_auto_send = QCheckBox("Auto sync after parameter changes")
        left_layout.addWidget(self.chk_auto_send)

        left_layout.addSpacing(16)
        left_layout.addWidget(QLabel("Log:"))
        self.txt_log = QTextEdit()
        self.txt_log.setReadOnly(True)
        left_layout.addWidget(self.txt_log)

        btn_clear = QPushButton("Clear log")
        btn_clear.clicked.connect(self.txt_log.clear)
        left_layout.addWidget(btn_clear)

        right_panel = QWidget()
        right_layout = QVBoxLayout(right_panel)
        right_layout.setContentsMargins(20, 20, 20, 20)

        grp_ctrl = QGroupBox("N32 focus and steering")
        ctrl_l = QVBoxLayout(grp_ctrl)

        h_amp = QHBoxLayout()
        h_amp.addWidget(QLabel("Amplitude mode:"))
        self.cmb_amp = QComboBox()
        self.cmb_amp.addItems(AMPLITUDE_MODES)
        self.cmb_amp.setCurrentText(MODE_FITTED)
        self.cmb_amp.currentIndexChanged.connect(self.calculate_array)
        h_amp.addWidget(self.cmb_amp)
        h_amp.addStretch()
        ctrl_l.addLayout(h_amp)

        self.lbl_az = QLabel()
        self.sld_az = self.make_slider(-int(MAX_STEER_DEG), int(MAX_STEER_DEG), 0)
        self.sld_az.valueChanged.connect(self.calculate_array)
        ctrl_l.addLayout(self.make_labeled_slider(self.lbl_az, self.sld_az))

        self.lbl_el = QLabel()
        self.sld_el = self.make_slider(-int(MAX_STEER_DEG), int(MAX_STEER_DEG), 0)
        self.sld_el.valueChanged.connect(self.calculate_array)
        ctrl_l.addLayout(self.make_labeled_slider(self.lbl_el, self.sld_el))

        h_dist = QHBoxLayout()
        self.lbl_dist = QLabel()
        self.lbl_dist.setFixedWidth(220)
        self.sld_dist = QSlider(Qt.Orientation.Horizontal)
        self.sld_dist.setRange(MIN_DISTANCE_MM, MAX_DISTANCE_MM)
        self.sld_dist.setValue(1000)
        self.sld_dist.valueChanged.connect(self.distance_slider_changed)
        self.spn_dist = QSpinBox()
        self.spn_dist.setRange(MIN_DISTANCE_MM, MAX_DISTANCE_MM)
        self.spn_dist.setValue(1000)
        self.spn_dist.setSuffix(" mm")
        self.spn_dist.valueChanged.connect(self.distance_spin_changed)
        h_dist.addWidget(self.lbl_dist)
        h_dist.addWidget(self.sld_dist)
        h_dist.addWidget(self.spn_dist)
        ctrl_l.addLayout(h_dist)

        self.lbl_limited = QLabel()
        self.lbl_limited.setStyleSheet("color: #FFC107;")
        ctrl_l.addWidget(self.lbl_limited)

        right_layout.addWidget(grp_ctrl)

        grp_array = QGroupBox("N32 array, E00..E31 maps to transducer_io[0]..[31]")
        array_l = QGridLayout(grp_array)
        array_l.setSpacing(6)
        positions = self.grid_positions()
        for ch, (row, col) in enumerate(positions):
            lbl = QLabel()
            lbl.setProperty("class", "transducer")
            lbl.setAlignment(Qt.AlignmentFlag.AlignCenter)
            lbl.setFixedSize(64, 48)
            array_l.addWidget(lbl, row, col, alignment=Qt.AlignmentFlag.AlignCenter)
            self.transducer_labels.append(lbl)
        right_layout.addWidget(grp_array, stretch=1)

        self.btn_send = QPushButton("Sync parameters to FPGA")
        self.btn_send.setObjectName("btn_action")
        self.btn_send.setMinimumHeight(55)
        self.btn_send.setFont(QFont("Microsoft YaHei UI", 14, QFont.Weight.Bold))
        self.btn_send.clicked.connect(self.send_to_fpga)
        right_layout.addWidget(self.btn_send)

        main_layout.addWidget(left_panel)
        main_layout.addWidget(right_panel)

    def make_slider(self, low, high, value):
        slider = QSlider(Qt.Orientation.Horizontal)
        slider.setRange(low, high)
        slider.setValue(value)
        return slider

    def make_labeled_slider(self, label, slider):
        layout = QHBoxLayout()
        label.setFixedWidth(220)
        layout.addWidget(label)
        layout.addWidget(slider)
        return layout

    def grid_positions(self):
        size = 13
        max_r = max(max(abs(x), abs(y)) for x, y in self.coords_mm)
        used = set()
        positions = []
        offsets = [(0, 0), (0, 1), (1, 0), (0, -1), (-1, 0), (1, 1), (-1, -1)]
        for x_mm, y_mm in self.coords_mm:
            col0 = int(round((x_mm + max_r) / (2.0 * max_r) * (size - 1)))
            row0 = int(round((max_r - y_mm) / (2.0 * max_r) * (size - 1)))
            chosen = (row0, col0)
            for dr, dc in offsets:
                row = max(0, min(size - 1, row0 + dr))
                col = max(0, min(size - 1, col0 + dc))
                if (row, col) not in used:
                    chosen = (row, col)
                    break
            used.add(chosen)
            positions.append(chosen)
        return positions

    def refresh_ports(self):
        self.cmb_port.clear()
        if serial is None:
            self.cmb_port.addItem("pyserial missing")
            return
        for port in serial.tools.list_ports.comports():
            self.cmb_port.addItem("{} ({})".format(port.device, port.description))

    def toggle_serial(self):
        if self.btn_open.text() == "Connect FPGA":
            ok, msg = self.logic.open(self.cmb_port.currentText(), self.cmb_baud.currentText())
            if ok:
                self.btn_open.setText("Disconnect FPGA")
                self.btn_open.setObjectName("btn_danger")
                self.btn_open.setStyleSheet("background-color: rgba(255, 87, 34, 0.6);")
                self.log("Connected to {}".format(self.cmb_port.currentText()), "green")
            else:
                QMessageBox.critical(self, "Serial error", "Open failed: {}".format(msg))
        else:
            self.logic.close()
            self.btn_open.setText("Connect FPGA")
            self.btn_open.setObjectName("btn_action")
            self.btn_open.setStyleSheet("")
            self.log("Disconnected", "red")

    def distance_slider_changed(self, value):
        self.spn_dist.blockSignals(True)
        self.spn_dist.setValue(value)
        self.spn_dist.blockSignals(False)
        self.calculate_array()

    def distance_spin_changed(self, value):
        self.sld_dist.blockSignals(True)
        self.sld_dist.setValue(value)
        self.sld_dist.blockSignals(False)
        self.calculate_array()

    def calculate_array(self):
        result = calculate_beam_params(
            distance_mm=self.sld_dist.value(),
            az_deg=self.sld_az.value(),
            el_deg=self.sld_el.value(),
            amp_mode=self.cmb_amp.currentText(),
            coords_mm=self.coords_mm,
        )
        self.current_amps = list(result.amplitudes)
        self.current_phases = list(result.phases)

        self.lbl_az.setText("Azimuth: {} deg".format(self.sld_az.value()))
        self.lbl_el.setText("Elevation: {} deg".format(self.sld_el.value()))
        self.lbl_dist.setText("Focus distance: {:.3f} m".format(result.distance_mm / 1000.0))

        limited = (
            abs(result.actual_az_deg - self.sld_az.value()) > 0.05
            or abs(result.actual_el_deg - self.sld_el.value()) > 0.05
        )
        if limited:
            self.lbl_limited.setText(
                "Vector limit active: az={:.2f} deg, el={:.2f} deg, combined={:.2f} deg".format(
                    result.actual_az_deg,
                    result.actual_el_deg,
                    result.combined_angle_deg,
                )
            )
        else:
            self.lbl_limited.setText("Combined steering angle: {:.2f} deg".format(result.combined_angle_deg))

        for ch, label in enumerate(self.transducer_labels):
            amp = self.current_amps[ch]
            phase = self.current_phases[ch]
            label.setText("E{:02d}\nA:{:03d}\nP:{:04d}".format(ch, amp, phase))
            intensity = max(20, min(220, amp))
            if phase == 0:
                color = "rgba(255, 87, 34, 0.68)"
                border = "#FF5722"
            else:
                color = "rgba(0, 229, 255, {:.2f})".format(0.10 + 0.45 * intensity / 255.0)
                border = "rgba(0, 229, 255, 0.45)"
            label.setStyleSheet(
                "background-color: {}; border: 1px solid {}; border-radius: 6px; color: #FFFFFF;".format(
                    color,
                    border,
                )
            )

        if self.chk_auto_send.isChecked():
            self.auto_timer.start(120)

    def make_packet(self):
        return build_packet(self.current_amps, self.current_phases)

    def send_to_fpga(self):
        packet = self.make_packet()
        ok, msg = self.logic.send(packet)
        if ok:
            self.log("Sent N32 params [{} bytes]<br>{}".format(len(packet), format_packet_hex(packet)), "#00E5FF")
        else:
            self.log("Send failed: {}".format(msg), "red")

    def log(self, text, color="#FFFFFF"):
        ts = datetime.datetime.now().strftime("%H:%M:%S.%f")[:-3]
        self.txt_log.append("<span style='color:{};'>[{}] {}</span>".format(color, ts, text))


if __name__ == "__main__":
    app = QApplication(sys.argv)
    bg_path = resource_path("bg.jpg").replace("\\", "/")
    app.setStyleSheet(STYLESHEET.replace("bg.jpg", bg_path))
    win = PhasedArrayWindow()
    win.show()
    sys.exit(app.exec())
