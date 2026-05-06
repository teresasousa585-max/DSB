//=============================================================================
// 5x5 ultrasonic array top level
//
// Control plane:
//   Host sends UART command frames at 921600 baud. Full 25-channel phase/amplitude
//   tables are received into a shadow buffer and committed on a 40 kHz carrier
//   boundary.
//
// Real-time plane:
//   ADC -> envelope sampler -> 25-channel PWM generator -> transducer outputs.
//=============================================================================

`timescale 1ns / 1ps

module top_ultrasound_array #(
    parameter integer CLK_FREQ       = 50_000_000,
    parameter integer UART_BAUD      = 921_600,
    parameter integer CARRIER_PERIOD = 1250
) (
    input  wire        clk,
    input  wire        rst_n,

    input  wire        uart_rx,
    output wire [7:0]  cmd_status,
    output wire [7:0]  last_cmd,
    output wire [7:0]  last_seq,
    output wire        frame_ok,
    output wire        frame_error,
    output wire        ultrasound_en,

    output wire        adc_cs,
    output wire        adc_clk,
    input  wire        adc_data,

    output wire [24:0] pwm_out
);

    localparam integer CARRIER_CNT_WIDTH = 11;
    localparam [CARRIER_CNT_WIDTH-1:0] CARRIER_PERIOD_REG = CARRIER_PERIOD;

    reg [CARRIER_CNT_WIDTH-1:0] carrier_cnt;
    wire carrier_tick;

    wire [11:0] adc_sample_data;
    wire        adc_sample_valid;
    wire        adc_sample_busy;

    wire [7:0]  envelope;
    wire        envelope_valid;

    wire [7:0] beam_amplitude [0:24];
    wire [9:0] beam_phase     [0:24];
    wire       beam_param_valid;
    wire       ultrasound_soft_rst_cmd;

    reg [3:0] soft_rst_cnt;
    wire      ultra_rst_n;
    wire      ultra_run;

    assign carrier_tick = (carrier_cnt == CARRIER_PERIOD_REG - 1'b1);
    assign ultra_rst_n  = rst_n && (soft_rst_cnt == 4'd0);
    assign ultra_run    = ultra_rst_n && ultrasound_en;

    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            soft_rst_cnt <= 4'd0;
        end else if (ultrasound_soft_rst_cmd) begin
            soft_rst_cnt <= 4'd8;
        end else if (soft_rst_cnt != 4'd0) begin
            soft_rst_cnt <= soft_rst_cnt - 1'b1;
        end
    end

    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            carrier_cnt <= {CARRIER_CNT_WIDTH{1'b0}};
        end else if (!ultra_rst_n) begin
            carrier_cnt <= {CARRIER_CNT_WIDTH{1'b0}};
        end else if (carrier_tick) begin
            carrier_cnt <= {CARRIER_CNT_WIDTH{1'b0}};
        end else begin
            carrier_cnt <= carrier_cnt + 1'b1;
        end
    end

    uart_beam_cmd #(
        .CLK_FREQ (CLK_FREQ),
        .BAUD     (UART_BAUD),
        .CH_NUM   (25)
    ) u_uart_beam_cmd (
        .clk                  (clk),
        .rst_n                (rst_n),
        .uart_rx              (uart_rx),
        .carrier_tick         (carrier_tick),
        .amplitude            (beam_amplitude),
        .phase                (beam_phase),
        .param_valid          (beam_param_valid),
        .ultrasound_en        (ultrasound_en),
        .ultrasound_soft_rst  (ultrasound_soft_rst_cmd),
        .frame_ok             (frame_ok),
        .frame_error          (frame_error),
        .last_cmd             (last_cmd),
        .last_seq             (last_seq),
        .status               (cmd_status)
    );

    mcp3201_driver u_mcp3201_driver (
        .clk           (clk),
        .rst_n         (rst_n),
        .adc_cs        (adc_cs),
        .adc_clk       (adc_clk),
        .adc_data      (adc_data),
        .sample_data   (adc_sample_data),
        .sample_valid  (adc_sample_valid),
        .sample_busy   (adc_sample_busy)
    );

    dsb_modulator u_dsb_modulator (
        .clk            (clk),
        .rst_n          (ultra_rst_n),
        .en             (ultra_run),
        .sample_data    (adc_sample_data),
        .sample_valid   (adc_sample_valid),
        .carrier_tick   (carrier_tick),
        .envelope       (envelope),
        .envelope_valid (envelope_valid)
    );

    pwm25_generator u_pwm25_generator (
        .clk             (clk),
        .rst_n           (ultra_rst_n),
        .en              (ultra_run),
        .carrier_tick    (carrier_tick),
        .carrier_period  (CARRIER_PERIOD_REG),
        .amplitude       (beam_amplitude),
        .phase           (beam_phase),
        .param_valid     (beam_param_valid),
        .envelope        (envelope),
        .envelope_valid  (envelope_valid),
        .pwm_out         (pwm_out)
    );

    // Keep debug-only status signals visible during simulation.
    // synopsys translate_off
    wire _unused_debug = &{1'b0, adc_sample_busy, 1'b0};
    // synopsys translate_on

endmodule
