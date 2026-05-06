//=============================================================================
// Module       : mcp3201_driver
// Description  : MCP3201 12-bit SPI ADC Driver
//                - 50MHz system clock input
//                - 1MHz SPI clock output (50x divider)
//                - 40kHz continuous sampling rate
//                - SPI Mode 0,0 (CPOL=0, CPHA=0)
//                - 17 SPI clock cycles per conversion
// Target       : Generic FPGA (Verilog-2001 synthesizable)
// Author       : FPGA Design Engineer
// Date         : 2024
//=============================================================================

module mcp3201_driver (
    input  wire       clk,           // 50MHz system clock
    input  wire       rst_n,         // Active-low asynchronous reset

    // ADC SPI interface
    output reg        adc_cs,        // ADC chip select, active low
    output reg        adc_clk,       // ADC SPI clock (1MHz max)
    input  wire       adc_data,      // ADC serial data output (MISO)

    // Sample data output
    output reg  [11:0] sample_data,  // 12-bit ADC conversion result
    output reg         sample_valid, // Data valid flag (1 clk pulse when new data ready)
    output reg         sample_busy   // Conversion in progress flag
);

    //=========================================================================
    // Parameter Definitions
    //=========================================================================
    // System and timing parameters (adjustable for different clock rates)
    parameter CLK_FREQ        = 50_000_000;   // System clock: 50 MHz
    parameter SPI_CLK_FREQ    = 1_000_000;    // SPI clock: 1 MHz (MCP3201 max 1.6MHz)
    parameter SAMPLE_RATE     = 40_000;       // Sampling rate: 40 kHz

    // Derived timing constants
    parameter SPI_HALF_PERIOD = CLK_FREQ / SPI_CLK_FREQ / 2;  // 25 (50MHz->1MHz)
    parameter SPI_FULL_PERIOD = SPI_HALF_PERIOD * 2;          // 50 clocks per SPI cycle
    parameter SAMPLE_INTERVAL = CLK_FREQ / SAMPLE_RATE;       // 1250 (40kHz period)

    //=========================================================================
    // Local State Definitions
    //=========================================================================
    // Three-state FSM: IDLE -> READING -> DONE -> IDLE
    localparam STATE_IDLE    = 2'b00;   // Waiting for next 40kHz tick
    localparam STATE_READING = 2'b01;   // SPI communication in progress
    localparam STATE_DONE    = 2'b10;   // Conversion complete, output data

    //=========================================================================
    // Internal Registers
    //=========================================================================
    reg [1:0]  state;              // Current FSM state
    reg [10:0] tick_cnt;           // 40kHz tick counter (0 ~ SAMPLE_INTERVAL-1)
                                   // Free-running at 50MHz, wraps every 1250 clocks
    reg [5:0]  spi_div_cnt;        // SPI clock divider (0 ~ SPI_FULL_PERIOD-1)
                                   // Divides 50MHz to 1MHz (50:1 ratio)
    reg [4:0]  bit_cnt;            // SPI bit counter (0 ~ 16, total 17 clock cycles)
                                   // MCP3201 requires 17 SPI clocks for full conversion:
                                   //   bit 0:  null bit (DOUT hi-Z during conversion)
                                   //   bit 1:  start bit (always logic '1')
                                   //   bit 2:  null bit
                                   //   bit 3~14: 12-bit ADC data B11~B0 (MSB first)
                                   //   bit 15~16: sub-LSB bits (optional, can be ignored)
    reg [11:0] data_buf;           // 12-bit data buffer for collecting ADC result

    //=========================================================================
    // FSM + Timing Control - Main Sequential Logic
    //=========================================================================
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            // Reset all registers to known state
            state       <= STATE_IDLE;
            tick_cnt    <= 11'd0;
            spi_div_cnt <= 6'd0;
            bit_cnt     <= 5'd0;
            data_buf    <= 12'd0;
            adc_cs      <= 1'b1;       // CS inactive (high)
            adc_clk     <= 1'b0;       // SPI clock idle low for Mode 0,0
            sample_data <= 12'd0;
            sample_valid<= 1'b0;
            sample_busy <= 1'b0;
        end else begin
            //---------------------------------------------------------------
            // Free-running 40kHz tick counter
            // Wraps every SAMPLE_INTERVAL (1250) clocks @ 50MHz
            // This ensures precise 40kHz sampling rate regardless of
            // conversion time (as long as conversion < 1250 clocks)
            // 17 SPI cycles x 50 = 850 clocks < 1250, so timing is safe
            //---------------------------------------------------------------
            tick_cnt <= (tick_cnt >= SAMPLE_INTERVAL - 1) ? 11'd0 : tick_cnt + 1'b1;

            //---------------------------------------------------------------
            // State Machine
            //---------------------------------------------------------------
            case (state)

                //===========================================================
                // STATE_IDLE: Wait for next 40kHz sampling tick
                // ADC_CS is high (inactive), ADC_CLK is low (idle)
                //===========================================================
                STATE_IDLE: begin
                    adc_cs       <= 1'b1;
                    adc_clk      <= 1'b0;
                    sample_valid <= 1'b0;
                    sample_busy  <= 1'b0;

                    // At tick_cnt==0, initiate a new ADC conversion
                    // This triggers exactly every 1250 clocks = 40kHz
                    if (tick_cnt == 11'd0) begin
                        state       <= STATE_READING;
                        adc_cs      <= 1'b0;         // Pull CS low to start conversion
                        adc_clk     <= 1'b1;         // Initialize CLK high so first
                                                     // transition is a falling edge
                        spi_div_cnt <= 6'd0;
                        bit_cnt     <= 5'd0;
                        data_buf    <= 12'd0;
                        sample_busy <= 1'b1;
                    end
                end

                //===========================================================
                // STATE_READING: SPI communication with MCP3201
                // Generate 1MHz SPI clock and sample data on falling edges
                //===========================================================
                STATE_READING: begin
                    //--- SPI Clock Generation (1MHz from 50MHz) ---
                    // adc_clk is high for first half (0~24), low for second half (25~49)
                    // Duty cycle ~50%, frequency = 50MHz/50 = 1MHz
                    if (spi_div_cnt < SPI_HALF_PERIOD[5:0])
                        adc_clk <= 1'b1;
                    else
                        adc_clk <= 1'b0;

                    //--- Data Sampling on Falling Edge ---
                    // In SPI Mode 0,0: data changes on rising edge, sampled on falling edge
                    // Falling edge occurs when spi_div_cnt == SPI_HALF_PERIOD (25)
                    if (spi_div_cnt == SPI_HALF_PERIOD[5:0]) begin
                        // MCP3201 output data format over 17 SPI clocks:
                        //  bit_cnt=0:   null bit (DOUT is hi-Z during internal conversion)
                        //  bit_cnt=1:   start bit = 1 (conversion started indicator)
                        //  bit_cnt=2:   null bit
                        //  bit_cnt=3:   B11 (MSB of 12-bit ADC result)
                        //  bit_cnt=4:   B10
                        //  ...
                        //  bit_cnt=14:  B0  (LSB of 12-bit ADC result)
                        //  bit_cnt=15:  sub-LSB (sub-resolution bit, may be ignored)
                        //  bit_cnt=16:  sub-LSB (sub-resolution bit, may be ignored)

                        // Capture only the 12 valid data bits (bit_cnt 3 through 14)
                        // Store B11 at data_buf[11], B10 at data_buf[10], ..., B0 at data_buf[0]
                        if (bit_cnt >= 5'd3 && bit_cnt <= 5'd14) begin
                            data_buf[14 - bit_cnt] <= adc_data;
                        end
                    end

                    //--- SPI Clock Divider Counter ---
                    // Count from 0 to SPI_FULL_PERIOD-1 (0 to 49)
                    if (spi_div_cnt >= SPI_FULL_PERIOD - 1) begin
                        // One complete SPI clock cycle done
                        spi_div_cnt <= 6'd0;

                        // After 17 SPI clock cycles (bit_cnt 0~16), conversion is complete
                        if (bit_cnt >= 5'd16) begin
                            state <= STATE_DONE;
                        end else begin
                            bit_cnt <= bit_cnt + 1'b1;
                        end
                    end else begin
                        spi_div_cnt <= spi_div_cnt + 1'b1;
                    end
                end

                //===========================================================
                // STATE_DONE: Conversion complete, latch output data
                // Deassert CS, pulse sample_valid, return to IDLE
                //===========================================================
                STATE_DONE: begin
                    adc_cs       <= 1'b1;          // Deactivate chip select
                    adc_clk      <= 1'b0;          // Return SPI clock to idle (low)
                    sample_data  <= data_buf;      // Output the 12-bit conversion result
                    sample_valid <= 1'b1;          // Pulse data valid flag (1 clock cycle)
                    sample_busy  <= 1'b0;          // Clear busy flag
                    state        <= STATE_IDLE;    // Return to IDLE for next cycle
                end

                //===========================================================
                // Default: Safety fallback to IDLE
                //===========================================================
                default: begin
                    state   <= STATE_IDLE;
                    adc_cs  <= 1'b1;
                    adc_clk <= 1'b0;
                end

            endcase
        end
    end

endmodule
