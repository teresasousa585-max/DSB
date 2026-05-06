//==============================================================================
// Module       : dsb_modulator
// Description  : Unsigned ADC envelope sampler for 40 kHz PWM drive
//
// The 40 kHz carrier is generated in pwm25_generator. This block only maps the
// latest unsigned ADC sample into an 8-bit envelope and publishes it on the next
// carrier boundary so PWM duty updates are synchronized to the carrier period.
//
// Input scale:
//   sample_data: 12-bit unsigned ADC value, 0..4095
//   envelope   : 8-bit unsigned envelope, 0..255
//==============================================================================

module dsb_modulator (
    input  wire        clk,
    input  wire        rst_n,
    input  wire        en,

    input  wire [11:0] sample_data,
    input  wire        sample_valid,

    input  wire        carrier_tick,

    output reg  [7:0]  envelope,
    output reg         envelope_valid
);

    reg [7:0] pending_envelope;
    reg       sample_seen;

    wire [7:0] sample_envelope;
    assign sample_envelope = sample_data[11:4];

    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            pending_envelope <= 8'd0;
            sample_seen      <= 1'b0;
            envelope         <= 8'd0;
            envelope_valid   <= 1'b0;
        end else if (!en) begin
            pending_envelope <= 8'd0;
            sample_seen      <= 1'b0;
            envelope         <= 8'd0;
            envelope_valid   <= 1'b0;
        end else begin
            envelope_valid <= 1'b0;

            if (sample_valid) begin
                pending_envelope <= sample_envelope;
                sample_seen      <= 1'b1;
            end

            if (carrier_tick) begin
                envelope       <= pending_envelope;
                envelope_valid <= sample_seen;
                sample_seen    <= 1'b0;
            end
        end
    end

endmodule
