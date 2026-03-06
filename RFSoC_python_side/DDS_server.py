
def hz_to_pinc48(freq_hz: float, fs_hz: float) -> int:
    f = ((float(freq_hz) + fs_hz/2.0) % fs_hz) - fs_hz/2.0
    return int(round((f / fs_hz) * (1 << 48))) & ((1 << 48) - 1)

def rad_to_poff48(rad: float) -> int:
    return int(round(((float(rad) % (2*np.pi)) / (2*np.pi)) * (1 << 48))) & ((1 << 48) - 1)

def deg_to_poff48(rad: float) -> int:
    return rad_to_poff48(rad)

def us_to_ticks(us: float, aclk_hz: float) -> int:
    return int(round(float(us) * 1e-6 * float(aclk_hz)))

def _u256_to_bytes_le(x: int) -> bytes:
    return int(x & ((1 << 256) - 1)).to_bytes(32, "little")

def pack_table_header(row_count: int) -> bytes:
    x  = (0xA6 & 0xFF) << 248
    x |= (int(row_count) & 0xFFFF) << 216
    return _u256_to_bytes_le(x)

def pack_row_header(hop_count: int, jump_flag: int) -> bytes:
    x  = (0xA5 & 0xFF) << 248
    x |= (int(hop_count) & 0xFFFF) << 216
    x |= (1 if jump_flag else 0) << 215
    return _u256_to_bytes_le(x)

def pack_hop(dwell_ticks: int, pinc48: int, poff48: int, resync: int) -> bytes:
    x  = (int(dwell_ticks) & ((1 << 32) - 1)) << 224
    x |= (int(pinc48)      & ((1 << 48) - 1)) << 144
    x |= (int(poff48)      & ((1 << 48) - 1)) << 96
    x |= (1 if resync else 0) << 95
    return _u256_to_bytes_le(x)

def build_dma_buffer_little_endian(descs):
    blob = b"".join(descs)
    words = np.frombuffer(blob, dtype="<u8")
    buf = allocate(shape=(len(words),), dtype=np.uint64)
    buf[:] = words
    return buf



def _rf_mhz_to_baseband_hz(freqs_mhz, center_mhz):
    out = []
    for f_mhz in freqs_mhz:
        if not (MIN_RF_MHZ <= f_mhz <= MAX_RF_MHZ):
            raise ValueError(f"RF {f_mhz} MHz out of range [{MIN_RF_MHZ}, {MAX_RF_MHZ}]")
        bb = (float(f_mhz) - float(center_mhz)) * 1e6
        if abs(bb) > 175e6 + 1:
            raise ValueError(f"|RF-center|={abs(bb)/1e6:.3f} MHz exceeds ±175 MHz")
        out.append(bb)
    return out

def build_table_descriptors_rf_mhz(rows_mhz, aclk_hz, fs_dac_hz, center_mhz=630.0):
    descs = []
    descs.append(pack_table_header(len(rows_mhz)))
    for _, freqs_mhz, phases, times_us, resyncs, jump_flag in rows_mhz:
        freqs_hz = _rf_mhz_to_baseband_hz(freqs_mhz, center_mhz)
        n = len(freqs_hz)
        descs.append(pack_row_header(n, jump_flag))
        phase_corr = 0.0
        sum_time = 0
        for i in range(n):
            
            dwell = max(1, us_to_ticks(times_us[i], aclk_hz))-5
            pinc  = hz_to_pinc48(freqs_hz[i]-160, fs_dac_hz)
            poff  = deg_to_poff48(phases[i] + phase_corr)
            descs.append(pack_hop(dwell, pinc, poff, int(resyncs[i])))
            if i < n-1:
                time_clk_cyl = (dwell+5)*1/(409.600097)
                # phase_corr = 2*np.pi*freqs_hz[i+1]*1e-6*(sum(times_us[:i+1]))
                print('sum us:',sum(times_us[:i+1]))
                sum_time += time_clk_cyl
                phase_corr = 2*math.pi*freqs_hz[i+1]*1e-6*(sum_time)
                print('sum clk cyl:',sum_time)

                print(phase_corr)
    return descs

def upload_table_rf_mhz(ol, rows_mhz, aclk_hz=409_600_097.0, fs_dac_hz=409_600_097.0, dma_name="axi_dma_0", center_mhz=630.0):
    dma = getattr(ol, dma_name)
    descs = build_table_descriptors_rf_mhz(rows_mhz, aclk_hz, fs_dac_hz, center_mhz)
    buf = build_dma_buffer_little_endian(descs)
    try:
        if not dma.sendchannel.running:
            dma.sendchannel.start()
        dma.sendchannel.transfer(buf)
        dma.sendchannel.wait()
    finally:
        buf.freebuffer()
    return len(descs)


def start_rows_server(host="0.0.0.0", port=8000):
    class RowsHandler(BaseHTTPRequestHandler):
        def do_POST(self):
            if self.path != "/upload_rows":
                self.send_response(404)
                self.end_headers()
                return

            length = int(self.headers.get("Content-Length", "0"))
            body = self.rfile.read(length)
            try:
                data = json.loads(body.decode("utf-8"))
            except Exception as e:
                self.send_response(400)
                self.send_header("Content-Type", "application/json")
                self.end_headers()
                self.wfile.write(json.dumps({"status": "error", "message": "invalid JSON"}).encode("utf-8"))
                return

            rows = data.get("rows")
            print(rows)
            if rows is None:
                self.send_response(400)
                self.send_header("Content-Type", "application/json")
                self.end_headers()
                self.wfile.write(json.dumps({"status": "error", "message": "missing 'rows'"}).encode("utf-8"))
                return

            aclk_hz = data.get("aclk_hz", 409_600_097.0)
            fs_dac_hz = data.get("fs_dac_hz", 409_600_097.0)
            center_mhz = data.get("center_mhz", 630.0)
            dma_name = data.get("dma_name", "axi_dma_0")
            
            try:
                count = upload_table_rf_mhz(
                    ol,
                    rows,
                    aclk_hz=aclk_hz,
                    fs_dac_hz=fs_dac_hz,
                    dma_name=dma_name,
                    center_mhz=center_mhz,
                )
                resp = {"status": "ok", "descriptors": count}
                self.send_response(200)
            except Exception as e:
                resp = {"status": "error", "message": str(e)}
                self.send_response(500)

            self.send_header("Content-Type", "application/json")
            self.end_headers()
            self.wfile.write(json.dumps(resp).encode("utf-8"))

        def log_message(self, format, *args):
            return

    server = HTTPServer((host, port), RowsHandler)
    print(f"Serving on {host}:{port}")
    server.serve_forever()

if __name__ == "__main__":
    import numpy as np
    from pynq import allocate
    import math
    from pynq import Overlay
    import xrfclk
    MIN_RF_MHZ = 455.0
    MAX_RF_MHZ = 805.0
    import json
    from http.server import BaseHTTPRequestHandler, HTTPServer

    # Load your bitstream & clocks (ZCU111 example)
    print("overlaying the bitstream")
    ol = Overlay("DDS_630BB.xsa", download=True)
    print("Bitstream uploaded")
    print("setting the clock frequencies")
    xrfclk.set_ref_clks(lmk_freq=122.88, lmx_freq=409.6)
    print("clock initialized")
    print("starting the server")
    start_rows_server(host="0.0.0.0", port=9009)
    print("server up and running")
    