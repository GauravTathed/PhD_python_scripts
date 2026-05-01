import serial
import requests
import time

# -----------------------------
# Serial settings
# -----------------------------
SERIAL_PORT = "COM14"     # Change this to your ESP32 port
BAUD_RATE = 115200

# -----------------------------
# InfluxDB settings
# -----------------------------
INFLUX_URL = "http://localhost:8086/api/v2/write"
INFLUX_ORG = "home"
INFLUX_BUCKET = "weather"
INFLUX_TOKEN = "sensor-token-12345"

DEVICE_NAME = "sparkfun_thing_plus_esp32c6_serial"

def send_to_influx(temperature_c, humidity_percent, pressure_hpa):
    url = (
        f"{INFLUX_URL}"
        f"?org={INFLUX_ORG}"
        f"&bucket={INFLUX_BUCKET}"
        f"&precision=s"
    )

    line_protocol = (
        f"bme280,device={DEVICE_NAME} "
        f"temperature_c={temperature_c},"
        f"humidity_percent={humidity_percent},"
        f"pressure_hpa={pressure_hpa}"
    )

    headers = {
        "Authorization": f"Token {INFLUX_TOKEN}",
        "Content-Type": "text/plain; charset=utf-8",
    }

    response = requests.post(url, headers=headers, data=line_protocol, timeout=5)

    if response.status_code == 204:
        print("Sent to InfluxDB:", line_protocol)
    else:
        print("InfluxDB error:", response.status_code, response.text)

def main():
    print(f"Opening serial port {SERIAL_PORT} at {BAUD_RATE} baud...")

    with serial.Serial(SERIAL_PORT, BAUD_RATE, timeout=2) as ser:
        time.sleep(2)
        print("Reading serial data...")

        while True:
            raw_line = ser.readline().decode("utf-8", errors="ignore").strip()

            if not raw_line:
                continue

            print("Serial:", raw_line)

            # Skip text lines and header lines
            if raw_line.startswith("ESP32"):
                continue
            if raw_line.startswith("Trying"):
                continue
            if raw_line.startswith("BME280"):
                continue
            if raw_line.startswith("temperature_c"):
                continue
            if raw_line.startswith("Stopping"):
                continue

            parts = raw_line.split(",")

            if len(parts) != 3:
                print("Skipping non-data line")
                continue

            try:
                temperature_c = float(parts[0])
                humidity_percent = float(parts[1])
                pressure_hpa = float(parts[2])
            except ValueError:
                print("Skipping invalid data line")
                continue

            send_to_influx(
                temperature_c=temperature_c,
                humidity_percent=humidity_percent,
                pressure_hpa=pressure_hpa,
            )

if __name__ == "__main__":
    main()