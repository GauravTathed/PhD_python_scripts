@echo off
title ESP32 BME280 Serial Logger

echo Starting ESP32 serial logger...
echo.

cd /d "C:\Users\iamga\OneDrive\Documents\GitHub\PhD_python_scripts\HPT_logger"

python serial_to_influx.py

echo.
echo Logger stopped.
pause