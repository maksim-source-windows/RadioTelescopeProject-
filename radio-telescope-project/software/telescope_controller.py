"""
УПРАВЛЕНИЕ ТЕЛЕСКОПОМ (Arduino, протокол: "SET <az> <alt>" -> "POS <az> <alt>")
Требования: pip install pyserial
"""

import serial
import time
import sys

SERIAL_PORT = 'COM3'
BAUD_RATE = 9600
AZ_LIMITS = (0, 180)    # обычный серво: только 180°
ALT_LIMITS = (10, 80)
STEP_DELAY_S = 0.02     # delay(20) на каждый градус в прошивке


def move_timeout(old, new):
    """Время движения в прошивке (макс. по осям) + запас"""
    steps = max(abs(new[0] - old[0]), abs(new[1] - old[1]))
    return steps * STEP_DELAY_S + 1.5


class TelescopeController:
    def __init__(self, port, baud=9600):
        self.pos = (90, 45)  # стартовая позиция прошивки
        try:
            self.ser = serial.Serial(port, baud, timeout=0.5)
            time.sleep(2)  # перезагрузка Arduino при открытии порта
            print(f"Подключено к {port}")
        except Exception as e:
            print(f"Ошибка подключения: {e}")
            self.ser = None

    def set_position(self, az, alt):
        """Навести и дождаться 'POS'. True только при подтверждении."""
        if not self.ser:
            return False

        az = int(max(AZ_LIMITS[0], min(AZ_LIMITS[1], az)))
        alt = int(max(ALT_LIMITS[0], min(ALT_LIMITS[1], alt)))

        self.ser.reset_input_buffer()  # выбросить READY и старые ответы
        self.ser.write(f"SET {az} {alt}\n".encode())

        deadline = time.time() + move_timeout(self.pos, (az, alt))
        while time.time() < deadline:
            line = self.ser.readline().decode(errors='ignore').strip()  # timeout 0.5 с
            if line.startswith("POS"):
                print(f"Позиция: {line}")
                self.pos = (az, alt)
                return True
            if line:
                print(f"Ответ: {line}")

        print("Таймаут ожидания POS от Arduino")
        return False

    def close(self):
        if self.ser:
            self.ser.close()


if __name__ == "__main__":
    telescope = TelescopeController(SERIAL_PORT, BAUD_RATE)

    if not telescope.ser:
        print("Не удалось подключиться к Arduino")
        sys.exit(1)

    print("\n=== ТЕСТОВАЯ ПОСЛЕДОВАТЕЛЬНОСТЬ ===")
    for az, alt in [(0, 45), (90, 45), (180, 45)]:
        print(f"Наведение: Az={az}°, Alt={alt}°")
        telescope.set_position(az, alt)

    telescope.close()
    print("\nГотово!")
