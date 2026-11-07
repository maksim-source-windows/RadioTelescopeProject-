#include <Servo.h>

Servo servoH, servoV;
int currentH = 90, currentV = 45;

void setup() {
  Serial.begin(9600);
  servoH.attach(9);
  servoV.attach(10);
  servoH.write(currentH);
  servoV.write(currentV);
  delay(1000);
  Serial.println("READY");
}

void moveTo(int az, int alt) {
  // Безопасные ограничения
  az = constrain(az, 0, 180);
  alt = constrain(alt, 10, 80);

  // Если уже на этой позиции, ничего не делаем
  if (az == currentH && alt == currentV) {
    Serial.print("POS ");
    Serial.print(currentH);
    Serial.print(" ");
    Serial.println(currentV);
    return;
  }

  int steps = max(abs(az - currentH), abs(alt - currentV));

  for (int i = 0; i <= steps; i++) {
    int h = currentH + (az - currentH) * i / steps;
    int v = currentV + (alt - currentV) * i / steps;
    servoH.write(h);
    servoV.write(v);
    delay(20);
  }

  currentH = az;
  currentV = alt;

  Serial.print("POS ");
  Serial.print(currentH);
  Serial.print(" ");
  Serial.println(currentV);
}

void loop() {
  if (Serial.available()) {
    String command = Serial.readStringUntil('\n');
    command.trim();

    if (command.startsWith("SET ")) {
      int space1 = command.indexOf(' ');
      int space2 = command.indexOf(' ', space1 + 1);

      if (space1 != -1 && space2 != -1) {
        int az = command.substring(space1 + 1, space2).toInt();
        int alt = command.substring(space2 + 1).toInt();

        moveTo(az, alt);
      }
    }
  }
}
