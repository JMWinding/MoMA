#include <Ticker.h>
#include <OneWire.h>
#include <SPI.h>
#include <SD.h>

// sudo chmod a+rw /dev/ttyACM0

//****************Setup Parameters************
const int Ts = 5; //Gap between sampling in ms
//********************************************

//*****************Setup Connections**********
int record = LOW;
const int chipSelect = 53;
const byte ECsensorPin = A1;  // EC Meter analog output,pin on analog 1
//********************************************

//*******************Variables****************
File myFile;

int fileIndex = 6;
char fileName[50];
int activated = 0;

unsigned int AnalogReading;
unsigned long AnalogSampleTime;

//*********************************************
void rx_receive (void);
Ticker timerRx (rx_receive, Ts, 0, MILLIS);

void setup ()
{  
  while (activated == 0)
  {
    Serial.begin (115200);
    Serial1.begin (115200);

    Serial.print ("\n\n\n\n\n");
    Serial.print ("Initializing SD card...");
    if (!SD.begin (chipSelect)) 
    {
      Serial.println ("initialization failed!");
    }
    else
    {
      Serial.println ("initialization done.");
      activated = 1;
    }
  }

  AnalogReading = 0;
  AnalogSampleTime = millis ();
  AnalogReading = analogRead (ECsensorPin);
  Serial1.write ('a');
}

void loop ()
{
  Serial1.write ('a');
  if (Serial1.available ())
  {
    switch (Serial1.read ())
    {
      case 'a':
        record = HIGH;
        break;
      case 'z':
        record = LOW;
        break;
      default:
        Serial.println ("Tx sent a wrong message");
        exit (0);
    }
  }

  switch (record)
  {
    case HIGH:
      if (!myFile)
      {
        Serial.print ("Opening File Number: ");
        Serial.println (fileIndex);
        sprintf (fileName, "%02d.txt", fileIndex);
        if (SD.exists (fileName)) {SD.remove (fileName);}
        myFile = SD.open (fileName, FILE_WRITE);
        if (myFile) {Serial.println ("File Created");}
        timerRx.start ();
      }
      if (myFile)
      {
        timerRx.update ();
      }
      break;
      
    case LOW:
      if (myFile)
      { 
        timerRx.stop ();
        Serial.print ("Closing File Number: ");
        Serial.println (fileIndex);
        Serial.print (AnalogSampleTime);
        Serial.print ('\t');
        Serial.println (AnalogReading);
        Serial.println ("");
        myFile.println ("END");
        myFile.close ();
        fileIndex++;
      }
      break;
  }
}

void rx_receive (void)
{
  AnalogSampleTime = millis ();
  AnalogReading = analogRead (ECsensorPin);
  myFile.print (AnalogSampleTime);
  myFile.print ('\t');
  myFile.println (AnalogReading);
  // Serial.println (AnalogReading);
}
