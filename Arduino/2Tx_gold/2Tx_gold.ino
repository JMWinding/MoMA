#include <Ticker.h>
#include <OneWire.h>
#include <SPI.h>
#include <SD.h>

#define debug false

// **** File Parameters ****
File recordFile;
char TxRecordFile[25];
File settingFile;
char TxSettingFile[25];
const int fileIndexStart = 36;
const int fileIndexEnd = 59;

// **** Protocol Parameters ****
const unsigned int prelen = 4;
const unsigned int ChipInterval = 150;
const unsigned int PulseOnTime = 45;
const unsigned int PulseOffTime = ChipInterval - PulseOnTime;
const unsigned int DisposedTime = 10000;

// **** Pump Parameters ****
const int nTx = 2;
const int motors[nTx] = {2, 3};
const int motorOffsets[nTx] = {20, 20+14* (1+prelen)};
const float motorCon[nTx] = {10, 10};
const unsigned int TotalSymbols = 100;

// **** CDMA Parameters ****
int cdmanum;
int goldlen;
int cdmalen;
int cdmacode[3][7];
int cdmaidx[nTx];
int pchip[nTx];
int pind[nTx];
int xbit[nTx];
int xind[nTx];
int totalChips;

// **** random ****
int randlist[7];
int motorOffsets2[nTx];

// **** TxRx Sync for Record ****
const int chipSelect = 53;
/* Arduino Uno
 * SD card attached to SPI bus as follows:
 ** CS   - pin 10 (Chip Select)
 ** MOSI - pin 11
 ** MISO - pin 12
 ** CLK  - pin 13
*/

// **** Tickers and Functions ****
void tx_transmit (int pump, int cbit);
void tx_transmit_over (int pump);
void txi_transmit (int i);

// pump reset function after each transmitted chip
void tx0_transmit_over (void) { tx_transmit_over (0); }
Ticker timerTx0Off (tx0_transmit_over, PulseOnTime, 1, MILLIS);
// periodic function to transmit chips
void tx0_transmit (void) { txi_transmit (0); timerTx0Off.start (); }
Ticker timerTx0On (tx0_transmit, ChipInterval, 0, MILLIS);
// function to control the start of the packet
int tx0OnInitCounter = 0;
void tx0OnInit (void) { 
  if (--tx0OnInitCounter == 0) {
    pchip[0] = 0; pind[0] = 0; xbit[0] = 1; xind[0] = 0; timerTx0On.start ();
  }
}
Ticker timerTx0OnInit (tx0OnInit, ChipInterval, 0, MILLIS);

// tx1
void tx1_transmit_over (void) { tx_transmit_over (1); }
Ticker timerTx1Off (tx1_transmit_over, PulseOnTime, 1, MILLIS);
void tx1_transmit (void) { txi_transmit (1); timerTx1Off.start (); }
Ticker timerTx1On (tx1_transmit, ChipInterval, 0, MILLIS);
int tx1OnInitCounter = 0;
void tx1OnInit (void) { 
  if (--tx1OnInitCounter == 0) {
    pchip[1] = 0; pind[1] = 0; xbit[1] = 1; xind[1] = 0; timerTx1On.start ();
  }
}
Ticker timerTx1OnInit (tx1OnInit, ChipInterval, 0, MILLIS);

// **** Setup ****
void setup () 
{
  // initialize serial communication with computer:
  Serial.begin (115200);
  Serial1.begin (115200);
  
  // setup pins
  for (int i = 0; i < nTx; i++)
  {
    pinMode (motors[i], OUTPUT);
  }

  // Initialize SD Card
  Serial.println ("Init SD");
  if (!SD.begin (chipSelect)) 
  {
    Serial.println ("Init fail");
    exit (0);
  }
  Serial.println ("Init done");

  // load CDMA settings
  sprintf (TxSettingFile, "gold.txt");
  if (!SD.exists (TxSettingFile))
  {
    Serial.println ("Setting not exist");
    Serial.println ("");
    exit (0);
  }
  Serial.print ("Setting file Name: ");
  Serial.println (TxSettingFile);
  settingFile = SD.open (TxSettingFile, FILE_READ);
  if (!settingFile)
  {
    Serial.println ("Setting open fail");
    Serial.println ("");
    exit (0);
  }
  cdmanum = settingFile.parseInt ();
  goldlen = settingFile.parseInt ();
  cdmalen = goldlen;
  totalChips = (prelen + TotalSymbols) * cdmalen;
  #if debug
  Serial.println (goldlen);
  #endif
  for (int i = 0; i < cdmanum; i++)
  {
    for (int j = 0; j < goldlen; j++)
    {
      cdmacode[i][j] = settingFile.parseInt ();
      #if debug
      Serial.print (cdmacode[i][j]);
      #endif
    }
    #if debug
    Serial.println ("");
    #endif
  }
  settingFile.close ();
  Serial.println ("Setting loaded");
  delay (500);

  // wait for rx to set up
  Serial.println ("Waiting for Rx to set up...");
  while (!Serial1.available ())
  {
    delay (1000);
    Serial.println ("Waiting...");
  }
  delay (200);
  if (Serial1.read () != 'a')
  {
    Serial.println ("Rx sent a wrong message");
    exit (0);
  }
  Serial.println ("Rx set up confirmed");
  
  // prepare the pumps
  Serial.println ("Set up");
  Serial.println ("Pipe in water");
  turn_off_all ();
  delay (DisposedTime);
  Serial.println ("Fill the tubes");
  turn_on_all ();
  delay (DisposedTime);
  Serial.println ("Clear the tubes");
  turn_off_all ();
  delay (DisposedTime);
  Serial.println ("Set up done");
  Serial.println ("");
  Serial.println ("");
  
  randomSeed (analogRead (0));
}

void loop () {
  delay (5000);

  for (int fileIndex = fileIndexStart; fileIndex <= fileIndexEnd; fileIndex++)
  {
    // check if write file exists not
    sprintf (TxRecordFile, "%02d.txt", fileIndex);
    if (SD.exists (TxRecordFile)) 
    {
      SD.remove (TxRecordFile);
      #if debug
      Serial.println ("File removed");
      #endif
    }
    
    // randomness
    for (int i = 0; i < cdmanum; i++)
    {
      randlist[i] = i;
    }
    for (int i = 0; i < nTx; i++)
    {
      int temp = random (cdmanum-i);
      cdmaidx[i] = randlist[temp];
      randlist[temp] = randlist[cdmanum-i-1];
      motorOffsets2[i] = random (cdmalen);
    }
    
    // create record file
    Serial.print ("File Name: ");
    Serial.println (TxRecordFile);
    recordFile = SD.open (TxRecordFile, FILE_WRITE);
    if (!recordFile)
    {
      Serial.println ("Record open fail");
      Serial.println ("");
      exit (0);
    }
    
    // start collecting data
    recordFile.print (ChipInterval);
    recordFile.print (" ");
    recordFile.print (TotalSymbols);
    recordFile.print (" ");
    recordFile.print (nTx);
    recordFile.print (" ");
    recordFile.print (goldlen);
    recordFile.print (" ");
    recordFile.print (prelen);
    recordFile.println ("");

    for (int i = 0; i < nTx; i++)
    {
      #if debug
      Serial.print ("(");
      Serial.print (motors[i]);
      Serial.print (", ");
      Serial.print (motorOffsets[i]);
      Serial.print (", ");
      Serial.print (motorOffsets2[i]);
      Serial.print (", ");
      Serial.print (motorCon[i]);
      Serial.print ("), ");

      Serial.print (" (code ");
      Serial.print (cdmaidx[i]);
      Serial.println (")");
      #endif
      
      recordFile.print ("(");
      recordFile.print (motors[i]);
      recordFile.print (", ");
      recordFile.print (motorOffsets[i]);
      recordFile.print (", ");
      recordFile.print (motorOffsets2[i]);
      recordFile.print (", ");
      recordFile.print (motorCon[i]);
      recordFile.print ("), ");

      recordFile.print (" (code ");
      recordFile.print (cdmaidx[i]);
      recordFile.println (")");
    }
    
    Serial.println ("Random bits");

    Serial.println ("START");
    recordFile.println ("START");

    // set counter
    tx0OnInitCounter = motorOffsets[0]+motorOffsets2[0];
    tx1OnInitCounter = motorOffsets[1]+motorOffsets2[1];

    // Notify Rx to start
    Serial1.write ('a');
    delay (200);

    timerTx0OnInit.start ();
    timerTx1OnInit.start ();
    while (1)
    {
      timerTx0OnInit.update ();
      if (timerTx0On.state () == status_t::RUNNING) { timerTx0OnInit.stop (); }
      timerTx0On.update ();
      if (timerTx0On.counter () >= totalChips) { timerTx0On.stop (); }
      timerTx0Off.update ();
      timerTx1OnInit.update ();
      if (timerTx1On.state () == status_t::RUNNING) { timerTx1OnInit.stop (); }
      timerTx1On.update ();
      if (timerTx1On.counter () >= totalChips) { timerTx1On.stop (); }
      timerTx1Off.update ();

      if (timerTx0OnInit.state () != status_t::STOPPED) { continue; }
      if (timerTx0On.state () != status_t::STOPPED) { continue; }
      if (timerTx0Off.state () != status_t::STOPPED) { continue; }
      if (timerTx1OnInit.state () != status_t::STOPPED) { continue; }
      if (timerTx1On.state () != status_t::STOPPED) { continue; }
      if (timerTx1Off.state () != status_t::STOPPED) { continue; }

      break;
    }

    delay (DisposedTime);

    // Notify Rx to end
    Serial1.write ('z');
    delay (200);
    Serial.println ("END");
    Serial.println ("");
    Serial.println ("");
    recordFile.close ();

  }
}

// Preparing pumps actions
void turn_off_all (void)
{
  for (int i = 0; i < nTx; i++)
  {
    digitalWrite (motors[i], LOW);
  }
}

void turn_on_all (void)
{
  for (int i = 0; i < nTx; i++)
  {
    digitalWrite (motors[i], HIGH);
  }
}

// General pump actions
void tx_transmit (int pump, int cbit)
{
  #if debug
  Serial.print (millis ());
  Serial.print (", ");
  Serial.print (pump);
  Serial.print (", ");
  Serial.println (cbit);
  #endif

  recordFile.print (millis ());
  recordFile.print (", ");
  recordFile.print (pump);
  recordFile.print (", ");
  recordFile.println (cbit);

  if (cbit)
  {
    digitalWrite (motors[pump], HIGH);
  }
}

void tx_transmit_over (int pump)
{
  digitalWrite (motors[pump], LOW);
}

void txi_transmit (int i)
{
  int tempind;
  int cidx; 
  int midx;
  if (pchip[i] == goldlen)
  {
    cidx = xind[i];
    midx = 0;
    tempind = (xbit[i] != midx) == cdmacode[cdmaidx[i]][cidx];

    if (xind[i] == 0)
    {
      #if debug
      Serial.print (millis ());
      Serial.print (", ");
      Serial.print (i);
      Serial.print (", bit ");
      Serial.println (xbit[i]);
      #endif
      
      recordFile.print (millis ());
      recordFile.print (", ");
      recordFile.print (i);
      recordFile.print (", bit ");
      recordFile.println (xbit[i]);
    }
    
    tx_transmit (i, tempind);
    
    xind[i]++;
    if (xind[i] == cdmalen)
    {
      xind[i] = 0;
      xbit[i] = random (2);
    }
  }
  else
  {
    cidx = pchip[i];
    
    if (pind[i] < prelen)
    {
      tempind = cdmacode[cdmaidx[i]][cidx];
    }
    else
    {
      tempind = 1 - cdmacode[cdmaidx[i]][cidx];
    }
      
    if (pind[i] == 0)
    {
      #if debug
      Serial.print (millis ());
      Serial.print (", ");
      Serial.print (i);
      Serial.print (", prem ");
      Serial.println (tempind);
      #endif
    
      recordFile.print (millis ());
      recordFile.print (", ");
      recordFile.print (i);
      recordFile.print (", prem ");
      recordFile.println (tempind);
    }
    
    tx_transmit (i, tempind);
    
    pind[i]++;
    if (pind[i] == prelen)
    {
      pind[i] = 0;
      pchip[i]++;
    }
  }
}
