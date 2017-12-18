
// This programm enables the fast reading the analogue input of
// an Arduino Due after a given delay of a trigger send to one
// of its digital I/O pins.
// The communication with a PC is done via the native USB port of an
// Arduino Due using the CmdMessenger v4.0 libary written by Neil Dudman,
// Dreamcat4, Thijs Elenbaas and Valeriy Kucherenko
// (https://github.com/thijse/Arduino-CmdMessenger).

#include <CmdMessenger.h>  // sending/receiving user-defined commands

// global variables
// ================
int input_res = 4096;
float input_max = 3.3;
float input_min = 0;

float read_delay = 38e-6;
bool record_intensity = true;
bool laser_is_on = false;

// get class instances for serial communication with PC
CmdMessenger cmd = CmdMessenger(SerialUSB);

// ==========================================================
// Initialisation of CmdMessenger and definition of commands
// ==========================================================

// definition of supported commands
enum
{
  check_arduino,
  ret_fast_intensity,
  sweep_delay,
  set_delay,
  get_delay,
  get_intensity_trace,
  oscilloscope_mode,
  set_laser_on,
  get_laser_on,
  ret_string,              // returns a string
  ret_int,                 // returns an integer
  ret_char,                // returns a character
  ret_float,               // returns a float value
  ret_double,              // returns a double value
  ret_bool,                // returns a boolean value
  ret_byte                 // returns a byte value
};

// callback methods for commands
void attachCommands()
{
  cmd.attach(UnknownCommand);
  cmd.attach(check_arduino, on_check_arduino);
  cmd.attach(sweep_delay, on_sweep_delay);
  cmd.attach(set_delay, on_set_delay);
  cmd.attach(get_delay, on_get_delay);
  cmd.attach(get_intensity_trace, on_get_intensity_trace);
  cmd.attach(oscilloscope_mode, on_oscilloscope_mode);
  cmd.attach(set_laser_on, on_set_laser_on);
  cmd.attach(get_laser_on, on_get_laser_on);
}

// ==============================
// CmdMessenger callback methods
// ==============================

void UnknownCommand() {
  cmd.sendCmd(ret_string, "Error: Command without attached callback!");
}

void on_check_arduino() {
  cmd.sendCmd(ret_string, "Arduino communication is fine.");
}

void on_set_delay() {  
  read_delay = cmd.readBinArg<float>();
}
void on_get_delay() {
  cmd.sendBinCmd(ret_float, read_delay);  
}

void on_sweep_delay() {  
  int data_points = cmd.readBinArg<int>();
  int num_pulses = cmd.readBinArg<int>();
  float delay_increment = cmd.readBinArg<float>();
  
  float pulse_off = 5e-3; // in seconds 
  float curr_sig1;
  float curr_delay;
  uint32_t RA;
  uint32_t RC;
  uint32_t clock_counts;

  for (int i = 0; i < data_points; i++) {
    curr_sig1 = 0;
    curr_delay = 0;    
    // set register A of channel 0 of Timer 1 to number of clock periods that correspond to current delay
    RA = (i+1)*delay_increment/23.8095e-9;
    TC_SetRA(TC1,0, RA); 
    // set register C to number of clock periods when (trigger) pulse should be off
    RC = pulse_off/23.8095e-9;
//    TC_SetRC(TC1,0, RC);

    for (int avg_no = 0; avg_no < num_pulses; avg_no++) {         
      curr_sig1 += (float(delayed_ADC_read(1)) / (input_res - 1)) * (input_max - input_min) / num_pulses;
      clock_counts = REG_TC1_CV0;         // read counter (was stopped after ADC conversion finished)
      curr_delay += clock_counts/num_pulses*23.8095e-9;           
//      while ((REG_TC1_SR0 & 0x10) == 0); // check TC status register for occurance of RC compare
      delay(4);  
    }
    // send (averaged) intensity to PC
    cmd.sendCmdStart(ret_fast_intensity);
    cmd.sendCmdBinArg(curr_sig1);
    cmd.sendCmdBinArg(curr_delay);
    cmd.sendCmdEnd();
  }
}

void on_get_intensity_trace() {

  int data_points = cmd.readBinArg<int>();
  int num_pulses = cmd.readBinArg<int>();
  int points_per_pulse = cmd.readBinArg<int>();
  float delay_increment = cmd.readBinArg<float>();  // in seconds
  
  float curr_sig1;
  float curr_delay;
  uint32_t RA;
  uint32_t RC;
  uint32_t clock_counts;

  for (int i = 0; i < data_points; i++) {
    curr_sig1 = 0;
    curr_delay = i*delay_increment;
    
    // set register A of channel 0 of Timer 1 to number of clock periods that correspond to current delay
    RA = read_delay/23.8095e-9;
    TC_SetRA(TC1,0, RA); 
    // set register C to number of clock periods when (trigger) pulse should be off
//    RC = pulse_off/23.8095e-9;
//    TC_SetRC(TC1,0, RC);

    for (int avg_no = 0; avg_no < num_pulses; avg_no++) {         
      curr_sig1 += (float(delayed_ADC_read(points_per_pulse)) / (input_res - 1)) * (input_max - input_min) / num_pulses;
//      clock_counts = REG_TC1_CV0;         // read counter (was stopped after ADC conversion finished)
//      curr_delay += clock_counts/avg_points*23.8095e-9;           
//      while ((REG_TC1_SR0 & 0x10) == 0); // check TC status register for occurance of RC compare
      delayMicroseconds(100);  
    }    
    // send (averaged) intensity to PC
    cmd.sendCmdStart(ret_fast_intensity);
    cmd.sendCmdBinArg(curr_sig1);
    cmd.sendCmdBinArg(curr_delay);
    cmd.sendCmdEnd();
    delayMicroseconds(delay_increment*1e6);    // wait until next intensity measurement
  } 
}

int delayed_ADC_read(int num_points) {

  int intensity1_level=0;
  uint32_t counter_status; 

  REG_PIOB_SODR |= (1<<25);           // emit trigger by set pin 2 HIGH (port B, pin 25)     
  counter_status = REG_TC1_SR0;       // clear status registers by reading them   
  counter_status = REG_TC1_CV0;        
  
  REG_TC1_CCR0 = 0x5;                 // reset and enable counter of channel 0 of timer 1      
  while ((REG_TC1_SR0 & 0x4) == 0);   // wait until RA compare by checking TC status register
  for (int i = 0; i < num_points; i++) {		//MESSPUNKTE PRO PULS = 20
    ADC->ADC_CR = 0x2;                  // start ADC conversion   
    while ((ADC->ADC_ISR & 0x80) == 0); // wait until AD conversion of channel 7 finished by checking Interrupts Status Register (ADC_ISC)
    intensity1_level += ADC->ADC_CDR[7]; //get values from Channel Data Register (ADC_CDR) of channel 7 
  } 
  REG_TC1_CCR0 = 0x2;                 // disable counter
  REG_PIOB_CODR |= (1<<25);           // set pin 2 LOW (port B, pin 25)
  return intensity1_level;    
}


/*
   reads the ADCs as fast as possible with the
   specified number of data points
*/
void on_oscilloscope_mode() {
  
  int data_points = cmd.readBinArg<int>();
  uint32_t clock_counts[data_points];
  float signal1[data_points];
  float time_delay;

  REG_TC1_CCR0 = 0x5;                     // reset and enable counter of channel 0 of timer 1   
  for (int i = 0; i < data_points; i++) {
    ADC->ADC_CR |= (1<<1);                // start ADC conversion    
    while ((ADC->ADC_ISR & 0x80) == 0);   // check Interrupts Status Register (ADC_ISC) for channel 7 (A0) and wait for conversion
    signal1[i] = ADC->ADC_CDR[7];         //get values from Channel Data Register (ADC_CDR) of channel 7
    clock_counts[i] = REG_TC1_CV0;        // read timer counter 
  }

  for (int i = 0; i < data_points; i++) {
    time_delay = clock_counts[i]*23.8095e-9;
    cmd.sendCmdStart(ret_fast_intensity);  
    cmd.sendCmdBinArg(signal1[i] / (input_res - 1) * (input_max - input_min));
    cmd.sendCmdBinArg(time_delay);
    cmd.sendCmdEnd();
  }  
}

void on_set_laser_on() {
  laser_is_on = cmd.readBinArg<bool>();
}
void on_get_laser_on() {
  cmd.sendBinCmd(set_laser_on, laser_is_on);
}



void const_laser_trigger() {  
//  digitalWrite(2,HIGH);
  REG_PIOB_SODR |= (1<<25);           // emit trigger by set pin 2 HIGH (port B, pin 25)
  delayMicroseconds(70);
//  digitalWrite(2,LOW);
  REG_PIOB_CODR |= (1<<25);           // set pin 2 LOW (port B, pin 25)  
//  cmd.sendCmd(ret_string, "laser trigger");
}




void setup()
{
  // setup pins
  pinMode(A0, INPUT); // ADC
  pinMode(2, OUTPUT); // trigger pin

  // set up Timer 1, channel 0 (ID_TC3) for delayed ADC measurement
  pmc_enable_periph_clk(TC3_IRQn); 
  REG_TC1_CMR0= 0x00008040; // set channel mode register (use waveform mode, MCL/2 clock (42 MHz), no effects on triggers, TIOA, and TIOB, timer stopped when it reaches RC)
  REG_TC1_IDR0= 0xFF; // disable all interrupts

  // set the ADC registers
  ADC->ADC_CHER = 0x80;       // enable ADC on pin A0 (corresponds to Ch7 on Atmel chip)
  //ADC->ADC_MR = 0x90330180;
  ADC->ADC_MR = 0x10000180;   // set the ADC Mode Register (ADC_MR)
  // default value: 0x9038ff00, binary: 0000 0000 1001 0000 0011 1000 1111 1111 0000 0000  
  // structure of register:
  //
  // 0000 0000    0  0    01      0000     0   0    00     0000  0000 0001   1     000 0000
  //        -- |USEQ|-|transfer|tracking|ANACH|-|settling|startup|prescale|freerun|..  
  //
  // this uses the following prameters:
  // 
  // * use free run mode, that is continuouse AD conversion without trigger
  // * prescale for ADC clock: 1 => 84 MHz/2/(1+1) = 21 MHz -> ADC clock period: 47.6 ns
  // * startup time: 0 ADC clock periods (time befor starting first conversion, only relevant when sleep mode enabled)
  // * settling time: 3 ADC clock periods (higher values only necessay when changing offset or gain between conversion of different channels)
  // * ANACH: 0 (no differnt analog settings for different ADC channels allowed)
  // * tracking time: 0 (set bits don't have any effect because of a hardware defect, tracking time is always 15 ACD clock periods, which is the time for one AD conversion)
  // * transfer time: 1 -> corresponds 5 to ADC clock periods, determines time between change of ADC channels and the corresponds to the time it takes to transfer the analog input from the sample to hold
  // * USEQ: 0 -> no user-defined sequence for ADC channels because just one channel 
  
  SerialUSB.begin(115200); // initialise USB port, actual baud rate irrelevant as always USB 2.0 speed
  attachCommands();       // attach CmdMessenger commands
  laser_is_on = false;
}

void loop()
{
  // process incoming serial data from PC and perform callbacks
  cmd.feedinSerialData();
  if(laser_is_on) {
    const_laser_trigger();
  }
}


