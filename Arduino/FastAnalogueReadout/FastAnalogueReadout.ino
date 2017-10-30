
// This programm enables the fast reading the analogue input of
// an Arduino Due after a given delay of a trigger send to one
// of its digital I/O pins.
// The communication with a PC is done via the native USB port of an
// Arduino Due using the CmdMessenger v4.0 libary written by Neil Dudman,
// Dreamcat4, Thijs Elenbaas and Valeriy Kucherenko
// (https://github.com/thijse/Arduino-CmdMessenger).
// The continuouse readout of the analogue-to-digital convertes (ADCs)
// is done by setting a timer interrupt using Ivan Seidel's
// DueTimer libary (https://github.com/ivanseidel/DueTimer)

#include <CmdMessenger.h>  // sending/receiving user-defined commands
#include <DueTimer.h>      // set user-defined timer and timer interrups

// pin definition
// ==============
const int input1_pin   = A0;
const int input2_pin   = A1;
const int trigger_pin = 2;

// global variables
// ================
float intensity1;
float intensity2;
int input_res = 4096;
float input_max = 3.3;
float input_min = 0;

float read_delay = 40e-6;
bool record_intensity = true;
bool cont_signal_readout = true;
int adc_readout_freq = 1e4;
int signal_avg = 10;     // number data point to average
int signal_avg_no = 0;
int sig1_level;
int sig2_level;
float signal1_temp = 0;
float signal2_temp = 0;

// get class instances for serial communication with PC
CmdMessenger cmd = CmdMessenger(SerialUSB, xmodem);

// ==========================================================
// Initialisation of CmdMessenger and definition of commands
// ==========================================================

// definition of supported commands
enum
{
  check_arduino,
  ret_fast_intensity,
  set_cont_signal_readout,
  get_cont_signal_readout,
  sweep_delay,
  set_delay,
  get_delay,
  fast_ADC_read,
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
  cmd.attach(set_cont_signal_readout, on_set_cont_signal_readout);
  cmd.attach(get_cont_signal_readout, on_get_cont_signal_readout);
  cmd.attach(sweep_delay, on_sweep_delay);
  cmd.attach(fast_ADC_read, on_fast_ADC_read);
  cmd.attach(set_delay, on_set_delay);
  cmd.attach(get_delay, on_get_delay);
}

// ==============================
// CmdMessenger callback methods
// ==============================

void UnknownCommand() {
  cmd.sendCmd(ret_string, "Error: Command without attached callback!");
}

void on_check_arduino() {
  cmd.sendCmd(ret_string, "Arduino communication fine.");
}


void on_set_cont_signal_readout() {
  cont_signal_readout = cmd.readBinArg<bool>();
}

void on_get_cont_signal_readout() {
  cmd.sendBinCmd(ret_bool, cont_signal_readout);
}

void on_sweep_delay() {
  
  int data_points = cmd.readBinArg<int>();
  int avg_points = cmd.readBinArg<int>();
  float delay_increment = cmd.readBinArg<float>();
  float pulse_off = 5e-3; // in seconds 
  int intensity1_level;
  int intensity2_level=0;
  float curr_sig1;
  float curr_sig2;
  float curr_delay=0;
  uint32_t counter_status; 
  uint32_t RA;
  uint32_t RC;
  uint32_t clock_counts;
  float clock_counts_avg;

  // set up Timer1, channel 0 (ID_TC3)
  pmc_enable_periph_clk(TC3_IRQn); 
  REG_TC1_CMR0= 0x00008040; // set channel mode register (use waveform mode, MCL/2 clock (42 MHz), no effects on triggers, TIOA, and TIOB, timer stopped when it reaches RC)
  REG_TC1_IDR0= 0xFF; // disable all interrupts

  for (int i = 0; i < data_points; i++) {

    curr_sig1 = 0;
    curr_sig2 = 0;
    curr_delay = 0;
    //    curr_delay = (i + 1) * delay_increment;

    // set register A to number of clock periods of delay
    RA = (i+1)*delay_increment/23.8095e-9;
    TC_SetRA(TC1,0, RA); 
    // set register C to number of clock periods when pulse should be off
    RC = pulse_off/23.8095e-9;
    TC_SetRC(TC1,0, RC);
    
    for (int avg_no = 0; avg_no < avg_points; avg_no++) {
      // emit trigger by set pin 2 HIGH (port B, pin 25)
      REG_PIOB_SODR |= (1<<25);
      //digitalWrite(trigger_pin, HIGH);      
      // wait for current delay time
      //delayMicroseconds(curr_delay * 1e6);
         
      counter_status = REG_TC1_SR0;
      counter_status = REG_TC1_CV0; // read registers to cler them
       
      REG_TC1_CCR0 = 0x5; // reset and enable counter      
      while ((REG_TC1_SR0 & 0x4) == 0); // check TC status register for occurance of RA compare
      // .. when RA compare has occured, start ADC conversion                      
      ADC->ADC_CR = 0x2;    
//      while ((ADC->ADC_ISR & 0x80) == 0); // check Interrupts Status Register (ADC_ISC) for channel 7 (A0) and wait for conversion
      while ((ADC->ADC_ISR & 0x40) == 0); //  wait for conversion by checking Interrupts Status Register (ADC_ISC) for channel 6 (A1), which is converted in second step
      intensity1_level = ADC->ADC_CDR[7]; //get values from Channel Data Register (ADC_CDR) of channel 7
      intensity2_level = ADC->ADC_CDR[6]; //get values from Channel Data Register (ADC_CDR) of channel 6 
      
//      REG_TC1_CCR0 = 0x2; // disable counter
      clock_counts = REG_TC1_CV0; // read counter
      curr_delay += clock_counts/avg_points*23.8095e-9;
          
      curr_sig1 += (float(intensity1_level) / (input_res - 1)) * (input_max - input_min) / avg_points;
      curr_sig2 += (float(intensity2_level) / (input_res - 1)) * (input_max - input_min) / avg_points;

//      cmd.sendBinCmd(ret_int, clock_counts);
//      delay(4);
      while ((REG_TC1_SR0 & 0x10) == 0); // check TC status register for occurance of RC compare
      // when RC compare has occured, set trigger pin to LOW and start new measurement
      // set pin 2 LOW (port B, pin 25)
      REG_PIOB_CODR |= (1<<25);
      //digitalWrite(trigger_pin, LOW);
      
    }
      
//    if (cont_signal_readout) {
      cmd.sendCmdStart(ret_fast_intensity);
      cmd.sendCmdBinArg(curr_sig1);
      cmd.sendCmdBinArg(curr_sig2);
      cmd.sendCmdBinArg(curr_delay);
//      cmd.sendCmdBinArg(curr_delay);
      cmd.sendCmdEnd();
//    }
  }
}

void on_set_delay() {  
  read_delay = cmd.readBinArg<float>();
}
void on_get_delay() {
  cmd.sendBinCmd(ret_float, read_delay);  
}

void on_intensity_trace() {

  int data_points = cmd.readBinArg<int>();
  int num_pulses = cmd.readBinArg<int>();
  float delay_increment = cmd.readBinArg<int>();  // in seconds
  float curr_sig1;
  float curr_sig2;
  float curr_delay;
  int intensity1_level;
  int intensity2_level;
  uint32_t counter_status; 
  uint32_t RA;
  uint32_t RC;

// set up Timer1, channel 0 (ID_TC3)
  pmc_enable_periph_clk(TC3_IRQn); 
  REG_TC1_CMR0= 0x00008040; // set channel mode register (use waveform mode, MCL/2 clock (42 MHz), no effects on triggers, TIOA, and TIOB, timer stopped when it reaches RC)
  REG_TC1_IDR0= 0xFF; // disable all interrupts

  for (int i = 0; i < data_points; i++) {
    curr_sig1 = 0;
    curr_sig2 = 0;
    curr_delay = i*delay_increment;
    
    for (int j = 0; j < num_pulses; j++) {      
      // set register A to number of clock periods of delay
      RA = read_delay/23.8095e-9;
      TC_SetRA(TC1,0, RA); 
      // set register C to number of clock periods when pulse should be off
      RC = 67e-6/23.8095e-9; // minimal time between two pulses (=15 kHz)
      TC_SetRC(TC1,0, RC);
    
      // emit trigger by set pin 2 HIGH (port B, pin 25)
      REG_PIOB_SODR |= (1<<25);         
      counter_status = REG_TC1_SR0;
      counter_status = REG_TC1_CV0; // read registers to cler them       
      REG_TC1_CCR0 = 0x5; // reset and enable counter      
      while ((REG_TC1_SR0 & 0x4) == 0); // check TC status register for occurance of RA compare
      // .. when RA compare has occured, start ADC conversion                      
      ADC->ADC_CR = 0x2;    
//    while ((ADC->ADC_ISR & 0x40) == 0); //  wait for conversion by checking Interrupts Status Register (ADC_ISC) for channel 6 (A1), which is converted in second step
      intensity1_level = ADC->ADC_CDR[7]; //get values from Channel Data Register (ADC_CDR) of channel 7
      intensity2_level = ADC->ADC_CDR[6]; //get values from Channel Data Register (ADC_CDR) of channel 6 
      
      curr_sig1 += (float(intensity1_level) / (input_res - 1)) * (input_max - input_min) / num_pulses;
      curr_sig2 += (float(intensity2_level) / (input_res - 1)) * (input_max - input_min) / num_pulses;

      while ((REG_TC1_SR0 & 0x10) == 0); // check TC status register for occurance of RC compare
      // when RC compare has occured, set trigger pin to LOW and start new measurement
      // set pin 2 LOW (port B, pin 25)
      REG_PIOB_CODR |= (1<<25); 
    }
    cmd.sendCmdStart(ret_fast_intensity);
    cmd.sendCmdBinArg(curr_sig1);
    cmd.sendCmdBinArg(curr_sig2);
    cmd.sendCmdBinArg(curr_delay);
    cmd.sendCmdEnd();
    delay(delay_increment*1000);    
  } 
}



/*
   reads the ADCs as fast as possible with the
   specified number of data points
*/
void on_fast_ADC_read() {
  
  int data_points = cmd.readBinArg<int>();
  float signal1[data_points];
  float signal2[data_points];
  for (int i = 0; i < data_points; i++) {
    ADC->ADC_CR |= (1<<1);  // start ADC conversion    
//    while ((ADC->ADC_ISR & 0x80) == 0); // check Interrupts Status Register (ADC_ISC) for channel 7 (A0) and wait for conversion
    while ((ADC->ADC_ISR & 0x40) == 0); //  wait for conversion by checking Interrupts Status Register (ADC_ISC) for channel 6 (A1), which is converted in second step
    signal1[i] = ADC->ADC_CDR[7]; //get values from Channel Data Register (ADC_CDR) of channel 7
    signal2[i] = ADC->ADC_CDR[6]; //get values from Channel Data Register (ADC_CDR) of channel 6 
  }
  cmd.sendCmdStart(ret_fast_intensity);
  for (int i = 0; i < data_points; i++) {
    cmd.sendCmdBinArg(signal1[i]);
  }
  cmd.sendCmdEnd();
}




// ===============
// Setup function
// ===============
void setup()
{
  /* mind default settings of ADC_MR:
   * ADC clock rate: MCL/(2*3) = 14 MHz -> 71.4 ns
   * startup time: 512 periods of ADC clock = 36.6 us
   * settling time: 17 periods of ADC clock = 1.2 us
   * tracking time: 1 ADC clock period = 71.4 ns => too fast (max. possible value: 15 x ADC clock period) => definitely to long for 1 us resolution
   * transfer time: 5 ADC clock periods = 357 ns
   */

  
  //default settings: ADC->ADC_MR= 0x9038ff00, binary: 0000 0000 1001 0000 0011 1000 1111 1111 0000 0000   
  // set the ADC registers
  ADC->ADC_CHER = 0xC0;       // enable ADC on pin A0 and A1 and that correspond to Ch7 and Ch6 on Atmel chip (for just A0: 0x80)
   
  // set all ADC_MR bits in once:
  //ADC->ADC_MR = 90000200; // binary: 1001 0000 0000 0000 0000 0010 0000 0000
  //ADC->ADC_MR = 90380280; // free running mode

  // free run mode, prescale: 1 (84/2/(1+1) = 21 MHz -> 47.6 ns),
  // startup: 24 ADC clock periods = 1.14 us, settling time: 5 ADC clock cycles = 238 ns,
  // tracking time: 5 ADC clock cycles = 238 ns,
  
// 0000 0000    1  0    01      0000     0   0    01     0011  0000 0001   1     000 0000
//        -- |USEQ|-|transfer|tracking|ANACH|-|settling|startup|prescale|freerun|..  
   ADC->ADC_MR = 0x90130180;     
   SerialUSB.println(REG_ADC_MR, HEX); 
   
// if calculated tracking time < 15 ADC clock periods: TRACKTIME = 0 and TRANSFER = 1
// => transfer periode = (1*2+3) ADC clock periods = 238 ns

// settling time only relevant when gain and offset of ADC channels changed => recommended settling time: 200 ns
// startup time only relevant when sleep mode on (and during first conversion?)
// max. ADC clock frequency: 22 MHz
// max. sampling frequency: 1 MHz (because conversion time take 20 ADC clock cycles)
   
  ADC->ADC_SEQR1= 0x67000000; // set the ADC read-out sequence, first measure A0 then A1, equals: (7<<24) | (6<<28) and 0000 0000 0110 0111 0000 0000 0000 0000 0000 0000)

  // setup pins
  //analogReadResolution(12);
  // analogReference(DEFAULT); // sets ADC reference to 3.3 V
  //pinMode(input1_pin, INPUT);
  //pinMode(input2_pin, INPUT);
  pinMode(trigger_pin, OUTPUT);
  pinMode(13, OUTPUT);

  // initialise USB port
  SerialUSB.begin(115200);
  delay(2000);
  
  SerialUSB.print("ADC Mode Register: ");
  SerialUSB.println(REG_ADC_MR, HEX); 

  // setup CmdMessenger
  attachCommands();
}


// Loop function
// =============
void loop()
{
  // process incoming serial data from PC and perform callbacks
  cmd.feedinSerialData();
}


