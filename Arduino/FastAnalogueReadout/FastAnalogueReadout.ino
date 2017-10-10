// *** AFM/STM Electronic ***

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
const int input1_pin   = A1;
const int input2_pin   = A2;
const int trigger_pin = 2;

// global variables
// ================
float intensity1;
float intensity2;
int input_res = 4096;
float input_max = 3.3;
float input_min = 0;

bool cont_signal_readout = true;
int adc_readout_freq = 1e4; 
int signal_avg = 10;     // number data point to average
int signal_avg_no = 0;
float signal1_temp = 0;
float signal2_temp = 0;

// get class instances for serial communication with PC
CmdMessenger cmd = CmdMessenger(SerialUSB);

// ==========================================================
// Initialisation of CmdMessenger and definition of commands
// ==========================================================

// definition of supported commands
enum
{
  check_arduino,
  get_intensity,
  ret_intensity,           // returns a float value of the last intensity measurements
  set_cont_signal_readout,
  get_cont_signal_readout,
  set_daq_freq,
  get_daq_freq,
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
  cmd.attach(get_intensity, on_get_intensity);
  cmd.attach(set_cont_signal_readout, on_set_cont_signal_readout);
  cmd.attach(get_cont_signal_readout, on_get_cont_signal_readout);
  cmd.attach(set_daq_freq, on_set_daq_freq);
  cmd.attach(get_daq_freq, on_get_daq_freq);
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

void on_get_intensity() {
  cmd.sendCmdStart(ret_intensity);
  cmd.sendCmdBinArg(intensity1);
  cmd.sendCmdBinArg(intensity2);
  cmd.sendCmdEnd();
}
  
void on_set_cont_signal_readout() {
  cont_signal_readout = cmd.readBinArg<bool>();
}

void on_get_cont_signal_readout() {
  cmd.sendBinCmd(ret_bool, cont_signal_readout);
}

void on_set_daq_freq() {
  adc_readout_freq = cmd.readBinArg<int>();
  Timer3.setFrequency(adc_readout_freq).start();
}

void on_get_daq_freq() {
  cmd.sendBinCmd(ret_int, adc_readout_freq);
}


// ======================================
// sleep functions to pause another timer
// ======================================

class Sleep {
  public:
    static void sleep(long delay_time) {
      Timer4.stop();
      Timer5.attachInterrupt(when_awake).setPeriod(delay_time).start();
    }
  private:    
    static void when_awake() {
      Timer5.stop();
      Timer4.start();  
    }
};

void readADCs()
{
  /* This method reads the ADCs.
     It is executed with a frequency 'adc_readout_freq' set by Timer3 and
     writes the averaged intensity into the global variables 'intensity1'
     and 'intensity2'.
  */
  int intensity1_level = analogRead(input1_pin);
  int intensity2_level = analogRead(input2_pin);
  float _intensity1 = (float(intensity1_level) / (input_res - 1)) * (input_max - input_min);
  float _intensity2 = (float(intensity2_level) / (input_res - 1)) * (input_max - input_min);
  
  // averaging of ADC signals  
  if (signal_avg_no < signal_avg) {
    signal1_temp += _intensity1 / signal_avg;
    signal2_temp += _intensity2 / signal_avg;
    signal_avg_no++;
  }
  else {
    intensity1 = signal1_temp;
    intensity2 = signal2_temp;    
    signal1_temp = _intensity1;
    signal2_temp = _intensity2;
    signal_avg_no = 1; 
    
    if (cont_signal_readout) {
      cmd.sendCmdStart(ret_intensity);
      cmd.sendCmdBinArg(intensity1);
      cmd.sendCmdBinArg(intensity2);
      cmd.sendCmdEnd();
    }  
  }  
}


// ===============
// Setup function
// ===============
void setup()
{
  // bits for fast ADC readout
  REG_ADC_MR = (REG_ADC_MR & 0xFFF0FFFF) | 0x00020000; // change ADC mode register
  REG_ADC_MR = (REG_ADC_MR & 0xFFFFFF0F) | 0x00000080; //enable FREERUN mode and change ADC to continuous mode
  
  // setup pins
  analogReadResolution(12);
  // analogReference(DEFAULT); // set ADC reference to 3.3 V
  pinMode(input1_pin, INPUT);
  pinMode(input2_pin, INPUT);
  pinMode(trigger_pin, OUTPUT);
  pinMode(13, OUTPUT);

  // initialise USB port
  SerialUSB.begin(115200);

  // setup CmdMessenger
  attachCommands();
  // setup a timer (interrupt) to constantly read the ADCs and
  // run the feedback if activated. The Arduino Due has 9
  // available Timer Couters (TC), with TC3, TC4, and TC5 not
  // beeing mapped to a physical pin. So TC3 is used here.
  Timer3.attachInterrupt(readADCs).setFrequency(adc_readout_freq).start();
}


// Loop function
// =============
void loop()
{
  // process incoming serial data from PC and perform callbacks
  cmd.feedinSerialData();
}


