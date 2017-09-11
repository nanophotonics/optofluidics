# -*- coding: utf-8 -*-
"""
Created on Tue Sep 05 17:18:10 2017

@author: Marius Weber
"""

import serial
import serial.tools
import numpy as np
import crcmod
import time

from NKTLaserExceptions import CorruptionException, UnexpectedResponseException
from NKTLaserExceptions import MsgTooLongException, MsgStartException

#global variable names (represent byte values to make more legible)
#start/end of telegram
SOT=0x0D
EOT=0x0A
#special character substitution byte
SOE=0x5E

#message types - list of all message types can be found in the manual
#acknowlege
ACK=3
READ=4
WRITE=5
#data response to read
DATAGRAM=8

#laser registers - all possible regisers may be found in the manual
emissionReg =  0x30
powerLevelReg = 0x37
pulsePickerRatioReg = 0x34
interlockReg = 0x32
#varia registers
shortWavelengthReg = 0x34
longWavelengthReg = 0x33
variaStatusReg = 0x66

class NKTLaserInterface():
    """A class for interfacing with an NKT Laser
    
    The class was developed for the NKT SuperK Extreme with attached Varia 
    filter box. However, it should work with any NKT laser though some 
    features may not apply. When using with a different laser verify the 
    register values are correct by checking with the manual.
    """
    def __init__(self, comPort='COM7', laserAdd=15, variaAdd=16, sourceAdd=161):
        """Initialise class
        
        Addresses default to values encountered when writing the class.
        SUBJECT TO CHANGE!!!
        Current addresses can be determined by looking at values displayed when
        connecting using NKT provided GUI.
        """
        #addresses
        self.laserAdd = laserAdd
        self.variaAdd = variaAdd
        self.sourceAdd = sourceAdd
        #connect to serial port with correct settings
        self.ser = serial.Serial(comPort, baudrate=115200, 
                            timeout=1, write_timeout=1)
        
    def __del__(self):
        """Close connection on destruct"""
        self.closeConnection()
        
    def closeConnection(self):
        """close connection to laser"""
        if self.ser.is_open:
            self.ser.close()
            
    def openConnection(self):
        """open connection to laser"""
        if not self.ser.is_open:
            self.ser.open()
    
    def isConnected(self):
        """returns true if connected to laser, else false"""
        return self.ser.is_open
                        
    def _checkCRC(self, data):
        """calculate CRC-CCITT (xModem) algorithm to verify message integrity
        
        returns bytearray [MSB, LSB]"""
        #CRC-CCITT (xModem) algorithm
        crc = crcmod.Crc(0x11021, rev=False, initCrc=0x0000, xorOut=0x0000)
        crc.update(bytes(data))
        return bytearray(crc.digest())
    
    #Telegram format: 
    #SOT DestinationAdd SourceAdd Type Register# dataBytes CRC-MSB CRC-LSB EOT
    def _send(self, telegram):
        """Final processing and sends message to device
        
        
        Calculates checksum, check special characters and attaches start and 
        stop bytes.
        """
        #calculate and append checksum
        CRCMSB, CRCLSB = self._checkCRC(telegram)
        telegram.append(CRCMSB)
        telegram.append(CRCLSB)
        #replace special characters (SOT, EOT, SOE)
        for i in range(len(telegram)):
            if telegram[i]==SOT or telegram[i]==EOT or telegram[i]==SOE:
                telegram[i] +=64
                #insert special charater flag
                telegram.insert(i,SOE)
        #start and end byte
        telegram.insert(0,SOT)
        telegram.append(EOT)
        
        self.ser.write(telegram)
    
    def _receive(self):
        """Receive answer from device"""
        #read byte for byte
        #message starte byte given by SOT
        #message end byte given by EOT
        telegram = bytearray(self.ser.read(1))
        if telegram[0] != SOT:
            raise MsgStartException(
            'Unexpected start of telegram!\n Received {}, but expect {}'.
            format(telegram[0], SOT))
        #telegram.append(self.ser.read(1)[0])
        #start at 1 as previously read start byte
        i=1
        while telegram[-1]!= EOT:
            telegram.append(self.ser.read(1)[0])
            if i >= 248:
                raise MsgTooLongException(
                        'The recieved message is too long!\n' + 
                        'Aborted after {} bytes. '.format(i+1) +
                        'Maximum permissible length is 248 bytes.')
            i +=1
        #remove start and stop byte
        telegram = telegram[1:-1]
        #find special characters
        for i in range(len(telegram)-telegram.count(bytearray([SOE]))):
            if telegram[i]==SOE:
                telegram[i+1] -=64 
                #remove special charater flag
                del telegram[i]
        
        #verify checksum
        CRCMSB, CRCLSB = self._checkCRC(telegram[0:-2])
        if CRCMSB != telegram[-2] or CRCLSB != telegram[-1]:
            raise CorruptionException('The recieved message is corrupted!\n' + 
                        'Checksum does not match message!')
        #remove checksum in answer
        return telegram[0:-2]
        
    
    def _read(self, destAdd, sourceAdd, register):
        """Send read request to device at destAdd for register
        
        all inputs are expected as uint8
        
        returns response data bytes
        """
        #write data in telegram
        telegram = bytearray([destAdd, sourceAdd, READ, register])
        resend = True
        i=0
        while resend:
            i += 1
            self._send(telegram)
            try:
                answer = self._receive()
                #answer=telegram
                resend = False
            except CorruptionException:
                #resend remains True if an exception occurs
                print('re-sending telegam due to corruption')
                #raise exception after 10 tries
                if i>10:
                    raise
            #except MsgTooLongException:
            #    pass
            #except MsgStartException:
            #    pass
        #unpack response
        ansDest = answer[0]
        ansSource = answer[1]
        ansType = answer[2]
        ansRegister = answer[3]
        
        #verify answer is as expected
        if (ansDest != sourceAdd or ansSource !=destAdd 
                or ansRegister != register or ansType != DATAGRAM ):
            raise UnexpectedResponseException(
            'The response did not match expectation!\n' + 
            'Received:\n' +
            'Destination: {}, Source: {}, Register: {}, Type: {}\n'.format(
            ansDest, ansSource, ansRegister, ansType) +
            'Expected:\n' +
            'Destination: {}, Source: {}, Register: {}, Type: {}\n'.format(
            sourceAdd, destAdd, register, DATAGRAM)      
            )
        
        #this data must still be decoded according to the expected format
        ansData = answer[4:]
        return ansData
        
    def _write(self, destAdd, sourceAdd, register, data):
        """
        Write data to register of device at destAdd 
        
        destAdd, sourceAdd and register are expected as uint8
        data is expected as a bytearray
        """
        #write data in telegram
        telegram = bytearray([destAdd, sourceAdd, WRITE, register])
        telegram.extend(data)
        resend = True
        i=0
        while resend:
            i += 1
            self._send(telegram)
            try:
                answer = self._receive()
                #don't resend if previous step succeded
                resend = False
            except CorruptionException:
                #resend remains True if an exception occurs
                print('re-sending telegam due to corruption')
                #raise exception after 10 tries
                if i>10:
                    raise
            #except MsgTooLongException:
            #    pass
            #except MsgStartException:
            #    pass
        
        #unpack response
        ansDest = answer[0]
        ansSource = answer[1]
        ansType = answer[2]
        ansRegister = answer[3]
        
        #check confirmation by device
        if (ansDest != sourceAdd or ansSource !=destAdd 
                or ansRegister != register or ansType != ACK ):
            raise UnexpectedResponseException(
            'The response did not match expectation!\n' + 
            'Received:\n' +
            'Destination: {}, Source: {}, Register: {}, Type: {}\n'.format(
            ansDest, ansSource, ansRegister, ansType) +
            'Expected:\n' +
            'Destination: {}, Source: {}, Register: {}, Type: {}\n'.format(
            sourceAdd, destAdd, register, DATAGRAM)      
            )
    def isGratingMoving(self):
        """Returns True if filter-box gratings are moving, else False"""
        status = self.readParam('variaStatus')
        if status >>12 & 1 or  status >>13 & 1 or  status >>14 & 1:
            return True
        else:
            return False
        
    def readParam(self, param):
        """
        Reads the specified parameter from the appropreate device
        
        The following is formated as follows
        parameter_string: describtion (retrun_format)
        
        Allowed laser parameters are
        emission: is emission on? (bool)
        powerLevel: power lever in permille (uint16)
        pulsePickerRatio: not implemented
        interlock: not implemented
        
        Allowed varia parameters are
        shortWavelength: lower wavelength limit /(0.1 nm) (uint16)
        longWavelength: upper wavelength limit /(0.1 nm) (uint16)
        variaStatus: varia status bits (uint16)
        """
        
        if param == 'emission':
            data = self._read(self.laserAdd, self.sourceAdd, emissionReg)
            if data[0] == 0:
                return False
            elif data[0] == 3:
                return True
            else:
                raise UnexpectedResponseException(
                'Emission read returned {}, but only 0 or 3 are allowed!'.
                format(data[0]))
        elif param == 'powerLevel':
            data = self._read(self.laserAdd, self.sourceAdd, powerLevelReg)
            if len(data) != 2:
                raise UnexpectedResponseException(
                'Power level read returned {} bytes, but expect 2 bytes!'.
                format(len(data)))
            #16 bit integer (LSB first)
            return data[0] + (data[1]<<8)
        elif param == 'pulsePickerRatioReg':
            raise NotImplementedError()
        elif param == 'interlockRe':
            raise NotImplementedError() 
        
        #varia from here
        elif param == 'shortWavelength':
            data = self._read(self.variaAdd, self.sourceAdd, shortWavelengthReg)
            if len(data) != 2:
                raise UnexpectedResponseException(
                'Short Wavelength read returned {} bytes, but expect 2 bytes!'.
                format(len(data)))
            #16 bit integer (LSB first)
            return data[0] + (data[1]<<8)
        elif param == 'longWavelength':
            data = self._read(self.variaAdd, self.sourceAdd, longWavelengthReg)
            if len(data) != 2:
                raise UnexpectedResponseException(
                'Long Wavelength read returned {} bytes, but expect 2 bytes!'.
                format(len(data)))
            #16 bit integer (LSB first)
            return data[0] + (data[1]<<8)
        elif param == 'variaStatus':
            data = self._read(self.variaAdd, self.sourceAdd, variaStatusReg)
            if len(data) != 2:
                raise UnexpectedResponseException(
                'Long Wavelength read returned {} bytes, but expect 2 bytes!'.
                format(len(data)))
            #16 bit integer (LSB first)
            return data[0] + (data[1]<<8)
        else:
            raise ValueError(
            "Read parameter '{}' not recognised.".format(param))
    
    def writeParam(self, param, value):
        """
        Writes the specified value to the appropreate device
        
        The following is formated as follows
        parameter_string: describtion (expected_value_format)
        
        Allowed laser parameters are
        emission: turn emission on? (bool)
        powerLevel: power lever in permille (uint16)
        pulsePickerRatio: not implemented
        interlock: not implemented
        
        Allowed varia parameters are
        shortWavelength: lower wavelength limit /(0.1 nm) (uint16)
        longWavelength: upper wavelength limit /(0.1 nm) (uint16)
        """
        if param == 'emission':
            if value:
                data = bytearray([3])
            else:
                data = bytearray([0])
            self._write(self.laserAdd, self.sourceAdd, emissionReg, data)
        elif param == 'powerLevel':
            #16 bit integer (LSB first)
            data = bytearray([np.uint8(value), np.uint8(value>>8)])
            self._write(self.laserAdd, self.sourceAdd, powerLevelReg, data)
        elif param == 'pulsePickerRatioReg':
            raise NotImplementedError()
        elif param == 'interlockRe':
            raise NotImplementedError() 
        
        #varia from here
        elif param == 'shortWavelength':
            #16 bit integer (LSB first)
            data = bytearray([np.uint8(value), np.uint8(value>>8)])
            self._write(self.variaAdd, self.sourceAdd, shortWavelengthReg, data)
        elif param == 'longWavelength':
            #16 bit integer (LSB first)
            data = bytearray([np.uint8(value), np.uint8(value>>8)])
            self._write(self.variaAdd, self.sourceAdd, longWavelengthReg, data)
        else:
            raise ValueError(
            "Read parameter '{}' not recognised.".format(param))


if  __name__ == '__main__':
    #print(ser.is_open)
    #ser.read(1)
    
    #print(readParam('shortWavelength'))
    #writeParam('shortWavelength',5800)
    #print(readParam('shortWavelength'))
    #writeParam('shortWavelength',5900)
    #print(readParam('shortWavelength'))
    
    laser = NKTLaserInterface()
    print(laser.readParam('emission'))
    
    print(laser.isGratingMoving())
    #laser.writeParam('longWavelength',6000)
    
    print(laser.isGratingMoving())
    time.sleep(2)
    print(laser.isGratingMoving())
    
    
    del laser
    
    
    print('\ndone\n')