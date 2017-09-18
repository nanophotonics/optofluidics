# -*- coding: utf-8 -*-
"""

@author: Richard

Requires hologram_engine to be running and simple_gratings_and_lenses.py must have been executed
(found on SLM lab computer desktop in python folder)
"""

import socket
import time
import numpy as np
#import cv2
import os
from cameraclass import camera
from NKTLaserInterface import NKTLaserInterface



def GetVertSumOld(img):
    """For coloured Camera"""
    height, width, channels = img.shape
    img_array=np.zeros((height,width))
    for w in range(0,width-1):
        for h in range(0,height-1):
            img_array[h][w]=img.item(h,w,0)+img.item(h,w,1)+img.item(h,w,2)
    values=np.sum(img_array,axis=0)
    return values
def GetVertSum(img):
    """For b&w camera"""
    values=np.sum(img,axis=0)
    return values

def UDP_send(message):
    sock.sendto(message,ENGINE)
    data, addr = sock.recvfrom(128)
def make_spots(spots):
    spots = np.array(spots)
    assert spots.shape[1]==4, "Spots are 4 elements long - your array must be (n,4)"
    UDP_send("""<data>
    <uniform id="0">
    %s
    </uniform>ge
    <uniform id="1">
    %d
    </uniform>
    </data>
    """ % (" ".join(np.char.mod("%f", np.reshape(spots,spots.shape[0]*spots.shape[1]))),spots.shape[0]))


laser = NKTLaserInterface()
laser.writeParam('emission',True)
#repeat for all desired wavelengths
#for wavelength in np.arange(485,840,10):
for wavelength in np.arange(815,840,10):
    #if for loop crashes do laser.ser.close(), cam.closecam()
    #and delete any data folders for runs you want to repeat before restart
    print('starting {} nm'.format(wavelength))
    print('adjusting laser')
    laser.writeParam('longWavelength', (wavelength+5)*10)
    laser.writeParam('shortWavelength', (wavelength-5)*10)
    while laser.isGratingMoving():
        time.sleep(0.1)
    
    
        
        
    #Created on Mon Dec 01 12:04:12 2014
    date='2017-09-11'
    #date = time.strftime("%Y-%m-%d %H_%M_%S")
    wav='{}nm+-5nm'.format(wavelength)
     
    path = 'R:\\3-Temporary\\maw95\\Measurements\\SLM_Calibration\\{}-{}\\'.format(date,wav)
    
    os.makedirs('R:\\3-Temporary\\maw95\\Measurements\\SLM_Calibration\\{}-{}\\'.format(date,wav))
    os.makedirs('R:\\3-Temporary\\maw95\\Measurements\\SLM_Calibration\\{}-{}\\data\\'.format(date,wav))
    
    
    UDP_IP = "127.0.0.1"
    UDP_PORT = 61557
    ENGINE = (UDP_IP, UDP_PORT)
    
    sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM) #make a UDP socket
    sock.settimeout(1.0)
    #for new camera
    cam = camera(10)
    #autoexposure is slow, so don't do in loop! (constant exposure should be fine)
    print("Autoexposing camera [exposure: {} ms]".format(cam.autoexposure()))
    
    u = -1
    maxu=1000
    
    
    while u < maxu:
        u = u + 1
        print "u is %d" % u
        str1 = ("""<data>
    <network_reply>
    </network_reply>
    <window_rect>
    1921,0,512,512
    </window_rect>
    <shader_source>
    uniform vec4 spots[1];
    uniform int n;
    float uu = """
    +str(u)+
    """;    
    const float k=15700; //units of mm-1
    const float f=500.0; //mm
    const vec2 slmsize=vec2(7.68,7.68);
    void main(){ // basic gratings and lenses for a single spot
      vec2 uv = gl_TexCoord[0].xy*slmsize;
      vec3 pos = vec3(k*uv/f, k*dot(uv,uv)/(2.0*f*f));
      float phase, real=0.0, imag=0.0;
      if((pos.x-120)<0)
      {
      phase = uu/"""
          +str(maxu)+
          """;
      }
      else
      {
      phase= 0;
      }
//     for(int i; i<n; i++){
//       float rr = sqrt(((pos.x-120)*(pos.x-120))+((pos.y-120)*(pos.y-120)));
//       phase=  2*atan((pos.y-120)/(pos.x-120));
//       phase =  6.28*uu/"""
        +str(maxu)+
        """+ 0.1*(pos.x-120);
//       phase = mod(phase,6.28);    
//       real += spots[i][3] * sin(phase);
//        imag += spots[i][3] * cos(phase);
//     }
// phase = atan(real, imag);
      float g = phase; // map -pi to pi onto zero to one
      gl_FragColor=vec4(g,g,g,1);
    }
    </shader_source>
    </data>
    """)
        UDP_send(str1)  
        UDP_send("""<data>
    <uniform id=0>
    1 1 0.0 1
    </uniform>
    <uniform id=1>
    1
    </uniform>
    </data>
    """)
        
            
        
        time.sleep(0.1);
        #time.sleep(2)
        #for old camera
        #cap = cv2.VideoCapture(2);
        #cap.set(3,1280);
        #cap.set(4,1024);
        
        #for new camera
        frame = cam.takeimagebw()
        
        #data_array=[]
        #ret, frame = cap.read()
        data_array = GetVertSum(frame)
        
        
        for i in range(0,9):
            #old camera
            #ret, frame = cap.read()
            #new camera
            frame = cam.takeimagebw()
            #print frame
            #print 'new'
            data_array=[x + y for x, y in zip(data_array, GetVertSum(frame))]
            #time.sleep(1)
        
        file = open(path+'data\\summed_data_averaged_'+str(i+2)+'_times_u_'+str(u+1000)+'.txt','w')
        
        for item in data_array:
            file.write(str(item))
            file.write('\n')
            
        file.close()
            
    #    # When everything done, release the capture
        #cap.release()    
        #cv2.destroyAllWindows()
        #cv2.imwrite(path+'800nm_SLM'+str(u)+'.jpg', frame, [int(cv2.IMWRITE_JPEG_QUALITY), 100])
        #gray = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY)
        ##cv2.imwrite(path+'610nm_SLM'+str(u)+'.bmp', frame)
        #cv2.imshow('image',frame)
    #    #cv2.waitKey(0)
        ##time.sleep(0.5)
    
    #for new camera (release camera)
    #old routine opens and closes in loop...
    cam.closecam()
    sock.close()
    # When everything done, release the capture
    #cap.release()
    #cv2.destroyAllWindows()
    
    #cv2.destroyAllWindows()

laser.writeParam('emission', False)
del laser