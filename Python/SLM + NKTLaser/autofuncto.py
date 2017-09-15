import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np
import time
import math
from scipy.ndimage.interpolation import zoom, rotate
import scipy
import cv2



def rotatept(origin, point, angle):
    """
    Rotate a point counterclockwise by a given angle around a given origin.

    The angle should be given in radians.
    """
    ox, oy = origin
    px, py = point

    qx = ox + math.cos(angle) * (px - ox) - math.sin(angle) * (py - oy)
    qy = oy + math.sin(angle) * (px - ox) + math.cos(angle) * (py - oy)
    return np.array([qx, qy])

def hexagon(origin, radius, angle):
    pointx=np.array([origin[0]])
    pointy=np.array([origin[1]])
    pointx=pointx-radius
    new=rotatept(origin, np.hstack([pointx[0],pointy[0]]),angle)
    pointx=[new[0]]
    pointy=[new[1]]
    for i in range(1,6):
        pointx=np.hstack([pointx, rotatept(origin, np.hstack([pointx[i-1],pointy[i-1]]),np.pi/3)[0]])
        pointy=np.hstack([pointy, rotatept(origin, np.hstack([pointx[i-1],pointy[i-1]]),np.pi/3)[1]])
    pointx=np.hstack([pointx,pointx[0]])
    pointy=np.hstack([pointy,pointy[0]])
    return pointx, pointy
        
class LineBuilder:
    def __init__(self, line):
        self.line = line
        self.counter=0
        self.xs = list(line.get_xdata())
        self.ys = list(line.get_ydata())
        self.cid = line.figure.canvas.mpl_connect('button_press_event', self)

    def __call__(self, event):
        if event.inaxes!=self.line.axes: return
        self.xs.append(event.xdata)
        self.ys.append(event.ydata)
        self.line.set_data(self.xs, self.ys)
        self.line.figure.canvas.draw()
        self.counter=self.counter+1
        if self.counter==6:
            self.line.set_data(np.concatenate((self.xs, [self.xs[0]])), np.concatenate((self.ys, [self.ys[0]])))
            self.line.figure.canvas.draw()
            
            
def areainterest(image):
    img=image
    fig = plt.figure()
    axes = fig.gca()
    axes.set_xlim([0,np.shape(img)[1]])
    axes.set_ylim([0,np.shape(img)[0]])
    
    plt.imshow(img,zorder=0,cmap='gray')
    plt.clim(0,80)
    ax = fig.add_subplot(111)
    ax.set_title('Click on hexagon corners in clockwise fashion.')
    line, = ax.plot([], [])  # empty line
    linebuilder = LineBuilder(line)
    plt.show()
    
    while linebuilder.counter<6:
        line.figure.canvas.start_event_loop(linebuilder.cid)
    fig.canvas.mpl_disconnect(linebuilder.cid)
    xs=linebuilder.xs
    ys=linebuilder.ys
    plt.close() 
    
    
    cent=np.array([np.mean(xs), np.mean(ys)])
    dist=0
    for i in range (0,np.shape(xs)[0]):
        dist=dist+np.linalg.norm(cent-np.array([xs[i],ys[i]]))/6
    rotx=np.array([xs[0]])
    roty=np.array([ys[0]])
    for i in range (1,np.shape(xs)[0]):
        rotx=np.hstack([rotx, rotatept(cent, np.hstack([xs[i],ys[i]]),np.pi/3*i)[0]])
        roty=np.hstack([roty, rotatept(cent, np.hstack([xs[i],ys[i]]),np.pi/3*i)[1]])
    angle=np.arctan((np.mean(roty)-cent[1])/(np.mean(rotx)-cent[0]))
    return cent, dist, angle

def maskgen(mask,currentimg):
    if np.array(np.shape(np.shape(mask)))[0]==3:
        mask=mask[:,:,0]
    mask=np.pad(mask,((np.shape(mask)[0],0),(int(np.shape(mask)[1]/2),int(np.shape(mask)[1]/2))),'constant')
    mask2=np.zeros(np.shape(mask))
    for i in range(0,6):
        mask2=mask2+rotate(mask,60*i,reshape=False)
    mask=mask2    
    ratio=(np.array(np.shape(currentimg))).astype(float)/(np.array(np.shape(mask))).astype(float)
    mask=zoom(mask,np.max(ratio))
    diff=np.array(np.shape(mask))-np.array(np.shape(currentimg))
    mask=mask[0:np.shape(mask)[0]-diff[0],0:np.shape(mask)[1]-diff[1]] #Note: here we assume mask produces a close to hexagon rectangle.
    return mask

def norm(mask,cutmask):
    mask=mask/np.sum(mask*cutmask)
    return mask

def cropcurrent(currentimg, cent, dist, angle):
    currentimg=currentimg[int(cent[1])-np.int(dist):int(cent[1])+np.int(dist), int(cent[0])-np.int(dist):int(cent[0])+np.int(dist)]
    currentimg=rotate(currentimg,angle/(2*np.pi)*360,reshape=False)
    currentimg=currentimg[int(np.shape(currentimg)[1]/2)-np.int(dist*np.sqrt(3)/2):int(np.shape(currentimg)[1]/2)+np.int(dist*np.sqrt(3)/2), int(np.shape(currentimg)[0]/2)-np.int(dist):int(np.shape(currentimg)[0]/2)+np.int(dist)]
    return currentimg

def optsetup(img, maskfile,cutmaskfile):
    cent, dist, angle=areainterest(img)
    currentimg=cropcurrent(img, cent, dist, angle)
    mask=maskgen(mpimg.imread(maskfile),currentimg)
    cutmask=maskgen(mpimg.imread(cutmaskfile),currentimg)
    cutmask=cutmask/np.amax(cutmask)
    cutmask=((cutmask.astype(float)/(np.amax(cutmask)))>0.5).astype(float)
    mask=norm(mask,cutmask)
    return mask, cutmask, np.array([cent, dist, angle])

def compare(img, mask,cutmask, hexcor):
    maxcorr=0
    for angcorr in np.arange(-0.02,0.03,0.01):
        currentimg=cropcurrent(img, hexcor[0], hexcor[1], hexcor[2]+angcorr)
        currentimg=norm(currentimg,cutmask)
        res=cv2.matchTemplate(np.pad(currentimg,((30,30),(30,30)),'constant').astype('float32'), mask.astype('float32'),cv2.TM_SQDIFF) 
        if (np.max(res)>maxcorr):
            maxcorr=np.max(res)   
    return maxcorr, currentimg
