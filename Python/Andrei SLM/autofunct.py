#this program provides a set of useful functions for automating SLM mode optimisation.
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np
import time
import math
from scipy.ndimage.interpolation import zoom, rotate
import scipy
import cv2



def rotatept(origin, point, angle):
	#Rotate a point counterclockwise by a given angle around a given origin.
    	#The angle should be given in radians.
    ox, oy = origin
    px, py = point

    qx = ox + math.cos(angle) * (px - ox) - math.sin(angle) * (py - oy)
    qy = oy + math.sin(angle) * (px - ox) + math.cos(angle) * (py - oy)
    return np.array([qx, qy])

def hexagon(origin, radius, angle):
	#Generates a series of x,y coordinate points that represent a hexagon.
    #Angle represents the orientation of the hexagon, radius its size and origin its centre.
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
	#This class enables the areainterest function to draw lines on a displayed image.
    def __init__(self, line, ngon):
        self.line = line
        self.ngon=ngon
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
        if self.counter==self.ngon:
            self.line.set_data(np.concatenate((self.xs, [self.xs[0]])), np.concatenate((self.ys, [self.ys[0]])))
            self.line.figure.canvas.draw()
            
            
def areainterest(image, ngon):
	#Lets the user identify an n-gonal region in an image interactively.
    image=image.astype(float)
    img=image/np.max(image)*255
    fig = plt.figure()
    axes = fig.gca()
    axes.set_xlim([0,np.shape(img)[1]])
    axes.set_ylim([0,np.shape(img)[0]])
    
    plt.imshow(img,zorder=0,cmap='gray')
    plt.clim(0,80)
    ax = fig.add_subplot(111)
    ax.set_title('Click on {}-gon corners in clockwise fashion.'.format(str(ngon)))
    line, = ax.plot([], [])
    linebuilder = LineBuilder(line, ngon)
    plt.show()
    
    while linebuilder.counter<ngon:
        line.figure.canvas.start_event_loop(linebuilder.cid)
    fig.canvas.mpl_disconnect(linebuilder.cid)
    xs=linebuilder.xs
    ys=linebuilder.ys
    plt.close() 
    
    
    cent=np.array([np.mean(xs), np.mean(ys)])
    dist=0
    for i in range (0,np.shape(xs)[0]):
        dist=dist+np.linalg.norm(cent-np.array([xs[i],ys[i]]))/ngon
    rotx=np.array([xs[0]])
    roty=np.array([ys[0]])
    for i in range (1,np.shape(xs)[0]):
        rotx=np.hstack([rotx, rotatept(cent, np.hstack([xs[i],ys[i]]),2*np.pi/ngon*i)[0]])
        roty=np.hstack([roty, rotatept(cent, np.hstack([xs[i],ys[i]]),2*np.pi/ngon*i)[1]])
    angle=np.arctan((np.mean(roty)-cent[1])/(np.mean(rotx)-cent[0]))
    angle=np.mod(angle,2*np.pi/ngon)
    return cent, dist, angle


def maskgenr(mask,currentimg,n):
	#Generates a mask by making n copies of an image (labelled i=1, 2, ... n)
    #it rotates each mask by 360/n*i degrees and combines the results.
    #Ensures that mask has same dimensions as currentimg.
    if np.array(np.shape(np.shape(mask)))[0]==3:
        mask=mask[:,:,0]
    mask2=np.zeros(np.shape(mask))
    for i in range(0,n):
        mask2=mask2+rotate(mask,360/n*i,reshape=False)
    mask=mask2    
    ratio=(np.array(np.shape(currentimg[:,:]))).astype(float)/(np.array(np.shape(mask))).astype(float)
    mask=zoom(mask,np.max(ratio))
    diff=np.array(np.shape(mask))-np.array(np.shape(currentimg[:,:]))
    mask=mask[0:np.shape(mask)[0]-diff[0],0:np.shape(mask)[1]-diff[1]] #Note: here we assume mask produces a close to hexagon rectangle.
    mask=mask.astype('float')/np.sum(mask.astype('float'))    
    return mask

def norm(mask,cutmask):
	#Normalisation function.
    mask=mask.astype('float32')/np.sum(mask.astype('float32')*cutmask.astype('float32'))
    return mask

def cropcurrent(currentimg, cent, dist, angle):
	#Crops an image (currentimg) to a given hexagonal region.
    currentimg=currentimg[int(cent[1])-np.int(dist):int(cent[1])+np.int(dist), int(cent[0])-np.int(dist):int(cent[0])+np.int(dist)]
    currentimg=rotate(currentimg,angle/(2*np.pi)*360,reshape=False)
    currentimg=currentimg[int(np.shape(currentimg)[1]/2)-np.int(dist*np.sqrt(3)/2):int(np.shape(currentimg)[1]/2)+np.int(dist*np.sqrt(3)/2), int(np.shape(currentimg)[0]/2)-np.int(dist):int(np.shape(currentimg)[0]/2)+np.int(dist)]
    return currentimg

def optsetup(img, maskfile,cutmaskfile, n):
	#Runs a selection of the predefined functions to start the optimisation process.
    cent, dist, angle=areainterest(img,6) #User input to identify hexagonal core.
    currentimg=cropcurrent(img, cent, dist, angle) #Crop image to get only the core.
    mask=maskgenr(mpimg.imread(maskfile),currentimg, n) #Generate a mask that has same size as cropped image (note that n is set to 1 when used)
    cutmask=maskgenr(mpimg.imread(cutmaskfile),currentimg, 6) #The cutmask is used to isolate the core: it is generated by rotating an equilateral triangle of 1s 6 times to generate the cutmask defined below.
    cutmask=cutmask/np.amax(cutmask) #This just tidies the result.
    cutmask=((cutmask.astype(float)/(np.amax(cutmask)))>0.5).astype(float)#This just tidies the result.
    return mask, cutmask, np.array([cent, dist, angle])

def compare(img, mask, cutmask, hexcor):
	#Calculates the fitness function.
    #img is the raw camera image
    #mask is the simulation result for comparison
    #cutmask is a matrix of the same size as mask but with 1 in the hexagonal region and 0 everywhere else
    #hexcor is a vector consisting of [cent, dist, angle] (i.e. coordinates to identify hexagonal core in raw image)
    mincorr=float(99999999)
    bestimg=[]
    for angcorr in np.arange(-0.02,0.03,0.01):
        currentimg=cropcurrent(img, hexcor[0], hexcor[1], hexcor[2]+angcorr)
        res=cv2.matchTemplate(np.pad(norm(currentimg, cutmask),((50,50),(50,50)),'constant').astype('float32'), norm(mask.astype('float32'), cutmask),cv2.TM_SQDIFF_NORMED) 
        if (np.amin(res)<mincorr):
            mincorr=np.amin(res)
            bestimg=currentimg
    return mincorr, bestimg