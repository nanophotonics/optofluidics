import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np
import time
import math
from scipy.ndimage.interpolation import zoom, rotate
import scipy.optimize
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
    """Allowes user to identify a regular n-gon shape"""
    image=image.astype(float)
    img=image/np.max(image)*255#*2
    fig = plt.figure()
    axes = fig.gca()
    axes.set_xlim([0,np.shape(img)[1]])
    axes.set_ylim([0,np.shape(img)[0]])
    
    plt.imshow(img,zorder=0,cmap='gray')
    plt.clim(0,80)
    ax = fig.add_subplot(111)
    ax.set_title('Click on {}-gon corners in clockwise fashion.'.format(str(ngon)))
    line, = ax.plot([], [])  # empty line
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
    
    
def bandgapinterest(image, ngon):
    """
    Allow user to identify bandgap fibre
    
    Inputs
    image: array - immage to show user
    pointNo: number of points to be specified (positive int)
    
    Outputs
    cent: vector (x,y) - centre position
    radius: scalar
    """
    image=image.astype(float)
    img=image/np.max(image)*255#*2
    fig = plt.figure()
    axes = fig.gca()
    axes.set_xlim([0,np.shape(img)[1]])
    axes.set_ylim([0,np.shape(img)[0]])
    
    plt.imshow(img,zorder=0,cmap='gray')
    plt.clim(0,80)
    ax = fig.add_subplot(111)
    ax.set_title('Click on the {} smaller corners  in a clockwise fashion.'.format(str(ngon)))
    line, = ax.plot([], [])  # empty line
    linebuilder = LineBuilder(line, ngon)
    plt.show()
    
    while linebuilder.counter<ngon:
        line.figure.canvas.start_event_loop(linebuilder.cid)
    fig.canvas.mpl_disconnect(linebuilder.cid)
    xs=linebuilder.xs
    ys=linebuilder.ys
    plt.close() 
    
    
    pointCoords = np.column_stack((xs,ys))
    #initail guess for least squares
    cent0=np.mean(pointCoords,axis=1)
    #rms distance
    dist0 = np.sum(np.sqrt( np.sum( (pointCoords-cent0)**2, axis=0) ))/ngon 
    
    #fit parameter vector
    x0=np.append(cent0,dist0)
    #now implement least squared fit
    xopt = scipy.optimize.fmin(chisquaredcircle, x0, args=(pointCoords))
    
    cent=xopt[0:1]
    radius = xopt[2]
    
    
    rotx=np.array([xs[0]])
    roty=np.array([ys[0]])
    for i in range (1,np.shape(xs)[0]):
        rotx=np.hstack([rotx, rotatept(cent, np.hstack([xs[i],ys[i]]),2*np.pi/ngon*i)[0]])
        roty=np.hstack([roty, rotatept(cent, np.hstack([xs[i],ys[i]]),2*np.pi/ngon*i)[1]])
    angle=np.arctan((np.mean(roty)-cent[1])/(np.mean(rotx)-cent[0]))
    angle=np.mod(angle,2*np.pi/ngon)
    
    return cent, radius, angle


def chisquaredcircle(param, points):
    """return sum of square deviation of points from circle
    
    Inputs
    param: vector (centrex, centrey, radius)
    points: 2D array [x/y, pointIndex]
    """
    cent=param[0:1]
    radius = param[2]
    return np.sum( np.sum((points-cent)**2 ,axis = 0) -radius**2 )




'''
def maskgen(mask,currentimg):
    if np.array(np.shape(np.shape(mask)))[0]==3:
        mask=mask[:,:,0]
    mask=np.pad(mask,((np.shape(mask)[0],0),(int(np.shape(mask)[1]/2),int(np.shape(mask)[1]/2))),'constant')
    mask2=np.zeros(np.shape(mask))
    for i in range(0,6):
        mask2=mask2+rotate(mask,60*i,reshape=False)
    mask=mask2    
    ratio=(np.array(np.shape(currentimg[:,:,0]))).astype(float)/(np.array(np.shape(mask))).astype(float)
    mask=zoom(mask,np.max(ratio))
    diff=np.array(np.shape(mask))-np.array(np.shape(currentimg[:,:,0]))
    mask=mask[0:np.shape(mask)[0]-diff[0],0:np.shape(mask)[1]-diff[1]] #Note: here we assume mask produces a close to hexagon rectangle.
    return mask    
'''

def maskgenr(mask,currentimg,n):#what's going on with n?
    if np.array(np.shape(np.shape(mask)))[0]==3:
        mask=mask[:,:,0]
    mask2=np.zeros(np.shape(mask))
    for i in range(0,n):
        mask2=mask2+rotate(mask,360/n*i,reshape=False)
    mask=mask2    
    #get zoom factors
    ratio=(np.array(np.shape(currentimg[:,:]))).astype(float)/(np.array(np.shape(mask))).astype(float)
    mask=zoom(mask,np.max(ratio))
    diff=np.array(np.shape(mask))-np.array(np.shape(currentimg[:,:]))
    mask=mask[0:np.shape(mask)[0]-diff[0],0:np.shape(mask)[1]-diff[1]] #Note: here we assume mask produces a close to hexagon rectangle.
    mask=mask.astype('float')/np.sum(mask.astype('float'))    
    return mask

def norm(mask,cutmask):
    mask=mask.astype('float32')/np.sum(mask.astype('float32')*cutmask.astype('float32'))
    return mask

def cropcurrent(currentimg, cent, dist, angle):
    currentimg=currentimg[int(cent[1])-np.int(dist):int(cent[1])+np.int(dist), int(cent[0])-np.int(dist):int(cent[0])+np.int(dist)]
    currentimg=rotate(currentimg,angle/(2*np.pi)*360,reshape=False)
    currentimg=currentimg[int(np.shape(currentimg)[1]/2)-np.int(dist*np.sqrt(3)/2):int(np.shape(currentimg)[1]/2)+np.int(dist*np.sqrt(3)/2), int(np.shape(currentimg)[0]/2)-np.int(dist):int(np.shape(currentimg)[0]/2)+np.int(dist)]
    return currentimg

def cropcurrentbandgap(currentimg, cent, dist):
    """Crop current image (assume circle)"""
    currentimg=currentimg[int(cent[1])-np.int(dist):int(cent[1])+np.int(dist),
                          int(cent[0])-np.int(dist):int(cent[0])+np.int(dist)]
    #currentimg=rotate(currentimg,angle/(2*np.pi)*360,reshape=False)
    #currentimg=currentimg[int(np.shape(currentimg)[1]/2)-
    #    np.int(dist*np.sqrt(3)/2):int(np.shape(currentimg)[1]/2)+np.int(dist*np.sqrt(3)/2),
    #    int(np.shape(currentimg)[0]/2)-np.int(dist):int(np.shape(currentimg)[0]/2)+np.int(dist)]
    return currentimg

def optsetup(img, maskfile,cutmaskfile, n):
    cent, dist, angle=areainterest(img,6)
    currentimg=cropcurrent(img, cent, dist, angle)
    mask=maskgenr(mpimg.imread(maskfile),currentimg, n)
    cutmask=maskgenr(mpimg.imread(cutmaskfile),currentimg, 6)
    cutmask=cutmask/np.amax(cutmask)
    cutmask=((cutmask.astype(float)/(np.amax(cutmask)))>0.5).astype(float)
    return mask, cutmask, np.array([cent, dist, angle])
    
def optsetupbandgap(img, maskfile,cutmaskfile, n):
    """Recognise fibre
    
    Input
    img: input image array
    maskfile: string (.png)
    cutmaskfile: string (.png)
    n: leave as 1
    (n is used if the mask file exploits symmetry (only stores nth of mode))
    
    Output
    mask, cutmask, np.array([cent, dist, angle])
    """
    cent, dist, angle=bandgapinterest(img,6)
    #is cropping to circle fine?
    currentimg=cropcurrentbandgap(img, cent, dist)
    mask=maskgenr(mpimg.imread(maskfile),currentimg, n)
    cutmask=maskgenr(mpimg.imread(cutmaskfile),currentimg, 6)
    cutmask=cutmask/np.amax(cutmask)
    cutmask=((cutmask.astype(float)/(np.amax(cutmask)))>0.5).astype(float)
    return mask, cutmask, np.array([cent, dist, angle])

def compare(img, mask,cutmask, hexcor):
    mincorr=float(99999999)
    bestimg=[]
    #try different angles
    for angcorr in np.arange(-0.02,0.03,0.01):
        #hexcor is np.array([cent, dist, angle])
        currentimg=cropcurrent(img, hexcor[0], hexcor[1], hexcor[2]+angcorr)
        #cv2.matchTemplate compares mask & image at different positions
        #amount of padding gives number of shifted positions tried
        res=cv2.matchTemplate(
                    np.pad( norm(currentimg, cutmask), ((50,50),(50,50)), 'constant').astype('float32'),
                            norm(mask.astype('float32'), cutmask),cv2.TM_SQDIFF_NORMED) 
        if (np.amin(res)<mincorr):
            mincorr=np.amin(res)
            bestimg=currentimg
    return mincorr, bestimg
