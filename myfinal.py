#!/usr/bin/env python 

# Some suitable functions and data structures for drawing a map and particles

import time
import random
import math

# -----------------------------------------------------------
# --------------------initialize -----------------------
# -----------------------------------------------------------
import brickpi
import os


sonar_port = 3 # port which ultrasoic sensor is plugged in to
interface=brickpi.Interface()
interface.initialize()
interface.sensorEnable(sonar_port, brickpi.SensorType.SENSOR_ULTRASONIC);

motors = [0,1]
sensor_motor = [3]

interface.motorEnable(motors[0])
interface.motorEnable(motors[1])

interface.motorEnable(sensor_motor[0])


motorParams0 = interface.MotorAngleControllerParameters()
motorParams0.maxRotationAcceleration = 11.0#6.0
motorParams0.maxRotationSpeed = 26.0#12.0
motorParams0.feedForwardGain =255/20.0
motorParams0.minPWM = 18.0
motorParams0.pidParameters.minOutput = -255
motorParams0.pidParameters.maxOutput = 255
motorParams0.pidParameters.k_p = 240.0
motorParams0.pidParameters.k_i = 600
motorParams0.pidParameters.k_d = 12


motorParams1 = interface.MotorAngleControllerParameters()
motorParams1.maxRotationAcceleration = 3
motorParams1.maxRotationSpeed = 3
motorParams1.feedForwardGain =255/20.0
motorParams1.minPWM = 18.0
motorParams1.pidParameters.minOutput = -255
motorParams1.pidParameters.maxOutput = 255
motorParams1.pidParameters.k_p = 240.0
motorParams1.pidParameters.k_i = 600
motorParams1.pidParameters.k_d = 12

motorParams2 = interface.MotorAngleControllerParameters()
motorParams2.maxRotationAcceleration = 18.0
motorParams2.maxRotationSpeed = 36.0
motorParams2.feedForwardGain =255/20.0
motorParams2.minPWM = 18.0
motorParams2.pidParameters.minOutput = -255
motorParams2.pidParameters.maxOutput = 255
motorParams2.pidParameters.k_p = 240.0
motorParams2.pidParameters.k_i = 600
motorParams2.pidParameters.k_d = 12



interface.setMotorAngleControllerParameters(motors[0],motorParams0)
interface.setMotorAngleControllerParameters(motors[1],motorParams0)

interface.setMotorAngleControllerParameters(sensor_motor[0],motorParams1)
# -----------------------------------------------------------
# --------------------global variables -----------------------
# -----------------------------------------------------------
preMotorL=0
preMotorR =0
current_position=[0,0,0]
#pX = 84.0
#pY = 30.0
#pT = 0.0

bias = 3.5#3.5


def init():
    global preMotorL
    global preMotorR 
    
    preMotorL =  (interface.getMotorAngles(motors))[0][0]
    preMotorR =  (interface.getMotorAngles(motors))[1][0]
    


# Functions to generate some dummy particles data:
sE=0.15#0.5
sF=0.02#0.01#0.1
sG=0.015#0.006#0.007

def getRandomX(tmpX,D,tmpT):
    e=random.gauss(0, sE)
    return tmpX+(D +e)*math.cos(tmpT)#*bias


def getRandomY(tmpY,D,tmpT):
    e=random.gauss(0, sE)
    return tmpY+(D +e)*math.sin(tmpT)#*bias


def getRandomF(T): 
    e=random.gauss(0, sF)#/32
    return T+e

def getRandomG(T): 
    e=random.gauss(0, sG)#/32
    return T+e

# A Canvas class for drawing a map and particles:
#     - it takes care of a proper scaling and coordinate transformation between
#      the map frame of reference (in cm) and the display (in pixels)
class Canvas:
    def __init__(self,map_size=210):
        self.map_size    = map_size;    # in cm;
        self.canvas_size = 768;         # in pixels;
        self.margin      = 0.05*map_size;
        self.scale       = self.canvas_size/(map_size+2*self.margin);

    def drawLine(self,line):
        x1 = self.__screenX(line[0]);
        y1 = self.__screenY(line[1]);
        x2 = self.__screenX(line[2]);
        y2 = self.__screenY(line[3]);
        print "drawLine:" + str((x1,y1,x2,y2))

    def drawParticles(self,data):
        display = [(self.__screenX(d[0]),self.__screenY(d[1])) + d[2:] for d in data];
        print "drawParticles:" + str(display);

    def __screenX(self,x):
        return (x + self.margin)*self.scale

    def __screenY(self,y):
        return (self.map_size + self.margin - y)*self.scale
    
    def targetsToSquare(self, data, amap):
        x = data[0]
        y = data[1]
        
        amap.add_wall((x - 1, y - 1, x - 1, y + 1))
        amap.add_wall((x - 1, y - 1, x + 1, y - 1))
        amap.add_wall((x + 1, y + 1, x - 1, y + 1))
        amap.add_wall((x + 1, y + 1, x + 1, y - 1))

# A Map class containing walls
class Map:
    def __init__(self):
        self.walls = [];

    def add_wall(self,wall):
        self.walls.append(wall);

    def clear(self):
        self.walls = [];

    def draw(self):
        for wall in self.walls:
            canvas.drawLine(wall);
            
    def getWalls(self):
        return self.walls;
    
# Simple Particles set
class Particles:
    def __init__(self):
        self.n = 100;    
        print "hhhhhhhasdfasdf: ", current_position
        self.data = [current_position for i in range(self.n) ];

    def updateD(self,D,i):

        (tmpX, tmpY, tmpT) = self.data[i]
        self.data[i] = (getRandomX(tmpX,D,tmpT),getRandomY(tmpY,D,tmpT), getRandomF(tmpT));
        
# print "updateD",tmpX,tmpY,tmpT

    def updateT(self,alpha,i):

        
        #print "length",len(self.data)
        (tmpX, tmpY, tmpT) = self.data[i]
        
        
        tmpPT = getRandomG(tmpT)+alpha
                
        if tmpPT> math.pi:
            
            tmpPT-= math.pi *2
        
        if tmpPT < -math.pi:
                
            tmpPT+= math.pi *2
        

        self.data[i] = (tmpX,tmpY,tmpPT);
        #print "error", self.data[i]
    def turnright90(self,rot,i):
        
        (tmpX, tmpY, tmpT) = self.data[i]
        
        tmpT-=rot
        
                
        if tmpT> math.pi:
            
            tmpT-= math.pi *2
        
        if tmpT < -math.pi:
                
            tmpT+= math.pi *2
        self.data[i] = (tmpX,tmpY,tmpT);
        
        
    
    def draw(self):
        canvas.drawParticles(self.data);
        
    def setParticles(self, newData):
        self.data = newData
    
    def getParticles(self):
        return self.data


def Rightdeg(rot):
    #global pT
    global preMotorL
    global preMotorR 
    global particles

    para=2.545#3.12
    rot=rot*para    
    interface.increaseMotorAngleReferences(motors,[-rot,rot])
    while not interface.motorAngleReferencesReached(motors) :
        motorAngles = interface.getMotorAngles(motors)
        newMotorL = interface.getMotorAngles(motors)[0][0]        
        newMotorR = interface.getMotorAngles(motors)[1][0]
        deltaL = newMotorL -preMotorL
        deltaR = newMotorR -preMotorR
                
        alpha = (( deltaL - deltaR)/16.0)*math.pi #*0.202798076678*1.2* math.pi 
        
        for i in range(100):
            #print "count  :",i
            particles.updateT(alpha,i)
        
        #pT += alpha
        #if pT> math.pi:
            #pT-= math.pi *2
        #if pT < -math.pi:
            #pT+= math.pi *2
        
        #if motorAngles :
            #print "Motor angles: ", motorAngles[0][0], ", ", motorAngles[1][0]
            #time.sleep(0.1)
        particles.draw()
        
        preMotorL = newMotorL
        preMotorR = newMotorR
        #time.sleep(0.0004)
    #print "partial destination reached"


def forward(angle):
    #global pX
    #global pY
    global preMotorL
    global preMotorR 
    global particles
    
    angle=angle*0.289#0.2855
    interface.increaseMotorAngleReferences(motors,[angle,angle])
    while not interface.motorAngleReferencesReached(motors) :
        motorAngles = interface.getMotorAngles(motors)
        newMotorL = interface.getMotorAngles(motors)[0][0]        
        newMotorR = interface.getMotorAngles(motors)[1][0]
        deltaL = newMotorL -preMotorL
        deltaR = newMotorR -preMotorR
        
        D = (deltaL+deltaR)/2*bias
        
        for i in range(100):                
                particles.updateD(D,i)
                
                
        #pX += D* math.cos(pT)
        #pY += D* math.sin(pT)
        
        #print pX,pY
        #if motorAngles:
# print "Motor angles: ", motorAngles[0][0], ", ", motorAngles[1][0]
            #time.sleep(0.1)
    
        particles.draw()
        
        preMotorL = newMotorL
        preMotorR = newMotorR
        #time.sleep(0.0004)
    
    #print "destination reached"




def getRobotInfoResampling():
    reading = interface.getSensorValue(sonar_port)
    if reading:
        z = reading[0]#interface.getSensorValue(sonar_port)[0]
    else:
        print "no reading"
    #z = interface.getSensorValue(sonar_port)[0]
    
    if not z or z< 5 or z > 155:
        #print "ultrosonic sensor location failed"
        #print "current position xyt:", mean()
        return mean()
    z=z+1
    print "z",z


            
    ZeroCount=0
    likelihoodArray = []
    
    for p in particles.getParticles():
        
        likelihood = calculate_likelihood(p[0],p[1],p[2],z)
        
        likelihoodArray.append(likelihood)
        if likelihood == 0:
            ZeroCount+= 1
            
    #print "zero",ZeroCount  
    
    
    if ZeroCount > 70:
        #print "Too many particles angle >40"        
        #print "current position xyt: ",mean()
        return mean()
    
    
    #print "mean",mean()  

    
    
    resampling(likelihoodArray) 
    particles.draw()
    return mean()
            
def mean():
    data = particles.getParticles()
    meanX = 0            
    meanY = 0
    meanT = 0
    
    tempMeanTsin = 0
    tempMeanTcos = 0
    
    for i in range(len(data)):
        tmpdata = data[i]
        
        meanX += tmpdata[0]
        meanY += tmpdata[1]
        #meanT += tmpdata[2]
        
        tempMeanTsin += math.sin(tmpdata[2])
        tempMeanTcos += math.cos(tmpdata[2])
    
    tempMeanTsin /= 100
    tempMeanTcos /= 100
    
    meanT = math.atan2(tempMeanTsin, tempMeanTcos)
    
    #print "meanT: ", meanT
    
    pX = meanX/100.0
    pY = meanY/100.0
    pT = meanT
    #pT = meanT
    
    #pT = normalRadianToSystemRadian(meanT)
    #print "mean",pX,pY,pT
    return (pX, pY, pT)
        


def calculate_likelihood( x,  y,  theta,  z):
    [m, beta] = Faced_Wall(x, y, theta)
    if m > 143: 
        return 0
    elif m > 140 and beta > (math.pi* 8/180):
        return 0
    elif beta > (math.pi* 15/180) and m > 100: 
        return 0
    elif  beta > (math.pi* 20/180) and m > 80:
        return 0
    elif  beta > (math.pi* 23/180) and m > 60:
        return 0
    elif  beta > (math.pi* 25/180) and m > 40:
        return 0
    elif  beta > (math.pi* 38/180) and m > 20:
        return 0
    elif  beta==-100:
        return 0
    else : 
        return Likehood(z,m)



    
    
def Faced_Wall(x0,  y0,  pT):
    
    z = interface.getSensorValue(sonar_port)[0]
    z=z+1

    #print "zzz",z
    
    m = 256
    beta=-100 #error
    for w in mymap.getWalls():
        x1 = w[0]
        y1 = w[1]
        x2 = w[2]
        y2 = w[3]
        if ((y2-y1)*math.cos(pT)-(x2-x1)*math.sin(pT))== 0:
            formula = -1#math.fabs((x2-x1)*(y1-y0)-(x1-x0)*(y2-y1))/math.sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1))
            
        else:
            formula = ((y2-y1)*(x1 - x0) - (x2 - x1)*(y1-y0))/((y2-y1)*math.cos(pT)-(x2-x1)*math.sin(pT))
        

        XX = x0 + formula*math.cos(pT)
        YY = y0 + formula*math.sin(pT)
        
        
        if XX <= max(x1,x2)+0.01 and XX>=min(x1,x2)-0.01 and YY<=max(y1,y2)+0.01 and YY>=min(y1,y2)-0.01 and formula< m and formula>0 and formula<z+15:       
            m = formula
            #beta = math.acos( (math.cos(pT)*(y1-y2) + math.sin(pT)*(x2-x1) ) / math.sqrt(pow(y1-y2,2) +pow(x2-x1,2) ))
            

            beta = math.acos( (math.cos(pT)*(y1-y2) + math.sin(pT)*(x2-x1) ) / math.sqrt(pow(y1-y2,2) +pow(x2-x1,2) ))
            
            #print x1, x2, y1, y2


    return m,beta


def Likehood(z,m):
    variance = 1.5#0.45    #system error
    k = 0.002#0.005          #garbage value fraction
    para=1
    gs = math.exp( -pow((z-m), 2)/(2*variance) ) + k
    '''
    if gs < k:
        return k
    else:
    '''
    return para*gs
        

def resampling(likelihoodArray):
    tmp_sum = sum(likelihoodArray)
    
    weight = []
    tmp_likelihood = 0
    for i in range(len(likelihoodArray)):
        tmp_likelihood+=likelihoodArray[i]
        tmp_likelihood=tmp_likelihood
        weight.append(tmp_likelihood)
        
    indexArray=[]



    for s in range(100):
    #while len(indexArray)!=100:  
        rd = random.random()
        rd *= tmp_sum
        for i in range(len(weight)):
            if i ==0:
                if rd < weight[i]:
                    indexArray.append(i);
                    break         
            else:
                if  weight[i-1] < rd and rd< weight[i]:
                    indexArray.append(i);
                    break  
                    
    

                
    # resample 100 particles with indexArray
    newData = []
    data = particles.getParticles()
    
    for i in range(len(indexArray)):
        
        newData.append(data[indexArray[i]])
    
    #print "newData.length : ", len(newData)
    particles.setParticles(newData)


    
    
def checkIfReached(point, target):
    x1 = point[0] - target[0]
    y1 = point[1] - target[1]
    
    distance = math.pow(math.pow(x1, 2) + math.pow(y1, 2), 0.5)
    
    print "current point and target DDD ", distance
    if distance<5:
        print "reached"
        
    return distance < 5

def navigateToWaypoint(x, y): 
    
    #global xnow,ynow,thetanow
    print "123current position: ",mean()

    
    (x_now,y_now,t_now)=mean()
    dx = x-x_now
    dy = y-y_now

    d = math.sqrt(math.pow(dx,2) + math.pow(dy,2))



    if dx == 0 and dy > 0:
        angle = math.pi/2
    elif dx == 0 and dy < 0:
        angle = -math.pi/2
    elif dx == 0 and dy == 0:
        angle = 0
    else:
        angle = math.atan(dy/float(dx))

    if dx<0:
        angle = angle + math.pi

    angle=angle-t_now
    
    if angle>math.pi:
        angle-=2*math.pi
    if angle<-math.pi:
        angle+=2*math.pi
    
    
    
    #print math.atan(1)                                                                              
    #print d                                                                                         
    #print angle 
    
    Rightdeg(-angle)

    if x==84 and y==30:
        if d>90:
            forward(90)
        else:
            forward(d)
    elif x==138 and y==168:
        if d>65:
            forward(65)
        else:
            forward(d)        
    else:
        if d>45:
            forward(45)
        else:
            forward(d)



    
    
    
    if (x==84 and y==30) or (x==138 and y==54):
        if x==84 and y==30:
        
            fast_Sensor_Rightdeg(1.1978)
            for i in range(100):
                    particles.turnright90(1.1978,i)
            getRobotInfoResampling()

            fast_Sensor_Rightdeg(-1.1978)
            for i in range(100):
                    particles.turnright90(-1.1978,i)
        else:
            
            fast_Sensor_Rightdeg(-0.5*math.pi)
            for i in range(100):
                    particles.turnright90(-0.5*math.pi,i)
            getRobotInfoResampling()

            fast_Sensor_Rightdeg(1*0.5*math.pi)
            for i in range(100):
                    particles.turnright90(1*0.5*math.pi,i)
            
    else:                
        fast_Sensor_Rightdeg(0.5*math.pi)
        for i in range(100):
                particles.turnright90(0.5*math.pi,i)
        getRobotInfoResampling()

        fast_Sensor_Rightdeg(-1*0.5*math.pi)
        for i in range(100):
                particles.turnright90(-1*0.5*math.pi,i)
    
    
    '''
    fast_Sensor_Rightdeg(-2*0.5*math.pi)
    for i in range(100):
            particles.turnright90(-2*0.5*math.pi,i)
    getRobotInfoResampling()
   
    
     
    fast_Sensor_Rightdeg(0.5*math.pi)
    for i in range(100):
            particles.turnright90(0.5*math.pi,i)
    getRobotInfoResampling()
    
    
        
    fast_Sensor_Rightdeg(1*0.5*math.pi)
    for i in range(100):
            particles.turnright90(1*0.5*math.pi,i)
    
    '''
    
    
    
    
    
    
    
    
    
    
    
    print "target position: ",x,y,d,angle


    print "current position: ",mean()
    
    #z = interface.getSensorValue(sonar_port)[0]
    #print "z",z
    #angle = float(input("Enter a angle to rotate (in radians): "))
    #time.sleep(10)
    #xnow=x
    #ynow=y
    #thetanow=angle+thetanow

    
    
    
#!/usr/bin/env python
# By Jacek Zienkiewicz and Andrew Davison, Imperial College London, 2014
# Based on original C code by Adrien Angeli, 2009



# -----------------------------------------------------------
# --------------------global variables -----------------------
# -----------------------------------------------------------
def fast_Sensor_Rightdeg(rot):
    interface.setMotorAngleControllerParameters(sensor_motor[0],motorParams2)

    para=1.06#3.12
    rot=rot*para    
    interface.increaseMotorAngleReferences(sensor_motor,[-rot])
    while not interface.motorAngleReferencesReached(sensor_motor) :
        time.sleep(0.00001)
        
def Sensor_Rightdeg(rot):
    #global pT
    #global preSensorMotor
    interface.setMotorAngleControllerParameters(sensor_motor[0],motorParams1)
    distanceArray = []
    count=0
    
    para=1.06#3.12
    rot=rot*para    
    interface.increaseMotorAngleReferences(sensor_motor,[-rot])
    
    preSensorMotor =  (interface.getMotorAngles(sensor_motor))[0][0]
    #print "pre", preSensorMotor
    
    #temp=0.0
    while not interface.motorAngleReferencesReached(sensor_motor) :
        newSensorMotor = interface.getMotorAngles(sensor_motor)[0][0]  
    
        delta=math.fabs(preSensorMotor-newSensorMotor)
        #print "delat",delta

        #print "delatdiff",delta-temp
        
        
        if delta>(math.pi/180+(count/180.0)*math.pi) and count<=359:
            reading = interface.getSensorValue(sonar_port)
            if reading:
                distance = reading[0]#interface.getSensorValue(sonar_port)[0]
            else:
                print "no reading"
                
            #print "distance: ", distance
            
            #print distance,(math.pi/180+(count/180.0)*math.pi)/math.pi*180,count
            
            distanceArray.append(distance)
            #print interface.getSensorValue(sonar_port)[0]
            
            count+=1
        
        '''
        if int(round(delta / (math.pi * 180))) - 1 < 360:

            print interface.getSensorValue(sonar_port)[0]
            distanceArray[int(round(delta / (math.pi * 180))) - 1] =interface.getSensorValue(sonar_port)[0]
            print distanceArray[int(round(delta / (math.pi * 180))) - 1]
            #time.sleep(1)
            
        time.sleep(0.001)  
        '''
    #print "distances", distanceArray
    #print "lendistances", len(distanceArray)
    return distanceArray

    
    
    
# Location signature class: stores a signature characterizing one location
class LocationSignature:
    def __init__(self, no_bins = 360):
        self.sig = [0] * no_bins
        
    def print_signature(self):
        for i in range(len(self.sig)):
            print self.sig[i]

# --------------------- File management class ---------------
class SignatureContainer():
    def __init__(self, size = 5):
        self.size      = size; # max number of signatures that can be stored
        self.filenames = [];
        
        # Fills the filenames variable with names like loc_%%.dat 
        # where %% are 2 digits (00, 01, 02...) indicating the location number. 
        for i in range(self.size):
            self.filenames.append('loc_{0:02d}.dat'.format(i))

    # Get the index of a filename for the new signature. If all filenames are 
    # used, it returns -1;
    def get_free_index(self):
        n = 0
        while n < self.size:
            if (os.path.isfile(self.filenames[n]) == False):
                break
            n += 1
            
        if (n >= self.size):
            return -1;
        else:    
            return n;

    # Delete all loc_%%.dat files
    def delete_loc_files(self):
        print "STATUS:  All signature files removed."
        for n in range(self.size):
            if os.path.isfile(self.filenames[n]):
                os.remove(self.filenames[n])
            
    # Writes the signature to the file identified by index (e.g, if index is 1
    # it will be file loc_01.dat). If file already exists, it will be replaced.
    def save(self, signature, index):
        filename = self.filenames[index]
        if os.path.isfile(filename):
            os.remove(filename)
            
        f = open(filename, 'w')

        for i in range(len(signature.sig)):
            s = str(signature.sig[i]) + "\n"
            f.write(s)
        f.close();

    # Read signature file identified by index. If the file doesn't exist
    # it returns an empty signature.
    def read(self, index):
        ls = LocationSignature()
        filename = self.filenames[index]
        if os.path.isfile(filename):
            f = open(filename, 'r')
            for i in range(len(ls.sig)):
                s = f.readline()
                if (s != ''):
                    ls.sig[i] = int(s)
            f.close();
        else:
            print "WARNING: Signature does not exist."
        
        return ls
        
# FILL IN: spin robot or sonar to capture a signature and store it in ls
def characterize_location(ls):
    print "TODO:    You should implement the function that captures a signature."

    alist = []
    alist.extend(Sensor_Rightdeg(2*math.pi))
    
    #print "alist", alist
    #print len(alist)
    fast_Sensor_Rightdeg(-2*math.pi)
    #Sensor_Rightdeg(-2*math.pi)
    #Sensor_Rightdeg(-math.pi)
    #alist.extend(Sensor_Rightdeg(math.pi))
    
    for i in range(len(ls.sig)):
        ls.sig[i] = alist[i]

        
def frequency(ls):
    epsilon = 10
    alist = [0]*int(300/epsilon)
    
    
    
    '''
    for i in ls.sig:
        alist.append(ls.sig.count(i))
    blist = list(set(zip(ls.sig,alist)))
    
    print blist
    
    clist = []
    '''
    for n in range(0,int(300/epsilon)):    
        for j in ls.sig:
            if j>n*epsilon and j<=(n+1)*epsilon:
                alist[n] +=1
    #print alist
    return alist


    #d, f = zip(*blist)
    #fL = list(f)
    #fsum = sum(fL)*1.0
    #print fsum
    #nfL = [round(e/fsum*100,2) for e in fL]
    #elist = zip(d,nfL)
    #return elist
    
    
            
        
# FILL IN: compare two signatures
def compare_signatures(ls1, ls2):
    #dist = 0
    dist1 = 0.0
    tempdist1 =0.0
    #Ls1 = ls1.sig#frequency(ls1)
    #Ls2 = ls2.sig#frequency(ls2)
    '''for i in range(len(Ls1)):
                #print i[0],j[0]
                #print i[1],j[1]
        dist += pow(Ls1[i]-Ls2[i],2)'''
    fls1 = frequency(ls1)
    fls2 = frequency(ls2)
    #return dist
    
    for i in range(len(fls1)):
        
        #if fls2[i]==0:
            #continue
        #if pow(fls1[i]-fls2[i],2)>10 :
        
        #dist1+=math.fabs(fls1[i]-fls2[i])
    
        dist1+=pow(fls1[i]-fls2[i],2)
    
    '''
    for i in fls1:
        for j in fls2:
            if i[0] == j[0]:
                print i[0],j[0]
                print i[1],j[1]
                dist1 += pow(i[1]-j[1],2)
                '''
    #print dist1
    return dist1

# This function characterizes the current location, and stores the obtained 
# signature into the next available file.
def learn_location():
    ls = LocationSignature()
    characterize_location(ls)
    idx = signatures.get_free_index();
    if (idx == -1): # run out of signature files
        print "\nWARNING:"
        print "No signature file is available. NOTHING NEW will be learned and stored."
        print "Please remove some loc_%%.dat files.\n"
        return
    
    signatures.save(ls,idx)
    print "STATUS:  Location " + str(idx) + " learned and saved."

# This function tries to recognize the current location.
# 1.   Characterize current location
# 2.   For every learned locations
# 2.1. Read signature of learned location from file
# 2.2. Compare signature to signature coming from actual characterization
# 3.   Retain the learned location whose minimum distance with
#      actual characterization is the smallest.
# 4.   Display the index of the recognized location on the screen
def recognize_location():
    ls_obs = LocationSignature();
    characterize_location(ls_obs);
    
    
    
    dist=0
    # FILL IN: COMPARE ls_read with ls_obs and find the best match
    alist = []
    for idx in range(signatures.size):
        print "STATUS:  Comparing signature " + str(idx) + " with the observed signature."
        ls_read = signatures.read(idx);
        #print "ls_read.sig: ",ls_read.sig
        dist    = compare_signatures(ls_obs, ls_read)
        
        
        print "dist" + str(idx) + ": " + str(dist)
        alist.append(dist)
        
    print "smalletst dist: " + str(min(alist)) + "index: " + str(alist.index(min(alist)))
    
    min_index=alist.index(min(alist))
    
    
    if min_index==1 or min_index==2:
        check=frequency(ls_obs)
        if check[6]>10:
            min_index=2
        else:
            min_index=1
    ls_read = signatures.read(alist.index(min(alist)));
    print "min_index", min_index

    #print "ls_obs.sig: ",ls_obs.sig
    
    #print "ls_read.sig: ",ls_read.sig


    
    #Sensor_Rightdeg(-0.5*math.pi)

    '''
    ls_obs1 = LocationSignature();
    characterize_location(ls_obs1);
    print "ls_obs.sig1",ls_obs1.sig
    
    
    print "111",frequency(ls_obs)
    print "222",frequency(ls_obs1)


    
    disttemp=0.0
    disttemp    = compare_signatures(ls_obs, ls_obs1)
    
    print "disttemp",disttemp
    
    #print "ls_read.sig",ls_read.sig
    
    
    
    ''' 
    #Sensor_Rightdeg(0.5*math.pi)
    #forward(2)
    '''
    ls_obs1 = LocationSignature();
    characterize_location(ls_obs1);
    print "ls_obs.sig1",ls_obs1.sig
    
    
    dist5=0
    templist = []
    
    new_obs_list = []
    temp=0
    for i in range(len(ls_obs.sig)):
        dist5=0
        new_obs_list= ls_obs1.sig[i:]+ls_obs1.sig[:i]
        
        #print "new_obs_list: ",new_obs_list
        
        
        for m in range(len(ls_obs1.sig)):
            temp=pow(new_obs_list[m]-ls_obs.sig[m],2)
            #if temp>5:
            dist5+=pow(new_obs_list[m]-ls_obs.sig[m],2)
        templist.append(dist5)
    print  "blist",templist
    print  "len",len(templist)
    print "min",min(templist)
    print "min_index",templist.index(min(templist))
    
    
    
    
    
    
    
    '''
    
    dist2=0
    blist = []
    
    new_obs_list = []
    temp=0
    for i in range(len(ls_obs.sig)):
        dist2=0
        new_obs_list= ls_obs.sig[i:]+ls_obs.sig[:i]
        
        #print "new_obs_list: ",new_obs_list
        
        
        for m in range(len(ls_obs.sig)):
            temp=pow(new_obs_list[m]-ls_read.sig[m],2)
            #if temp>10:
            dist2+=temp
        blist.append(dist2)
    print  "blist: ",blist
    #print  "len: ",len(blist)
    print "min: ",min(blist)
    print "min_index: ",blist.index(min(blist))
    
    print "ls_read",ls_read.sig
    return [min_index,blist.index(min(blist))]

    #return [0,0]
# Prior to starting learning the locations, it should delete files from previous
# learning either manually or by calling signatures.delete_loc_files(). 
# Then, either learn a location, until all the locations are learned, or try to
# recognize one of them, if locations have already been learned.
#Sensor_Rightdeg(-math.pi)
#Sensor_Rightdeg(math.pi)

#a = [1,2,3,1,2,3,1,2,3,4,4,4,4,4,4,5,5,5,5,53,42,34,3,3,3,3,32,2,3,4,23,4,1]
#frequency1(a)

#signatures = SignatureContainer();

#signatures.delete_loc_files()

#learn_location();
#recognize_location();


    
    
    
####################################main###########################################


init()

canvas = Canvas();    # global canvas we are going to draw on

mymap = Map();
mypoint= Map();
# Definitions of walls
# a: O to A
# b: A to B
# c: C to D
# d: D to E
# e: E to F
# f: F to G
# g: G to H
# h: H to O
mymap.add_wall((0,0,0,168));        # a
mymap.add_wall((0,168,84,168));     # b
mymap.add_wall((84,126,84,210));    # c
mymap.add_wall((84,210,168,210));   # d
mymap.add_wall((168,210,168,84));   # e
mymap.add_wall((168,84,210,84));    # f
mymap.add_wall((210,84,210,0));     # g
mymap.add_wall((210,0,0,0));        # h



t = 0;

targetPoints = [ (84, 30),
                 (180, 30),
                 (180, 54),
                 (138, 54),
                 (138, 168)#,
                 #(138,30)
                  ]







for point in targetPoints:
    canvas.targetsToSquare(point, mypoint)
mymap.draw();
mypoint.draw();




signatures = SignatureContainer();

#signatures.delete_loc_files()

#learn_location();
startposition=[0,0]
startposition=recognize_location();
print "start",startposition
point_index=startposition[0]

start_point_x= targetPoints[startposition[0]][0]
start_point_y= targetPoints[startposition[0]][1]
start_point_th= startposition[1]/180.0*math.pi


print "th",start_point_th
targetPoints=targetPoints[point_index+1:]+targetPoints[:point_index+1]


current_position=(start_point_x,start_point_y,start_point_th)




particles = Particles();

print targetPoints

for (x, y) in targetPoints:
    

    #print "target: ",x,y
    #print "mean position: ",mean()
    (curX, curY, theta) = mean()#getRobotInfoResampling()
    
    
    
    while not checkIfReached((curX, curY), (x, y)):
        
        #print "target: ",x,y
        #print "mean position: ",mean()
        navigateToWaypoint(x, y)
        (curX, curY, theta) = getRobotInfoResampling()
    
      
    
    
    
    
    
    
    
    
    

    
    
    
    
    
    
    

    
    
    
    
    
    
    
    
    
    
    
    