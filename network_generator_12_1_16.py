#1d/2d percolation
#This script randomly generates a network of nanowires and stores the data representing it in a few files. The basic idea is to generate each segment with 3 random variables - 2 for the position, one for the orientation. The user specifies the space size, piece size, and number of pieces. The generator also has several modes for modifying nanowire orientation, which the user can choose from and specify a degree of strength for. Lastly, the model checks if a percolating path exists across the network - that is, if it is possible to move across the network in a continuous path of nanowires. See my research papers and presentations for more info.
#Milind Jagota

from scipy import random as rdom
import random
import math
import pylab
import json

print("Welcome to nanowire network generator. \nEnter prompted information and network will be generated and tested for percolation. Image of the network with percolating cluster in red will be saved to project folder. \n")

#stores all the segments, sorted into "zones" - zones are used to shorten the calculation of finding all intersections - we only need to check segments in adjacent zones
y=[]

#stores all the segments
x=[]

#stores all the intersections, where each is a pair of the index of the first segment and the index of the second segments
intersects=[]

#stores all the intersection locations, where each is an x-y coordinate pair
intersect_locations=[]

#use to find and store clusters of indirectly connected segments
cluster=[]

print("Modes: \n 0: Random position, random orientation \n 1: Random position, bimodal orientation \n 2: Random position, restricted uniform orientation (horizontal) \n 3: Random position, restricted uniform orientation (verticacal) \n 4: Random position, Gaussian orientation \n")
dropmode = -1
ofac = -1
while (dropmode < 0 or dropmode > 4):
      dropmode=int(input("Choose mode: "))
pf=0
if (dropmode == 1 or dropmode == 2 or dropmode == 3 or dropmode == 4):
    while (ofac < 0 or ofac > 1):
        ofac = float(input("Enter factor of orienation: "))
sidelengthx=float(input("Box width? "))*(3.5e-05)
sidelengthy=float(input("Box height? "))*(3.5e-05)
numpieces=int(input("How many segments? "))
piecelength=float(input("Segment length? "))*(3.5e-05)


if math.fmod(sidelengthx,piecelength)==0:
    zonesrow=sidelengthx/piecelength
else:
    zonesrow=math.floor(sidelengthx/piecelength)+1
if math.fmod(sidelengthy,piecelength)==0:
    zonescol=sidelengthy/piecelength
else:
    zonescol=math.floor(sidelengthy/piecelength)+1
zonestotal=int(zonesrow*zonescol)
for a in range (0, zonestotal):
    y.append([])
    
percolation=False

#generates a single random segment and adds its data to the list of segments, with the orienation distribution based on the user input.
def drop_segment(x,y):
    if dropmode==0 or dropmode==6:
        c=(random.uniform(0,1)* math.pi)-((math.pi)/2)
        a=random.uniform(0,sidelengthx)
        b=random.uniform(0,sidelengthy)
    if dropmode==1:
        a=random.uniform(0,sidelengthx)
        b=random.uniform(0,sidelengthy)
        e=random.randint(0,1)
        if e==0:
            c=math.pi*ofac
        elif e==1:
            c=-math.pi*ofac
    if dropmode==2:
        a=random.uniform(0,sidelengthx)
        b=random.uniform(0,sidelengthy)
        c=((random.uniform(0,1)*math.pi)-((math.pi)/2))*ofac
    if dropmode==3:
        a=random.uniform(0,sidelengthx)
        b=random.uniform(0,sidelengthy)
        c=(random.uniform(0,1)*math.pi*ofac)+(math.pi*(1-ofac)/2.0)
    if dropmode==4:
        a=random.uniform(0,sidelengthx)
        b=random.uniform(0,sidelengthy)
        c=random.gauss(0,ofac*math.pi)
    d=math.tan(c)
    distance=((a**2)+(b**2))**0.5
    zone=int(math.floor(a/piecelength)+(zonesrow*math.floor(b/piecelength)))
    if dropmode==6:
        g=piecelength*random.uniform(1-pf,1+pf)
    else:
        g=piecelength
    x.append([distance,a,b,d,zone,0,c,g])
    return x

#determines if an intersection exists between two segments
def intersect(a,b):
    if a[2]==b[2]:
        return(False)
    else:
        z=((a[2]*a[0])-a[1]-(b[2]*b[0])+b[1])/(
            a[2]-b[2])
        y=(a[2]*z)-(a[2]*a[0])+a[1]
        if ((z-a[0])**2+(y-a[1])**2)**.5<=.5*a[6] and ((z-b[0])**2+(y-b[1])**2)**.5<=.5*b[6]:
            if z>sidelengthx or z<0 or y>sidelengthy or y<0:
                return(False)
            else:
                intersect_locations.append([z,y])
                return(True)
        else:
            return(False)

#finds and stores all intersections between segments in the network
def find_intersects(x,y):
    for a in range (0,len(y)):
        for b in range (0,len(y[a])-1):
            for c in range(b+1,len(y[a])):
                if intersect(x[y[a][b]],x[y[a][c]])==True:
                    intersects.append([y[a][b],y[a][c],1,1])
    for a in range (0,len(y)-1):
        for b in range (a,len(y)):
            if math.fabs(a-b)==1 or math.fabs(a-b)==zonesrow or math.fabs(a-b)==zonesrow-1 or math.fabs(a-b)==zonesrow+1:
                for c in range(0,len(y[a])):
                    for d in range(0,len(y[b])):
                        if intersect(x[y[a][c]],x[y[b][d]])==True:
                            intersects.append([y[a][c],y[b][d],1,1])

#Using the cluster variable, which is always initialized to have one segment, finds all segments that are indirectly connected to that cluster, using x, which here is the list of intersections
def find_cluster(a,x):
    for b in cluster:
        for c in x:
            if c[3]==1:
                if c[0]==b[0] or c[0]==b[1] or c[1]==b[1] or c[1]==b[0]:
                    cluster.append(c)
                    c[2]=0
                    c[3]=0

#this function plots the entire network, with the percolating cluster in red if there is one, and all others in blue
def endpoints(x):
    for a in range(0,len(x)):
        plotx=[]
        ploty=[]
        f=-((math.cos(math.atan(x[a][2])))*(x[a][6]/2))+(x[a][0])
        g=((math.cos(math.atan(x[a][2])))*(x[a][6]/2))+(x[a][0])
        plotx.append(f)
        plotx.append(g)
        h=((math.sin((math.atan(x[a][2]))))*(x[a][6]/2))+(x[a][1])
        i=-((math.sin((math.atan(x[a][2]))))*(x[a][6]/2))+(x[a][1])
        ploty.append(i)
        ploty.append(h)
        if x[a][4]==1:
            pylab.plot(plotx,ploty,'r')
        else:
            pylab.plot(plotx,ploty,'b')
        pylab.ylim([0-((2*piecelength)/3),sidelengthy+((2*piecelength)/3)])
        pylab.xlim([0-((2*piecelength)/3),sidelengthx+((2*piecelength)/3)])
        pylab.axes().set_aspect('equal')


#The main drops the correct number of segments, finds all intersections, then generates every cluster and checks if any of them percolate. It returns true if any do
for a in range (0,numpieces):
    drop_segment(x,y)
x.sort()
for a in range(0,len(x)):
    del x[a][0]
    y[x[a][3]].append(a)
f=open('segments.txt','w')
json.dump(x,f)
find_intersects(x,y)

#generate every possible cluster
for a in range(0,len(intersects)):
    if percolation==True:
        break
    cluster=[intersects[a]]
    if intersects[a][2]==1:
        #build the cluster
        find_cluster(a,intersects)
        #chuck if the cluster crosses both the left and right borders or both the top and bottom borders
        for b in cluster:
            if percolation==True:
                break
            if sidelengthx>=sidelengthy:
                #check if any segment in the cluster crosses the left border
                if -((math.cos((math.atan(x[b[0]][2]))))*(x[b[0]][6]/2))+(x[b[0]][0])<=0 or -((math.cos((math.atan(x[b[1]][2]))))*(x[b[1]][6]/2))+(x[b[1]][0])<=0:
                    for c in cluster:
                        #check if any segment in the cluster crosses the right border
                        if ((math.cos((math.atan(x[c[0]][2]))))*(x[c[0]][6]/2))+(x[c[0]][0])>=sidelengthx or ((math.cos((math.atan(x[c[1]][2]))))*(x[c[1]][6]/2))+(x[c[1]][0])>=sidelengthx:
                            percolation=True
                            for d in cluster:
                                x[d[0]][4]=1
                                x[d[1]][4]=1
                            break
            if sidelengthy>=sidelengthx:
                #check if any segment in the cluster crosses the bottom border
                if ((math.sin((math.atan(x[b[0]][2]))))*(x[b[0]][6]/2))+(x[b[0]][1])<=0 or -((math.sin((math.atan(x[b[0]][2]))))*(x[b[0]][6]/2))+(x[b[0]][1])<=0 or ((math.sin((math.atan(x[b[1]][2]))))*(x[b[1]][6]/2))+(x[b[1]][1])<=0 or -((math.sin((math.atan(x[b[1]][2]))))*(x[b[1]][6]/2))+(x[b[1]][1])<=0:
                    for c in cluster:
                        #check if any segment in the cluster crosses the top border
                        if ((math.sin((math.atan(x[c[0]][2]))))*(x[c[0]][6]/2))+(x[c[0]][1])>=sidelengthy or -((math.sin((math.atan(x[c[0]][2]))))*(x[c[0]][6]/2))+(x[c[0]][1])>=sidelengthy or ((math.sin((math.atan(x[c[1]][2]))))*(x[c[1]][6]/2))+(x[c[1]][1])>=sidelengthy or -((math.sin((math.atan(x[c[1]][2]))))*(x[c[1]][6]/2))+(x[c[1]][1])>=sidelengthy:
                            percolation=True
                            for d in cluster:
                                x[d[0]][4]=1
                                x[d[1]][4]=1
                            break

if (percolation):
    print("Network percolates")
else:
    print("Network doesn't percolate")

g=open('intersects.txt','w')
json.dump(intersects,g)
h=open('intersect_locations.txt','w')
json.dump(intersect_locations,h)
v=open('parameters.txt','w')
v.write(str(piecelength)+'\n')
v.write(str(sidelengthx))
endpoints(x)
pylab.plot([0,sidelengthx,sidelengthx,0,0],[0,0,sidelengthy,sidelengthy,0],"k")
pylab.savefig("plot.pdf")
#pylab.show()
f.close()
g.close()
h.close()
v.close()
