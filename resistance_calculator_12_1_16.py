#resistance calculator for random networks
#This script analyzes a network of nanowires, stored in the files generated by the network generator, and calculates its conductivity. The basic idea is to write an electrical equation at each node using Kirchhoff's Law, then solve this large matrix equation to get voltages across the matrix and determine network conductivity. Implementing this is pretty complicated, especially since contact resistance has to be incorporated. See my research papers and presentations for more info.
#Milind Jagota

import math
import sys
import pylab
import matplotlib.cm as cm
import time
from scipy import linalg
import json
import os

print("Welcome to nanowire network resistance calculator. \nMake sure to run the network generator so that this script has a network to analyze. Do not change any of the created text files from the network generator. The calculator will fail or return non-sensical values if there is no x-direction percolation in the passed network. \n")

execute = ""
while (execute != "Y"):
    execute = str(raw_input("Enter Y to continue, N to quit: "))
    if execute == "N":
        print ("Exiting...")
        sys.exit()

if (not os.path.isfile('./segments.txt')):
    print("Segments file is missing. Run the network generator and try again. ")
    sys.exit();
if (not os.path.isfile('./intersects.txt')):
    print("Intersects file is missing. Run the network generator and try again. ")
    sys.exit();
if (not os.path.isfile('./intersect_locations.txt')):
    print("Intersect locations file is missing. Run the network generator and try again. ")
    sys.exit();
if (not os.path.isfile('./parameters.txt')):
    print("Parameters file is missing. Run the network generator and try again. ")
    sys.exit();

print("Choose diagram to generate: \n 0: No diagram \n 1: Map of voltages at each node (10V max, 0V min) \n 2: Map of current traveling through each wire (Red is highest, then blue, then yellow) \n 3: Map of voltage drop across the junction at each node \n 4: Voltage of node vs x-coordinate scatterplot \n")
diagramChoice = -1
while (diagramChoice < 0 or diagramChoice > 4):
    diagramChoice=int(input("Choose diagram: "))

e=open('segments.txt','r')
asdfgh=open('intersects.txt','r')
g=open('intersect_locations.txt','r')
v=open('parameters.txt','r')
segments=json.load(e)
intersects=json.load(asdfgh)
intersect_locations=json.load(g)
piecelength1=v.readline()
sidelength1=v.readline()
piecelength=float(piecelength1[0:-1])
sidelength=float(sidelength1)
plotx=[]
ploty=[]
upper=[]
app_voltage = 10
junctionresistance=1500

#will hold the voltage at every junction
voltages=[]
#will hold the junctions each segment connects to
junctions=[[] for x in range(0,len(segments))]
#will hold the indices of junctions that are to the left of the left border
zero_crossers=[]
#will hold the indices of junctions that are to the right of the right border
far_side_crossers=[]

#Calculates distance between two points
def distance_formula(a,b,c,d):
    #a and b are x and y coordinates of first point, c and d are x and y coordinates of second point
    distance = ((float((c-a))**2)+(float((d-b))**2))**.5
    return distance

def find_endpoints(a,b,c,d):
    #returns endpoints of segement given x of center(a), y of center(b), theta(c), and length(d)
    alpha=(math.cos(c)*d)/2
    beta=(math.sin(c)*d)/2
    endpoints=[(a-alpha),(b-beta),(a+alpha),(b+beta)]
    return endpoints

#Takes list of intersects and segments and splits them into new list of segments. Each new segment starts
#at an intersection and ends at either an intersection or the end of a rod. The multidimensional array is
#stored as each node along with start and end points of each of four segments going into the node. It also contains the index of the junction each segment comes from, for a total of 20 values per junction
def split_segs(a,b,c,d):
    #a is list of intersects, b is list of intersect locations, c is list of segments, d is sidelength of box

    #build junction list as described at top
    for x in range(0,len(c)):
        for y in a:
            if y[0]==x:
                junctions[x].append(a.index(y))
            elif y[1]==x:
                junctions[x].append(a.index(y))

    #holds data for the newly formatted segments. Each element is a junction. Contains the start and end points of the four segments going into the junction, and the index of the junction the segment comes from
    split_segs=[[] for i in range(0,len(a)+(2*len(c)))]
    counter=0
    for i in range(0,len(junctions)):
        #for each segment
        if junctions[i]==[]:
            #if connected to no junctions
            counter+=1
        
        #get endpoints of segment i
        endpoints=find_endpoints(c[i][0],c[i][1],c[i][5],c[i][6])

        #This number assigns a unique junction number to each segment endpoint, all of which are after the numberings for the normal intersections
        endpoint_number=len(a)+(2*(i))-(2*counter)

        #this will store the intersection coordinate of a segment with each of the junctions it hits, along with the index of the intersection
        cross_spots=[[] for h in range(0,len(junctions[i]))]

        #fill the cross spots with the appropriate coordinates and indices
        for h in range(0,len(junctions[i])):
            cross_spots[h]=[b[junctions[i][h]][0],b[junctions[i][h]][1],junctions[i][h]]

        if cross_spots!=[]:
            #add "cross spots" for the endpoints using the endpoint number
            cross_spots.append([endpoints[0],endpoints[1],endpoint_number])
            cross_spots.append([endpoints[2],endpoints[3],endpoint_number+1])
            cross_spots.sort()
            if len(cross_spots)>1:
                for x in range(0,len(cross_spots)-1):
                    if x==0:
                        split_segs[cross_spots[x][2]].append([cross_spots[x+1][0],cross_spots[x+1][1],cross_spots[x][0],cross_spots[x][1],cross_spots[x+1][2]])
                    else:
                        split_segs[cross_spots[x][2]].append([cross_spots[x][0],cross_spots[x][1],cross_spots[x+1][0],cross_spots[x+1][1],cross_spots[x+1][2]])
                    split_segs[cross_spots[x+1][2]].append([cross_spots[x][0],cross_spots[x][1],cross_spots[x+1][0],cross_spots[x+1][1],cross_spots[x][2]])
    length=len(split_segs)
    counter=0
    while counter<length:
        if split_segs[counter]==[]:
            del split_segs[counter]
            length-=1
        else:
            counter+=1
    return(split_segs)

#This function creates the upper triangle of a conductance matrix. The element at row i, column j is -1/R where R is the resistance between node i and j (0 if no connection). If i and j are equal, then the quantity is 1/(sum(R)). The triangle is set up like this so that if we perform matrix multiplication with the row column voltage vector, we get a row column current vector representing the net current flowing into each node. See slide 6 of my Siemens Competition slideshow (posted in this repo) for a graphical depiction of what this looks like.
def gen_upper(a,b,c,d,e,f):
    #a is list of intersects, b is list of intersect locations, c is list of segments, d is sidelength of box, e is resistivity of material, f is radius of segments
    
    #We want to make this matrix banded to reduce the computation time of the matrix algebra as much as possible. This first chunk of code renumbers all the junctions in order of increasing x-coordinate, which significantly reduces the bandwidth of the matrix. maxdif tracks the bandwidth of the matrix
    maxdif=0        
    for a in range(0,len(new_segs)):
        intersect_renumber[a]=[new_segs[a][0][2],a]
    intersect_renumber.sort()
    for b in range(0,len(intersect_renumber)):
        intersect_renumber[b]=intersect_renumber[b][1]
    for i in range(0,len(new_segs)):
        new_number=intersect_renumber.index(i)
        for j in range(0,len(new_segs[i])):
            if new_segs[i][j]!=junctionresistance:
                segnumber=intersect_renumber.index(new_segs[i][j][4])
                if segnumber>new_number:
                    if (segnumber-new_number)+1>maxdif:
                        maxdif=segnumber-new_number+1

    #Declare the banded upper triangle matrix with the bandwidth just calculated
    global upper
    upper=[[0 for g in range(0,len(new_segs))] for w in range(0,maxdif)]

    #for each junction
    for i in range(0,len(new_segs)):
        new_number=intersect_renumber.index(i)
        #for each resistor coming into junction i
        for j in range(0,len(new_segs[i])):
            #if it's not the contact resistor
            if new_segs[i][j]!=junctionresistance:
                #get the index of the junction it comes from
                segnumber=intersect_renumber.index(new_segs[i][j][4])
                
                #calculate the resistance of the resistor
                distance=distance_formula(new_segs[i][j][0],new_segs[i][j][1],new_segs[i][j][2],new_segs[i][j][3])
                resistance=(distance*e)/(math.pi*(f**2))
                
                #Add the appropriate conductance to the diagnol element
                upper[-1][new_number]+=(1/resistance)
                
                #Add the appropriate conductance to the non-diagnol element, but only half the time, since we only store half the matrix
                if segnumber>new_number:
                    upper[new_number-segnumber-1][segnumber]=(-1/resistance)
        
            #if it's the junction resistor
            elif new_segs[i][j]==junctionresistance:
                #add the appropriate conductance to the diagnol element
                upper[-1][new_number]+=(1.0/float(junctionresistance))
                
                #Add the appropriate conductance to the non-diagnol element, but only half the time, since we only store half the matrix
                if 1-(2*(i%2))==1:
                    upper[-2][new_number+1]=(-1.0/float(junctionresistance))

#Prepares the banded upper matrix for solving by applying 2 numerical methods. First, we connect every node to "ground" (at 0V) via a very high resistance resistor. This prevents unsolvable networks when there are disconnected clusters without appreciably affecting the resistance numbers. Second, we find all endpoint junctions to the left of the space and apply a given voltage to them, and do the same for endpoint junctions to the right of the space. The numerical implementations of this are explained below, and are also explained in my Intel Science Talent Search Paper, which is in this repo
def strip_voltage_2(a,b,c,d):
    #a is left side voltage, b is right side voltage, c is split segs list, d is side length
    big_n=10**6
    
    #for each endpoint junction
    for x in range(2*len(intersects),len(c)):
        new_num=intersect_renumber.index(x)
        #if it is to the left of the space border
        if c[x][0][0]<0 or c[x][0][2]<0:
            #Set the voltage to the parameter for left voltage. This is done by connecting the endpoint node to a node of the desired voltage through a very low resistance resistor. Practically, this means that the diagnol element for this junction becomes much larger. The current element for this junction then gains a value equal to the desired voltage times that diagnol element.
            upper[-1][new_num]=upper[-1][new_num]*(big_n)
            currents[new_num]=upper[-1][new_num]*(a)
            zero_crossers.append(c[x])
        #if it is to the right of the space border
        if c[x][0][0]>d or c[x][0][2]>d:
            #set the voltage to the parameter for right voltage. See above.
            upper[-1][new_num]=upper[-1][new_num]*(big_n)
            currents[new_num]=upper[-1][new_num]*(b)
            far_side_crossers.append(c[x])
    #This for loop "grounds" every node as described above. Practically, since the ground is 0 V and connected through a huge resistance, this corresponds to adding a tiny number to every diagnol element.
    for a in range(0,len(upper[-1])):
        upper[-1][a]+=0.000000001


#rebuilds split segs to add contact resistance. each intersection junction becomes two junctions each with 3 segments going into it (3 from the original split segs, one new one from the contact resistance. The endpoint junctions at the end are unchanged
def junction_resistance(a,b):
    #Takes split segments and incorporates junction resistance
    #a is new_segs, b is junction resistance (Ohms)
    
    #insert the required new space
    for i in range(0, 2*len(intersects), 2):
        a.insert(i+1,[[],[],[]])
    
    #copy the data into the newly created junction
    for h in range(0, 2*len(intersects), 2):
        a[h+1][0]=a[h][2]
        a[h+1][1]=a[h][3]
        a[h+1][2]=float(b)
        del a[h][3]
        a[h][2]=float(b)

    #for each new junction
    for q in range(0, len(a)):
        #store the index of the first junction connected to the current junction
        reference1=a[q][0][4]
        
        #this if and the following elif adjust the indices of the junctions that each junction refers to based on the insertion of new junctions.
        if reference1<len(intersects):
            for t in range(0,2):
                if (a[2*reference1][t][0:2]==a[q][0][2:4] and a[2*reference1][t][2:4]==a[q][0][0:2]) or (a[2*reference1][t][2:4]==a[q][0][2:4] and a[2*reference1][t][0:2]==a[q][0][0:2]):
                    a[q][0][4]=2*reference1
                elif (a[(2*reference1)+1][t][0:2]==a[q][0][2:4] and a[(2*reference1)+1][t][2:4]==a[q][0][0:2]) or (a[(2*reference1)+1][t][2:4]==a[q][0][2:4] and a[(2*reference1)+1][t][0:2]==a[q][0][0:2]):
                    a[q][0][4]=(2*reference1)+1
        elif reference1>=len(intersects):
            a[q][0][4]=reference1+len(intersects)

        #This if and elif do the same thing for all junctions that are not endpoint junctions (the initial if checks for this
        if len(a[q])>1:
            reference2=a[q][1][4]
            if reference2<len(intersects):
                for t in range(0,2):
                    if (a[2*reference2][t][0:2]==a[q][1][2:4] and a[2*reference2][t][2:4]==a[q][1][0:2]) or (a[2*reference2][t][2:4]==a[q][1][2:4] and a[2*reference2][t][0:2]==a[q][1][0:2]):
                        a[q][1][4]=2*reference2
                    elif (a[(2*reference2)+1][t][0:2]==a[q][1][2:4] and a[(2*reference2)+1][t][2:4]==a[q][1][0:2]) or (a[(2*reference2)+1][t][2:4]==a[q][1][2:4] and a[(2*reference2)+1][t][0:2]==a[q][1][0:2]):
                        a[q][1][4]=(2*reference2)+1
            elif reference2>=len(intersects):
                a[q][1][4]=reference2+len(intersects)

#Given the voltage row vector, edge crossers, and parameters of the network, finds the network conductivity. We do this by summing the current in all the initial segments at the left side of the network. Then we use this total current and the applied voltage to find the resistance using Ohm's Law
def resistance_calculator(a,b,c,d,e,f):
    #calculates sheet resistance from voltages
    #a is list of voltages, b is zero_crossers, c is new_segs, d is resistivity of metal, e is radius of segments, f is applied voltage, g is far_side_crossers
    
    #This will hold the voltage changes over each initial segment
    one_step_voltages1=[]
    
    #this will hold the length of each initial segment
    distances1=[]
    
    #this will hold the resistance of each initial segment, calculated using resistivity
    resistances1=[]
    
    #this will hold the current through each initial segment
    initial_currents1=[]
    
    #this will hold the total current, the sum of the above vectors. The two should be approximately equal
    total_current1=0
    
    #for each junction left of the border
    for i in b:
        if i[0][2]<0:
            #Find the junction it is connected to
            for h in range(0,2*len(intersects)):
                for j in range(0,2):
                    if (c[h][j][0]==i[0][2] and c[h][j][1]==i[0][3]):
                        #Get the voltage difference and distance
                        one_step_voltages1.append(a[intersect_renumber.index(h)])
                        distance=distance_formula(c[h][j][0],c[h][j][1],i[0][0],i[0][1])
                        distances1.append(distance)

    #get the resistances, get the currents, sum them
    for v in distances1:
        resistances1.append((v*d)/(math.pi*(e**2)))
    for l in range(0,len(resistances1)):
        initial_currents1.append((app_voltage-one_step_voltages1[l])/resistances1[l])
    for h in initial_currents1:
        total_current1+=h
    print(str(total_current1/f)+" Siemens")
    print(str(f/total_current1)+" Ohms")
                   
def voltage_plotter(a,b):
#a is voltage list, b is intersect location list
    for i in range(0,len(b)):
        plotx.append(intersect_locations[i][0]*(app_voltage/sidelength))
        plotx.append(intersect_locations[i][0]*(app_voltage/sidelength))
        ploty.append(voltages[intersect_renumber.index(2*i)])
        ploty.append(voltages[intersect_renumber.index(2*i+1)])
    pylab.plot(plotx,ploty,'bo')
    pylab.ylim([0,app_voltage])
    pylab.xlim([0,sidelength*(app_voltage/sidelength)])
    pylab.axes().set_aspect('equal')

def heatmap_plotter(x,y,z,z1,diagramChoice):
    volt_diff=[]
    colors=[]
    scatterx=[]
    scattery=[]
    current_mapper_list=[]
    current_mapper_max=-10
    for a in range(0,2*len(intersects),2):
        if a%2==0:
            volt_diff.append(abs(y[intersect_renumber.index(a)]-y[intersect_renumber.index(a+1)]))
    volt_diff.sort()
    for a in range(0,2*len(intersects)):
        for b in range (0,2):
            voltage1=-10
            voltage2=-10
            if b==0:
                voltage2=voltages[intersect_renumber.index(a)]
                if [[x[a][b][2],x[a][b][3],x[a][b][0],x[a][b][1],a]] in zero_crossers:
                    voltage1=10.0
                else:
                    intersect1=[x[a][b][0],x[a][b][1]]
                    if intersect1 in intersect_locations:
                        voltage1=voltages[intersect_renumber.index(new_segs[a][b][4])]
                    else:
                        voltage1=voltage2
            if b==1:
                voltage1=voltages[intersect_renumber.index(a)]
                if [[x[a][b][0],x[a][b][1],x[a][b][2],x[a][b][3],a]] in far_side_crossers:
                    voltage2=0.0
                else:
                    intersect2=[x[a][b][2],x[a][b][3]]
                    if intersect2 in intersect_locations:
                        voltage2=voltages[intersect_renumber.index(new_segs[a][b][4])]
                    else:
                        voltage2=voltage1
            distance=distance_formula(x[a][b][0],x[a][b][1],x[a][b][2],x[a][b][3])
            resistance=(z*(distance))/((math.pi)*((z1)**2))
            current_mapper=((voltage1-voltage2)/resistance)
            if current_mapper<0:
                current_mapper=0
            if current_mapper>current_mapper_max:
                current_mapper_max=current_mapper
            current_mapper_list.append(current_mapper)    
    for a in range(0,2*len(intersects)):
        for b in range(0,2):
            plotx=[]
            ploty=[]
            plotx.append(x[a][b][0])
            plotx.append(x[a][b][2])
            ploty.append(x[a][b][1])
            ploty.append(x[a][b][3])
            current_color=(current_mapper_list[(2*a)+b]/current_mapper_max)
            width_value=(current_mapper_list[(2*a)+b]/current_mapper_max)
            if width_value<(0.25):
                width_value=1
            elif width_value<0.75:
                width_value=3
            else:
                width_value=4
            if (diagramChoice == 2):
                pylab.plot(plotx,ploty,color=(4*(current_color-0.5)**2,(1-current_color),4*(current_color*(1-current_color))),lw=width_value,zorder=width_value)
            else:
                pylab.plot(plotx,ploty,color=(0,0,1),zorder=1)
    if (diagramChoice!=2):
        for i in range(0,2*len(intersects),2):
            voltdrop=(abs(y[intersect_renumber.index(i)]-y[intersect_renumber.index(i+1)]))
            voltage=(y[intersect_renumber.index(i)]+y[intersect_renumber.index(i+1)])/2
            if (diagramChoice == 1):
                colors.append(voltage)
            elif (diagramChoice == 3):
                colors.append(voltdrop)
            scatterx.append(x[i][0][2])
            scattery.append(x[i][0][3])
        pylab.scatter(scatterx,scattery,s=20,c=colors,cmap=cm.bwr,zorder=2)
    pylab.ylim([0-piecelength/2,sidelength+piecelength/2])
    pylab.xlim([0-piecelength/2,sidelength+piecelength/2])
    pylab.axes().set_aspect('equal')
                                                                                                                                                                                                                                                                                      
new_segs=split_segs(intersects,intersect_locations,segments,sidelength)
junction_resistance(new_segs,junctionresistance)

#This will be a row vector of currents where the ith element is the net current entering junction i
currents=[0 for efg in range(0, len(new_segs))]

intersect_renumber=[0 for x in range(0,len(new_segs))]
gen_upper(intersects,intersect_locations,segments,sidelength,1.59*(10**(-8)),6.5*(10**(-8)))
strip_voltage_2(app_voltage,0,new_segs,sidelength)

#Use a cholesky decomposition of upper to solve for the voltages at every node in the network. The network resistance can than be calculated from these voltages
chol_f=linalg.cholesky_banded(upper)
voltages=linalg.cho_solve_banded([chol_f,False],currents)

resistance_calculator(voltages,zero_crossers,new_segs,1.59*(10**(-8)),6.5*(10**(-8)),app_voltage)

#Plot the chosen diagram
if (diagramChoice == 4):
    voltage_plotter(voltages,intersect_locations)
elif (diagramChoice != 0):
    heatmap_plotter(new_segs,voltages,1.59*(10**(-8)),6.5*(10**(-8)),diagramChoice)
pylab.plot([0,sidelength,sidelength,0,0],[0,0,sidelength,sidelength,0],"k")
if (diagramChoice == 1 or diagramChoice == 3):
    pylab.colorbar()
pylab.savefig('heatmap.png', bbox_inches='tight')
#pylab.show()

e.close()
asdfgh.close()
g.close()
v.close()
