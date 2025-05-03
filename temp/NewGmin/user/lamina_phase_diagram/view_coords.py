#This script will read in the coords files to visuliase the elastic filaments
import numpy as np
import matplotlib.pyplot as plt

#1) Choose the file we want to visualise
path = "~/Desktop/NewGmin/user/lamina_phase_diagram"

filename="/Users/katemcgarva/Desktop/NewGmin/user/lamina_phase_diagram/lamina_NROD10_NSEG_200_RAT5/0.0267F_S0K30.001BsRat1.0Obc0.1_1"

mycoords=(filename + "/lowest")
#X=(filename + "/coords.perturb")
nheader=4
myinits=(filename + "/initpos")
nheaderinit=0

#2) User inputs
N_ROD=10
N_SEG=200
L=5

#3) Read in the coordinates
fid = open('/Users/katemcgarva/Desktop/NewGmin/user/lamina_phase_diagram/lamina_NROD10_NSEG_200_RAT5/0.0267F_S0K30.001BsRat1.0Obc0.1_1/lowest', 'r')

angles=np.loadtxt(fid,skiprows=nheader)
fid.close
print(len(angles))

fid=open(myinits,'r')
inits=np.loadtxt(fid,skiprows=nheaderinit)
fid.close

#4) Restructure the data and convert to xy
LS=L/N_SEG
coords=np.zeros([N_ROD,N_SEG,2])
coords[:,0,0]=inits[:,0]
coords[:,0,1]=inits[:,1]
s=0
for i in range(0,N_ROD):
	for j in range(1,N_SEG):
		coords[i,j,0]=coords[i,j-1,0]+angles[s+1990]*np.cos(angles[s]) #1990
		coords[i,j,1]=coords[i,j-1,1]+angles[s+1990]*np.sin(angles[s])
		s=s+1

#5) View the results
fig1=plt.figure()
ax1 = plt.axes()

for i in range(0,N_ROD):
	plt.plot(coords[i,:,0],coords[i,:,1])

plt.show()

#6) print the length of each rod
s=0
for i in range(0,N_ROD):
		lengths = np.zeros(N_ROD)
		lengths[i]=np.sum(angles[199*i+199*N_ROD:199*i+199*(N_ROD+1)])+LS

		print(f"length {i+1}={lengths[i]}")
