#Efficient updating strategy uses vertex-vertex index to update the vertex coverage table

import csv
import ast
import hdbscan
import numpy as np
import pandas as pd
with open('trajectories.csv') as csvfile: 
	readcsv=csv.reader(csvfile,delimiter=',')
	row1=0
	points=[]
	#locs_in_traj is the trajectory vertex index i.e., the locations each trajectory passes through
	locs_in_traj={}
	for row in readcsv:
		if row1!=0 and row1<50:
			x=ast.literal_eval(row[8])
			key=row1
			list=[]
			for z in x:
				obj={}
				obj['x']=z[0]
				obj['y']=z[1]
				flg=0
				list.append(obj)
				for point in points:
					if point['x']==obj['x'] and point['y']==obj['y']:
						flg=1
				if flg==0:
					points.append(obj)
			locs_in_traj[key]=list
			row1=row1+1
		else:
			row1=row1+1

print("-------Efficient updating strategy-------")
no_of_points=len(points)
print("No of points identified in the region are:-",end=" ")
print(no_of_points)

#inv_traj_index is the vertex-trajectory index and it has the trajectories passing throughit
inv_traj_index = {}

for point in points:
	k=hash(str((point))) #to generate a unique hash code for each point of the form (x,y)
	inv_traj_index[k]=[]

i=1
while(i<50):
	for z in locs_in_traj[i]:
		k=hash(str((z)))
		inv_traj_index[k].append(i)
	i=i+1

#print("vertex-trajectory index is as follows:")
#for point in points:
#	k=hash(str((point)))
	#print(point,end=":- ")
	#print(inv_traj_index[k],end="||")
	
#ver_ver_index is the vertex-vertex index and it has the number of trajectories passing through both of them
ver_ver_index={}

#Initialising vertex-vertex index for each pair
i=0
while(i<no_of_points):
	j=i+1
	while(j<no_of_points):	
		if(i!=j):
			k1=hash(str(points[i])) #hashes for the two points picked now
			k2=hash(str(points[j]))
			ver_ver_index[k1,k2]=0 
			ver_ver_index[k2,k1]=0
		j=j+1
	i=i+1

#computing the ver_ver_index
i=0
while(i<no_of_points):
	j=i+1
	while(j<no_of_points):	
		if(i!=j):
			k1=hash(str(points[i]))
			k2=hash(str(points[j]))
			for x in inv_traj_index[k1]:
				for y in inv_traj_index[k2]:
					if(x==y):
						ver_ver_index[k1,k2]=ver_ver_index[k1,k2]+1
						ver_ver_index[k2,k1]=ver_ver_index[k2,k1]+1

		j=j+1
	i=i+1

#print("vertex-vertex index is as follows:-")
i=0
while(i<no_of_points):
	j=i+1
	while(j<no_of_points):	
		if(i!=j):
			k1=hash(str(points[i]))
			k2=hash(str(points[j]))
			#print(points[i],end=" and ")
			#print(points[j],end=" :- ")
			#print(ver_ver_index[k1,k2],end="||")
		j=j+1
	i=i+1
	
most_inf_loc=4

#v_final is list consisting of k-most influential locations
v_final=[]

#it contains the trajectories that will be covered if this location is selected
vertex_coverage_table={}

#initial coverage of each vertex is its inv_traj_index value
for point in points:
	k=hash(str((point)))
	vertex_coverage_table[k]=[]
	#looping is done to eliminate duplicates
	#print(inv_traj_index[k],end="--")
	for z in inv_traj_index[k]: 
		flag=0
		for y in vertex_coverage_table[k]:
			if(z==y):
				flag=1
		if(flag==0):
			vertex_coverage_table[k].append(z)
#to store only the number of trajectories that will be covered			
reduced_vertex_coverage_table={}
for point in points:
	k=hash(str((point)))
	reduced_vertex_coverage_table[k]=len(vertex_coverage_table[k])
	#print(reduced_vertex_coverage_table[k],end="--")
		
traj_final=[]
#returns the location that covers maximum no.of trajectories
def maxim():
	z=[item for item in points if item not in v_final]
	max_so_far=0
	choosen_vertex=z[0]
	for x in z:
		y=reduced_vertex_coverage_table[hash(str(x))]
		if(max_so_far<y):
			max_so_far=y
			choosen_vertex=x
	return choosen_vertex

i=0
no_of_basic=0
no_of_efficient=0
#according to the greedy heuristic,loop runs for most_inf_loc number of times
while(i<most_inf_loc):
	local_max=maxim()
	#print("The local maximum is :- ",end=" ")
	#print(local_max)
	traj_com=[]
	traj_new=[]
	for x in inv_traj_index[hash(str(local_max))]:
		flag=0
		for y in traj_final:
			if(x==y):
				flag=1
		if(flag==1):
			traj_com.append(x)
		else:
			traj_new.append(x)
	#print(inv_traj_index[hash(str(local_max))])
	#print(traj_final)
	#print(traj_new)
	#print(traj_com)
	#perform basic updating strategy
	if(len(traj_new)<=len(traj_com)):
		no_of_basic=no_of_basic+1
		for z in inv_traj_index[hash(str(local_max))]:
			for q in locs_in_traj[z]:
				vertex_coverage_table[hash(str(q))].remove(z)
		for point in points:
			k=hash(str((point)))
			reduced_vertex_coverage_table[k]=len(vertex_coverage_table[k])
	
	#perform efficient updating strategy
	else:
		no_of_efficient=no_of_efficient+1
		#locations that are not yet selected into the final group
		z=[item for item in points if item not in v_final]
		
		#subtract the no of shared trajectories from the ver_ver_index for each remaining vertices
		for x in z:
			if(x!=local_max):
				reduced_vertex_coverage_table[hash(str(x))]-=ver_ver_index[hash(str(x)),hash(str(local_max))]
		
		#adding the coverage values that are subtracted in previous iterations
		for y in traj_com:
			for p in locs_in_traj[y]:
				flag=0
				for q in z:
					if(p==q and q!=local_max):
						flag=1
				if(flag==1):
					reduced_vertex_coverage_table[hash(str(p))]+=1
	
	#ading the location to the final set
	v_final.append(local_max)
	
	#updating the covered trajectories list with those covered by the new location
	for x in inv_traj_index[hash(str(local_max))]:
		flag=0
		for y in traj_final:
			if(x==y):
				flag=1
		if(flag==0):
			traj_final.append(x)
	i=i+1
print()
print("The most influential k-location set in the region is:")
print(v_final)	
