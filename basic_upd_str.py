import csv
import ast
import hdbscan
import numpy as np
import pandas as pd
with open('trajectories.csv') as csvfile: 
	readcsv=csv.reader(csvfile,delimiter=',')
	row1=0
	points=[]
	locs_in_traj={}
	for row in readcsv:
		if row1!=0 and row1<50:
			x=ast.literal_eval(row[8])
			key=row1
			#print(key)
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

#print(locs_in_traj[320])
no_of_points=len(points)
inv_traj_index = {}

for point in points:
	k=hash(str((point))) #to generate a unique hash code for each point of the form (x,y)
	inv_traj_index[k]=[]

i=1
while(i<50): #because 10 trajectories are considered initially
	#print(locs_in_traj[i])
	for z in locs_in_traj[i]:
		k=hash(str((z)))
		#print(k)
		inv_traj_index[k].append(i)
	i=i+1


for point in points:
		k=hash(str(point))
		#print(inv_traj_index[k],end="--")	

most_inf_loc=6
#v_final is list consisting of k-most influential locations
v_final=[]

vertex_coverage_table={}
#initial coverage of each vertex is its inv_traj_index value
for point in points:
	k=hash(str((point)))
	vertex_coverage_table[k]=[]
	for z in inv_traj_index[k]:
		flag=0
		for y in vertex_coverage_table[k]:
			if(z==y):
				flag=1
		if(flag==0):
			vertex_coverage_table[k].append(z)

def maxim():
	z=[item for item in points if item not in v_final]
	#print(len(z))
	max_so_far=0
	choosen_vertex=z[0]
	for x in z:
		#print(vertex_coverage_table[hash(str(x))])
		y=len(vertex_coverage_table[hash(str(x))])
		if(max_so_far<y):
			max_so_far=y
			choosen_vertex=x
			#print(y)
	return choosen_vertex
	
#most_inf_loc loop in greedy heuristic	
i=0
while(i<most_inf_loc):
	local_max=maxim()
	v_final.append(local_max)
	loc=hash(str(local_max))
	#print(inv_traj_index[loc])
	for z in inv_traj_index[loc]:
		for q in locs_in_traj[z]:
			if z in vertex_coverage_table[hash(str(q))]: vertex_coverage_table[hash(str(q))].remove(z)
	i=i+1
print(v_final)	


	
	
	
	
	
	
	
	
	