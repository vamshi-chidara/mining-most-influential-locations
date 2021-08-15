#I have used hdbscan library for clustering
import itertools
import csv
import ast
import hdbscan
import numpy as np
import pandas as pd
#trajectories.csv file has the trajectory data.
with open('trajectories.csv') as csvfile: 
	readcsv=csv.reader(csvfile,delimiter=',')
	row1=0
	points=[]
	#locs_in_traj is used as trajectory-vertex index.i.e.,locations in each trajectory.
	locs_in_traj={}
	#this for loop is used to identify unique locations in the spatial region using trajectory data.
	for row in readcsv:
		if row1!=0 and row1<10:
			x=ast.literal_eval(row[8])
			key=row1
			#print(key)
			list1=[]
			for z in x:
				obj={}
				obj['x']=z[0]
				obj['y']=z[1]
				flg=0
				list1.append(obj)
				#to check if the point(location) is previously found in earlier trajectories.
				for point in points:
					if point['x']==obj['x'] and point['y']==obj['y']:
						flg=1 #flag=1 indices that it is previously found.
				if flg==0:
					points.append(obj)
				#locs_in_traj[k] gives the locations that are involved in trajectory number k.
			locs_in_traj[key]=list1
			row1=row1+1
		else:
			row1=row1+1
df=pd.DataFrame(points)
print("-------Group pruning optimal algorithm--------")
#Total number of locations identified
count=df.shape[0]
print("The No of Locations in the spatial region are:",end=" ")
print(count)
#print("The Locations in the spatial region are:")
#print(df)
#hdbscan cluster object with atleast 10 locations in each cluster. This variable can be changed as per requirements.
clusterer = hdbscan.HDBSCAN(min_cluster_size=10)
#fit the locations identified above to the object.
cluster_labels = clusterer.fit(df)
print("clusters formed are:")
#locations belonging to the same cluster have the same index in the cluster.labels_
print(clusterer.labels_)

inv_traj_index = {}
for point in points:
	k=hash(str(point)) #to generate a unique hash code for each point of the form (x,y)
	inv_traj_index[k]=[]

i=1
while(i<10): #because 10 trajectories are considered initially
	for z in locs_in_traj[i]:
		k=hash(str(z))
		inv_traj_index[k].append(i)
	i=i+1

hdb=clusterer.labels_
no_of_clusters=clusterer.labels_.max()+1
print()
print("no of clusters formed - ",end="")
print(no_of_clusters)
#the locations in each cluster
clusters = np.empty((no_of_clusters, 0)).tolist()
i=0
for x in hdb:
	if(x!=-1):
		clusters[x].append(points[i])
		i=i+1
	else:
		i=i+1

i=0
#while(i<no_of_clusters):
	#print(" ")
	#print(i)
	#print(" ")
	#print(clusters[i])
	#i=i+1
#contains the trajectories covered by the locations in that cluster
cluster_coverage = np.empty((no_of_clusters, 0)).tolist()

#to get the coverage of each cluster
j=0
while(j<no_of_clusters): #for each cluster j out of "no_of_clusters" clusters
	for x in clusters[j]: #for each point x in cluster number j
		k=hash(str(x))
		for y in inv_traj_index[k]: #for each trajectory y that has passed through the point x
			flag=0
			for z in cluster_coverage[j]: #check if the trajectory y has already been counted in the coverage of the cluster j
				if(z==y):
					flag=1 #flag=1 indicates that the trajectory y is already added to the coverage of cluster j
			if(flag==0):
				cluster_coverage[j].append(y)
	j=j+1

#K value
most_inf_k=2
cluster_list=[] 
i=0
while(i<no_of_clusters):
	cluster_list.append(i)
	i=i+1
#print(cluster_list)

#all combinations of size most_inf_k in cluster_list without repititions
combin = list(itertools.combinations(cluster_list, most_inf_k)) 
comb_list=[list(t) for t in combin]
print(" ")
print("All k grp combinations:")
print(comb_list)
comb_list_len=len(comb_list)

#coverage of each combination
k_grp_set_coverage = np.empty((comb_list_len, 0)).tolist()
i=0
while(i<comb_list_len):
	temp_list=comb_list[i]
	j=0
	while(j<most_inf_k):
		z=temp_list[j]
		for a in cluster_coverage[z]:
			flag=0
			for b in k_grp_set_coverage[i]:
				if(a==b):
					flag=1
			if(flag==0):
				k_grp_set_coverage[i].append(a)
		j=j+1
	i=i+1

#count of the coverage of each combination.
count_of_combs=[]
i=0
while(i < comb_list_len):
	count_of_combs.append(len(k_grp_set_coverage[i]))
	i=i+1

#sort combinations(i.e., comb_list) in ascending order bases on their coverage(i.e., count_of_combs)
comb_list = [comb_list for _,comb_list in sorted(zip(count_of_combs,comb_list))]

#sort k_grp_set_coverage in descending order based on their cardinality
k_grp_set_coverage = [k_grp_set_coverage for _,k_grp_set_coverage in sorted(zip(count_of_combs,k_grp_set_coverage))]
k_grp_set_coverage.reverse()
#print(k_grp_set_coverage)

#For Best-first-pruning (considering combinaation with highest coverage first)
comb_list.reverse()
print()
print("Combinations in sorted order for best first pruning:")
print(comb_list)

#Sort count_of_combs in descending order
count_of_combs.sort(reverse=True)
#print(count_of_combs)

#maximum coverage so far
max_so_far=0
final_k_loc_set=[]
#considering each k-cluster group in the decreasing order of their coverages
i=0
for k_grp in comb_list:
	if(count_of_combs[i]>max_so_far and i<6 and i>2):
		locs_in_comb=[]
		for j in comb_list[i]:
			for loc in clusters[j]:
				locs_in_comb.append(loc)
		#print(len(locs_in_comb))
		loc_sets = list(itertools.combinations(locs_in_comb, most_inf_k)) 
		k_locs_list=[list(t) for t in loc_sets]
		#print(len(k_locs_list))
		for k_locs_set in k_locs_list:
			k_locs_set_coverage=[]
			for temp in k_locs_set:
				var=hash(str(temp))
				for tr in inv_traj_index[var]:
					flag=0
					for z in k_locs_set_coverage:
						if(tr==z):
							flag=1
					if(flag==0):
						k_locs_set_coverage.append(tr)
			if(max_so_far<len(k_locs_set_coverage)):
				#print(len(k_locs_set_coverage))
				max_so_far=len(k_locs_set_coverage)
				final_k_loc_set=list(k_locs_set)
	i=i+1
print(" ")
print("Most influential k location set:")
print(final_k_loc_set)



