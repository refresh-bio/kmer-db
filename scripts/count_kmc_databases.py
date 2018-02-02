#!/usr/bin/python3

import os
import shutil
import subprocess
import queue
import threading 
kmc="../../kmc/bin/kmc"
kmc_tools="../../kmc/bin/kmc_tools"
kmer_len = 18
num_thread=32

input_files = 'files.txt'

def count_file(src, tmp):
	dest=src.split('/')[-1].split(".fna.gz")[0]
	if os.path.exists(tmp):
		print("Error: {0} already exists!!!".format(tmp))
	else:
		os.makedirs(tmp)
	count_command="{0} -r -t1 -k{1} -ci1 -cs1 -fm -p7 {2} {3} {4}".format(kmc, kmer_len, src, dest, tmp)
	sort_command="{0} -t1 transform {1} sort {2}".format(kmc_tools, dest, dest+'_s')
	rm_command="rm {0}.kmc_*".format(dest)
	
	FNULL = open(os.devnull, 'w')
	subprocess.call([count_command], shell=True,stdout=FNULL)
	subprocess.call([sort_command], shell=True,stdout=FNULL)
	subprocess.call([rm_command], shell=True,stdout=FNULL)
	shutil.rmtree(tmp)

q=queue.Queue()

def worker():
	while True:
		(src,tmp) = q.get()
		count_file(src, tmp)
		q.task_done()


for i in range(0,num_thread):
	t=threading.Thread(target=worker)
	t.daemon=True
	t.start()

f=open(input_files)
lines = f.readlines()
nr=0
for l in lines:	
	q.put((l[:-1], 'tmp' + str(nr)))
	nr+=1

q.join()

