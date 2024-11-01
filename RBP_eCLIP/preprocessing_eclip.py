import glob
import os
import time
import subprocess
from queue import Queue
import threading
import sys

# python preprocessing_eclip.py K562

def run_preprocessing_py(queue):
	while not queue.empty():
		j = str(queue.get(True,1))
		sample = j.rsplit("/",1)[1].split(".")[0]
		Out1 = bam_r2_dir + sample + ".r2.bam"
		Out2 = bed_r2_dir + sample + ".r2.bed"
		Out3 = shifted_bed_r2_dir + sample + ".r2.shifted.bed"
		Out4 = plus_minus_bedgraph_dir + sample + ".r2.plus.bedgraph"
		Out5 = plus_minus_bedgraph_dir + sample + ".r2.minus.bedgraph"
		Out6 = coverage_dir + sample + ".plus_minus.coverage.bed"
		os.system("python " +  "calculate_coverage.py " + sample + " " + cell)
		print("Done: " + sample)
		queue.task_done()
#
cell = str(sys.argv[1])

bam_dir = ""
files = glob.glob(bam_dir + "*.bam")
bam_r2_dir = ""
bed_r2_dir = ""
shifted_bed_r2_dir = ""
plus_minus_bedgraph_dir = ""
coverage_dir = ""

queue=Queue(maxsize=300000)
for i in files:
	queue.put(i)
threads = []
for i in range(10):
	t=threading.Thread(target=run_preprocessing_py,args=(queue,))
	threads.append(t)
for i in range(10):
	threads[i].start()
queue.join()
print("ALL DONE at " + time.asctime( time.localtime(time.time()) ))

#samtools view -hb -f 128  IGF2BP1_1.bam  >  IGF2BP1_1.r2.bam
#bedtools bamtobed -i IGF2BP1_1.r2.bam > IGF2BP1_1.r2.bed &
#bedtools shift -m 1 -p -1 -i IGF2BP1_1.r2.bed -g sizes.genome > IGF2BP1_1.r2.shifted.bed &
#bedtools genomecov -bg -strand + -5 -i IGF2BP1_1.r2.shifted.bed -g sizes.genome > IGF2BP1_1.r2.plus.bedgraph &
#bedtools genomecov -bg -strand - -5 -i IGF2BP1_1.r2.shifted.bed -g sizes.genome > IGF2BP1_1.r2.minus.bedgraph &
#perl make_coverage.pl IGF2BP1_1 &