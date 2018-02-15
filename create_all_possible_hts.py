#!/usr/bin/python

import sys
import itertools

n = int(sys.argv[1])
outfile = 'hts.txt'

if n > 10:
	sys.exit("Currently only supporting a maximum of 10 loci")

global_list = []
combs = list(set(list(itertools.combinations('abab',2))))
hts = []
#global_list.append(combs)
#combs = list(set(list(itertools.combinations('2323',2))))
#global_list.append(combs)
#combs = list(set(list(itertools.combinations('4545',2))))
#global_list.append(combs)

for i in range(n):
	global_list.append(combs)

#	hts=[]
#	for element in list(itertools.product(*global_list)):
#		first = ''
#		second = ''
#		for some in element:
#			first += some[0]
#			second += some[1]
#		hts.extend([first,second])


#print len(global_list)
#for locus in global_list:
#	print locus

#for j in range(len(global_list)):
#	for i in range(4):
#		for k in range(4):
#			string = ''
#			for z in range(len(global_list)):
##				string += "".join(global_list[z][i])
#				if z == j:
#					string += "".join(global_list[j][k])
#				else:
#					string += "".join(global_list[z][i])
##			string += "".join(global_list[j][k])
##			print "locus %s state %s" %
##			print len(locus),locus
#			print string

#global_list=list(set(global_list))
fh = open('gts.txt', 'w')
for gt in list(itertools.product(*global_list)):
	string = ''
	for g in gt:
		string += "".join(g)
	first=''
	second=''
	for i in xrange(0, len(string), 2):
		first += string[i]
	
	for i in xrange(1, len(string), 2):
		second += string[i]

	hts.extend([first,second])
	fh.write(string+"\n")

fh.close()

hts = list(set(hts))
fh = open('hts.txt', 'w')
for ht in hts:
	fh.write(ht+'\n')
	


#testlist = [[0,1],[2,3],[4,5]]

#testlist2=['01','23','45']
##print list(itertools.product(*testlist2))

#test3 = ['0','1']

#test4 = '0101'

#for j in range(1,11):
#	global_list = []
#	for i in range(j):
#		global_list.append(list(set(list(itertools.combinations('0101',2)))))
#
#	print j,len(list(itertools.product(*global_list)))
#
#
#	hts = []
#	for element in list(itertools.product(*global_list)):
#		first = ''
#		second = ''
#		for some in element:
#			first += some[0]
#			second += some[1]
#		hts.extend([first,second])
#
#	print j,len(list(set(hts)))
##	print j,list(set(hts))
#	print
