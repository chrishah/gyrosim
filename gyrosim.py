#!/usr/bin/python

import sys
from datetime import datetime
import random
from collections import defaultdict
import numpy as np

#seed = 12345
#random.seed(seed)

fish = "".join(datetime.utcnow().strftime('%Y%m%d%H%M%S.%f').split('.'))

start_pressure = 2 # each tick each individual is challenged with a probability of 0.05. if challenged the individul dies with probability of 0.05
pressure = start_pressure

#For density dependent regime
density_dependence_interval = 10 # in hours at which population size will be checked and compared against the specified
density = 50 # desired population size

time_dependence_interval = 96
incremental_pressure_increase = 1

last_only = False
genome_size = 100

if not len(sys.argv) > 2:
	sys.exit("\nPlease provide the desired number of hours for which to run the experiment and the mode of reproduction for second and subsequent births (plus size of genome optionally).\n")
ticks = int(sys.argv[1])+1
mode = str(sys.argv[2])
if len(sys.argv) == 4:
	genome_size = int(sys.argv[3])



#'sexual', 'apomixis', 'agd', 'atf', 'acf', 'arf'
#mode='agd' 


def create_automictic_gt(A,prob_hom):#, seed=seed):
	
#	random.seed(seed)

	if random.randint(0,100) < prob_hom:
		if random.randint(0,100) < 50:
                	offspring = A[0]+A[0]
                else:
                	offspring = A[1]+A[1]
	else:
		offspring = A
	
	return offspring


def reproduce(mode, parentA, parentB=[], hom_prob=[]):#, seed=seed):

#	random.seed(seed)

	offspring = []
	if mode == 'sexual':
		if not parentB:
			parentB = parentA

#		print "Reproductive mode: Sexual"
		if len(parentA) == len(parentB):
			for i in range(len(parentA)):
#				string = parentA[i]+parentA[i]
#				offspring.append(string[random.randint(0,3)])
#				string = parentB[i]+parentB[i]
#				offspring[i] += string[random.randint(0,3)]
				offspring.append(parentA[i][random.randint(0,1)]+parentB[i][random.randint(0,1)])
#				offspring[-1] = "".join(sorted(offspring[-1]))
		else:
			sys.exit("Sexual selected - need valid genotpyes from 2 individuals")
	else:
#		print "Parthenogenesis"
		if not hom_prob:
			sys.exit("To simulate parthenogenetic reproduction I need the expected probabilities for transition to homozygosity")
		if parentB:
			print "Parthenogenesis selected - Genotypes from parent B will be ignored!"

		if mode == 'apomixis':
#			print "Reproductive mode: Apomixis"
			for i in range(len(parentA)):
				offspring.append(parentA[i])
		elif mode == 'agd': #'auto_gamete_duplication':
#			print "Reproductive mode: Automictic parthenogenesis - gamete duplication"
			for i in range(len(parentA)):
				if random.randint(0,100) < 50:
					offspring.append(parentA[i][0]+parentA[i][0])
				else:
					offspring.append(parentA[i][1]+parentA[i][1])
#				offspring[-1] = "".join(sorted(offspring[-1]))

		elif mode == 'atf': #'auto_terminal_fusion':
#			print "Reproductive mode: Automictic parthenogenesis - terminal fusion"
			for i in range(len(parentA)):
                                offspring.append(create_automictic_gt(parentA[i],hom_prob[i]))
				offspring[-1] = "".join(sorted(offspring[-1]))
				
		elif mode == 'acf': #'auto_central_fusion':
#			print "Reproductive mode: Automictic parthenogenesis - central fusion"
			for i in range(len(parentA)):
                                offspring.append(create_automictic_gt(parentA[i],hom_prob[i]))
#				offspring[-1] = "".join(sorted(offspring[-1]))
		
		elif mode == 'arf': #'auto_random_fusion':
#			print "Reproductive mode: Automictic parthenogenesis - random fusion"
			for i in range(len(parentA)):
                                offspring.append(create_automictic_gt(parentA[i],hom_prob[i]))
#				offspring[-1] = "".join(sorted(offspring[-1]))

	return offspring

def count_gts(gts):
	
	from collections import defaultdict
	counts=defaultdict(int)
	for i in range(len(gts)):
		counts["".join(sorted(gts[i]))] += 1

	return counts

def generate_genome(n):
	#n is the number of heterozygous loci in the genome
	genome = []	

	for i in range(n):
		genome.append('ab')

	return genome

def assign_homozygous_transition_prob(genotypes, model):#, seed=seed):

#	random.seed(seed)

	#genotypes is a list of genotypes
	probabilities = []

	for i in range(len(genotypes)):
		if model == 'apomixis':
			probabilities.append(0)

		elif model == 'agd': #'auto_gamete_duplication':
			probabilities.append(100)
		elif model == 'atf': #'auto_terminal_fusion':
			probabilities.append(random.randint(33,100))
		elif model == 'acf': #'auto_central_fusion':
			probabilities.append(random.randint(0,33))
		elif model == 'arf': #'auto_random_fusion':
			probabilities.append(33)
		elif model == 'sexual': #sexual
			probabilities.append(50)
		else:
			sys.exit('You need to specifiy a valid model of parthenogenesis')	

	return probabilities

#genome = generate_genome(genome_size)

#probs = assign_homozygous_transition_prob(genotypes=genome, model=mode)
#offspring = reproduce(mode=mode, parentA=genome, hom_prob=probs, parentB=genome)

#for i in range(len(genome)):
#	print "\t%s\t%s\t%s" %(probs[i], genome[i], offspring[i])

#freqs = count_gts(offspring)
#print "\taa: %s\tab: %s\tbb: %s" %(freqs['aa'], freqs['ab'], freqs['bb'])

def make_prob(maxi, mu=0, sigma=5):#, seed=seed):

	import numpy as np

#	random.seed(seed)
	
	s = np.random.normal(mu, sigma, 100000)
	count = defaultdict(int)
	for n in s:
		out = int(round(n,0))+maxi
		if out <= maxi:
			count[out] += 1

#	print count
	divide_by=count[maxi]
#	print divide_by

	for key in count:
		count[key] = int(round(count[key]/float(divide_by)*100))

	return count

def add_minimum(minimum, dictionary):

	outdict = dictionary.copy()
	for i in range(sorted(outdict)[-1]+1):
		if not i in outdict:
			outdict[i] = minimum
		else:
			if outdict[i] < minimum:
				outdict[i] = minimum

	return outdict

first_birth_prob=make_prob(maxi=50)
sec_and_sub_prob=make_prob(maxi=100)
initial_death_prob=make_prob(maxi=340, sigma=10)

#print "Death probability FIRST:"
#for hour in sorted(death_prob):
#	print hour,death_prob[hour]

first_birth_prob = add_minimum(minimum=0, dictionary=first_birth_prob)
sec_and_sub_prob = add_minimum(minimum=0, dictionary=sec_and_sub_prob)

death_prob = add_minimum(minimum=pressure, dictionary=initial_death_prob)

pressure_regime = {}#{1:1, 240:5, 400:11}
#increase_interval = 10000 #240
#up = 1 #increase the pressure by that many percent

genome = generate_genome(genome_size)
apomixis_probs = assign_homozygous_transition_prob(genotypes=genome, model='apomixis')
probs = assign_homozygous_transition_prob(genotypes=genome, model=mode)
#agd_probs = assign_homozygous_transition_prob(genotypes=genome, model='agd') 
#acf_probs = assign_homozygous_transition_prob(genotypes=genome, model='acf') 
#atf_probs = assign_homozygous_transition_prob(genotypes=genome, model='atf') 
#arf_probs = assign_homozygous_transition_prob(genotypes=genome, model='arf') 

founder = "".join(datetime.utcnow().strftime('%Y%m%d%H%M%S.%f').split('.'))
daughter = ''

population = {founder: {'mother':'origin', 'births': int(0), 'sex': 'female', 'born':0, 'gave_birth':0, 'genome':genome}}
individuals = []
sex_distribution = defaultdict(int)
alleles = []
av_heterozygosity = []
all_homo = False
final_out = ''
for i in range(len(population[founder]['genome'])):
	alleles.append('')
	av_heterozygosity.append('')	

birth_counts = {1:0, 2:0, 3:0, 4:0}

prev_sizes = []

for i in range(1,ticks):

	#at regular intervals
	ok = True
	if i % density_dependence_interval == 0:
		if density > 0:
			std = np.std(prev_sizes,keepdims=True)[0]
			mean = np.mean(prev_sizes)
#			print prev_sizes,np.mean(prev_sizes),std,int(std)
			prev_sizes = []
			if len(population) > mean:
				if len(population) >= density*0.9:
					if len(population) < density:
						print "slight adjust - increase (above) pressure by %i" %int(round(std,0)/2)
						pressure += int(round(std,0)/2)
					elif len(population) >= density:
						print "adjust - increase (above) pressure by %i" %int(round(std,0))
						pressure += int(round(std,0))
					ok = False
			elif len(population) < mean:
				if len(population) <= density*1.1:
					if len(population) > density:
						print "slight adjust - decrease (below) pressure by %i" %int(round(std,0)/2)
						pressure -= int(round(std,0)/2)
					elif len(population) <= density:
						print "adjust - decrease (below) pressure by %i" %int(round(std,0))
						pressure -= int(round(std,0))
					ok = False
	prev_sizes.append(len(population))

	if i % time_dependence_interval == 0:
		print "time dependent increase by %i" %incremental_pressure_increase
		start_pressure += incremental_pressure_increase
		if start_pressure > pressure:
			pressure = start_pressure
		ok = False

	if not ok:
		if pressure < start_pressure:
			pressure = start_pressure
#			print "new pressure: %s" %pressure
		death_prob = add_minimum(minimum=pressure, dictionary=initial_death_prob)


	if pressure_regime and i in pressure_regime:
#	if i % increase_interval == 0:
#		pressure += up
#		print "# DAY: %.1f - INCREASING THE PRESSURE TO: %s" %((float(i)/24), pressure)
		pressure = pressure_regime[i]
		death_prob = add_minimum(minimum=pressure_regime[i], dictionary=initial_death_prob)


	for j in range(len(genome)):
		alleles[j] = ''
		av_heterozygosity[j] = float(0)

	ages = []
	born=0
	died=0
	sex_distribution['female'] = 0
	sex_distribution['male'] = 0
	individuals = population.keys()
	for ind in individuals:
		age = i-population[ind]['born']
		if population[ind]['births'] == 0:
			if random.randint(0,100) < first_birth_prob[age]:
				daughter = "".join(datetime.utcnow().strftime('%Y%m%d%H%M%S.%f').split('.'))
#				print "# DAY: %.1f - %s (age: %s) gives birth (%s) to %s" %((float(i)/24), ind,age,population[ind]['births']+1, daughter)
				population[ind]['births'] += 1
				population[ind]['gave_birth'] = i
				born+=1
				population[daughter] = {'mother':ind, 'births': int(0), 'sex':'female', 'born':i, 'gave_birth':0, 'genome':reproduce(mode='apomixis', parentA=population[ind]['genome'], hom_prob=apomixis_probs)}
				sex_distribution[population[daughter]['sex']] += 1
				ages.append(0)
				birth_counts[population[ind]['births']] += 1

#				print "Add new genome - first born"
				for j in range(len(genome)):
					alleles[j] += population[daughter]['genome'][j]
					if population[daughter]['genome'][j] == 'ab':
						av_heterozygosity[j] += 1
#				print ind,population[ind]

				
		
		else:
			if (random.randint(0,100) < sec_and_sub_prob[i-population[ind]['gave_birth']]):
				daughter = "".join(datetime.utcnow().strftime('%Y%m%d%H%M%S.%f').split('.'))
#				print "# DAY: %.1f - %s (age: %s) gives birth (%s) to %s" %((float(i)/24), ind,age,population[ind]['births']+1, daughter)
	                        population[ind]['births'] += 1
				population[ind]['gave_birth'] = i
				born+=1
#				population[daughter] = {'mother':ind, 'births': int(0), 'sex':'female', 'born':i, 'gave_birth':0, 'genome':reproduce(mode='sexual', parentA=population[ind]['genome'], hom_prob=[], parentB=population[ind]['genome'])}#arf_probs)}
				population[daughter] = {'mother':ind, 'births': int(0), 'sex':'female', 'born':i, 'gave_birth':0, 'genome':reproduce(mode=mode, parentA=population[ind]['genome'], hom_prob=probs)}
#				if population[ind]['sex'] == 'female':
				population[ind]['sex'] = 'male'
#				elif population[ind]['sex'] == 'male':
#					population[ind]['sex'] = 'sex'
				sex_distribution[population[daughter]['sex']] += 1
				ages.append(0)
				birth_counts[population[ind]['births']] += 1


#				print "Add new genome - second born"
				for j in range(len(genome)):
					alleles[j] += population[daughter]['genome'][j]
					if population[daughter]['genome'][j] == 'ab':
						av_heterozygosity[j] += 1
#				print ind,population[ind]


		if random.randint(0,100) < death_prob[age]:
#			print "# DAY: %.1f - CHALLENGE: %s\tage: %s\tsex: %s" %((float(i)/24), ind, age, population[ind]['sex'])
			
			r = random.randint(1,100)
			if r < death_prob[age]:
#				print "# DAY %.1f - %s dies (age: %s) (prop = %s; random number: %s)" %((float(i)/24), ind, age, death_prob[age], r)
#				print ind,population[ind]
				del population[ind]
				died+=1
				continue
			else:
#				print "Phew!"
				pass

		sex_distribution[population[ind]['sex']] += 1
		ages.append(age)
#		print "Add genome"
		for j in range(len(genome)):
			alleles[j] += population[ind]['genome'][j]
			if population[ind]['genome'][j] == 'ab':
				av_heterozygosity[j] += 1

	if len(population) == 0:
		print "[ %s - POP. EXTINCTION ]\tdays: %.1f" %(fish, float(i)/24)
		ages.append(0)
		fema_perc = float(0)
		male_perc = float(0)
		break
	else:
		frequencies = []
		percent_homozygous = float(0)
		
		for j in range(len(genome)):
			counts = defaultdict(int)
			counts['a'] = 0
			counts['b'] = 0
			for allele in alleles[j]:
				for a in allele.split():
					counts[a] += 1

			tot_num_alleles = float(0)
			for key in counts:
				tot_num_alleles += counts[key]
#			tot_num_alleles = float(counts[sorted(counts)[0]]+counts[sorted(counts)[1]])
#			print "%s/%s = %.2f/%.2f" %(sorted(counts)[0],sorted(counts)[1],counts[sorted(counts)[0]]/tot_num_alleles*100,counts[sorted(counts)[1]]/tot_num_alleles*100)
#			print tot_num_alleles,counts
			freq = frequencies.append("%.2f / %.2f" %((counts['a']/tot_num_alleles*100),(counts['b']/tot_num_alleles*100)))
			if counts['a'] == 0 or counts['b'] == 0:
				percent_homozygous += 1 

			av_heterozygosity[j] = av_heterozygosity[j]/len(population)*100

		percent_homozygous = percent_homozygous/len(genome)*100		

#		print "###POPULATION SIZE AFTER %.1f DAYS ( %s hours): %s\t\tborn: %s\t\tdied: %s\t\tmales: %.2f %% \t\tfemales: %.2f %%\t\tavg.age: %.0f" %((float(i)/24), i, len(population), born, died, float(sex_distribution['male'])/len(population)*100, float(sex_distribution['female'])/len(population)*100, np.mean(ages)) 

		if not len(population) == 0:
			male_perc = float(sex_distribution['male'])/len(population)*100
			fema_perc = float(sex_distribution['female'])/len(population)*100
		else:
			male_perc = float(0)
			fema_perc = float(0)

		if not last_only:
			print "[ %s - VERBOSE OUTPUT ]\tdays: %.1f\thours: %s hours\tpopulation size: %s\tborn: %s\tdied: %s\tmales: %.2f %%\tfemales: %.2f %%\tavg.age: %.0f\tpressure: %s\t" %(fish, (float(i)/24), i, len(population), born, died, male_perc, fema_perc, np.mean(ages), pressure), 
			print "N loci: %s\tpercent fixed: %.2f\taverage heterozygosity: %.2f\tfirst births: %s\tsecond births: %s\tthird births: %s\tfourth birth: %s" %(len(genome),percent_homozygous, np.mean(av_heterozygosity), birth_counts[1], birth_counts[2], birth_counts[3], birth_counts[4])
		
#			print "percent homozygous: %.2f\t%s" %(percent_homozygous, frequencies)
#			print "percent fixed: %.2f" %(percent_homozygous)
#			print av_heterozygosity
#			print np.mean(av_heterozygosity)
		else:
			if np.mean(av_heterozygosity) == 0:
				if not all_homo:
					final_out = "[ %s - ALL HOMOZYGOUS ]\tdays: %.1f\thours: %s hours\tpopulation size: %s\tmales: %.2f %%\tfemales: %.2f %%\tavg.age: %.0f\tpressure: %s\t" %(fish, (float(i)/24), i, len(population), male_perc, fema_perc, np.mean(ages), pressure)
	                        	final_out += "N loci: %s\tpercent fixed: %.2f\taverage heterozygosity: %.2f\tfirst births: %s\tsecond births: %s\tthird births: %s\tfourth birth: %s" %(len(genome),percent_homozygous, np.mean(av_heterozygosity), birth_counts[1], birth_counts[2], birth_counts[3], birth_counts[4])
					all_homo = True
			

if last_only:
	if final_out:
		print final_out
	print "[ %s - SIMULATION DONE ]\tdays: %.1f\thours: %s \tpopulation size: %s\tmales: %.2f %%\tfemales: %.2f %%\tavg.age: %.0f h\tpressure: %s\t" %(fish, (float(i)/24), i, len(population), male_perc, fema_perc, np.mean(ages), pressure), 
		
#	print "percent homozygous: %.2f\t%s" %(percent_homozygous, frequencies)
	print "N loci: %s\tpercent fixed: %.2f\taverage heterozygosity: %.2f\tfirst births: %s\tsecond births: %s\tthird births: %s\tfourth birth: %s" %(len(genome),percent_homozygous, np.mean(av_heterozygosity), birth_counts[1], birth_counts[2], birth_counts[3], birth_counts[4])

#print "\n###POPULATION SIZE AFTER %.1f DAYS: %s" %((float(i)/24), len(population))









import base64
def encode(key, clear):
    enc = []
    for i in range(len(clear)):
        key_c = key[i % len(key)]
        enc_c = chr((ord(clear[i]) + ord(key_c)) % 256)
        enc.append(enc_c)
    return base64.urlsafe_b64encode("".join(enc))

def decode(key, enc):
    dec = []
    enc = base64.urlsafe_b64decode(enc)
    for i in range(len(enc)):
        key_c = key[i % len(key)]
        dec_c = chr((256 + ord(enc[i]) - ord(key_c)) % 256)
        dec.append(dec_c)
    return "".join(dec)

sponge = []
for i in range(100):
	dt = "".join(datetime.utcnow().strftime('%Y%m%d%H%M%S.%f').split('.'))
	encoded = encode('pass',dt)
	decoded = decode('pass',encoded)
#	print dt,encoded,decoded
	sponge.append(encoded)


