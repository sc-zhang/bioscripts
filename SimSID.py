#!/usr/bin/env python
import os
import sys
import argparse
import random


def ArgParser():
	group = argparse.ArgumentParser()
	group.add_argument('-s', '--snp', type=float, help='snp ratio of whole genome, percentage, default: 0.01', default=0.01)
	group.add_argument('-i', '--insertion', type=float, help='insertion ratio of whole genome, percentage, default: 0.01', default=0.01)
	group.add_argument('--insert_length', type=int, help='max length of insertion, default: 10', default=10)
	group.add_argument('-d', '--deletion', type=float, help='delection ratio of whole genome, percentage, default: 0.01', default=0.01)
	group.add_argument('--delete_length', type=int, help='max length of deletion, default: 10', default=10)
	group.add_argument('--random_length', action="store_true", help='use this argument for generate random length of indels', default=False)
	group.add_argument('-v', '--verbose', action="store_true", help='print detail information', default=False)
	group.add_argument('-r', '--ref', help='origin fasta file of genome', required=True)
	group.add_argument('-o', '--out', help='prefix of simulated data', required=True)

	return group.parse_args()


def ReadFASTA(inputFASTA):
	fastaDB = {}
	posDB = {}
	with open(inputFASTA, 'r') as fIN:
		id = ''
		seq = ''
		for line in fIN:
			if line[0] == '>':
				if seq != '':
					fastaDB[id] = seq
				id = line.strip()[1:]
				seq = ''
			else:
				seq += line.strip()
		fastaDB[id] = seq
	for chrn in fastaDB:
		posDB[chrn] = [1]*len(fastaDB[chrn])
	return fastaDB, posDB


def IsInRegions(regionDB, queryPos):
	s = 0
	e = len(regionDB)-1
	refRegions = sorted(regionDB)
	if len(refRegions) == 0:
		return False
	while s<=e:
		mid = int((s+e)/2)
		if refRegions[mid][0] < queryPos:
			s = mid+1
		elif refRegions[mid][0] > queryPos:
			e = mid-1
		else:
			return True
	if refRegions[e][1] >= queryPos:
		return True
	else:
		return False


def GenDelRegions(fastaDB, posDB, delRatio, delLength, isRandom, isVerbose):
	delRegions = {}
	for chrn in fastaDB:
		chrLen = len(fastaDB[chrn])
		if isRandom:
			avgDelLen = int(delLength/2)
		else:
			avgDelLen = delLength
		cntDel = int(chrLen*delRatio/avgDelLen)
		print("%s\tdelections count: %d"%(chrn, cntDel))
		delRegions[chrn] = []
		for i in range(0, cntDel):
			if isVerbose:
				print("Generating: %d"%(i+1))
			if isRandom:
				curDelLen = random.randint(1, delLength)
			else:
				curDelLen = delLength
			sp = random.randint(0, chrLen-curDelLen)
			ep = sp+curDelLen
			while posDB[chrn][sp] == 0 or posDB[chrn][ep-1] == 0:
				sp = random.randint(0, chrLen-curDelLen)
				ep = sp+curDelLen
			delRegions[chrn].append([sp, ep])
			for i in range(sp, ep):
				posDB[chrn][i] = 0
	for chrn in delRegions:
		delRegions[chrn] = sorted(delRegions[chrn])
	
	return delRegions


def GenSeq(seqLen):
	nucType = ['A', 'T', 'G', 'C']
	seq = ''
	for i in range(0, seqLen):
		seq += nucType[random.randint(0, 3)]
	return seq


def GenInsPosSeqs(delRegions, fastaDB, posDB, insRatio, insLength, isRandom, isVerbose):
	insList = {}
	insSeqs = {}
	for chrn in fastaDB:
		chrLen = len(fastaDB[chrn])
		if isRandom:
			avgInsLen = int(insLength/2)
		else:
			avgInsLen = insLength
		cntIns = int(chrLen*insRatio/avgInsLen)
		print("%s\tinsertions count: %d"%(chrn, cntIns))
		insList[chrn] = []
		insSeqs[chrn] = []
		for i in range(0, cntIns):
			if isVerbose:
				print("Generating: %d"%(i+1))
			if isRandom:
				curInsLen = random.randint(1, insLength)
			else:
				curInsLen = insLength
			pos = random.randint(0, chrLen-1)
			while posDB[chrn][pos] == 0:
				pos = random.randint(0, chrLen-1)
			insList[chrn].append(pos)
			posDB[chrn][pos] = 0
			insSeqs[chrn].append(GenSeq(curInsLen))
	return insList, insSeqs


def GenSNPPos(delRegions, fastaDB, posDB, snpRatio, insPos, isVerbose):
	snpSeq = {}
	snpPos = {}
	nucType = ['A', 'T', 'G', 'C']
	for chrn in fastaDB:
		snpSeq[chrn] = []
		snpPos[chrn] = []
		chrLen = len(fastaDB[chrn])
		cntSNP = int(chrLen*snpRatio)
		print("%s\tSNPs count: %d"%(chrn, cntSNP))
		for i in range(cntSNP):
			if isVerbose:
				print("Generating: %d"%(i+1))
			pos = random.randint(0, chrLen-1)
			while posDB[chrn][pos] == 0:
				pos = random.randint(0, chrLen-1)
			SNP = nucType[random.randint(0, 3)]
			while SNP == fastaDB[chrn][pos]:
				SNP = nucType[random.randint(0, 3)]
			snpSeq[chrn].append(SNP)
			snpPos[chrn].append(pos)
			posDB[chrn][pos] = 0
	return snpPos, snpSeq


def SimSID(snpRatio, insRatio, delRatio, insLength, delLength, isRandom, isVerbose, inputFASTA, outPrefix):
	print("SNP Ratio = %.2f%%\nINS Ratio = %.2f%%\nDEL Ratio = %.2f%%\nINS Length = %d\nDEL Length = %d\nRandom: %s\nVerbose: %s\nInput file: %s\nOut prefix: %s"%(snpRatio*100, insRatio*100, delRatio*100, insLength, delLength, isRandom, isVerbose, inputFASTA, outPrefix))
	random.seed()
	print("Reading fasta")
	fastaDB, posDB = ReadFASTA(inputFASTA)
	print("Generating deletions")
	delRegions = GenDelRegions(fastaDB, posDB, delRatio, delLength, isRandom, isVerbose)

	print("Generating insertions")
	insPos, insSeq = GenInsPosSeqs(delRegions, fastaDB, posDB, insRatio, insLength, isRandom, isVerbose)

	print("Generating SNPs")
	snpPos, snpSeq = GenSNPPos(delRegions, fastaDB, posDB, snpRatio, insPos, isVerbose)
	
	print("Writing infomations")
	with open(outPrefix+"_snps.txt", 'w') as fSNP:
		fSNP.write("Chromosome\tPosition\tOrigin\tNew\n")
		writeStrings = []
		for chrn in sorted(fastaDB):
			for i in range(0, len(snpPos[chrn])):
				pos = snpPos[chrn][i]
				writeStrings.append([chrn, pos+1, fastaDB[chrn][pos], snpSeq[chrn][i]])
		for wString in sorted(writeStrings):
			fSNP.write('\t'.join(list(map(str, wString)))+'\n')
	
	with open(outPrefix+"_indel.txt", 'w') as fIndel:
		fIndel.write("Chromosome\tPosition\tOrigin\tNew\n")
		writeStrings = []
		for chrn in sorted(fastaDB):
			for i in range(0, len(delRegions[chrn])):
				sp = delRegions[chrn][i][0]
				ep = delRegions[chrn][i][1]
				writeStrings.append([chrn, sp+1, fastaDB[chrn][sp: ep], '[]'])
		for chrn in sorted(fastaDB):
			for i in range(0, len(insPos[chrn])):
				pos = insPos[chrn][i]
				writeStrings.append([chrn, pos+1, '[]', insSeq[chrn][i]])
		for wString in sorted(writeStrings):
			fIndel.write('\t'.join(list(map(str, wString)))+'\n')
	
	print("Writing fasta")
	with open(outPrefix+"_sim.fasta", 'w') as fSim:
		for chrn in sorted(fastaDB):
			newSeq = list(fastaDB[chrn])
			for i in range(0, len(snpPos[chrn])):
				newSeq[snpPos[chrn][i]] = snpSeq[chrn][i]
			for i in range(0, len(insPos[chrn])):
				newSeq[insPos[chrn][i]] = insSeq[chrn][i] + newSeq[insPos[chrn][i]]
			for i in range(0, len(delRegions[chrn])):
				for j in range(delRegions[chrn][i][0], delRegions[chrn][i][1]):
					newSeq[j] = ''
			fSim.write(">%s\n%s\n"%(chrn, ''.join(newSeq)))
	print("Success")


if __name__ == "__main__":
	opts = ArgParser()
	snpRatio = opts.snp/100.0
	insRatio = opts.insertion/100.0
	delRatio = opts.deletion/100.0
	inputFASTA = opts.ref
	outPrefix = opts.out
	insLength = opts.insert_length
	delLength = opts.delete_length
	isRandom = opts.random_length
	isVerbose = opts.verbose
	SimSID(snpRatio, insRatio, delRatio, insLength, delLength, isRandom, isVerbose, inputFASTA, outPrefix)
