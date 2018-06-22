#!/usr/bin/env python
#coding: utf8


####import module
from __future__ import division
import sys, os
import re
import glob
import gzip

#####Document Description
''' This script is created for calculating duplicated reads Rate that occurred in chip 715 and 800
    the script won't work if certain format is not met.
'''
##########log#############
'''
    0.2version -- bugs fixed for block counting; efficiency improved
'''
#####Acknowledgement
''' I would like express my deepest apprecitation to my colleague, Dongyang Xu, whose contribution means a lot to my work.
'''
Date = '2018-04-03'
Version='0.2.0'

usage = '''
      Version %s Ya Ding %s

      Usage: %s <fastq.gz> <chip type> >STDOUT
''' %(Version, Date, os.path.basename(sys.argv[0]))




def getfovIdx(fq):

    idline = gzip.open(fq,'r').readline().decode()

    m = re.search('C\d{3}R\d{3}', idline)
    return m.start()

def generateFovDict(fq,idx):
    fovDict = {}

    fh = gzip.open(fq)
    rec = [fh.readline().decode() for i in range(4)]
    fov = rec[0][idx:idx+8]
    rid = int(rec[0].split('/')[0][idx+8:])
    fovDict[rid] = rec[1].strip()
    while True:
        rec = [fh.readline().decode() for i in range(4)]
        if not rec[0]:
            break
        newFov = rec[0][idx:idx+8]
        rid = int(rec[0].split('/')[0][idx+8:])
        if newFov != fov:
            yield fovDict
            ## reinit fovDict
            fovDict = {}
            fov = newFov
            fovDict[rid] = rec[1].strip()
        else:
            fovDict[rid] = rec[1].strip()
    fh.close()
    yield fovDict

def block_find(ID):
    a=[]
    b=[]

    for i in range(0,8):
        for j in range(0,8):
            a.append(y[i]*x[j])

    for i in range(1,65):
            b.append(sum(a[:i])-1)

    if ID <= x[0]*y[0]-1:
        blk= 1
    else:
        for q in range(0,63):
            if ID<=b[q+1] and ID>b[q]:
                blk= q+2
    return blk,b

def concordinate(blk):

    if blk%8 ==0:
        Xpo=7
    else:
        Xpo=blk%8-1

    Ypo=(blk-1)//8
    return Xpo, Ypo

def neighbor_find(ID):
    blk, b=block_find(ID)
    Xpo, Ypo=concordinate(blk)
    b.insert(0,-1)

    L=[]
    R=[]
    for i in range(1,y[Ypo]-1):
      L.append(b[blk-1]+1+i*x[Xpo])
      R.append(b[blk-1]+(i+1)*x[Xpo])

    if ID==b[blk-1]+1:
      return b[blk-1]+2, b[blk-1]+1+x[Xpo]
    elif ID==b[blk-1]+x[Xpo]:
        return b[blk-1]+x[Xpo]-1, b[blk-1]+2*(x[Xpo])
    elif ID==b[blk]-x[Xpo]+1:
        return b[blk]-x[Xpo]+2, b[blk]-2*(x[Xpo])+1
    elif ID ==b[blk]:
        return b[blk]-1, b[blk]-x[Xpo]
    elif ID in range(b[blk-1]+2,b[blk-1]+x[Xpo]):
        return ID-1, ID+1, ID+x[Xpo]
    elif ID in range(b[blk]-x[Xpo]+2, b[blk]):
        return ID-1, ID+1, ID-x[Xpo]
    elif ID in L:
        return ID-x[Xpo], ID+1, ID+x[Xpo]
    elif ID in R:
        return ID-x[Xpo], ID-1, ID+x[Xpo]
    else:
        return ID-x[Xpo], ID+1, ID+x[Xpo],ID-1

def Dup_Rate(read_ID_dict):
    count=0
    for ID in read_ID_dict.keys():
        subread=read_ID_dict[ID][:9]
        sub=neighbor_find(ID)
        for i in sub:
            if i in read_ID_dict.keys():
                if read_ID_dict[i][:9] == subread:
                    count=count+1
                else:
                    continue
            else:
                continue
    return count/2/TotalDNB


def getFovDicts_AvgDup(fq):
    ## get fov idx
    idx = getfovIdx(fq)
    Avg_Dup=[]
    for dic in generateFovDict(fq, idx):
        # print (len(dic.keys()))
        Avg_Dup.append(Dup_Rate(dic))
    return sum(Avg_Dup)/len(Avg_Dup)




def main():
    if len(sys.argv) !=3:
        print(usage)
        sys.exit(1)
    fq=sys.argv[1]
    global x,y,TotalDNB
    if sys.argv[2] == '715':
        x=[63,85,195,217,217,195,85,63]
        y=[51,69,141,195,195,177,33,51]
        TotalDNB=1021440
    else:
        x=[67,109,165,193,193,165,109,67]
        y=[45,61,125,173,173,125,61,45]
        TotalDNB=862944 

    print('The duplicated reads rate is  ' + str(getFovDicts_AvgDup(fq)))

if __name__ == '__main__':
   main()
