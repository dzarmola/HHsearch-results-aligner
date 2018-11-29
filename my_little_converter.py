#!/usr/bin/env python2

import sys
import glob
import re

def from_save(filename):
    with open(filename) as save:
#        clusters = eval(save.readline())
#        labels = eval(save.readline())
#        data = eval(save.readline())
#        lens = eval(save.readline())
        clusters,labels,data,lens = map(eval,save.readlines())
    return data,labels,clusters

def read_hhm(fname):
    with open(fname) as input:
        data = input.read()
    data = re.split("#",data)[1]
    seq = re.findall("([A-Z-])\s[0-9]+",data)
    return seq

def get(dict,key):
    for k,v in dict.items():
        if key in k or k in key:
            return v
    raise IndexError("No such key:{}".format(key))

def main(sfile,dir):
    data,labels,clusters = from_save(sfile)
    profiles = {}
    for file in glob.glob(dir+"/*.hhm"):
        name = file.split("/")[-1].strip('.hhm')
        name = name.strip('_core')
        profiles[name] = read_hhm(file)
    for l,d in zip(labels,data):
        print ">"+l
        row = []
        for _ in d:
            if _ is None:
                row.append("-")
            else:
                row.append(get(profiles,str(l))[_-1])
        print "".join(row)

if __name__ == "__main__":
    main(sys.argv[1],sys.argv[2]) 
