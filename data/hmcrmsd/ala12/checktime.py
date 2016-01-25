#!/usr/bin/env python



# check the continuity of rmsd.log


import os, sys, glob




fnrmsdlog = "rmsd.log"

s = open(fnrmsdlog).readlines()
print "checking %s" % fnrmsdlog
dt = 1000
tm = 0
for ln in s:
  newtm = int( ln.split()[0] )
  if newtm - tm != dt:
    print "time corruption at time t %s -> %s" % (tm, newtm)
  tm = newtm


fnmdlog = "md.log"
s = open(fnmdlog).readlines()
print "checking %s" % fnmdlog
dt = 10000
tm = -dt
i = 0
n = len(s)
err = 0
while i < n - 1:
  ln = s[i].strip()
  if not ln.startswith("Step   "):
    i += 1
    continue

  newtm = int( s[i+1].strip().split()[0] )

  if newtm != tm + dt:
    print "time %s -> %s" % (tm, newtm)
    err += 1

  tm = newtm
  i += n

if not err:
  print "%s is fine!" % fnmdlog
