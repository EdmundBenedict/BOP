#!/usr/bin/env python

from numpy import *
from small_utils import *





nb = 7
nt = 3


def symm(a):
    for i in range(a.shape[0]):
        for j in range(i,a.shape[1]):
            a[i,j] = a[j,i]
    return a

def pa(n,m,f=' %6.3f'):
    print n+':\n', sprint_nar(m, frm=f)

                
nn = nb+nt

random.seed(64562)

a11 = random.random((nb,nb)) + 0j
a12 = random.random((nb,nt)) + 0j
a22 = random.random((nt,nt)) + 0j

a11 = symm(a11)
a22 = symm(a22)

a = empty((nn,nn))

a[:nb, :nb] = a11
a[nb:, :nb] = a12.transpose()
a[:nb, nb:] = a12
a[nb:, nb:] = a22

pa('a',a)

nia = linalg.inv(a)
pa('nia',nia)


ia11 = linalg.inv(a11)

b = dot(ia11,a12)

ch = a22 - dot(a12.transpose(),b)
w,v = linalg.eigh(ch)

iw = 1/w

c = empty((nt,nt))
for i in range(nt):
    for j in range(i,nt):
        c[i,j] = sum(v[i,:]*iw[:]*v[j,:])
        c[j,i] = c[i,j]

iws = sqrt(iw)        
cs = empty((nt,nt))
for i in range(nt):
    for j in range(i,nt):
        cs[i,j] = sum(v[i,:]*iws[:]*v[j,:])
        cs[j,i] = cs[i,j]

pa('cs',cs)
bcs = dot(b,cs)
        
bcb = dot(bcs,bcs.transpose())

bc = dot(b,c)

pia = empty((nn,nn))

pia[:nb, :nb] = ia11 + bcb
pia[:nb, nb:] = -bc
pia[nb:, :nb] = pia[:nb, nb:].transpose()
pia[nb:, nb:] = c

pa('pia',pia)

diff = pia-nia

pa('diff',diff)

print sqrt(sum(diff*diff))












