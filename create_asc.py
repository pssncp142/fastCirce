from numpy import *

f_name = 'asc/2015-06-17-%04d_asc.txt'

arr = arange(55, 65)

shfx, shfy = 100, 6
xx = zeros(6) + 100
yy = zeros(6) + 6



for i in arr:

    ref = loadtxt('asc_92.txt')
    nd = 62 - i
    ref[:,1] = -xx*nd + ref[:,1]
    ref[:,2] = -yy*nd + ref[:,2]
    savetxt(f_name % i, ref)
