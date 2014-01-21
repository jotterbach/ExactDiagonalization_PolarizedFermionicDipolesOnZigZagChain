    # -*- coding: utf-8 -*-
"""
Created on Mon Sep 16 11:26:39 2013

@author: jotterbach
"""

from multiprocessing import Process
from multiprocessing import *
from multiprocessing import Pool
import os
from time import *
from numpy import *
from resource import *


def sleeping(secs,x):
    #info('function sleeping')
   # print secs, ' sleeping', x
    sleep(secs)
   # print 'awake ', secs
    return secs, x*secs

def info(title):
    print title
    #print 'module name:', __name__
    if hasattr(os, 'getppid'):  # only available on Unix
        print 'parent process:', os.getppid()
    print 'process id:', os.getpid()

def f(name):
    #info('function f')
    print 'hello', name

if __name__ == '__main__':
    #info('main line')
    #print multiprocessing.cpu_count()
    rsrc = RLIMIT_DATA
    soft, hard = getrlimit(rsrc)
    print soft, hard
    setrlimit(rsrc, (1024, hard))
    rsrc = RLIMIT_DATA
    print getrlimit(rsrc)    
    
    v = zeros([8,1])
    print 'cpu_count() = %d\n' % cpu_count()
    pool = Pool()
    it = [pool.apply_async(sleeping, (i,2)) for i in range(2,6)]
    for r in it:
        print '\n', r.get(), r.get()[0], r.get()[1]
        v[r.get()[0]] = r.get()[1]

    print v
    
    
 
    
    
    
    