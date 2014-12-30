#!/bin/python
from numpy import *;
import random;

def roulette_pick(p):
    v = random.random();
    s = 0;
    for i in xrange(len(p)):
        s += p[i];
        if s >= v:
            return i;
    return len(p)-1;

