#!/usr/bin/python2

import argparse
import os
import string
from inspect import getsourcefile

mydir=os.path.dirname(os.path.abspath(getsourcefile(lambda:0)))


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("type",type=int)
    parser.add_argument("lat",type=float)
    parser.add_argument("lon",type=float)
    parser.add_argument("depth",type=float)
    parser.add_argument("mag",type=float)
    parser.add_argument("stk",type=float)
    parser.add_argument("dip",type=float)
    parser.add_argument("rak",type=float)
    parser.add_argument("dt",type=float)

    args=parser.parse_args()

    with open(os.path.join(mydir,"setSrfParams.py.template"),'r') as fr:
       lines = fr.readlines()
       str_lines=''.join(lines)
       str_lines2= string.Template(str_lines).substitute({'type':args.type,'lat':args.lat,'lon':args.lon,'depth':args.depth,'mag':args.mag, 'stk':args.stk, 'dip':args.dip,'rak':args.rak, 'dt':args.dt})

       with open(os.path.join(os.path.curdir,'setSrfParams.py'),'w') as fw:
           fw.write(str_lines2)

if __name__=='__main__':
    main()
