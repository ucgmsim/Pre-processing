#!/usr/bin/python2

import argparse
import os
import string
from inspect import getsourcefile
import sys
mydir=os.path.dirname(os.path.abspath(getsourcefile(lambda:0)))
template_dir = os.path.join(mydir,'setSrfParams.py.template')

default_prefix = 'Srf/source'
default_stoch  = 'Stoch'
default_SEED_INC = 1

def gen_params(params_all,f):
    src_type = int(params_all['TYPE'])
    #print type(src_type)
    order = get_order(src_type)
    
    ##Generic params
    ##type,prefix,stoch
    f.write("## Sets Variables for SRF/Stoch Generation\n")
    f.write("#"*20+"\n")
    f.write("#!!!Important!!\n")
    f.write("#"*20+"\n")
    f.write("#Please change the values accordingly if '= None'\n\n")
    f.write("# TYPE:\n# 1: point source to point source srf\n")
    f.write("# 2: point source to finite fault srf\n")
    f.write("# 3: finite fault to finite fault srf\n")
    f.write("# 4: multi-segment finite fault srf\n")
    f.write("TYPE = %s\n"%params_all['TYPE'])
    f.write("# specify basename for gsf, srf and stoch file created\n")
    f.write("# PREFIX for gsf/srf/stoch files\n")
    f.write("# if prefix ends with '_', automatic naming follows\n") 
    if params_all['PREFIX'] == None:
        params_all['PREFIX'] = default_prefix #set prefix to default value
    f.write("PREFIX = '%s'\n"%params_all['PREFIX'])
    f.write("# directory for Stoch file(s)\n")
    f.write("# set to None to not produce the stoch file\n")
    if params_all['STOCH'] == None:
        params_all['STOCH'] = default_stoch
    f.write("STOCH = '%s'\n"%params_all['STOCH'])

    ##params for 1,2,3
    if src_type in [1,2,3]:
        f.write("# latitude (float)\n")
        f.write("LAT = %s\n"%params_all['LAT'])
        f.write("# longitude (float)\n")
        f.write("LON = %s\n"%params_all['LON']) 
        f.write("# depth (float)\n")
        f.write("DEPTH = %s\n"%params_all['DEPTH'])
        f.write("# magnitude (float)\n")
        f.write("MAG = %s\n"%params_all['MAG'])
        f.write("# strike (int)\n")
        f.write("STK = %s\n"%params_all['STK'])
        f.write("# dip (int)\n")
        f.write("DIP = %s\n"%params_all['DIP'])
        f.write("# rake (int)\n")
        f.write("RAK = %s\n"%params_all['RAK'])
        f.write("# rupture timestep (float)\n")
        f.write("DT = %s\n"%params_all['DT'])
    #type 1 specific 
    if src_type == 1:
        f.write("# specify seismic moment directly (-1 to use magnitude)\n")
        f.write("MOM = %s\n"%params_all['MOM'])
    
    #params for 2,3,4
    if src_type in [2,3,4]:
        f.write("# GENSLIP VERSION:\n# '3.3', '5.2.3a'\n")
        f.write("GENSLIP = '%s'\n"%params_all['GENSLIP'])
        f.write("# roughness of fault only for genslip 5.0+\n# 0.1 is a good value to use - Rob Graves\n# 0.0050119 = (10^(-2.3)) - Shi & Day 2014. Used as default in 5.2.3a\n")
        f.write("ROUGH = %s\n"%params_all['ROUGH'])
        f.write("# (float), ex. RVFRAC = 0.8\n")
        f.write("RVFRAC = %s\n"%params_all['RVFRAC'])
        f.write("# (float), ex. SLIP_COV = 0.85\n")
        f.write("SLIP_COV = %s\n"%params_all['SLIP_COV'])
        f.write("# seed for stoch\n")
        f.write("SEED = %s\n"%params_all["SEED"])
    if src_type == 4:
        if params_all['SEED_INC'] == None:
            params_all['SEED_INC'] = default_SEED_INC
        f.write("# only used in type 4 to go over multiple seeds\nSEED_INC =%s\n"%params_all['SEED_INC'])
    if src_type == 2:
        f.write("# Mw Scaling Relation (string), one of:\n")
        f.write("# HanksBakun2002, BerrymanEtAl2002, VillamorEtAl2001\n")
        f.write("MWSR = '%s\n'"%params_all[MWSR])
    if src_type == 3:
        f.write("# (float), ex. 45.7. Should be the same as FWID\n")
        f.write("FLEN = %s\n"%params_all['FLEN'])
        f.write("# (float), ex. 0.10. Should be the same as DWID\n")
        f.write("DLEN = %s\n"%params_all['DLEN'])
        f.write("# (float), ex. 45.7")
        f.write("FWID = %s\n"%params_all['FWID'])
        f.write("# (float), ")
        f.write("DWID = %s\n"%params_all['DWID'])
        f.write("# (float), ex. 0.0")
        f.write("DTOP = %s\n"%params_all['DTOP'])
        f.write("# km east\nSHYPO = %s\n"%params_all['SHYPO'])
        f.write("# km down\nDHYPO = %s\n"%params_all['DHYPO'])
    if src_type == 4:
        f.write("# how many scenarios to run. different scenarios can have randomised parameters (randomisation amount set below)\n")
        f.write("# N_SCENARIOS = 5")
        f.write("N_SCENARIOS = %s\n"%params_all['N_SCENARIOS'])
        f.write("# how many scenarios to run before incrementing the seed value.\n")
        f.write("# this will randomise parameters but keep the seed the same.\n")
        f.write("# e.g N_SEED_INC = 20\n")
        f.write("N_SEED_INC = %s\n"%params_all['N_SEED_INC'])

        #array type params
        f.write("\n\n# the array/list size of each params below must be identical\n")
        f.write("# description of main segments\n")
        f.write("# e.g, len(CASES) == len(M_MAG) && len(M_MAG) == len(M_NSEG) && ...\n")
        f.write("# example: CASES = ['NNWEpic', 'CharX', 'GrNW', 'GrEW', 'GrEast', 'Horotar', 'NEStepOv']\n")
        f.write("CASES = None\n")
        f.write("# master segments can be made up of multiple sub-segments\n")
        f.write("# example: M_MAG = [5.9, 6.5, 6.4, 6.7, 6.1, 5.9, 3.35]\n")
        f.write("M_MAG = None\n")
        f.write("# as above, if MOM invalid ( <= 0), calculated from MAG\n")
        f.write("# eg. M_MOM = [-1, -1, -1, -1, -1, -1, -1]\n")
        f.write("M_MOM = None\n")
        f.write("# how many segments each master is made of\n")
        f.write("# eg. M_NSEG = [1, 1, 1, 1, 1, 1, 1]\n")
        f.write("M_NSEG = None\n")
        f.write("# e.g M_SEG_DELAY = [\"0\", \"1\", \"1\", \"1\", \"1\", \"1\", \"1\"]\n")
        f.write("M_SEG_DELAY = None\n")
        f.write("# e.g M_RVFAC_SEG = [\"-1\", \"0.5\", \"0.5\", \"0.5\", \"0.5\", \"0.5\", \"0.5\"]\n")
        f.write("M_RVFAC_SEG = None\n")
        f.write("# e.g M_GWID = [\"-1\", \"8.0\", \"8.0\", \"8.0\", \"8.0\", \"8.0\", \"8.0\"]")
        f.write("M_GWID = None\n")
        f.write("#e.g M_RUP_DELAY = [0, 5, 11.0, 7.0, 18.0, 18.0, 17.0]\n")
        f.write("M_RUP_DELAY = None\n")
        
        # 2 dimention arrays
        f.write("\n\n# 2 dimentional arrays for each set of segments\n# demensions should match M_NSEG\n")
        f.write("# if values are same, shortcut is [value] * repetitions, e.g M_FLEN = [[8]]*6\n")
        f.write("# example values:\n")
        f.write("# M_FLEN = [[8], [10], [16], [20], [14], [7], [11]]\n")
        f.write("# M_DLEN = [ [0.2], [0.2], [0.2], [0.2], [0.2], [0.2], [0.2]]\n")
        f.write("# M_FWID = [[8], [15], [18], [18], [18], [8], [15]]\n")
        f.write("# M_DWID = [ [0.2], [0.2], [0.2], [0.2], [0.2], [0.2], [0.2]]\n")
        f.write("# M_DTOP = [ [2.0], [0.0], [0.0], [0.0], [0.0], [0.0], [0.0]]\n")
        f.write("# M_STK  = [[150], [35], [121], [87], [87], [216], [40]]\n")
        f.write("# M_RAK  = [[6], [90], [-180], [180], [180], [78], [171]]\n")
        f.write("# M_DIP  = [[54], [70], [105], [80], [78], [50], [80]]\n")
        f.write("# M_ELON = [[172.1811], [172.125], [172.035], [172.195], [172.380], [171.9434], [172.3115]]\n")
        f.write("# M_ELAT = [[-43.5095], [-43.57], [-43.581], [-43.590], [-43.573], [-43.5779], [-43.5509]]\n")
        f.write("# M_SHYPO= [[1], [5], [4], [-1], [-7], [-3.5], [-5.5]]\n")
        f.write("# M_DHYPO= [[4], [10], [6],[6], [6],[4], [6]]\n")



        f.write("M_FLEN = None\n")
        f.write("M_DLEN = None\n")
        f.write("M_FWID = None\n")
        f.write("M_DWID = None\n")
        f.write("M_DTOP = None\n")
        f.write("M_STK  = None\n")
        f.write("M_RAK  = None\n")
        f.write("M_DIP  = None\n")
        f.write("M_ELON = None\n")
        f.write("M_ELAT = None\n")
        f.write("M_SHYPO= None\n")
        f.write("M_DHYPO= None\n")
    ## Stochastic params that is used for type 3 and 4 only.
    # if type 3, len(var) == 1.
    # if type ==4 , len(var) == len(CASES)
    if src_type == 3 or src_type == 4:
        f.write("###\n### RELATING TO STOCHASTIC GENERATION\n###\n")
        f.write("### set variability, it will change the relative values between scenarios\n")
        f.write("\n# used for calculating M0total in magnitude variability, e.g MW_TOTAL = 7.2\n")
        f.write("MW_TOTAL = %s\n"%params_all['MW_TOTAL'])
        f.write("### for type 3 srf, variables below will have only one value and will be use as absolute variability\n")
        f.write("### for type 4 srf, variables should have the same length as CASES and will be use as relative moment variability.\n")
        f.write("# V_MAG = [0.0, 0.0, 0.0, 0.4, 0.0, 0.0, 0.0]\n")
        f.write("V_MAG = None\n")
        f.write("# fault width/length variability (proportion, eg. 0.1 = 10%)\n")
        f.write("# e.g V_FWID = [0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00]\n")
        f.write("# e.g V_FLEN = [0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00]")
        f.write("V_FWID = None\n")
        f.write("V_FLEN = None\n")
        f.write("# rupture time delay variability, requires dependency specification below\n")
        f.write("# absolute value\n")
        f.write("# V_RDELAY = [0, 1, 3, 0.5, 0, 0, 0]\n")
        f.write("V_RDELAY = None\n")
        f.write("# rupture dependency\n")
        f.write("# None means it is a starting point\n")
        f.write("# otherwise it shoould be the segment that triggered this segment\n")
        f.write("# e.g [None, 0,1,1,3] means the seg.0 is the starting point, seg.1 is triggered by seg.0. seg.2 seg.3 is triggered by seg.1, seg.4 is triggered by seg.3 \n")
        f.write("D_RDELAY = %s\n"%None)
        f.write("# hypocentre variabiity - absolute for first segment, [SHYPO, DHYPO]\n")
        f.write("V_HYPO = %s\n"%None)
        
    
    #f.write()
    #f.write()
    #f.write()

def get_default(param, file_path):
    with open(file_path,'r') as fr:
        #lines = fr.readlines;
        for line in fr:
            #only process lines without '#'
            if line[0] != '#':
                #filter out space and \n, then split the line with '='
                line_split = line.replace(' ','').replace('\n','').rsplit('=',1)
                param_in_temp = line_split[0]
                if param_in_temp == param:
                    if line_split[1][0] == '$':
                        #print "No default value for %s in template"%param
                        return None
                    param_default = line_split[1]
                    return param_default
    print "failed to find %s in template"%param
    return None

def get_order(src_type):
    if src_type == 1:
        return ['LAT', 'LON', 'DEPTH', 'MAG', 'MOM', 'STK', 'RAK', 'DIP', 'DT','PREFIX','STOCH']
    if src_type == 2:
        return ['LAT', 'LON', 'MAG', 'STK', 'RAK', 'DIP', 'DT', 'PREFIX', 'SEED', 'RVFRAC', 'ROUGH', 'SLIP_COV', 'DEPTH', 'MWSR', 'STOCH', 'GENSLIP']
    if src_type == 3:
        return ['LAT', 'LON', 'MAG', 'STK', 'RAK', 'DIP', 'DT', 'PREFIX','SEED', 'RVFRAC', 'ROUGH', 'SLIP_COV', 'FLEN', 'DLEN', 'FWID', 'DWID', 'DTOP','SHYPO', 'DHYPO', 'STOCH', 'GENSLIP']
    if src_type == 4:
        return ['M_NSEG', 'M_SEG_DELAY', 'M_MAG', 'M_MOM','M_RVFAC_SEG', 'M_GWID', 'M_RUP_DELAY', 'M_FLEN', 'M_DLEN', 'M_FWID', 'M_DWID', 'M_DTOP', 'M_STK', 'M_RAK', 'M_DIP', 'M_ELON', 'M_ELAT', 'M_SHYPO', 'M_DHYPO', 'DT', 'SEED', 'SEED_INC', 'RVFRAC', 'ROUGH', 'SLIP_COV', 'PREFIX', 'N_SCENARIOS', 'N_SEED_INC', 'CASES', 'GENSLIP', 'STOCH']

def get_order_from_temp(file_path):
    order = []
    with open(file_path,'r') as fr:
        for line in fr:
            if line[0] != '#' and line != '\n':
                line_split = line.replace(' ','').replace('\n','').rsplit('=',1)
                order.append(line_split[0])
    return order


def check_params(params_all):
    check_error = False
    for i in params_all:
        if params_all[i] == None:
            #print '%s not specified, attemp to read default value'%i
            i_default = get_default(i,template_dir)
            if i_default == None:
                print "Error: %s is missing and have no default value in template"%i
                check_error = True
            params_all[i] = i_default
    return check_error
        
def confirm_params(params_all):
    #print params_all['TYPE']
    order = get_order(int(params_all['TYPE']))
    if int(params_all['TYPE']) == 4:
        for i in order:
            if i in params_all:
                if params_all[i] != None:
                    print i,': ',params_all[i]
        
        print "Only params above are set.\nPlease manually modify the rest of params by editing setSrfParams.py.\n"
        return 
    #print params_all['TYPE']
    for i in order:
        if i in params_all:
            print i,': ',params_all[i]
    #for i in params_all:
    #    print i,': ',params_all[i]
    


def main():

    parser = argparse.ArgumentParser()

    #parser.add_argument("type",type=int, help= "the type of required Source file, e.g: 1~4")
    subparser = parser.add_subparsers(help='Choose the type of srf you want to generate, usage of each type use: gen_src_params.py type -h ', dest="TYPE")

    parser_type1 = subparser.add_parser('1', help='Generate type 1 srf (Point source srf).')
    parser_type2 = subparser.add_parser('2', help='Generate type 2 srf (Point source to Finite Fault srf).')
    parser_type3 = subparser.add_parser('3', help='Generate type 3 srf (Finite Fault descriptor to Finite Fault srf).')
    parser_type4 = subparser.add_parser('4', help='Generate type 4 srf (Multiple Segment Finite Fault srf). ')

    #all the params needed for type 1 (Point Source)
    parser_type1.add_argument('--LAT',type=float,nargs='?',default=None)
    parser_type1.add_argument('--LON',type=float,nargs='?',default=None)
    parser_type1.add_argument('--DEPTH',type=float,nargs='?',default=None)
    parser_type1.add_argument('--MAG',type=float,nargs='?',default=None)
    parser_type1.add_argument('--MOM',type=float, nargs='?',default=None)
    parser_type1.add_argument('--STK', type=float,nargs='?',default=None)
    parser_type1.add_argument('--RAK', type=float,nargs='?',default=None)
    parser_type1.add_argument('--DIP', type=float,nargs='?',default=None)
    parser_type1.add_argument('--DT', type=float, nargs='?',default=None)
    #can be read from template
    parser_type1.add_argument('--PREFIX', type=str,nargs='?', default=default_prefix, help='The path/prefix you want for SRF files. Ending with _ will auto append extra info to file name. Defualt value will be used if not specified.')
    parser_type1.add_argument('--STOCH', type=str, nargs='?', default=default_stoch, help='String to set stoch generation. Default vaule will be used if not specified')

    #all the params needed for type 2 (Point Source to Finite Fault)
    parser_type2.add_argument('--LAT',type=float,nargs='?',default=None)
    parser_type2.add_argument('--LON',type=float,nargs='?',default=None)
    parser_type2.add_argument('--MAG',type=float,nargs='?',default=None)
    parser_type2.add_argument('--STK',type=float,nargs='?',default=None)
    parser_type2.add_argument('--RAK',type=float,nargs='?',default=None)
    parser_type2.add_argument('--DIP',type=float,nargs='?',default=None)
    parser_type2.add_argument('--DT',type=float,nargs='?',default=None)

    parser_type2.add_argument('--DEPTH',type=float,nargs='?',default=None)
    

    #can be read from template
    parser_type2.add_argument('--SEED',type=float,nargs='?',default=None)
    parser_type2.add_argument('--RVFRAC',type=float,nargs='?',default=None)
    parser_type2.add_argument('--ROUGH',type=float,nargs='?',default=None)
    parser_type2.add_argument('--SLIP_COV',type=float,nargs='?',default=None) 
    parser_type2.add_argument('--MWSR',type=str,nargs='?',default=None,help='Mw Scaling Relation (string), one of:HanksBakun2002,BerrymanEtAl2002,VillamorEtAl2001')
    parser_type2.add_argument('--STOCH',type=str,nargs='?',default=None,help='String to set stoch gener    ation. Default vaule will be used if not specified')
    parser_type2.add_argument('--GENSLIP',type=str,nargs='?',default=None,help='Version of genslip')

    parser_type2.add_argument('--PREFIX', type=str, nargs='?', default=None ,help='The prefix you want     for SRF files. defualt= source. Ending with _ will auto append extra info to file name. Defualt value will be used if not specified.')
    
    
    #all params needed for type 3 (ffd to ff)
    parser_type3.add_argument('--LAT',type=float,nargs='?',default=None)
    parser_type3.add_argument('--LON',type=float,nargs='?',default=None)
    parser_type3.add_argument('--MAG',type=float,nargs='?',default=None)
    parser_type3.add_argument('--STK',type=float,nargs='?',default=None)
    parser_type3.add_argument('--RAK',type=float,nargs='?',default=None)
    parser_type3.add_argument('--DIP',type=float,nargs='?',default=None)
    parser_type3.add_argument('--DT',type=float,nargs='?',default=None)
    #can be read from template
    parser_type3.add_argument('--SEED',type=float,nargs='?',default=None)
    parser_type3.add_argument('--RVFRAC',type=float,nargs='?',default=None)
    parser_type3.add_argument('--ROUGH',type=float,nargs='?',default=None)
    parser_type3.add_argument('--SLIP_COV',type=float,nargs='?',default=None)
    parser_type3.add_argument('--FLEN',type=float,nargs='?',default=None)
    parser_type3.add_argument('--DLEN',type=float,nargs='?',default=None)
    parser_type3.add_argument('--FWID',type=float,nargs='?',default=None)
    parser_type3.add_argument('--DWID',type=float,nargs='?',default=None)
    parser_type3.add_argument('--DTOP',type=float,nargs='?',default=None)
    parser_type3.add_argument('--SHYPO',type=float,nargs='?',default=None)
    parser_type3.add_argument('--DHYPO',type=float,nargs='?',default=None)
    parser_type3.add_argument('--STOCH',type=str,nargs='?',default=None)
    parser_type3.add_argument('--PREFIX',type=str,nargs='?',default=None)
    parser_type3.add_argument('--GENSLIP',type=str,nargs='?',default=None)
    #parser_type3.add_argument('',type=float,nargs='?',default=None)
     
    #all params needed for type 4 (multi-seg) 
    #parser_type4.add_argument('--M_NSEG',type=float,nargs='?',default=None)
    #parser_type4.add_argument('--M_SEG_DELAY',type=float,nargs='?',default=None)
    #parser_type4.add_argument('--M_MAG',type=float,nargs='?',default=None)
    #parser_type4.add_argument('--M_MOM',type=float,nargs='?',default=None)
    #parser_type4.add_argument('--M_RVFAC_SEG',type=float,nargs='?',default=None)
    #parser_type4.add_argument('--M_GWID',type=float,nargs='?',default=None)
    #parser_type4.add_argument('--M_RUP_DELAY',type=float,nargs='?',default=None)
    #parser_type4.add_argument('--M_FLEN',type=float,nargs='?',default=None)
    #parser_type4.add_argument('--M_DLEN',type=float,nargs='?',default=None)
    #parser_type4.add_argument('--M_FWID',type=float,nargs='?',default=None)
    #parser_type4.add_argument('--M_DWID',type=float,nargs='?',default=None)
    #parser_type4.add_argument('--M_DTOP',type=float,nargs='?',default=None)
    #parser_type4.add_argument('--M_STK',type=float,nargs='?',default=None)
    #parser_type4.add_argument('--M_RAK',type=float,nargs='?',default=None)
    #parser_type4.add_argument('--M_DIP',type=float,nargs='?',default=None)
    #parser_type4.add_argument('--M_ELON',type=float,nargs='?',default=None)
    #parser_type4.add_argument('--M_ELAT',type=float,nargs='?',default=None)
    #parser_type4.add_argument('--M_SHYPO',type=float,nargs='?',default=None)
    #parser_type4.add_argument('--M_DHYPO',type=float,nargs='?',default=None)
    parser_type4.add_argument('--DT',type=float,nargs='?',default=None)
    parser_type4.add_argument('--SEED',type=int,nargs='+',default=None)
    parser_type4.add_argument('--SEED_INC',type=int,default=None)
    parser_type4.add_argument('--RVFRAC',type=float,nargs='?',default=None)
    parser_type4.add_argument('--ROUGH',type=float,nargs='?',default=None)
    parser_type4.add_argument('--SLIP_COV',type=float,nargs='?',default=None)
    parser_type4.add_argument('--PREFIX',type=str,nargs='?',default=None)
    parser_type4.add_argument('--N_SCENARIOS',type=int,nargs='?',default=None)
    parser_type4.add_argument('--N_SEED_INC',type=int,nargs='?',default=None)
    #parser_type4.add_argument('--CASES',type=str,nargs='+',default=None)
    parser_type4.add_argument('--GENSLIP',type=str,nargs='?',default=None) 
    parser_type4.add_argument('--STOCH',type=str,nargs='?',default=None)
    parser_type4.add_argument('--MW_TOTAL',type=float,nargs='?',default=None) 
    #testdir = os.path.join(mydir,'setSrfParams.py.template')
    #print  get_default("PREFIX",testdir
     
    #import runpy
    #test=runpy.run_path(testdir)
    #teststr = ("t={PREFIX}").format(**test)
    #print teststr 
    #import imp

    #testimp = imp.load_source("a.py.tmp",mydir)
    #print "test:",testimp
    
    #test2 = imp.find_module("a", [os.path.join(mydir,'')] )
    #print test2
   
    #import importlib
    #testdir = os.path.join(mydir,'a.py.tmp')
    #test3 = importlib.import_module(testdir) 
    
    
    
    #parser_type1.add_argument('', type=float)
    #parser_type1.add_argument('', type=float)
    #parser_type1.add_argument('', type=float)

    #group2 = parser.add_argument_group('type 2', 'Additional arguments for type 2 srf')
    #group3 = parser.add_argument_group('type 3', 'Additional arguments for type 3 srf')
    #group4 = parser.add_argument_group('type 4', 'Additional arguments for type 4 srf')
    
    #arguments needed for type 1
    #group2.add_argument("lat",type=float,help="latitude")
    #group2.add_argument("--stoch", help="Set 1 to use stochastic generation")
    #parser.add_argument("lat",type=float,help= "latitue")
    #parser.add_argument("lon",type=float)
    #parser.add_argument("depth",type=float)
    #parser.add_argument("mag",type=float)
    #parser.add_argument("stk",type=float)
    #parser.add_argument("dip",type=float)
    #parser.add_argument("rak",type=float)
    #parser.add_argument("dt",type=float)

    args=parser.parse_args()
    #print args 
    params_all = vars(args)
    
    #if check_params(params_all):
    #    print "Input error, please use -h/--help for usage"
    #    sys.exit(0)
   
    confirm_params(params_all) 
    #print get_order(template_dir)
        #print i,': ',params_all[i]
    #for i in :
    #    print i,': ',a[i]
    f = open(os.path.join(mydir,"setSrfParams.py"),'w')
    gen_params(params_all,f) 
    #with open(os.path.join(mydir,"setSrfParams.py.template"),'r') as fr:
    #   lines = fr.readlines()
    #   str_lines=''.join(lines)
    #   str_lines2= string.Template(str_lines).substitute({'type':args.type,'lat':args.lat,'lon':args.lon,'depth':args.depth,'mag':args.mag, 'stk':args.stk, 'dip':args.dip,'rak':args.rak, 'dt':args.dt})
    #order = get_order(template_dir)
    #with open(os.path.join(os.path.curdir,'setSrfParams.py'),'w') as fw:
    #    for i in order:
    #        if i in params_all:
    #            fw.write("%s = %s"%
    #      fw.write(str_lines2)

if __name__=='__main__':
    main()
