#!/usr/bin/python2

import argparse
import os
import string
from inspect import getsourcefile
import sys
mydir=os.path.dirname(os.path.abspath(getsourcefile(lambda:0)))
template_dir = os.path.join(mydir,'setSrfParams.py.template')

default_params= {}
default_params['PREFIX'] = 'Srf/source'
default_params['STOCH']  = 'Stoch'
default_params['SEED_INC'] = 1
msg_specify = "#!!!PLEASE SPECIFY!!!#"
msg_not_assigned = "#This param was not specified, please edit setSrfParams.py manually."



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
    f.write("#Please change the values accordingly if '%s'\n\n"%msg_specify)
    f.write("# TYPE:\n# 1: point source to point source srf\n")
    f.write("# 2: point source to finite fault srf\n")
    f.write("# 3: finite fault to finite fault srf\n")
    f.write("# 4: multi-segment finite fault srf\n")
    f.write("TYPE = %s\n"%params_all['TYPE'])
    f.write("# specify basename for gsf, srf and stoch file created\n")
    f.write("# PREFIX for gsf/srf/stoch files\n")
    f.write("# if prefix ends with '_', automatic naming follows\n") 
    f.write("PREFIX = '%s'\n"%params_all['PREFIX'])
    f.write("# directory for Stoch file(s)\n")
    f.write("# set to None to not produce the stoch file\n")
    f.write("STOCH = '%s'\n"%params_all['STOCH'])
    f.write("# rupture timestep (float), e.g. DT = 0.025\n")
    f.write("DT = %s\n"%params_all['DT'])
    
    ##params for 1,2,3
    if src_type in [1,2,3]:
        f.write("# latitude (float), e.g. LAT = -43.5095\n")
        f.write("LAT = %s\n"%params_all['LAT'])
        f.write("# longitude (float), e.g. LON = 172.1811\n")
        f.write("LON = %s\n"%params_all['LON']) 
        f.write("# depth (float), ex. DEPTH = 2.0\n")
        f.write("DEPTH = %s\n"%params_all['DEPTH'])
        f.write("# magnitude (float), e.g. MAG = 5.5\n")
        f.write("MAG = %s\n"%params_all['MAG'])
        f.write("# strike (int), e.g. STK = 150\n")
        f.write("STK = %s\n"%params_all['STK'])
        f.write("# dip (int), e.g DIP = 54\n")
        f.write("DIP = %s\n"%params_all['DIP'])
        f.write("# rake (int), e.g RAK = 6\n")
        f.write("RAK = %s\n"%params_all['RAK'])
    #type 1 specific 
    if src_type == 1:
        f.write("# specify seismic moment directly (e.g. -1 to use magnitude)\n")
        f.write("MOM = %s\n"%params_all['MOM'])
    
    #params for 2,3,4
    if src_type in [2,3,4]:
        f.write("# GENSLIP VERSION:\n# '3.3', '5.2.3a'\n# e.g. GENSLIP = '5.2.3a'\n")
        if params_all['GENSLIP'] == msg_specify:
            f.write("GENSLIP = ''%s\n"%msg_specify)
        else:
            f.write("GENSLIP = '%s'\n"%params_all['GENSLIP'])
        f.write("# roughness of fault only for genslip 5.0+\n# 0.1 is a good valu   e to use - Rob Graves\n# 0.0050119 = (10^(-2.3)) - Shi & Day 2014. Used as default in 5.2.3a\n")
        f.write("ROUGH = %s\n"%params_all['ROUGH'])
        f.write("# (float), e.g. RVFRAC = 0.8\n")
        f.write("RVFRAC = %s\n"%params_all['RVFRAC'])
        f.write("# (float), e.g. SLIP_COV = 0.85\n")
        f.write("SLIP_COV = %s\n"%params_all['SLIP_COV'])
        f.write("# seed for stoch, e.g. SEED = 103245\n")
        f.write("SEED = %s\n"%params_all["SEED"])
    if src_type == 4:
        f.write("# only used in type 4 to go over multiple seeds\n")
        f.write("# e.g. SEED_INC = 1\n")
        f.write("SEED_INC =%s\n"%params_all['SEED_INC'])
    if src_type == 2:
        f.write("# Mw Scaling Relation (string), one of:\n")
        f.write("# HanksBakun2002, BerrymanEtAl2002, VillamorEtAl2001\n")
        if params_all['MWSR'] == msg_specify:
            f.write("MWSR = ''%s\n"%params_all['MWSR'])
        else:
            f.write("MWSR = '%s'\n"%params_all['MWSR'])
    if src_type == 3:
        f.write("# (float), ex. FLEN = 45.7. Should be the same as FWID\n")
        f.write("FLEN = %s\n"%params_all['FLEN'])
        f.write("# (float), ex. DLEN = 0.10. Should be the same as DWID\n")
        f.write("DLEN = %s\n"%params_all['DLEN'])
        f.write("# (float), ex. 45.7\n")
        f.write("FWID = %s\n"%params_all['FWID'])
        f.write("# (float), \n")
        f.write("DWID = %s\n"%params_all['DWID'])
        f.write("# (float), ex. DTOP = 0.0\n")
        f.write("DTOP = %s\n"%params_all['DTOP'])
        f.write("# km east\n# e.g. SHYPO = 0.0\n")
        f.write("SHYPO = %s\n"%params_all['SHYPO'])
        f.write("# km down\n# e.g. DHYPO = 22.8544094807\n")
        f.write("DHYPO = %s\n"%params_all['DHYPO'])
    if src_type == 4:
        f.write("# how many scenarios to run. different scenarios can have randomised parameters (randomisation amount set below)\n")
        f.write("# e.g. N_SCENARIOS = 5\n")
        f.write("N_SCENARIOS = %s\n"%params_all['N_SCENARIOS'])
        f.write("# how many scenarios to run before incrementing the seed value.\n")
        f.write("# this will randomise parameters but keep the seed the same.\n")
        f.write("# e.g. N_SEED_INC = 20\n")
        f.write("N_SEED_INC = %s\n"%params_all['N_SEED_INC'])

        #array type params
        f.write("\n\n# the array/list size of each params below must be identical\n")
        f.write("# description of main segments\n")
        f.write("# len(CASES) == len(M_MAG) && len(M_MAG) == len(M_NSEG) && ...\n")
        f.write("# e.g. CASES = ['NNWEpic', 'CharX', 'GrNW', 'GrEW', 'GrEast', 'Horotar', 'NEStepOv']\n")
        f.write("CASES = %s\n"%params_all['CASES'])
        f.write("# master segments can be made up of multiple sub-segments\n")
        f.write("# example: M_MAG = [5.9, 6.5, 6.4, 6.7, 6.1, 5.9, 3.35]\n")
        f.write("M_MAG = %s\n"%params_all['M_MAG'])
        f.write("# as above, if MOM invalid ( <= 0), calculated from MAG\n")
        f.write("# eg. M_MOM = [-1, -1, -1, -1, -1, -1, -1]\n")
        f.write("M_MOM = %s\n"%params_all['M_MOM'])
        f.write("# how many segments each master is made of\n")
        f.write("# eg. M_NSEG = [1, 1, 1, 1, 1, 1, 1]\n")
        f.write("M_NSEG = %s\n"%params_all['M_NSEG'])
        f.write("# e.g M_SEG_DELAY = [\"0\", \"1\", \"1\", \"1\", \"1\", \"1\", \"1\"]\n")
        f.write("M_SEG_DELAY = %s\n"%params_all['M_SEG_DELAY'])
        f.write("# e.g M_RVFAC_SEG = [\"-1\", \"0.5\", \"0.5\", \"0.5\", \"0.5\", \"0.5\", \"0.5\"]\n")
        f.write("M_RVFAC_SEG = %s\n"%params_all['M_RVFAC_SEG'])
        f.write("# e.g M_GWID = [\"-1\", \"8.0\", \"8.0\", \"8.0\", \"8.0\", \"8.0\", \"8.0\"]\n")
        f.write("M_GWID = %s\n"%params_all['M_GWID'])
        f.write("#e.g M_RUP_DELAY = [0, 5, 11.0, 7.0, 18.0, 18.0, 17.0]\n")
        f.write("M_RUP_DELAY = %s\n"%params_all['M_RUP_DELAY'])
        
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



        f.write("M_FLEN = %s\n"%params_all['M_FLEN'])
        f.write("M_DLEN = %s\n"%params_all['M_DLEN'])
        f.write("M_FWID = %s\n"%params_all['M_FWID'])
        f.write("M_DWID = %s\n"%params_all['M_DWID'])
        f.write("M_DTOP = %s\n"%params_all['M_DTOP'])
        f.write("M_STK  = %s\n"%params_all['M_STK'])
        f.write("M_RAK  = %s\n"%params_all['M_RAK'])
        f.write("M_DIP  = %s\n"%params_all['M_DIP'])
        f.write("M_ELON = %s\n"%params_all['M_ELON'])
        f.write("M_ELAT = %s\n"%params_all['M_ELAT'])
        f.write("M_SHYPO= %s\n"%params_all['M_SHYPO'])
        f.write("M_DHYPO= %s\n"%params_all['M_DHYPO'])
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
        f.write("V_MAG = %s\n"%params_all['V_MAG'])
        f.write("# fault width/length variability (proportion, eg. 0.1 = 10%)\n")
        f.write("# e.g V_FWID = [0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00]\n")
        f.write("# e.g V_FLEN = [0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00]\n")
        f.write("V_FWID = %s\n"%params_all['V_FWID'])
        f.write("V_FLEN = %s\n"%params_all['V_FLEN'])
        f.write("# rupture time delay variability, requires dependency specification below\n")
        f.write("# absolute value\n")
        f.write("# e.g.V_RDELAY = [0, 1, 3, 0.5, 0, 0, 0]\n")
        f.write("V_RDELAY = %s\n"%params_all['V_RDELAY'])
        f.write("# rupture dependency\n")
        f.write("# None means it is a starting point\n")
        f.write("# otherwise it shoould be the segment that triggered this segment\n")
        f.write("# e.g. [None, 0,1,1,3] means the seg.0 is the starting point, seg.1 is triggered by seg.0. seg.2 seg.3 is triggered by seg.1, seg.4 is triggered by seg.3 \n")
        f.write("D_RDELAY = %s\n"%params_all['D_RDELAY'])
        f.write("# hypocentre variabiity - absolute for first segment, [SHYPO, DHYPO]\n")
        f.write("V_HYPO = %s\n"%params_all['V_HYPO'])
        
    
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
        return ['TYPE','LAT', 'LON', 'DEPTH', 'MAG', 'MOM', 'STK', 'RAK', 'DIP', 'DT','PREFIX','STOCH']
    if src_type == 2:
        return ['TYPE','LAT', 'LON', 'MAG', 'STK', 'RAK', 'DIP', 'DT', 'PREFIX', 'SEED', 'RVFRAC', 'ROUGH', 'SLIP_COV', 'DEPTH', 'MWSR', 'STOCH', 'GENSLIP']
    if src_type == 3:
        return ['TYPE','LAT', 'LON', 'DEPTH','MAG', 'STK', 'RAK', 'DIP', 'DT', 'PREFIX','SEED', 'RVFRAC', 'ROUGH', 'SLIP_COV', 'FLEN', 'DLEN', 'FWID', 'DWID', 'DTOP','SHYPO', 'DHYPO', 'STOCH', 'GENSLIP','MW_TOTAL', 'V_MAG', 'V_FWID', 'V_FLEN', 'V_RDELAY', 'D_RDELAY', 'V_HYPO']
    if src_type == 4:
        return ['TYPE','M_NSEG', 'M_SEG_DELAY', 'M_MAG', 'M_MOM','M_RVFAC_SEG', 'M_GWID', 'M_RUP_DELAY', 'M_FLEN', 'M_DLEN', 'M_FWID', 'M_DWID', 'M_DTOP', 'M_STK', 'M_RAK', 'M_DIP', 'M_ELON', 'M_ELAT', 'M_SHYPO', 'M_DHYPO', 'DT', 'SEED', 'SEED_INC', 'RVFRAC', 'ROUGH', 'SLIP_COV', 'PREFIX', 'N_SCENARIOS', 'N_SEED_INC', 'CASES', 'GENSLIP', 'STOCH','MW_TOTAL', 'V_MAG', 'V_FWID', 'V_FLEN', 'V_RDELAY', 'D_RDELAY', 'V_HYPO']

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

def print_not_assigned(list_missing_params,params_all):
    if len(list_missing_params) == 0:
        return
    print "#"*20
    print "the params below are not specified, please manually edit setSrfParams.py to specify the value."
    print "#"*20
    for i in list_missing_params:
        #print i,': ',params_all[i]
        print i,'='


        
def confirm_params(params_all):
    #print params_all['TYPE']
    src_type = int(params_all['TYPE'])
    order = get_order(src_type)
    list_missing_params = [] 
    #for i in order :
    #    try:
    #        params_all[i] 
    #    except:
    #        print i," not assigned"
    '''
    if src_type == 4:
        for i in order:
            if i in params_all:
                if params_all[i] != None:
                    print i,': ',params_all[i]
            else:
                params_all[i] = None
                print i,': ',params_all[i],' ',msg_not_assigned
        
        print "Only params above are set.\nPlease manually modify the rest of params by editing setSrfParams.py.\n"
        return 
    '''
    #print params_all['TYPE']
    for i in order:
        if i in params_all:
            if params_all[i] != None:
                print i,'= ',params_all[i]
            else:
                if i in default_params:
                    params_all[i] = default_params[i]
                    print i,'= ',params_all[i]
                else:
                    #print "cannot find",i
                    params_all[i] = msg_specify
                    list_missing_params.append(i)
                    
        else:
            params_all[i] = msg_specify
            #print i,': ',params_all[i]
            list_missing_params.append(i)
    #print list of missing params
    print_not_assigned(list_missing_params,params_all)
        

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
    parser_type1.add_argument('--LAT',type=float,nargs='?',default=None,help='latitude (float)')
    parser_type1.add_argument('--LON',type=float,nargs='?',default=None,help='longitude (float)')
    parser_type1.add_argument('--DEPTH',type=float,nargs='?',default=None,help='depth (float)')
    parser_type1.add_argument('--MAG',type=float,nargs='?',default=None,help='magnitude (float)')
    parser_type1.add_argument('--MOM',type=float, nargs='?',default=None,help='seismic moment,-1 to use magnitude')
    parser_type1.add_argument('--STK', type=int,nargs='?',default=None,help='strike (int)')
    parser_type1.add_argument('--RAK', type=int,nargs='?',default=None,help='rake (int)')
    parser_type1.add_argument('--DIP', type=int,nargs='?',default=None,help='dip (int)')
    parser_type1.add_argument('--DT', type=float, nargs='?',default=None,help='rupture timestep')
    #can be read from template
    parser_type1.add_argument('--PREFIX', type=str,nargs='?', default=default_params['PREFIX'], help='The path/prefix you want for SRF files. Ending with _ will auto append extra info to file name. Defualt value will be used if not specified.')
    parser_type1.add_argument('--STOCH', type=str, nargs='?', default=default_params['STOCH'], help='String to set stoch generation. Default vaule will be used if not specified')

    #all the params needed for type 2 (Point Source to Finite Fault)
    parser_type2.add_argument('--LAT',type=float,nargs='?',default=None,help='latitude (float)')
    parser_type2.add_argument('--LON',type=float,nargs='?',default=None,help='longitude (float)')
    parser_type2.add_argument('--MAG',type=float,nargs='?',default=None,help='magnitude (float)')
    parser_type2.add_argument('--STK',type=float,nargs='?',default=None,help='strike (int)')
    parser_type2.add_argument('--RAK',type=float,nargs='?',default=None,help='rake (int)')
    parser_type2.add_argument('--DIP',type=float,nargs='?',default=None,help='dip (int)')
    parser_type2.add_argument('--DT',type=float,nargs='?',default=None,help='rupture timestep')
    parser_type2.add_argument('--DEPTH',type=float,nargs='?',default=None)
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
    parser_type3.add_argument('--DEPTH',type=float,nargs='?',default=None)
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
    parser_type3.add_argument('--MW_TOTAL',type=float,nargs='?',default=None)
    #parser_type3.add_argument('--V_MAG',type=str,nargs='?',default=None)
    #parser_type3.add_argument('--V_FLEN',type=str,nargs='?',default=None)   
    #parser_type3.add_argument('--V_RDELAY',type=str,nargs='?',default=None)
    #parser_type3.add_argument('--D_RDELAY',type=str,nargs='?',default=None)
    #parser_type3.add_argument('--V_HYPO',type=str,nargs='?',default=None)
     
    #all params needed for type 4 (multi-seg) 
    #parser_type4.add_argument('--M_NSEG',type=float,nargs='?',default=None)
    #parser_type4.add_argument('--M_SEG_DELAY',type=float,nargs='?',default=None)
    #parser_type4.add_argument('--M_MAG',type=float,nargs='?',default=None)
    #parser_type4.add_argument('--M_MOM',type=float,nargs='?',default=None)
    #parser_type4.add_argument('--M_RVFAC_SEG',type=float,nargs='?',default=None)
    #parser_type4.add_argument('--M_GWID',type=float,nargs='?',default=None)
    #parser_type4.add_argument('--M_RUP_DELAY',type=float,nargs='?',default=None)
    #parser_type4.add_argument('--M_FLEN',type=float,nargs='+',default=None)
    #parser_type4.add_argument('--M_DLEN',type=float,nargs='+',default=None)
    #parser_type4.add_argument('--M_FWID',type=float,nargs='+',default=None)
    #parser_type4.add_argument('--M_DWID',type=float,nargs='+',default=None)
    #parser_type4.add_argument('--M_DTOP',type=float,nargs='+',default=None)
    #parser_type4.add_argument('--M_STK',type=float,nargs='?',default=None)
    #parser_type4.add_argument('--M_RAK',type=float,nargs='?',default=None)
    #parser_type4.add_argument('--M_DIP',type=float,nargs='?',default=None)
    #parser_type4.add_argument('--M_ELON',type=float,nargs='?',default=None)
    #parser_type4.add_argument('--M_ELAT',type=float,nargs='?',default=None)
    #parser_type4.add_argument('--M_SHYPO',type=float,nargs='?',default=None)
    #parser_type4.add_argument('--M_DHYPO',type=float,nargs='?',default=None)
    parser_type4.add_argument('--DT',type=float,nargs='?',default=None)
    parser_type4.add_argument('--SEED',type=int,nargs='?',default=None)
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
    print "#"*20,"\n","setSrfParams.py is located at: %s\n"%os.path.join(os.getcwd(),"setSrfParams.py"),"#"*20

    f = open(os.path.join(os.getcwd(),"setSrfParams.py"),'w')
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
