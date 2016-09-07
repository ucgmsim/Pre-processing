import os
import os.path
import sys
from numpy import testing
import numpy as np
import subprocess

vel_files_path='/nesi/projects/nesi00213/RunFolder/baes/AlpineFault400m_m7p90-411p0x17p3_s1129570_s2n_VMSI_v1p64_400m-h0p400_EMODv3p0p4_160823/LF/Vel'


eps=1e-4 #tolerance
sigfig=4

def execute_cmd(cmd):
    p=subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    outputList = p.stdout.readlines()
    output = ''.join(outputList)
    return output


def values_in_line(line):
    values = line.strip('\n').split(' ')
    values = [float(v) for v in values if v!='']
    return values

def fail(msg, *argv):
    fail_soft(msg,argv)
    sys.exit()


def fail_soft(msg, *argv):
    print("FAIL: %s" %msg)
    for arg in argv:
        print(arg)


def succeed(msg, *argv):
    print("PASS: %s" %msg)
    for arg in argv:
        print(arg)
    print("\n\n\n")

def compare(v1,v2,sigfig=sigfig):
    try:
        testing.assert_approx_equal(v1, v2, significant=sigfig)
    except AssertionError:
        return abs(v1-v2)<=eps

    #    return False
    else:
        return True



if __name__ == '__main__':
    exts=['ver','000','090']

    if len(sys.argv)<3:
        print("Usage:%s VIRTUAL ACTUAL" %os.path.basename(sys.argv[0]))
        print(" eg. VIRTUAL can be 123739F for 123739F.zip")
        print(" eg. ACTUAL can be DSLC for DSLC.ver DSLC.000 DSLC.090")
        sys.exit()


   # fprefix_v='0C10859'
   # fprefix_a= vel_files_path+'/DSLC'
    fprefix_v = sys.argv[1]
    fprefix_a = vel_files_path+'/'+sys.argv[2]
    if os.path.exists(fprefix_v+'.zip'):
        execute_cmd('unzip -o %s.zip'%fprefix_v)

    for ext in exts:
        fname_v=fprefix_v+'.'+ext
        fname_a=fprefix_a+'.'+ext


        f_v=open(fname_v,'r')
        f_a=open(fname_a,'r')

        lines_v=f_v.readlines()
        lines_a=f_a.readlines()

        header_v = lines_v[0].strip('\n')
        header_a = lines_a[0].strip('\n')


        print("Virtual station seismogram : %s" %os.path.abspath(fname_v)) 
        print("   Header : %s" %header_v) 
        print("   Num. lines : %d" %len(lines_v))

        print("Actual station seismogram : %s" %os.path.abspath(fname_a)) 
        print("   Header : %s" %header_a) 
        print("   Num. lines : %d" %len(lines_a))

        total_fails = 0
        if len(lines_v) != len(lines_a):
            fail("The length of two files are different")

        diffs=[]
        for i in range(1,len(lines_v)):
            line_v=lines_v[i]
            line_a=lines_a[i]
            values_v = values_in_line(line_v)
            values_a = values_in_line(line_a)

            if len(values_v)!=len(values_a):
                fail("The number of fields in two lines are different",line_v,line_a)

            for j in range(len(values_v)):
                diffs.append(abs(values_v[j]-values_a[j]))
                if compare(values_v[j],values_a[j]) != True:
                    fail_soft("The difference between two values is greater than %d significant figures and eps=%f"%(sigfig,eps), "Line number:%d Column number:%d" %(i+1,j+1), "%s vs %s" %(str(values_v[j]), str(values_a[j])), "From:", line_v, line_a)
                    total_fails += 1


        diffs_array=np.array(diffs)
        print("Total of %d entries. Diff. mean = %f  std = %f  max = %f min = %f" %(len(diffs),diffs_array.mean(), diffs_array.std(), diffs_array.max(), diffs_array.min()))

        if total_fails==0:
            succeed("Comparison of %s and %s - All good" %(os.path.basename(fname_v), os.path.basename(fname_a)))
        else:
            fail_soft("Check the stats")

            

        
