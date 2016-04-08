

# reads a parameter from the parameters file
def par_value (variable):
    result = ''
    par_handle = open('e3d.par', 'r')
    for line in par_handle:
        if line.startswith(variable + '='):
            # keep going and get the last result
            result = line
    par_handle.close()
    return ''.join(result.split('=')[1:]).rstrip('\n')

