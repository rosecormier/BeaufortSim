f = open('input.txt', 'r')

import sys

for line in f:

    var_len = 0
    line_len = len(line)

    if (line_len > 0) and (line[0] != "#"):

        for char in line:
            if char == '=':
                break
            else:
                var_len = var_len + 1

        var = line[0:var_len].strip()

        try:
            val = float(line[var_len+1:len_len-1].strip())
        except:
            val = line[var_len+1:line_len].strip()

        print(var, val)
        #setattr(param_data, var, val)
        sys.exit()
f.close
