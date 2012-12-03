import operator
import bisect

def sort_table(table, col=0):
    return sorted(table, key=operator.itemgetter(col))


# check if two floats are close to each other
def float_eq(a,b):
    return abs(a-b)<1e-10*max(1e-5,abs(a),abs(b))


def linear_interpol(x,x1,x2,y1,y2):
    assert(x1<=x)
    assert(x2>=x)
    assert(x1<=x2)

    alpha = (x - x1) / (x2-x1)
    return (1.-alpha)*y1 + alpha*y2


def read_table(filename):
    table=[]
    for line in open(filename).readlines():
        if (line[0]!='#'):
            numbers = map(float, line.split())
            table.append(numbers)
    return table

def cut_table(table,min, max):
    tablen=[]
    for i in range(min,max,1):
        tablen.append(table[i,:])
    return tablen

def lookup_and_interpolate(table_x, table_y, x_value):	
    idx = bisect.bisect_left(table_x, x_value) - 1
    if (idx < 0):
        return table_y[0]
    elif (idx < len(table_x)-1):
        return linear_interpol(x_value, table_x[idx], table_x[idx+1], \
					     table_y[idx], table_y[idx+1])
    else:
        return table_y[idx]
