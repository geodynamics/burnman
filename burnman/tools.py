"""
    BurnMan- a lower mantle toolkit
    Copyright (C) 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
"""


import operator
import bisect
import os

def pretty_print_table(table,use_tabs=False):
    """
    Takes a 2d table and prints it in a nice text based format. If 
    use_tabs=True then only \t is used as a separator. This is useful for
    importing the data into other apps (Excel, ...). The default is to pad
    the columns with spaces to make them look neat. The first column is
    left aligned, while the remainder is right aligned.
    """
    if use_tabs:
        for r in table:
            print "\t".join(r)
        return        
    
    def col_width(table, colidx):
        return max([len(str(row[colidx])) for row in table])

    # create a format string with the first column left aligned, the others right  
    # example:   {:<27}{:>11}{:>6}{:>8}
    frmt = "".join([ ('{:<' if i==0 else '{:>')+str(1+col_width(table,i))+'}' for i in range(len(table[0])) ])
    for r in table:
        print frmt.format(*r)

def sort_table(table, col=0):
    return sorted(table, key=operator.itemgetter(col))


# check if two floats are close to each other
def float_eq(a,b):
    return abs(a-b)<1e-10*max(1e-5,abs(a),abs(b))


def linear_interpol(x, x1, x2, y1, y2):
    assert(x1<=x)
    assert(x2>=x)
    assert(x1<=x2)

    alpha = (x - x1) / (x2-x1)
    return (1.-alpha)*y1 + alpha*y2

def open_burnman_file(filename):
    path = os.path.split(os.path.dirname(os.path.abspath(__file__)))[0]
    fullname = os.path.join(path,filename)
    return open(fullname)

def read_table(filename):
    table=[]
    
    for line in open_burnman_file(filename).readlines():
        if (line[0]!='#'):
            numbers = map(float, line.split())
            table.append(numbers)
    return table

def cut_table(table, min, max):
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
