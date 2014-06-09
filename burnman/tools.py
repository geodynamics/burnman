# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

import operator
import bisect
import os
import pkgutil
import numpy as np

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
            print "\t".join(r).replace("_","\_")
        return

    def col_width(table, colidx):
        return max([len(str(row[colidx])) for row in table])

    # create a format string with the first column left aligned, the others right
    # example:   {:<27}{:>11}{:>6}{:>8}
    frmt = "".join([ ('{:<' if i==0 else '{:>')+str(1+col_width(table,i))+'}' for i in range(len(table[0])) ])
    for r in table:
        print frmt.format(*r)

def sort_table(table, col=0):
    """
    Sort the table according to the column number
    """
    return sorted(table, key=operator.itemgetter(col))


def float_eq(a,b):
    """
    Test if two floats are almost equal to each other
    """
    return abs(a-b)<1e-10*max(1e-5,abs(a),abs(b))


def linear_interpol(x, x1, x2, y1, y2):
    """
    Linearly interpolate to point x, between
    the points (x1,y1), (x2,y2)
    """
    assert(x1<=x)
    assert(x2>=x)
    assert(x1<=x2)

    alpha = (x - x1) / (x2-x1)
    return (1.-alpha)*y1 + alpha*y2

def read_table(filename):
    datastream = pkgutil.get_data('burnman', 'data/'+filename)
    datalines = [ line.strip() for line in datastream.split('\n') if line.strip() ]
    table=[]

    for line in datalines:
        if (line[0]!='#'):
            numbers = np.fromstring( line , sep =' ')
            table.append(numbers)
    return np.array(table)

def cut_table(table, min_value, max_value):
    tablen=[]
    for i in range(min_value,max_value,1):
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

def molar_volume_from_unit_cell_volume(unit_cell_v, z):
    """
    takes unit cell volume in Angstroms^3, as is often reported,
    and the z number for the cell (number of atoms per unit cell,
    NOT number of atoms per molecular formula), and calculates
    the molar volume, as expected by the equations of state.
    """
    N_a = 6.0221415e23
    return  unit_cell_v*N_a/1e30/z

