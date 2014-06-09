'''This module provides function for working with significant
figures.
from Chris Burn's website  http://users.obs.carnegiescience.edu/~cburns/site/?p=22
'''
import types,re,string

epat = re.compile(r'^([^e]+)e(.+)$')

def round_sig(x, n):
   '''round floating point x to n significant figures'''
   if type(n) is not types.IntType:
      raise TypeError, "n must be an integer"
   try:
      x = float(x)
   except:
      raise TypeError, "x must be a floating point object"
   form = "%0." + str(n-1) + "e"
   st = form % x
   num,expo = epat.findall(st)[0]
   expo = int(expo)
   fs = string.split(num,'.')
   if len(fs) < 2:
      fs = [fs[0],""]
   if expo == 0:
      return num
   elif expo > 0:
      if len(fs[1]) < expo:
         fs[1] += "0"*(expo-len(fs[1]))
      st = fs[0]+fs[1][0:expo]
      if len(fs[1][expo:]) > 0:
         st += '.'+fs[1][expo:]
      return st
   else:
      expo = -expo
      if fs[0][0] == '-':
         fs[0] = fs[0][1:]
         sign = "-"
      else:
         sign = ""
      return sign+"0."+"0"*(expo-1)+fs[0]+fs[1]
      
def round_sig_error(x, ex, n, paren=False):
   '''Find ex rounded to n sig-figs and make the floating point x
   match the number of decimals.  If [paren], the string is
   returned as quantity(error) format'''
   stex = round_sig(ex,n)
   if stex.find('.') < 0:
      extra_zeros = len(stex) - n
      sigfigs = len(str(int(x))) - extra_zeros
      stx = round_sig(x,sigfigs)
   else:
      num_after_dec = len(string.split(stex,'.')[1])
      stx = ("%%.%df" % num_after_dec) % (x)
   if paren:
      if stex.find('.') >= 0:
         stex = stex[stex.find('.')+1:]
      return "%s(%s)" % (stx,stex)
   return stx,stex

def format_table(cols, errors, n, labels=None, headers=None, latex=False):
   '''Format a table such that the errors have n significant
   figures.  [cols] and [errors] should be a list of 1D arrays
   that correspond to data and errors in columns.  [n] is the number of
   significant figures to keep in the errors.  [labels] is an optional
   column of strings that will be in the first column.  [headers] is
   an optional list of column headers.  If [latex] is true, format
   the table so that it can be included in a LaTeX table '''
   if len(cols) != len(errors):
      raise ValueError, "Error:  cols and errors must have same length"

   ncols = len(cols)
   nrows = len(cols[0])

   if headers is not None:
      if labels is not None:
         if len(headers) == ncols:
            headers = [""] + headers
         elif len(headers) == ncols+1:
            pass
         else:
            raise ValueError, "length of headers should be %d" % (ncols+1)
      else:
         if len(headers) != ncols:
            raise ValueError, "length of headers should be %d" % (ncols)

   if labels is not None:
      if len(labels) != nrows:
         raise ValueError, "length of labels should be %d" % (nrows)

   strcols = []
   for col,error in zip(cols,errors):
      strcols.append([])
      strcols.append([])
      for i in range(nrows):
         val,err = round_sig_error(col[i], error[i], n)
         strcols[-2].append(val)
         strcols[-1].append(err)

   lengths = [max([len(item) for item in strcol]) for strcol in strcols]
   format = ""
   if labels is not None:
      format += "%%%ds " % (max(map(len, labels)))
      if latex:
         format += "& "
   for length in lengths: 
      format += "%%%ds " % (length)
      if latex:
         format += "& "
   if latex:
      format = format[:-2] + " \\\\"
   output = []
   if headers:
      if labels:
         hs = [headers[0]]
         for head in headers[1:]:
            hs.append(head)
            hs.append('+/-')
      else:
         hs = []
         for head in headers:
            hs.append(head)
            hs.append('+/-')
      output.append(format % tuple(hs))
   for i in range(nrows):
      if labels is not None:
         output.append(format % tuple([labels[i]] + [strcol[i] for strcol in strcols]))
      else:
         output.append(format % tuple([strcol[i] for strcol in strcols]))
   return output

def round_sig_error2(x, ex1, ex2, n):
   '''Find min(ex1,ex2) rounded to n sig-figs and make the floating point x
   and max(ex,ex2) match the number of decimals.'''
   minerr = min(ex1,ex2)
   minstex = round_sig(minerr,n)
   if minstex.find('.') < 0:
      extra_zeros = len(minstex) - n
      sigfigs = len(str(int(x))) - extra_zeros
      stx = round_sig(x,sigfigs)
      maxstex = round_sig(max(ex1,ex2),sigfigs)
   else:
      num_after_dec = len(string.split(minstex,'.')[1])
      stx = ("%%.%df" % num_after_dec) % (x)
      maxstex = ("%%.%df" % num_after_dec) % (max(ex1,ex2))
   if ex1 < ex2:
      return stx,minstex,maxstex
   else:
      return stx,maxstex,minstex

