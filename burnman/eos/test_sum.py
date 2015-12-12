excesses = {'G': 100., 'dGdP': 10.}

excesses_2 = {'G': 100., 'dGdP': 10.}

for key, value in excesses_2.iteritems():
    excesses[key] += value


print excesses
