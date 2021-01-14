from glob import glob
from os import path
import csv

def collect(location, param, dedated):
    from snif import _parse_basename
    tag = '*_settings.json'
    mask = path.join(location, tag)
    filenames = glob(mask)
    csv_basenames = [fn[:1 - len(tag)] for fn in filenames]
    collection_fn = path.join(location, 'collected ' + param + '.csv')

    first_table = read_snif_csv(csv_basenames[0] + '.csv')
    rows = len(first_table) - 1
    names = [first_table[i][-1] for i in range (1, rows + 1)]
    names = get_unique_parts(names)
    description = _parse_basename(path.basename(csv_basenames[0]), not dedated)
    collection = [[desc[0] for desc in description] + names]

    for csvbn in csv_basenames:
        table = read_snif_csv(csvbn + '.csv')
        rows = len(table) - 1
        param_index = table[0].index(param)
        values = [table[i][param_index] for i in range (1, rows + 1)]
        description = _parse_basename(path.basename(csvbn), not dedated)
        collection.append([desc[1] for desc in description]  + values)
        csv.writer(open(collection_fn, "w", newline='')).writerows(collection)

def read_snif_csv(filename):
    table = csv.reader(open(filename), delimiter = ',')
    next(table) # skip the delimiter specification line
    return list(table)

def get_unique_parts(collection):

    if len(collection) == 1:
        fn = path.basename(collection[0])
        fn_no_extension = path.splitext(fn)[0]
        return [fn_no_extension]

    commonprefix = path.commonprefix(collection)
    commonpostfix = path.commonprefix([''.join(reversed(x)) for x in collection])
    p = len(commonprefix)
    q = len(commonpostfix)
    unique = [x[p:-q] for x in collection]
    return unique

if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser()
    parser.add_argument('directory', help = 'path with the SNIF output to process')
    parser.add_argument('parameter', help = 'name of the parameter to collect')
    parser.add_argument('-d', '--dedated', help = 'specify that the SNIF output has been de-dated (date-time tag removed from basenames).', action = "store_true", default = False)
    args = parser.parse_args()

    collect(args.directory, args.parameter, args.dedated)