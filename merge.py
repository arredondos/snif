from os import path, mkdir, rename
from glob import glob

def _create(fn):
    if path.isfile(fn):
        raise RuntimeError("the file {} already exists in the target directory".format(fn))
    else:
        return open(fn, 'w')

def _read_csv(filename):
    import csv
    return csv.reader(open(filename), delimiter = ',')

def _get_line(row):
    return ','.join(row) + '\n'

def merge(directory, auto_name):
    tag = '_settings.json'
    pattern = path.join(directory, '*' + tag)
    filenames = glob(pattern)
    basenames = [fn[:-len(tag)] for fn in filenames]

    params_fn = path.join(directory, "merged.csv")
    params = _create(params_fn)
    
    curves_fn = path.join(directory, "merged_curves.csv")
    curves = _create(curves_fn)
    names = [
        'time-',
        'source-iicr-',
        'source-cdf-',
        'source-pdf-',
        'inferred-iicr-',
        'inferred-cdf-',
        'inferred-pdf-',
        'model-iicr-',
        'model-cdf-',
        'model-pdf-',
    ]
    k = len(names)

    p = _read_csv(basenames[0] + '.csv')
    params.write(_get_line(next(p)))
    params.write(_get_line(next(p)))
    
    i = 0
    for bn in basenames:
        p = _read_csv(bn + '.csv')
        next(p) # skip delimiter line
        next(p) # skip header        
        c = list(_read_csv(bn + '_curves.csv'))
        for row in p:
            new_row = [str(i + 1)] + row[1:]
            params.write(_get_line(new_row))

            n = 0
            for j in range(k):
                new_row = [names[j] + str(i + 1)] + c[k * n + j][1:]
                curves.write(_get_line(new_row))

            n += 1
            i += 1
    
    print("merged {} tests!".format(i))

    if auto_name:
        from shutil import copy2

        k = len('yymmdd-hhmmss_')
        basename_tags = path.basename(basenames[0]).split('_')
        og_tags = 5
        if basename_tags[1][0] == 's':
            og_tags = 8
            basename_tags[3] = 's{:03d}'.format(i)
        
        new_basename = '_'.join(basename_tags[1 : og_tags])
        rename(params_fn, path.join(directory, new_basename + ".csv"))
        rename(curves_fn, path.join(directory, new_basename + "_curves.csv"))
        copy2(basenames[0] + tag, path.join(directory, new_basename + tag))


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('directory', help = 'Directory containing the SNIF outputs to be merged')
    parser.add_argument('-a', '--autoName', help = 'Specify wether to attempt to automatically name the merged file', action = "store_true")
    args = parser.parse_args()
    merge(args.directory, args.autoName)