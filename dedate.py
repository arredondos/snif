from os import path, rename
from glob import glob

def dedate(directory):
    k = len('yymmdd-hhmmss_')
    filenames = glob(path.join(directory, '*.*'))

    for fn in filenames:
        bn = path.basename(fn)
        new_bn = bn[k:]
        new_fn = path.join(directory, new_bn)
        rename(fn, new_fn)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('directory', help = 'Directory containing the files to be de-dated')
    args = parser.parse_args()
    dedate(args.directory)