from os import path, mkdir, rename
from glob import glob

def clean(directory):
    outdir = path.join(directory, 'completed')
    if not path.isdir(outdir):
        mkdir(outdir)

    tag = '_settings.json'
    pattern = path.join(directory, '*' + tag)
    json_filenames = glob(pattern)

    def move(fn):
        rename(fn, path.join(outdir, path.basename(fn)))

    def check_and_move(fn):
        if path.isfile(fn):
            move(fn)

    for fn in json_filenames:
        basename = fn[:-len(tag)]

        move(fn)
        move(basename + ".csv")
        move(basename + "_curves.csv")
        check_and_move(basename + ".msout")
        check_and_move(basename + ".psmcfa")
        check_and_move(basename + ".psmc")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('directory', help = 'Directory containing the SNIF outputs to clean')
    args = parser.parse_args()
    clean(args.directory)