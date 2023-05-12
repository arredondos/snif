from os import path

def convert(directory, g, mu):
    from glob import glob

    filenames = glob(path.join(directory, '*.iicr'))
    (theta, rho, s) = (0.05, 0.01, 100)
    N0 = theta / (4 * mu * s)

    for fn in  filenames:
        with open(fn, 'r') as f:
            lines = f.readlines()
        
        time = [float(years) / (2 * N0 * g) for years in lines[0].strip().split('\t')]
        iicr = [float(psmc) / N0 for psmc in lines[1].strip().split('\t')]

        if len(time) != len(iicr):
            raise Exception('invalid .IICR file')
        
        lines = ['//\n']
        lines.append(f'TR\t{theta}\t{rho}\n')
        for i, (x, y) in enumerate(zip(time, iicr)):
            lines.append(f'RS\t{i}\t{x}\t{y}\n')
        lines.append('//\n')
        
        with open(fn + '.psmc', 'w') as f:
            f.writelines(lines)

if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser()
    parser.add_argument('directory', help = 'Directory containing the .IICR files to be converted to .PSMC files.')
    parser.add_argument('g', help = 'Specify the generation time in years.')
    parser.add_argument('mu', help = 'Specify the per-site per-generation mutation rate.')
    args = parser.parse_args()
    
    convert(args.directory, float(args.g), float(args.mu))