from math import log10
from glob import glob
from os import path
import csv

def _read_csv(filename):
    return csv.reader(open(filename, encoding='utf-8-sig'), delimiter = ',')

def _write_csv(table, filename, delimiter='\t'):
    csv.writer(open(filename, 'w', newline=''), delimiter=delimiter).writerows(table)

def _parse_snif_output(filename):
    snif_data = []

    table = _read_csv(filename)
    first_line = next(table)
    if first_line[0] == 'SEP=':
        headers = next(table)
    else:
        headers = first_line

    cols = len(headers)
    c = int((cols - 10) / 3)
    
    for line in table:
        current_problem = {}
        for j in range(cols - 1):
            current_problem[headers[j]] = float(line[j])        
        current_problem[headers[cols - 1]] = line[cols - 1]
        
        snif_data.append(current_problem)
    
    return (c, snif_data)

def _get_dates(directory):
    filenames = glob(path.join(directory, '*.csv'))
    dates = []

    for fn in filenames:
        (c, data) = _parse_snif_output(fn)        
        for problem in data:
            for k in range(1, c):
                dates.append(problem[f'inf. t{k}'])
    
    return dates

def _distance(dates, events, g):
    log_dates = [log10(g * t) for t in dates]
    log_events = [log10(t) for t in events]
    d = 0

    for t in log_dates:
        d += min(abs(t - event) for event in log_events)
    
    return d

def process(directory, events_fn):
    from numpy import logspace
    from pathlib import Path

    dates = _get_dates(directory)
    min_gt = 0.1
    max_gt = 50

    events = next(_read_csv(events_fn))
    events = [float(t) for t in events]

    generation_time = list(logspace(log10(min_gt), log10(max_gt), 500))

    parent_dir = Path(directory).parent
    
    distance = [_distance(dates, events, g) for g in generation_time]

    species_name = path.basename(path.dirname(directory)) 
    out_filename = path.join(parent_dir, species_name + '.csv')
    out_table = [['gen. time'] + generation_time, [species_name] + distance]
    out_table = list(map(list, zip(*out_table)))
    
    _write_csv(out_table, out_filename, ',')

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('directory', help = 'Directory containing the SNIF csv files to consider')
    parser.add_argument('events', help = 'File containing the times in years of the geological events of interest')
    args = parser.parse_args()
    process(args.directory, args.events)