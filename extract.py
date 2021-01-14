import json
import csv

def extract(basename, test_numbers, get_inferred = True, get_simulated = False, ms_command = False):
    with open(basename + '_settings.json', 'r') as f:
        settings = json.load(f)

    p = settings["inference_parameters"]["number_of_components"]    
    table = list(csv.reader(open(basename + '.csv'), delimiter=','))
    ext = 'ms' if ms_command else 'json'
    
    for i in test_numbers:
        row = table[i + 1]
        if get_inferred:
            scenario = {
                "n" : int(row[7]),
                "t" : _f(row[8 : 8 + p]),
                "M" : _f(row[8 + p : 8 + 2 * p]),
                "s" : _f(row[8 + 2 * p : 8 + 3 * p])
            }
            if settings["inference_parameters"]["infer_scale"]:
                scenario.update({"Nref" : float(row[8 + 3 * p])})
            
            with open(f'{basename}_inferred-model-{i:03d}.{ext}', 'w') as f:
                if ms_command:
                    f.write(_to_ms_command(scenario))
                else:
                    json.dump({"PSNIC" : scenario}, f, indent = 4)

        if get_simulated:
            q = settings["simulation_parameters"]["number_of_components"]
            scenario = {
                "n" : int(row[9 + 3 * p]),
                "t" : _f(row[10 + 3 * p : 10 + 3 * p + q]),
                "M" : _f(row[10 + 3 * p + q : 10 + 3 * p + 2 * q]),
                "s" : _f(row[10 + 3 * p + 2 * q : 10 + 3 * p + 3 * q])
            }
            if settings["simulation_parameters"]["simulate_scale"]:
                scenario.update({"Nref" : float(row[10 + 3 * p + 3 * q])})

            with open(f'{basename}_simulated-model-{i:03d}.{ext}', 'w') as f:
                if ms_command:
                    f.write(_to_ms_command(scenario))
                else:
                    json.dump({"PSNIC" : scenario}, f, indent = 4)
        


def _f(iterator):
    return [float(x) for x in iterator]

def _to_ms_command(scenario):
    reps = int(1e5)
    n = scenario['n']
    t = scenario['t']
    m = scenario['M']
    s = scenario['s']
    c = len(t)
    if 'Nref' in scenario:
        t = [ti / (2 * scenario['Nref']) for ti in t]
    
    command = f'ms 2 {reps} -T -L -I {n} 2 '
    command += '0 ' * (n - 1) + str(m[0])
    command += f' -eN 0 {s[0]}'    
    for i in range(1, c):
        command += f' -eM {0.5 * t[i]} {m[i]}'
        command += f' -eN {0.5 * t[i]} {s[i]}'

    return command

if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser()
    parser.add_argument('filename', help = 'Base name (without the extension) of the main csv file containing the SNIF output.')
    parser.add_argument('-i', '--inferred', action = 'store_true', help = 'Specifies that the inferred scenarios should be extracted')
    parser.add_argument('-s', '--simulated', action = 'store_true', help = 'Specifies that the simulated scenarios should be extracted')
    parser.add_argument('-ms', '--ms', action = 'store_true', help = 'Outputs the scenarios to ms commands instead of JSON files')
    parser.add_argument('tests', metavar = 'N', type = int, nargs='+', help = 'The test numbers to extract')
    args = parser.parse_args()
    
    ms_command = args.ms
    get_inferred = args.inferred
    get_simulated = args.simulated
    if not(get_inferred or get_simulated):
        get_inferred = True

    extract(args.filename, args.tests, get_inferred, get_simulated, ms_command)