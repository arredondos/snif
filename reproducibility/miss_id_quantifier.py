from math import log10
import numpy as np
import random

import snif2tex

basenames = [
    # './.results/validation_main/XQ/NN/sc01NN_C_s400_c01NN_dA_w100_XQ',
    # './.results/validation_main/XQ/NN/sc02NN_C_s400_c02NN_dA_w100_XQ',
    # './.results/validation_main/XQ/NN/sc03NN_C_s400_c03NN_dA_w100_XQ',
    # './.results/validation_main/XQ/NN/sc04NN_C_s487_c04NN_dA_w100_XQ',
    './.results/validation_main/XQ/NN/sc05NN_C_s445_c05NN_dA_w100_XQ',
    # './.results/validation_main/XQ/NN/sc06NN_C_s445_c06NN_dA_w100_XQ',

    # './.results/validation_main/XQ/NY/sc01NY_C_s400_c01NY_dA_w100_XQ',
    # './.results/validation_main/XQ/NY/sc02NY_C_s400_c02NY_dA_w100_XQ',
    # './.results/validation_main/XQ/NY/sc03NY_C_s391_c03NY_dA_w100_XQ',
    # './.results/validation_main/XQ/NY/sc04NY_C_s430_c04NY_dA_w100_XQ',
    # './.results/validation_main/XQ/NY/sc05NY_C_s438_c05NY_dA_w100_XQ',
    # './.results/validation_main/XQ/NY/sc06NY_C_s416_c06NY_dA_w100_XQ',

    # './.results/validation_main/XQ/YN/sc01YN_C_s400_c01YN_dA_w100_XQ',
    # './.results/validation_main/XQ/YN/sc02YN_C_s400_c02YN_dA_w100_XQ',
    # './.results/validation_main/XQ/YN/sc03YN_C_s391_c03YN_dA_w100_XQ',
    # './.results/validation_main/XQ/YN/sc04YN_C_s445_c04YN_dA_w100_XQ',
    # './.results/validation_main/XQ/YN/sc05YN_C_s439_c05YN_dA_w100_XQ',
    # './.results/validation_main/XQ/YN/sc06YN_C_s435_c06YN_dA_w100_XQ',

    # './.results/validation_main/XQ/YY/sc02YY_C_s400_c02YY_dA_w100_XQ',
    # './.results/validation_main/XQ/YY/sc03YY_C_s391_c03YY_dA_w100_XQ',
    # './.results/validation_main/XQ/YY/sc04YY_C_s405_c04YY_dA_w100_XQ',
    # './.results/validation_main/XQ/YY/sc05YY_C_s411_c05YY_dA_w100_XQ',
    # './.results/validation_main/XQ/YY/sc06YY_C_s405_c06YY_dA_w100_XQ',
]

def any():
    return random.choice([0, 1, 3, 4])

def far():
    return random.choice([0, 4])

def random_indexes(tests, samples):
    return random.choices(tests, k=samples)    

def bad_indexes(tests, samples, errors, threshold):
    filtered = [test for test in tests if errors[test] > threshold]
    return random.choices(filtered, k=samples)

for basename in basenames:
    data = snif2tex._parse_snif_output(basename)
    snif = snif2tex._read_json(basename + '_settings.json')

    widths_2 = [log10(test["sim. t3"] / test["sim. t2"]) for test in data]
    err_2 = [100 * abs(test["inf. M2"] - test["sim. M2"]) / test["sim. M2"] for test in data]

    sim = []
    inf = []
    for i in range(5):
        sim.append([test["sim. M" + str(i)] for test in data])
        inf.append([test["inf. M" + str(i)] for test in data])

    all_tests = list(range(len(sim[0])))
    samples = 200
    threshold = 10

    p33 = np.percentile(widths_2, 33)
    p66 = np.percentile(widths_2, 67)
    partitioned_tests = [
        [test for test in all_tests if widths_2[test] < p33],
        [test for test in all_tests if p33 <= widths_2[test] < p66],
        [test for test in all_tests if p66 <= widths_2[test]]
    ]

    repetitions = 100

    from statistics import mean, stdev
    print("stat|mean|std_dev")
    for filter in [True, False]:

        global_r = []
        far_r = []
        far_best_r = []
        near_best_r = []
        for rep in range(repetitions):
            for i in range(3):
                tests = partitioned_tests[i]
                # print(f"\nBatch {i+1} with {len(tests)} tests:")

                # all
                if filter:
                    indexes = bad_indexes(tests, samples, err_2, threshold)
                else:
                    indexes = random.choices(tests, k=samples)
                x = [inf[2][i] for i in indexes]
                y = [sim[any()][i] for i in indexes]
                r2 = snif2tex._r_squared(x, y)
                global_r.append(r2)
                # print(f"global r2:\t{r2}")

                # any far
                if filter:
                    indexes = bad_indexes(tests, samples, err_2, threshold)
                else:
                    indexes = random.choices(tests, k=samples)
                x = [inf[2][i] for i in indexes]
                y = [sim[far()][i] for i in indexes]
                r2 = snif2tex._r_squared(x, y)
                far_r.append(r2)
                # print(f"global far r2:\t{r2}")

                # best far
                if filter:
                    indexes = bad_indexes(tests, samples, err_2, threshold)
                else:
                    indexes = random.choices(tests, k=samples)
                x = [inf[2][i] for i in indexes]
                y = [sim[0][i] if abs(sim[0][i]-inf[2][i]) < abs(sim[4][i]-inf[2][i]) else sim[4][i] for i in indexes]
                r2 = snif2tex._r_squared(x, y)
                far_best_r.append(r2)
                # print(f"best far r2:\t{r2}")

                # best far
                if filter:
                    indexes = bad_indexes(tests, samples, err_2, threshold)
                else:
                    indexes = random.choices(tests, k=samples)
                x = [inf[2][i] for i in indexes]
                y = [sim[1][i] if abs(sim[1][i]-inf[2][i]) < abs(sim[3][i]-inf[2][i]) else sim[3][i] for i in indexes]
                r2 = snif2tex._r_squared(x, y)
                near_best_r.append(r2)
                # print(f"best near r2:\t{r2}")
            
        print(f"random_global_{filter}|{mean(global_r)}|{stdev(global_r)}")
        print(f"random_far_{filter}|{mean(far_r)}|{stdev(far_r)}")
        print(f"best_far_{filter}|{mean(far_best_r)}|{stdev(far_best_r)}")
        print(f"best_near_{filter}|{mean(near_best_r)}|{stdev(near_best_r)}")
