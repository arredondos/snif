import os
import sys
import argparse
import subprocess
from enum import Enum
from collections import namedtuple
from math import log10, floor, ceil

from collect import get_unique_parts

class TimeScale(Enum):
	Coalescent, Generations, Years = range(3)

class OutputStyle(Enum):
    Full, Minimal, Excluded = range(3)

Configuration = namedtuple('Configuration', [
    'SNIF_basename',
    'plot_width_cm',
    'plot_height_cm',
	'IICR_plots_style',
	'PDF_plots_style',
	'CDF_plots_style',
    'islands_plot_style',
    'Nref_plot_style',
    'test_numbers',
    'one_file_per_test',
    'versus_plot_style',
	'CG_style',
    'CG_size_history',
	'CG_source',
	'CG_source_legend',
	'Nref_histograms_bins',
	'islands_histograms_bins',
	'time_histograms_bins',
	'migration_histograms_bins',
	'size_histograms_bins',
    'scaling_units',
    'generation_time'
])

default_configuration = Configuration(
    SNIF_basename = '.results/gadma/200327-034722_c01YN_dA_w100_MT_high_res',
    plot_width_cm = 8,
    plot_height_cm = 4,
	IICR_plots_style  = OutputStyle.Full,
	PDF_plots_style = OutputStyle.Excluded,
	CDF_plots_style = OutputStyle.Excluded,
    islands_plot_style = OutputStyle.Full,
    Nref_plot_style = OutputStyle.Full,
    test_numbers = "all",
    one_file_per_test = False,
    versus_plot_style = OutputStyle.Excluded,
	CG_style = OutputStyle.Full,
    CG_size_history = False,
	CG_source = '',
	CG_source_legend = '',
	Nref_histograms_bins = 100,
	islands_histograms_bins = 100,
	time_histograms_bins = 100,
	migration_histograms_bins = 100,
	size_histograms_bins = 100,
    scaling_units = TimeScale.Years,
    generation_time = 25,
)


# # # # # # # # # # # # # #
#                         #
#    PUBLIC  INTERFACE    #
#                         #
# # # # # # # # # # # # # #


def TeXify(config):
    basename = config.SNIF_basename
    if not os.path.exists(basename):
        os.mkdir(basename)
        
    data = _parse_snif_output(basename)
    snif = _read_json(basename + '_settings.json')
    tex_filenames = []

    plot_types = []
    if config.IICR_plots_style != OutputStyle.Excluded: plot_types.append("iicr")
    if config.CDF_plots_style != OutputStyle.Excluded: plot_types.append("cdf")
    if config.PDF_plots_style != OutputStyle.Excluded: plot_types.append("pdf")
    if len(plot_types) > 0:
        curves_data_fn = _generate_curves_data(data, config, snif)
        for plot_type in plot_types:
            tex_filenames.extend(_generate_tests_plots(
                data, 
                config, 
                snif, 
                plot_type, 
                curves_data_fn
            ))
    
    plot_types = []
    if config.islands_plot_style != OutputStyle.Excluded: plot_types.append("n")
    if config.Nref_plot_style != OutputStyle.Excluded: plot_types.append("N_ref")

    simulating = "simulation_parameters" in snif
    if len(plot_types) > 0 or (simulating and config.versus_plot_style != OutputStyle.Excluded):
        params2_data_fn = _generate_params2_data(data, config, snif)
        fitness_stats_fn = _generate_fitness_stats(data, config, snif)

    if len(plot_types) > 0 or (config.CG_style != OutputStyle.Excluded):
        hists_data_fn = _generate_hists_data(data, config, snif)

    if config.CG_style != OutputStyle.Excluded:
        params_data_fn = _generate_params_data(data, config, snif)
        tex_filenames.extend(_generate_cg_plot(
            data,
            config,
            snif,
            params_data_fn,
            hists_data_fn
        ))

    if len(plot_types) > 0:
        for plot_type in plot_types:
            tex_filenames.extend(_generate_1D_plot(
                data, 
                config, 
                snif, 
                plot_type,
                params2_data_fn,
                hists_data_fn
            ))

    if config.versus_plot_style != OutputStyle.Excluded:
        if simulating:
            tex_filenames.extend(_generate_versus_plots(
                data,
                config,
                snif,
                params2_data_fn,
                fitness_stats_fn
            ))
        else:
            Warning("Cannot generate versus plots for a SNIF output with no simulations!")

    if len(tex_filenames) > 0:
        _compile_tex_files(tex_filenames)
    else:
        Warning("No plots were specified. Please check the SNIF2TeX configuration.")


# # # # # # # # # # # # #
#                       #
#   PRIVATE INTERFACE   #
#                       #
# # # # # # # # # # # # #


_curve_names = [
    'time',
    'source-iicr',
    'source-cdf',
    'source-pdf',
    'inferred-iicr',
    'inferred-cdf',
    'inferred-pdf',
    'model-iicr',
    'model-cdf',
    'model-pdf',
]

def _parse_snif_output(basename):
    snif_data = []

    table = _read_csv(basename + '.csv')
    first_line = next(table)
    if first_line[0] == 'SEP=':
        headers = next(table)
    else:
        headers = first_line

    cols = len(headers)
    headers.extend(_curve_names)
    k = len(_curve_names)
    
    curves = list(_read_csv(basename + '_curves.csv'))
    
    i = 0
    for line in table:
        current_problem = {}
        for j in range(cols - 1):
            current_problem[headers[j]] = float(line[j])        
        current_problem[headers[cols - 1]] = line[cols - 1]

        for j in range(k):
            curve = curves[k * i + j]
            current_problem[headers[j - k]] = [float(x) for x in curve[1:]]
        
        i += 1
        snif_data.append(current_problem)
    
    return snif_data
def _generate_versus_plots(data, config, snif, params2_data_fn, fitness_stats_fn):
    from shutil import copy2
    
    n_limits = _get_inferred_param_bounds(snif, "bounds_islands")
    if n_limits[0] == 2:
        n_limits[0] = 0
    
    t_limits = _get_inferred_param_bounds(snif, "bounds_event_times")
    m_limits = _get_inferred_param_bounds(snif, "bounds_migrations_rates")
    s_limits = _get_inferred_param_bounds(snif, "bounds_deme_sizes")
    nref_limits = _get_inferred_param_bounds(snif, "bounds_effective_size")
    
    vc = min(snif["inference_parameters"]["number_of_components"], 
             snif["simulation_parameters"]["number_of_components"])

    subs = {
        'PARAMETERS2-DATA-FILENAME': os.path.basename(params2_data_fn),
        'FITNESS-STATS-FILENAME': os.path.basename(fitness_stats_fn),
        'PLOT-WIDTH': str(config.plot_width_cm) + 'cm',
        'N-LOWER-LIMIT': n_limits[0],
        'N-UPPER-LIMIT': n_limits[1],
        'T-LOWER-LIMIT': t_limits[0],
        'T-UPPER-LIMIT': t_limits[1],
        'M-LOWER-LIMIT': m_limits[0],
        'M-UPPER-LIMIT': m_limits[1],
        'S-LOWER-LIMIT': s_limits[0],
        'S-UPPER-LIMIT': s_limits[1],
        'NREF-LOWER-LIMIT': nref_limits[0],
        'NREF-UPPER-LIMIT': nref_limits[1],
        'SIMULATED-NREF': snif["simulation_parameters"]["effective_size"][0],
        'NUMBER-OF-TESTS' : len(data),
        'MIDDLE-COMPONENT-INDEXES': _join(range(1, vc - 1)),
        'LAST-COMPONENT-INDEX': vc - 1
    }

    if vc > 1:
        subs['%KEEP-IF-1C%'] = ''
    
    scaled = snif["inference_parameters"]["infer_scale"]
    sizes = not (s_limits[0] == 1.0 and s_limits[1] == 1.0)
    template_fn = './templates/versus-{}{}.tex'.format(
        'y' if scaled else 'n',
        'y' if sizes  else 'n'
    )
    
    basename = config.SNIF_basename
    for tag in ['middle', 'bottom']:
        support_fn = './templates/vrow-{}-tm{}.tex'.format(tag, 's' if sizes else '')
        copy2(support_fn, basename)

    target_fn = os.path.join(basename, os.path.basename(basename) + '-params-versus.tex')
    _replace_and_save(template_fn, target_fn, subs)
    return [target_fn]
def _generate_tests_plots(data, config, snif, plot_type, curves_data_fn):
    if plot_type == "iicr":
        style = config.IICR_plots_style
    elif plot_type == "cdf":
        style = config.CDF_plots_style
    elif plot_type == "pdf":
        style = config.PDF_plots_style
    else:
        raise ValueError("Invalid plot type string")
    
    if style == OutputStyle.Excluded:
        return []
    elif style == OutputStyle.Minimal:
        template_fn = './templates/{}-minimal.tex'.format(plot_type)
    else:
        template_fn = './templates/{}-full.tex'.format(plot_type)
    
    inf_c = snif["inference_parameters"]["number_of_components"]
    if "simulation_parameters" in snif:
        simulating = True
        sim_c = snif["simulation_parameters"]["number_of_components"]
    else:
        simulating = False

    (x_scale, x_limits, x_fit_limits, x_label) = _get_x_data(config, snif)
    y_max = _get_curve_upper_bounds(data, plot_type)

    basename = config.SNIF_basename
    generated_filenames = []
    
    tests = len(data)
    if isinstance(config.test_numbers, list):
        indexes = list(set(range(tests)) & set([x - 1 for x in config.test_numbers]))
    else:
        indexes = range(tests)

    subs = {
        'CURVES-DATA-FILENAME': os.path.basename(curves_data_fn),
        'PLOT-WIDTH': str(config.plot_width_cm) + 'cm',
        'PLOT-HEIGHT': str(config.plot_height_cm) + 'cm',
        'X-AXIS-LABEL': x_label,
        'X-LOWER-LIMIT': x_limits[0],
        'X-UPPER-LIMIT': x_limits[1],
        'X-LOWER-DIST-LIMIT': x_fit_limits[0],
        'X-UPPER-DIST-LIMIT': x_fit_limits[1],
    }

    if config.one_file_per_test:
        for i in indexes:
            target_fn = os.path.join(basename, os.path.basename(basename) + '-IICR-plot-{}.tex'.format(i + 1))
            generated_filenames.append(target_fn)
            inferred_times = (x_scale * data[i]["inf. t" + str(j)] for j in range(1, inf_c))
            simulated_times = []
            if simulating:
                simulated_times = (x_scale * data[i]["sim. t" + str(j)] for j in range(1, sim_c))
            else:
                try:
                    (_, t, _, _, _) = _read_source_scenario(data[i]["source filename"], config)
                    simulated_times = t
                except Exception:
                    pass
                    
            subs.update({
                'Y-UPPER-LIMIT': y_max[i],
                'INFERRED-EVENT-TIMES': _join(inferred_times, True),
                'SIMULATED-EVENT-TIMES': _join(simulated_times, True),
                'TEST-INDEXES': i + 1,
            })
            _replace_and_save(template_fn, target_fn, subs)
    else:
        target_fn = os.path.join(basename, os.path.basename(basename) + '-IICR.tex')
        generated_filenames.append(target_fn)
        inferred_times = (x_scale * data[i]["inf. t" + str(j)] for j in range(1, inf_c) for i in indexes)

        simulated_times = []
        if simulating:
            simulated_times = (x_scale * data[i]["sim. t" + str(j)] for j in range(1, sim_c) for i in indexes)
        else:
            try:
                for i in indexes:
                    (_, t, _, _, _) = _read_source_scenario(data[i]["source filename"], config)
                    simulated_times.extend(t)
            except Exception:
                pass        

        subs.update({
            'Y-UPPER-LIMIT': max([y_max[i] for i in indexes]),
            'INFERRED-EVENT-TIMES': _join(inferred_times, True),
            'SIMULATED-EVENT-TIMES': _join(simulated_times, True),
            'TEST-INDEXES': _join(i + 1 for i in indexes),
        })
        _replace_and_save(template_fn, target_fn, subs)
    
    return generated_filenames
def _generate_1D_plot(data, config, snif, plot_type, params2_data_fn, hists_data_fn):
    if plot_type == 'n':
        (bound_name, template_tag, label) = ("bounds_islands", "islands", "Islands")
        style = config.islands_plot_style
    elif plot_type == 'N_ref':
        (bound_name, template_tag, label) = ("bounds_effective_size", "ref-sizes", "$N_{\\mathrm{ref}}$")
        style = config.Nref_plot_style
    else:
        raise ValueError("Invalid plot type string")

    limits = _get_inferred_param_bounds(snif, bound_name)
    limits[0] = 0
    # if limits[0] == 2:
    #     limits[0] = 0

    simulated_vals = []
    if "simulation_parameters" in snif:
        simulated_vals.extend(test["sim. " + plot_type] for test in data)
    
    source_fn = config.CG_source
    if source_fn != '':
        (n, _, _, _, N_ref) = _read_source_scenario(source_fn, config)
        simulated_vals.append(eval(plot_type))

    subs = {
        'PLOT-HEIGHT' : str(config.plot_height_cm) + "cm",
        'PARAMETERS2-DATA-FILENAME': os.path.basename(params2_data_fn),
        'Y-AXIS-LABEL': label,
        'Y-LOWER-LIMIT': limits[0],
        'Y-UPPER-LIMIT': limits[1],
        'PARAMETER-NAME': plot_type,
        'SIMULATED-VALUES' : '\n'.join(str(x) for x in simulated_vals)
    }

    if style == OutputStyle.Full:        
        subs['HISTOGRAMS-DATA-FILENAME'] = os.path.basename(hists_data_fn)
        hists = 'y'
    else:
        hists = 'n'
    
    template_fn = './templates/1d-{}.tex'.format(hists)
    basename = config.SNIF_basename
    target_fn = os.path.join(basename, os.path.basename(basename) + '-{}.tex'.format(template_tag))
    _replace_and_save(template_fn, target_fn, subs)
    return [target_fn]
def _generate_cg_plot(data, config, snif, params_data_fn, hists_data_fn):
    style = config.CG_style
    if style == OutputStyle.Excluded:
        return []

    (x_scale, x_limits, _, x_label) = _get_x_data(config, snif)
    M_limits = _get_inferred_param_bounds(snif, 'bounds_migrations_rates')
    tests = len(data)
    subs = {
        'PLOT-WIDTH' : str(config.plot_width_cm) + "cm",
        'PLOT-HEIGHT' : str(config.plot_height_cm) + "cm",
        'PARAMETERS-DATA-FILENAME': os.path.basename(params_data_fn),
        'X-AXIS-LABEL': x_label,
        'X-LOWER-LIMIT': x_limits[0],
        'X-UPPER-LIMIT': x_limits[1],
        'M-LOWER-LIMIT': M_limits[0],
        'M-UPPER-LIMIT': M_limits[1],
    }
    
    sizes = config.CG_size_history
    simulating = "simulation_parameters" in snif
    histograms = (simulating) or (style == OutputStyle.Full)
    source_fn = config.CG_source
    
    if sizes:
        s_limits = _get_inferred_param_bounds(snif, 'bounds_deme_sizes')
        subs.update({
            'S-LOWER-LIMIT': s_limits[0],
            'S-UPPER-LIMIT': s_limits[1], 
        })
    
    if simulating:
        sim_m = snif["simulation_parameters"]["migrations_rates"]
        sim_t = snif["simulation_parameters"]["event_times"]
        if snif["simulation_parameters"]["simulate_scale"]:
            sim_n0 = snif["simulation_parameters"]["effective_size"]
            scaled_sim_t = (x_scale * 2 * n0 * t for t in sim_t for n0 in sim_n0)
        else:
            scaled_sim_t = sim_t
        
        subs.update({
            'NUMBER-OF-TESTS': tests,
            'SIMULATED-TIMES' : _join(scaled_sim_t),
            'SIMULATED-MIG-RATES' : _join(sim_m),
        })
        if sizes:
            sim_s = snif["simulation_parameters"]["deme_sizes"]
            subs.update({
                'SIMULATED-SIZES': _join(sim_s),
            })            
    else:
        subs.update({
            'LEGEND-ENTRIES': _get_legend_entries(data),
            'LIST-OF-INDEXES': _join(range(1, tests + 1)),
        })

    if histograms:
        subs.update({
            'HISTOGRAMS-DATA-FILENAME': os.path.basename(hists_data_fn),
        })

    if source_fn != '':
        (_, t, m, s, _) = _read_source_scenario(source_fn, config)
        m_coordinates = list(zip(t, m))
        subs.update({
            '%REMOVE-IF-SOURCE%' : '',
            'SOURCE-M-COORDINATES': ' '.join(str(x) for x in m_coordinates),
            'SOURCE-LEGEND-ENTRY' : config.CG_source_legend,
            'SIMULATED-TIMES' : _join(t),
            'SIMULATED-MIG-RATES' : _join(m),
        })
        if sizes:
            s_coordinates = list(zip(t, s))
            subs.update({
                'SOURCE-S-COORDINATES': ' '.join(str(x) for x in s_coordinates),
                'SIMULATED-SIZES': _join(s),
            })

    basename = config.SNIF_basename
    template_fn = _get_cg_template_fn(sizes, simulating, histograms)
    target_fn = os.path.join(basename, os.path.basename(basename) + '-connectivity_graph.tex')
    _replace_and_save(template_fn, target_fn, subs)
    return [target_fn]
def _compile_tex_files(tex_filenames):
    current_dir = os.getcwd()
    os.chdir(os.path.dirname(tex_filenames[0]))
    for i in range(len(tex_filenames)):
        output_basename = os.path.basename(tex_filenames[i])
        output_basename = os.path.splitext(output_basename)[0]
        print('\nAttempting to compile LaTeX file ' + os.path.basename(tex_filenames[i]))
        os.system('pdflatex -interaction=batchmode -halt-on-error ' + os.path.basename(tex_filenames[i]))
        if os.path.exists(output_basename + '.pdf'):
            print('compilation successfull!')
            if os.path.exists(output_basename + '.log'):
                os.remove(output_basename + '.log')
            if os.path.exists(output_basename + '.aux'):
                os.remove(output_basename + '.aux')
            if os.path.exists(output_basename + '.out'):
                os.remove(output_basename + '.out')
        else:
            print('error in compilation. Please review the LaTeX .log file')
            print('file {} could not be located!'.format(output_basename + '.pdf'))    
    
    os.chdir(current_dir)

def _generate_params_data(data, config, snif):
    table = []
    c = snif["inference_parameters"]["number_of_components"]
    for i in range(len(data)):
        if config.scaling_units == TimeScale.Years:
            time_scale = config.generation_time
        else:
            time_scale = 1.0
        
        t = ['1e-10'] + [str(time_scale * data[i]['inf. t' + str(j)]) for j in range(1, c)]
        t.append('1e10')
        table.append(['t_' + str(i + 1)] + t)
        
        m = [str(data[i]['inf. M' + str(j)]) for j in range(c)]
        m.append(m[-1])
        table.append(['M_' + str(i + 1)] + m)

        s = [str(data[i]['inf. s' + str(j)]) for j in range(c)]
        s.append(s[-1])
        table.append(['s_' + str(i + 1)] + s)

        if "simulation_parameters" in snif:
            sc = snif["simulation_parameters"]["number_of_components"]

            t = ['1e-10'] + [str(time_scale * data[i]['sim. t' + str(j)]) for j in range(1, sc)]
            t.append('1e10')
            table.append(['s_t_' + str(i + 1)] + t)
            
            m = [str(data[i]['sim. M' + str(j)]) for j in range(sc)]
            m.append(m[-1])
            table.append(['s_M_' + str(i + 1)] + m)

            s = [str(data[i]['sim. s' + str(j)]) for j in range(sc)]
            s.append(s[-1])
            table.append(['s_s_' + str(i + 1)] + s)        
        
    table = zip(*(table))
    basename = config.SNIF_basename
    table_filename = os.path.join(basename, os.path.basename(basename) + '_parameters.dat')
    _write_csv(table, table_filename)
    return table_filename
def _generate_params2_data(data, config, snif):    
    simulating = "simulation_parameters" in snif
    ic = snif["inference_parameters"]["number_of_components"]

    if simulating:
        sc = snif["simulation_parameters"]["number_of_components"]

    head = ['id']
    head.extend(
        ['n']
      + ['t_' + str(i) for i in range(1, ic)]
      + ['M_' + str(i) for i in range(ic)]
      + ['s_' + str(i) for i in range(ic)]
      + ['N_ref']
    )

    if simulating:
        head.extend(
            ['sn'] 
          + ['st_' + str(i) for i in range(1, sc)]
          + ['sM_' + str(i) for i in range(sc)]
          + ['ss_' + str(i) for i in range(sc)]
          + ['sN_ref']
        )
    
    table = [head]
    for test in data:
        line = [str(test['id'])]
        line.extend(
            [str(test['inf. n'])]
          + [str(test['inf. t' + str(i)]) for i in range(1, ic)]
          + [str(test['inf. M' + str(i)]) for i in range(ic)]
          + [str(test['inf. s' + str(i)]) for i in range(ic)]
          + [str(test['inf. N_ref'])]
        )

        if simulating:
            line.extend(
                [str(test['sim. n'])] 
              + [str(test['sim. t' + str(i)]) for i in range(1, sc)]
              + [str(test['sim. M' + str(i)]) for i in range(sc)]
              + [str(test['sim. s' + str(i)]) for i in range(sc)] 
              + [str(test['sim. N_ref'])]              
            )

        table.append(line)
    
    basename = config.SNIF_basename
    table_filename = os.path.join(basename, os.path.basename(basename) + '_parameters2.dat')
    _write_csv(table, table_filename)
    return table_filename
def _generate_curves_data(data, config, snif):
    curves = []
    (x_scale, _, _, _) = _get_x_data(config, snif)
    for i in range(len(data)):
        raw_t = data[i][_curve_names[0]]
        curves.append(["{}-{}".format(_curve_names[0], i + 1)] + [x_scale * t for t in raw_t])
        curves.extend(["{}-{}".format(cn, i + 1)] + data[i][cn] for cn in _curve_names[1:])
    
    longest = max(len(line) for line in curves)
    for i in range(len(curves)):
        dif = longest - len(curves[i])
        if dif > 0:
            curves[i].extend([curves[i][-1]] * dif)

    # if config.scaling_units == TimeScale.Years:
    #     for i in range(len(data)):
    #         curves[i][1:] = [curves[i][j] * config.generation_time for j in range(1, longest)]

    curves = list(zip(*(curves)))
    
    basename = config.SNIF_basename
    curves_filename = os.path.join(basename, os.path.basename(basename) + '_curves.dat')
    _write_csv(curves, curves_filename)
    return curves_filename
def _generate_hists_data(data, config, snif):
    from numpy import histogram
    from math import log10
    
    ####  t  ####
    (x_scale, x_limits, _, _) = _get_x_data(config, snif)
    c = snif["inference_parameters"]["number_of_components"]
    all_inferred_times = []
    for test in data:
        all_inferred_times.extend((x_scale * test["inf. t" + str(i)] for i in range(1, c)))

    t_hist = histogram(
        [log10(x) for x in all_inferred_times], 
        bins = config.time_histograms_bins, 
        range = (log10(x_limits[0]), log10(x_limits[1]))
    )
    t_bins = [str(10**x) for x in t_hist[1]]
    t_vals = [str(x) for x in t_hist[0]] + ['0']

    ####  M  ####
    simulating = "simulation_parameters" in snif
    m_limits = _get_inferred_param_bounds(snif, 'bounds_migrations_rates')
    special_m = [1]#0.1, 1, 10, 50]
    all_inferred_mrates = []
    all_inferred_special_mrates = []
    for test in data:
        if simulating:
            for i in range(c):
                sim_m = test["sim. M" + str(i)]
                inf_m = test["inf. M" + str(i)]
                if sim_m in special_m:
                    all_inferred_special_mrates.append(inf_m)
                else:
                    all_inferred_mrates.append(inf_m)
        else:
            all_inferred_mrates.extend((test["inf. M" + str(i)] for i in range(c)))
    m_hist = histogram(
        [log10(x) for x in all_inferred_mrates], 
        bins = config.migration_histograms_bins, 
        range = (log10(m_limits[0]), log10(m_limits[1]))
    )
    m_bins = [str(10**x) for x in m_hist[1]]
    m_vals = [str(x) for x in m_hist[0]] + ['0']

    special_m_hist = histogram(
        [log10(x) for x in all_inferred_special_mrates], 
        bins = config.migration_histograms_bins, 
        range = (log10(m_limits[0]), log10(m_limits[1]))
    )
    special_m_bins = [str(10**x) for x in special_m_hist[1]]
    special_m_vals = [str(x) for x in special_m_hist[0]] + ['0']

    ####  s  ####
    s_limits = _get_inferred_param_bounds(snif, 'bounds_deme_sizes')
    all_inferred_sizes = []
    for test in data:
        all_inferred_sizes.extend((test["inf. s" + str(i)] for i in range(c)))
    s_hist = histogram(
        all_inferred_sizes, 
        bins = config.size_histograms_bins, 
        range = s_limits
    )
    s_bins = [str(x) for x in s_hist[1]]
    s_vals = [str(x) for x in s_hist[0]] + ['0']

    ####  n  ####
    n_limits = _get_inferred_param_bounds(snif, 'bounds_islands')
    all_inferred_demes = [test["inf. n"] for test in data]
    n_hist = histogram(
        all_inferred_demes, 
        bins = config.islands_histograms_bins, 
        range = n_limits
    )
    n_bins = [str(x) for x in n_hist[1]]
    n_vals = [str(x) for x in n_hist[0]] + ['0']

    ####  N_ref  ####
    Nref_limits = _get_inferred_param_bounds(snif, 'bounds_effective_size')
    all_inferred_Nrefs = [test["inf. N_ref"] for test in data]
    Nref_hist = histogram(
        all_inferred_Nrefs, 
        bins = config.Nref_histograms_bins, 
        range = Nref_limits
    )
    Nref_bins = [str(x) for x in Nref_hist[1]]
    Nref_vals = [str(x) for x in Nref_hist[0]] + ['0']

    def _pad(v, min_len):
        diff = min_len - len(v)
        if diff > 0:
            v += [v[-1]] * diff

    all_vals = [t_vals, m_vals, special_m_vals, s_vals, n_vals, Nref_vals]
    all_bins = [t_bins, m_bins, special_m_bins, s_bins, n_bins, Nref_bins]
    max_len = max(*[len(v) for v in all_vals])
    for v in all_vals + all_bins:
        _pad(v, max_len)

    histograms_table = [['t-bins'] + t_bins, ['t-vals'] + t_vals]
    histograms_table.extend([['m-bins'] + m_bins, ['m-vals'] + m_vals])
    histograms_table.extend([['special-m-bins'] + special_m_bins, ['special-m-vals'] + special_m_vals])
    histograms_table.extend([['s-bins'] + s_bins, ['s-vals'] + s_vals])
    histograms_table.extend([['n-bins'] + n_bins, ['n-vals'] + n_vals])
    histograms_table.extend([['N_ref-bins'] + Nref_bins, ['N_ref-vals'] + Nref_vals])

    histograms_table = zip(*(histograms_table))
    basename = config.SNIF_basename
    table_filename = os.path.join(basename, os.path.basename(basename) + '_histograms.dat')
    _write_csv(histograms_table, table_filename)
    return table_filename
def _generate_fitness_stats(data, config, snif):
    ic = snif["inference_parameters"]["number_of_components"]
    sc = snif["simulation_parameters"]["number_of_components"]    
    c = min(ic, sc)
    n_tests = len(data)
    scaled = snif["inference_parameters"]["infer_scale"]

    head = ['c', 'n', 'M0']
    if scaled:
        head.insert(2, 'N_ref')
    for i in range(1, c):
        head.extend(['t' + str(i), 'M' + str(i)])

    r2_line = ['r2']
    for param in head[1:]:
        if False:#param.startswith('t'):
            x = [log10(test['sim. ' + param]) for test in data]
            y = [log10(test['inf. ' + param]) for test in data]
        else:
            x = [test['sim. ' + param] for test in data]
            y = [test['inf. ' + param] for test in data]
        # if param == 'N_ref':
        #     from statistics import stdev, mean
        #     r2 = stdev(y) / mean(y)
        # else:
        #     x = [test['sim. ' + param] for test in data]
        #     r2 = _r_squared(x, y)
        r2 = _nrmsd(x, y)
        r2_s = str(round(r2, 3))
        # if(int(r2_s[0]) == 0 and float(r2_s) > 0):
        #     r2_s = r2_s[1:]

        r2_line.append(r2_s)
    
    table = [head]#, r2_line]
    
    p_values = list(range(1, 101)) 
    # p_values += [10 * i for i in range(1, 10)] 
    # p_values += [100 * i for i in range(1, 11)]
    # p_values = list(range(1, 1001))
    for p in p_values:
        line = [str(p)]
        for param in head[1:]:
            p_n = 0
            for test in data:
                x = test['sim. '+param]
                x0 = test['inf. '+param]
                if - p * x0 < 100 * (x0 - x) < p * x: #100 * abs(x - x0) / x < p:
                    p_n += 1
            ratio = 100 * p_n / n_tests
            ratio_s = str(round(ratio, 3))
            if(int(ratio_s[0]) == 0 and float(ratio_s) > 0):
                ratio_s = ratio_s[1:]
            line.append(ratio_s)

        table.append(line)
    
    basename = config.SNIF_basename
    table_filename = os.path.join(basename, os.path.basename(basename) + '_fitness_stats.dat')
    _write_csv(table, table_filename)
    return table_filename

def _r_squared(x, y):
    from scipy.stats import linregress
    slope, intercept, r_value, p_value, std_err = linregress(x, y)
    return r_value ** 2
def _nrmsd(x, y):
    from statistics import mean
    from math import sqrt
    n = len(x)
    rmsd = sqrt(sum((x[i] - y[i]) ** 2 for i in range(n)) / n)
    return rmsd / mean(x)


def _replace_and_save(source_fn, target_fn, substitutions):
    source_file = open(source_fn, 'r')
    target_file = open(target_fn, 'w')
    for line in source_file:
        new_line = line
        for key, val in substitutions.items():
            new_line = new_line.replace(key, str(val), 1)

        target_file.write(new_line)    
def _get_x_data(config, snif):
    limits = snif["inference_parameters"]["data_cutoff_bounds"]
    fit_limits = snif["inference_parameters"]["distance_computation_interval"]

    scale = 1.0
    if snif["inference_parameters"]["infer_scale"]:
        timeScale = config.scaling_units
        if timeScale == TimeScale.Generations:
            label = 'Time in generations'
        elif timeScale == TimeScale.Years:
            scale = config.generation_time
            label = 'Time in years ({} y/gen)'.format(scale)
            limits = [scale * t for t in limits]
            fit_limits = [scale * t for t in fit_limits]
            # limits[0] *= scale
            # limits[1] *= scale            
            # fit_limits[0] *= scale
            # fit_limits[1] *= scale
        else:
            raise ValueError("Invalid scaling units")    
    else:
        label = 'Coalescent time'
    
    return (scale, limits, fit_limits, label)
def _get_curve_upper_bounds(data, curve_type):
    tests = len(data)
    y_max = [1] * tests
    if curve_type == "cdf":
        return y_max
    
    from math import log10, floor, ceil
    for i in range(tests):
        source_curve = _validate_curve(data[i]["source-" + curve_type])
        inferred_curve = _validate_curve(data[i]["inferred-" + curve_type])
        y_max[i] = max(1, *source_curve, *inferred_curve)

        magnitude = 10 ** floor(log10(y_max[i]))
        y_max[i] = magnitude * ceil(y_max[i] / magnitude) 

    return y_max
def _validate_curve(curve):
    from math import inf
    return [x for x in curve if 0 < x and x < inf]
def _get_inferred_param_bounds(snif, parameter):
    limits = snif["inference_parameters"][parameter]
    if len(limits) == 0:
        return [1, 1]

    if not(isinstance(limits[0], int) or isinstance(limits[0], float)):
        low = min(l[0] for l in limits)
        high = max(l[1] for l in limits)
        limits = [low, high]
    
    return limits
def _get_legend_entries(data):
    from collect import get_unique_parts

    source_filenames = [test['source filename'] for test in data]
    if len(source_filenames) == 1:
        return _sanitize_for_latex(os.path.basename(source_filenames[0]))
    
    parts = get_unique_parts(source_filenames)
    if parts[0] == '':
        parts = ['Test #{}'.format(int(test['id'])) for test in data]

    return ','.join(_sanitize_for_latex(x) for x in parts)
def _sanitize_for_latex(line):
    special = ['&', '$', '%', '#', '_', '{', '}']
    for x in special:
        line = line.replace(x, '\\' + x)
        
    return line
def _get_cg_template_fn(sizes, simulating, histograms):
    if sizes:
        if simulating:
            return './templates/cg-sim-ms.tex'
        elif histograms:
            return './templates/cg-yy.tex'
        else:
            return './templates/cg-ny.tex'
    else:
        if simulating:
            return './templates/cg-sim-m.tex'
        elif histograms:
            return './templates/cg-yn.tex'
        else:
            return './templates/cg-nn.tex'
def _read_source_scenario(filename, config):

    timeScale = config.scaling_units
    
    source = _read_json(filename)["PSNIC"]
    n = source["n"]    

    t = source["t"]
    t.append(1e10)
    t[0] = 1e-10
        
    if "Nref" in source:
        nref = source["Nref"]
        if timeScale == TimeScale.Coalescent:
            t = [ti / (2 * nref) for ti in t]
        elif timeScale == TimeScale.Years:
            t = [ti * config.generation_time for ti in t]
    else:
        nref = 0
        if timeScale != TimeScale.Coalescent:
            Warning("The demographic model in the specified JSON description cannot be scaled because there is no 'Nref' information.")

    m = source["M"]
    m.append(m[-1])    

    s = source["s"]
    s.append(s[-1])

    return (n, t, m, s, nref)
def _join(iterator, braces = False):
    string = ','.join(str(x) for x in iterator)
    if braces:
        return '{' + string + '}'
    
    return string
def _sanitize_for_latex(line):
    special = ['&', '$', '%', '#', '_', '{', '}']
    for x in special:
        line = line.replace(x, '\\' + x)
        
    return line

def _read_json(filename):
    import json
    with open(filename) as f:
        output_dict = json.load(f)
    
    return output_dict
def _read_csv(filename):
    import csv
    return csv.reader(open(filename), delimiter = ',')
def _write_csv(table, filename):
    import csv
    csv.writer(open(filename, 'w', newline=''), delimiter='\t').writerows(table)


# # # # # # # # 
#             # 
#   M A I N   # 
#             # 
# # # # # # # # 


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--settings', help = 'A JSON file containing a configuration for SNIF2TeX.')
    args = parser.parse_args()

    if args.settings is not None:
        json_config = _read_json(args.settings)
        current_configuration = Configuration(**json_config)
        TeXify(current_configuration)
    else:
        TeXify(default_configuration)