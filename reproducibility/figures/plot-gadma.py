from snif2tex import *

for c in ['02', '03', '04', '05']:
    for r in ['', '_restricted']:
        for w in ['025', '050', '100']:

            fn = './.results/gadma/c{}YN_dA_w{}_MT{}'.format(c, w, r)

            config = Configuration(
                SNIF_basename = fn,
                plot_width_cm = 13,
                plot_height_cm = 6,
                IICR_plots_style  = OutputStyle.Full,
                PDF_plots_style = OutputStyle.Excluded,
                CDF_plots_style = OutputStyle.Excluded,
                islands_plot_style = OutputStyle.Excluded,
                Nref_plot_style = OutputStyle.Excluded,
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

            TeXify(config)