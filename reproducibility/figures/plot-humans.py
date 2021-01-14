from snif2tex import *

filenames = [
    './.results/humans_4/c02YN_dA_w100_SQ__yr',
    './.results/humans_4/c03YN_dA_w100_SQ__yr',
    './.results/humans_4/c04YN_dA_w100_SQ__yr',
    './.results/humans_4/c05YN_dA_w100_SQ__yr',
    './.results/humans_4/c02YN_dA_w050_SQ__yr',
    './.results/humans_4/c03YN_dA_w050_SQ__yr',
    './.results/humans_4/c04YN_dA_w050_SQ__yr',
    './.results/humans_4/c05YN_dA_w050_SQ__yr',
    './.results/humans_4/c02YN_dA_w020_SQ__yr',
    './.results/humans_4/c03YN_dA_w020_SQ__yr',
    './.results/humans_4/c04YN_dA_w020_SQ__yr',
    './.results/humans_4/c05YN_dA_w020_SQ__yr',
    './.results/humans_4/c02YN_dA_w100_SQ__nr',
    './.results/humans_4/c03YN_dA_w100_SQ__nr',
    './.results/humans_4/c04YN_dA_w100_SQ__nr',
    './.results/humans_4/c05YN_dA_w100_SQ__nr',
    './.results/humans_4/c02YN_dA_w050_SQ__nr',
    './.results/humans_4/c03YN_dA_w050_SQ__nr',
    './.results/humans_4/c04YN_dA_w050_SQ__nr',
    './.results/humans_4/c05YN_dA_w050_SQ__nr',
    './.results/humans_4/c02YN_dA_w020_SQ__nr',
    './.results/humans_4/c03YN_dA_w020_SQ__nr',
    './.results/humans_4/c04YN_dA_w020_SQ__nr',
    './.results/humans_4/c05YN_dA_w020_SQ__nr',
]

for fn in filenames:
    config = Configuration(
        SNIF_basename = fn,
        plot_width_cm = 13,
        plot_height_cm = 6,
        IICR_plots_style  = OutputStyle.Full,
        PDF_plots_style = OutputStyle.Excluded,
        CDF_plots_style = OutputStyle.Excluded,
        islands_plot_style = OutputStyle.Full,
        Nref_plot_style = OutputStyle.Full,
        test_numbers = "all",
        one_file_per_test = True,
        versus_plot_style = OutputStyle.Excluded,
        CG_style = OutputStyle.Full,
        CG_size_history = False,
        CG_source = './.data/humans-rodriguez2018.json',
        CG_source_legend = "Rodr\\'iguez 2018",
        Nref_histograms_bins = 100,
        islands_histograms_bins = 100,
        time_histograms_bins = 100,
        migration_histograms_bins = 100,
        size_histograms_bins = 100,
        scaling_units = TimeScale.Years,
        generation_time = 25,
    )

    TeXify(config)