from snif2tex import *

filenames = [
    './.results/humans_4/seqsim-validations/yoruba/snif/200312-040005_c05YN_dA_w020_SS_yb',
]

for fn in filenames:
    config = Configuration(
        SNIF_basename = fn,
        plot_width_cm = 10,
        plot_height_cm = 5,
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
        CG_source = './.results/humans_4/selected/c05YN_dA_w020_SQ__yr_inferred-model-005-Yoruba.json',
        CG_source_legend = "Inferred history from Yoruba PSMC",
        Nref_histograms_bins = 100,
        islands_histograms_bins = 100,
        time_histograms_bins = 100,
        migration_histograms_bins = 100,
        size_histograms_bins = 100,
        scaling_units = TimeScale.Years,
        generation_time = 25,
    )

    TeXify(config)