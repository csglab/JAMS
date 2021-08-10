#!/bin/env python3
import base64
from IPython.display import Image
import os
import sys

## Usage: create_html.py experiment_name path_to_images

experiment_name = sys.argv[1]
in_path = sys.argv[2]
# in_path = "./"

out_html = in_path + experiment_name + ".html"

# Read log file
log_path = in_path + "/log.txt"
log_f = open(log_path, "r")
log = "".join(log_f.readlines()).replace("\n", "<BR>")

#  plots
img_path = in_path + "methylation_correlation_with_binding.png"
methyl_corr_img = base64.b64encode(open(img_path, 'rb').read()).decode('utf-8').replace('\n', '')

img_path = in_path + "methylation_stats.png"
methyl_stats_img = base64.b64encode(open(img_path, 'rb').read()).decode('utf-8').replace('\n', '')

img_path = in_path + "plot_DNA_acc_around_ptfbs.png"
dna_acc_lin_img = base64.b64encode(open(img_path, 'rb').read()).decode('utf-8').replace('\n', '')

img_path = in_path + "acc_boxplot.jpg"
acc_log_img = base64.b64encode(open(img_path, 'rb').read()).decode('utf-8').replace('\n', '')

img_path = in_path + "plot_length_vs_tags.jpg"
len_tags_img = base64.b64encode(open(img_path, 'rb').read()).decode('utf-8').replace('\n', '')

img_path = in_path + "mCpG_only_model.motif_logo.png"
logo_mCpG_img = base64.b64encode(open(img_path, 'rb').read()).decode('utf-8').replace('\n', '')

img_path = in_path + "CpG_only_coefficients.png"
CpG_coef_img = base64.b64encode(open(img_path, 'rb').read()).decode('utf-8').replace('\n', '')

img_path = in_path + "mCpG_only_model_total_cv_scatterplot_cv_poisson.jpg"
scatter_mCpG_img = base64.b64encode(open(img_path, 'rb').read()).decode('utf-8').replace('\n', '')

img_path = in_path + "All_mC_model.motif_logo.png"
logo_mC_img = base64.b64encode(open(img_path, 'rb').read()).decode('utf-8').replace('\n', '')

img_path = in_path + "all_mC_coefficients.png"
mC_coef_img = base64.b64encode(open(img_path, 'rb').read()).decode('utf-8').replace('\n', '')

img_path = in_path + "All_mC_model_total_cv_scatterplot_cv_poisson.jpg"
scatter_mC_img = base64.b64encode(open(img_path, 'rb').read()).decode('utf-8').replace('\n', '')

img_path = in_path + "Merged_model.motif_logo.png"
logo_merge_img = base64.b64encode(open(img_path, 'rb').read()).decode('utf-8').replace('\n', '')

img_path = in_path + "mCpG_only_model_scatterplot_hek293_model_this_tissue.jpg"
mCpG_scatter_other_img = base64.b64encode(open(img_path, 'rb').read()).decode('utf-8').replace('\n', '')

img_path = in_path + "all_mC_model_scatterplot_hek293_model_this_tissue.jpg"
mC_scatter_other_img = base64.b64encode(open(img_path, 'rb').read()).decode('utf-8').replace('\n', '')


html_page = """
<HEAD>
<TITLE>Methyl-ChIP report</TITLE>
</HEAD>
<BODY BGCOLOR="WHITE">
<CENTER>
<H1>LOG</H1>
{log}
<H1>Descriptive plots</H1>
<H1>Methylation correlation with binding</H1>
<img src="data:image/png;base64,{methyl_corr}" alt="Image" width="{width}">

<H1>Methylation stats</H1>
<img src="data:image/png;base64,{methyl_stats}" alt="Image" width="{width}">

<H1>Average DNA accessibility in the TFBS region</H1>
<img src="data:image/png;base64,{dna_acc_lin}" alt="Image" width="{width}">

<H1>Log of DNA accessibility in the TFBS region</H1>
<img src="data:image/png;base64,{acc_log}" alt="Image" width="{width}">

<H1>Lenght vs tags</H1>
<img src="data:image/png;base64,{len_tags}" alt="Image" width="{width}">

<H1>CpG methylation only model</H1>
<img src="data:image/png;base64,{logo_mCpG}" alt="Image" width="{width}">

<H1>DNA acc. coefficients (mCpG)</H1>
<img src="data:image/png;base64,{CpG_coef}" alt="Image" width="{width}">

<H1>Binding affinity: predictions vs observed (mCpG)</H1>
<img src="data:image/png;base64,{scatter_mCpG}" alt="Image" width="{width}">

<H1>All Cytosine methylation mode</H1>
<img src="data:image/png;base64,{logo_mC}" alt="Image" width="{width}">

<H1>DNA acc. coefficients (all mC)</H1>
<img src="data:image/png;base64,{mC_coef}" alt="Image" width="{width}">

<H1>Binding affinity: predictions vs observed (All mC)</H1>
<img src="data:image/png;base64,{scatter_mC}" alt="Image" width="{width}">

<H1>Merge model </H1>
<img src="data:image/png;base64,{logo_merge}" alt="Image" width="{width}">

<H1>Using Hek293 trained model in this tissue: CpG model</H1>
<img src="data:image/png;base64,{mCpG_scatter_other}" alt="Image" width="{width}">

<H1>Using Hek293 trained model in this tissue: All mC  model</H1>
<img src="data:image/png;base64,{mC_scatter_other}" alt="Image" width="{width}">
</body>
""".format(log = log, \
               methyl_corr = methyl_corr_img, \
               methyl_stats = methyl_stats_img, \
               dna_acc_lin = dna_acc_lin_img, \
               acc_log = acc_log_img, \
               len_tags = len_tags_img, \
               logo_mCpG = logo_mCpG_img, \
               CpG_coef = CpG_coef_img, \
               scatter_mCpG = scatter_mCpG_img, \
               logo_mC = logo_mC_img, \
               mC_coef = mC_coef_img, \
               scatter_mC = scatter_mC_img, \
               logo_merge = logo_merge_img, \
               mCpG_scatter_other = mCpG_scatter_other_img , \
               mC_scatter_other = mC_scatter_other_img, \
               width = "1500")

with open(out_html, "w") as html_f:

    html_f.write(html_page)










