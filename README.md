# 20230316_normalization_of_CCLE_miRNA_expression_data.R

# 解析の目的
CCLEのmiRNA発現データを様々な方法でノーマライズし、meanSdPlotでそれらのノーマライズがうまくいっているか評価する

# ファイルの説明
・20230316_normalization_of_CCLE_miRNA_expression_data.R
<br>CCLEのmiRNA発現データを様々な方法でノーマライズし、meanSdPlotでそれらのノーマライズがうまくいっているか評価するためのスクリプト。スクリプトを走らせるにあたって必要なファイルとその場所はスクリプト中に記載されている。CCLEで始まる以下４つのtxtファイルを生成する。
<br>
<br>・CCLE_miRNA_quantile_normalization.txt
<br>四分位数がすべてのサンプルで同じになるようにノーマライズ
<br>
<br>・CCLE_miRNA_median_normalization.txt
<br>中央値がすべてのサンプルで同じになるようにノーマライズ
<br>
<br>・CCLE_miRNA_variance_normalization_.txt
<br>すべてのサンプルで分散が安定するようにノーマライズ
<br>
<br>・CCLE_miRNA_mean_normalization_.txt
<br>ある遺伝子発現量をすべてのサンプル間における平均値で割ってノーマライズ
