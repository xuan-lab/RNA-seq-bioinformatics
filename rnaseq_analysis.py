# Filename: RNA-seq.py
# Author: Xiexuan <xiexuan@kernel-dev.com>

import os
import subprocess
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from rpy2.robjects import r, pandas2ri
import logging

# 启用 R-DataFrame 到 pandas 的转换
pandas2ri.activate()

# 设置日志记录
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# 全局配置
FASTQ_DIR = "./data/fastq"
RESULTS_DIR = "./results"
GENOME_INDEX = "./genome/hisat2_index/genome"  # Hisat2 索引文件的前缀
GTF_FILE = "./genome/annotations.gtf"  # 基因组注释文件

# 创建目录
os.makedirs(RESULTS_DIR, exist_ok=True)


# Step 1: 数据质量控制 (FastQC)
def run_fastqc(fastq_files, output_dir):
    logging.info("Running FastQC...")
    os.makedirs(output_dir, exist_ok=True)
    for fastq in fastq_files:
        cmd = f"fastqc {fastq} -o {output_dir}"
        subprocess.run(cmd, shell=True, check=True)
    logging.info("FastQC completed.")


# Step 2: 数据修剪 (Trim Galore)
def run_trim_galore(fastq_files, output_dir):
    logging.info("Running Trim Galore...")
    os.makedirs(output_dir, exist_ok=True)
    for fastq in fastq_files:
        cmd = f"trim_galore {fastq} --output_dir {output_dir}"
        subprocess.run(cmd, shell=True, check=True)
    logging.info("Trim Galore completed.")


# Step 3: 比对到参考基因组 (Hisat2)
def run_hisat2(fastq_files, genome_index, output_dir):
    logging.info("Running Hisat2...")
    os.makedirs(output_dir, exist_ok=True)
    for fastq in fastq_files:
        sam_file = os.path.join(output_dir, os.path.basename(fastq).replace(".fastq", ".sam"))
        cmd = f"hisat2 -x {genome_index} -U {fastq} -S {sam_file}"
        subprocess.run(cmd, shell=True, check=True)
    logging.info("Hisat2 completed.")


# Step 4: 基因计数 (FeatureCounts)
def run_featurecounts(sam_files, gtf_file, output_file):
    logging.info("Running FeatureCounts...")
    sam_files_str = " ".join(sam_files)
    cmd = f"featureCounts -a {gtf_file} -o {output_file} {sam_files_str}"
    subprocess.run(cmd, shell=True, check=True)
    logging.info("FeatureCounts completed.")


# Step 5: 差异表达分析 (用 DESeq2 in R)
def run_deseq2(count_matrix_file, condition_file):
    logging.info("Running DESeq2 in R for differential expression analysis...")
    r_code = f"""
    library("DESeq2")
    count_data <- read.csv("{count_matrix_file}", row.names=1)
    coldata <- read.csv("{condition_file}", row.names=1)

    dds <- DESeqDataSetFromMatrix(countData = count_data, colData = coldata, design = ~ condition)
    dds <- DESeq(dds)
    res <- results(dds)

    write.csv(as.data.frame(res), file="{RESULTS_DIR}/deseq2_results.csv")
    """

    r(r_code)  # 在 Python 中执行 R 代码
    logging.info("DESeq2 analysis completed and results saved.")


# Step 6: 结果可视化 (差异表达基因的火山图)
def plot_volcano(deseq2_results_file):
    logging.info("Plotting volcano plot...")
    df = pd.read_csv(deseq2_results_file)

    # 火山图
    plt.figure(figsize=(10, 6))
    sns.scatterplot(x='log2FoldChange', y='-log10(padj)', data=df)
    plt.xlabel('Log2 Fold Change')
    plt.ylabel('-Log10 Adjusted P-value')
    plt.title('Volcano Plot of Differentially Expressed Genes')
    plt.savefig(os.path.join(RESULTS_DIR, 'volcano_plot.png'))
    plt.show()
    logging.info("Volcano plot saved.")


# 主函数
if __name__ == "__main__":
    # 获取所有 FASTQ 文件
    fastq_files = [os.path.join(FASTQ_DIR, f) for f in os.listdir(FASTQ_DIR) if f.endswith(".fastq")]

    # 1. FastQC - 质量控制
    fastqc_output_dir = os.path.join(RESULTS_DIR, "fastqc")
    run_fastqc(fastq_files, fastqc_output_dir)

    # 2. Trim Galore - 数据修剪
    trimmed_output_dir = os.path.join(RESULTS_DIR, "trimmed")
    run_trim_galore(fastq_files, trimmed_output_dir)

    # 获取修剪后的 FASTQ 文件
    trimmed_fastq_files = [os.path.join(trimmed_output_dir, f) for f in os.listdir(trimmed_output_dir) if
                           f.endswith(".fq")]

    # 3. Hisat2 - 比对到参考基因组
    hisat2_output_dir = os.path.join(RESULTS_DIR, "hisat2")
    run_hisat2(trimmed_fastq_files, GENOME_INDEX, hisat2_output_dir)

    # 4. FeatureCounts - 基因表达计数
    sam_files = [os.path.join(hisat2_output_dir, f) for f in os.listdir(hisat2_output_dir) if f.endswith(".sam")]
    count_matrix_file = os.path.join(RESULTS_DIR, "counts.txt")
    run_featurecounts(sam_files, GTF_FILE, count_matrix_file)

    # 5. DESeq2 - 差异表达分析
    condition_file = "./conditions.csv"  # 条件文件，包含样本名称和条件信息
    run_deseq2(count_matrix_file, condition_file)

    # 6. Volcano Plot - 差异表达分析火山图
    deseq2_results_file = os.path.join(RESULTS_DIR, "deseq2_results.csv")
    plot_volcano(deseq2_results_file)
