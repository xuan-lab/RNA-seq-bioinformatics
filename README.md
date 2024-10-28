# RNA-seq 数据分析流程

本项目提供了一个基于 Python 实现的 RNA-seq 数据分析流程，涵盖从原始数据的质量控制到差异表达分析的所有关键步骤。该流程集成了多个主流的生物信息学工具，如 **FastQC**、**Trim Galore**、**Hisat2** 和 **FeatureCounts**，并使用 R 语言中的 **DESeq2** 进行差异表达分析。

## 目录

- [环境要求](#环境要求)
- [工具安装](#工具安装)
- [Python 库安装](#Python-库安装)
- [项目结构](#项目结构)
- [使用方法](#使用方法)
- [文件说明](#文件说明)
- [结果展示](#结果展示)
- [参考文献](#参考文献)

## 环境要求

要在本项目中运行 RNA-seq 数据分析流程，请确保在系统上满足以下要求：

- **操作系统**：Linux 或 macOS（Windows 用户建议使用 Windows 子系统 Linux (WSL)）
- **Python 版本**：3.6 或更高
- **R 版本**：3.5 或更高（用于差异表达分析所需的 DESeq2）
- **硬件要求**：充足的存储和内存资源，以处理较大规模的 RNA-seq 数据集

## 工具安装

请确保安装以下生物信息学工具，所有命令应通过命令行执行：

1. **FastQC**：用于评估原始 RNA-seq 数据的质量。
    ```bash
    sudo apt-get install fastqc
    ```

2. **Trim Galore**：用于去除低质量碱基和接头序列。
    ```bash
    sudo apt-get install trim-galore
    ```

3. **Hisat2**：用于将 RNA-seq 读段高效比对到参考基因组。
    ```bash
    sudo apt-get install hisat2
    ```

4. **FeatureCounts**（Subread 包的一部分）：用于从比对结果生成基因表达计数。
    ```bash
    sudo apt-get install subread
    ```

5. **DESeq2**（R 包）：用于差异基因表达分析。
    在 R 控制台中运行以下命令安装 DESeq2：
    ```R
    install.packages("BiocManager")
    BiocManager::install("DESeq2")
    ```

## Python 库安装

本项目依赖几个 Python 库。可以通过以下命令安装这些库：

```bash
pip install pandas matplotlib seaborn rpy2
```

- `pandas`：用于处理表格数据（如基因计数矩阵）。
- `matplotlib` & `seaborn`：用于生成图表和可视化结果。
- `rpy2`：用于在 Python 环境下调用 R 语言功能，包括 DESeq2。

## 项目结构

该项目的目录结构如下所示：

```text
.
├── data/
│   ├── fastq/                # 存放原始 FASTQ 数据
├── genome/
│   ├── hisat2_index/         # 存放 Hisat2 索引文件
│   ├── annotations.gtf       # 基因组注释文件（GTF 格式）
├── results/                  # 存储分析输出结果
│   ├── fastqc/               # FastQC 质量控制结果
│   ├── trimmed/              # 修剪后的 FASTQ 文件
│   ├── hisat2/               # 比对结果（SAM/BAM 文件）
│   ├── counts.txt            # FeatureCounts 生成的计数矩阵
│   ├── deseq2_results.csv    # 差异表达分析结果
│   ├── volcano_plot.png      # 差异表达基因的火山图
├── conditions.csv            # 样本条件文件，定义差异表达分析中的实验组别
├── rnaseq_analysis.py        # 实现 RNA-seq 全流程分析的 Python 脚本
└── README.md                 # 项目说明文档（即本文件）
```

## 使用方法

### 1. 准备数据：

- 将 RNA-seq 的原始 FASTQ 文件放置在 `./data/fastq/` 目录下。
- 将 Hisat2 的参考基因组索引文件放置在 `./genome/hisat2_index/` 目录下。
- 确保基因组注释文件（`annotations.gtf`）位于 `./genome/` 目录中。

### 2. 定义样本条件：

在项目根目录下创建文件 `conditions.csv`，用于定义各个样本对应的实验条件。文件格式如下：

```csv
sample,condition
sample1,treatment
sample2,control
sample3,treatment
sample4,control
```

- **sample**：样本名称，需与 FASTQ 文件名对应。
- **condition**：样本所属的实验条件（如处理组或对照组）。

### 3. 执行 RNA-seq 分析脚本：

通过命令行在项目根目录下运行以下命令：

```bash
python rnaseq_analysis.py
```

该脚本将依次执行以下步骤：

- **数据质量评估**：使用 FastQC 对原始数据进行初步质量控制。
- **数据修剪**：使用 Trim Galore 去除接头和低质量序列。
- **比对到参考基因组**：使用 Hisat2 将修剪后的读段比对到参考基因组。
- **读段计数**：使用 FeatureCounts 生成每个基因的表达计数矩阵。
- **差异表达分析**：使用 DESeq2 进行差异基因表达分析。
- **可视化**：生成火山图，用于展示差异表达基因的显著性和倍数变化。

### 4. 查看结果：

所有分析结果将保存在 `./results/` 目录中。主要结果包括：
- **差异表达分析结果**：`deseq2_results.csv` 文件包含差异基因表达分析的详细结果。
- **火山图**：`volcano_plot.png` 文件展示了差异表达基因的可视化结果。

## 文件说明

- **`rnaseq_analysis.py`**：主分析脚本，涵盖从数据质量控制到差异表达分析的整个流程。
- **`conditions.csv`**：定义样本与实验条件关系的文件，用于差异表达分析。
- **`counts.txt`**：由 FeatureCounts 生成的基因表达计数矩阵。
- **`deseq2_results.csv`**：DESeq2 生成的差异表达分析结果，包括基因的倍数变化、p 值和调整后的 p 值。
- **`volcano_plot.png`**：火山图，展示差异表达基因的显著性与倍数变化。

## 结果展示

### 1. 火山图 (Volcano Plot)：

火山图展示了所有差异表达基因的显著性（p 值）与倍数变化（log2FoldChange）的关系。下图为火山图的示例：

![Volcano Plot](./results/volcano_plot.png)

### 2. DESeq2 分析结果：

DESeq2 生成的差异表达结果保存在 `deseq2_results.csv` 文件中，主要包括以下字段：

- **log2FoldChange**：基因表达的对数倍数变化，表示在实验条件间的差异表达。
- **pvalue**：差异表达的显著性 p 值。
- **padj**：使用 Benjamini-Hochberg 方法校正的 p 值，用于控制多重检验的假阳性率。

## 作者
- 姓名：Xiexuan
- 邮箱：xiexuan@njfu.edu.cn
- 单位：南京林业大学

## 参考文献

- **FastQC**：Andrews, S. (2010). FastQC: A Quality Control Tool for High Throughput Sequence Data.
- **Trim Galore**：Krueger, F. (2012). Trim Galore: A Wrapper Tool Around Cutadapt and FastQC to Consistently Apply Quality and Adapter Trimming to FastQ Files.
- **Hisat2**：Kim, D., et al. (2015). HISAT: A fast spliced aligner with low memory requirements. *Nature Methods*, 12(4), 357-360.
- **FeatureCounts**：Liao, Y., et al. (2014). featureCounts: an efficient general-purpose program for assigning sequence reads to genomic features. *Bioinformatics*, 30(7), 923-930.
- **DESeq2**：Love, M.I., et al. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology*, 15(12), 550.
