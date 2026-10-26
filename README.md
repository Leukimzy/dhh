# 人类泛基因组论文三幅图的 R 复现

这是报告配套的代码与数据。选取 Liao 等人在 *Nature* 2023 年发表的《A draft human pangenome reference》中的 Fig.3h、Fig.3b、Fig.6c。三张图都从论文补充表或公开发布的数据文件整理而来，没有从图片上估读数值。

## 如何运行

安装 R 4.2 或更新版本，在 R 中执行：

```r
install.packages(c("ggplot2", "hexbin"))
```

然后在本仓库根目录运行：

```bash
Rscript reproduce.R
```

运行后会在 `results/` 中生成三张 PNG，并在终端打印 Fig.6c 三组位点数和相关系数。Fig.3h、Fig.3b 只用 R 自带函数；Fig.6c 使用 `ggplot2` 和 `hexbin`。

## 文件说明

| 路径 | 内容 |
| --- | --- |
| `reproduce.R` | 三张图的完整作图代码，按 Fig.3h、Fig.3b、Fig.6c 排列 |
| `data/fig3h_five_classes.csv` | 补充表 12 计数换算出的五类比例 |
| `data/fig3h_six_class_counts.csv` | 补充表 12 的六类原始计数，供核对 |
| `data/fig3b_counts_wide.csv` | 两个 Year1 v1.0 VCF 的 18 项计数整理成六行 |
| `data/fig3b_vcf_counts.csv` | VCF 逐条统计的原始分类结果 |
| `data/fig6c_autosomal_sites.csv` | Zenodo 两表匹配后用于画图的逐位点频率 |
| `figures/` | 报告采用的三张复现图，运行代码可重新生成 |
| `raw_workflow/` | 从原始大文件整理 CSV 的辅助脚本，供检查处理步骤 |
| `数据来源与可访问性.md` | 下载入口、文件大小、核对值和获取困难 |

## 三张图分别回答什么问题

**Fig.3h** 比较四种变异检测方法得到的小变异，在常见、低频、罕见等类别中的比例。图中每一横条加起来是 100%。这里把补充表里占比很小的 Low confidence 合并进 Missing，以匹配论文图中的五类图例；原始六类计数仍保留在 `data/` 中。

**Fig.3b** 比较 MC 和 PGGB 在 SNP、MNP、INDEL 三种小变异类型上的位点数，并把每根柱子按等位基因数分层。这里统计的是 HPRC 公开的 Year1 v1.0 **后处理 VCF**；论文 Fig.3b 的图分解流程与它并不完全相同，所以柱高不能当作论文原图数字的逐项精确复算。

**Fig.6c** 对照结构变异在组装面板和 PanGenie 推断中的等位基因频率。每个六边形聚合一片区域中的真实位点，颜色越深表示该区域位点越多。只保留常染色体后，删除、插入、其他结构变异的相关系数分别约为 0.965、0.930、0.899。

## 数据取得与核对

本仓库直接附带作图所需的整理后 CSV，因此下载仓库即可重画。原始 VCF 和 Zenodo 压缩文件较大，不放入仓库。下载地址、MD5 和从原始文件处理成 CSV 的步骤见[数据来源与可访问性](数据来源与可访问性.md)。

论文：Liao W-W et al. *A draft human pangenome reference*. Nature 617, 312–324 (2023). DOI: [10.1038/s41586-023-05896-x](https://doi.org/10.1038/s41586-023-05896-x)。
