# SCENIC 简单使用指南

## 文件位置
所有SCENIC数据库文件都放在工作目录下的 `scenic_extdata/` 文件夹中。

## 下载文件
```bash
# 下载所有SCENIC数据库文件
./inst/scripts/download_scenic_files.sh
```

## 检查文件
```bash
# 检查文件是否完整
./inst/scripts/step4_scenic.sh --check-only
```

## 运行SCENIC分析
```bash
# 运行完整的SCENIC分析
./inst/scripts/step4_scenic.sh
```

## 文件说明
`scenic_extdata/` 文件夹包含以下文件：
1. `hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather` - SCENIC rankings数据库
2. `motifs-v9-nr.hgnc-m0.001-o0.0.tbl` - SCENIC motifs数据库  
3. `hs_hgnc_tfs.txt` - 人类转录因子列表

## 手动下载（如果自动下载失败）
如果自动下载失败，可以手动从以下链接下载：
1. https://zenodo.org/records/18901480/files/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather
2. https://zenodo.org/records/18901480/files/motifs-v9-nr.hgnc-m0.001-o0.0.tbl
3. https://zenodo.org/records/18901480/files/hs_hgnc_tfs.txt

下载后放入 `scenic_extdata/` 文件夹即可。