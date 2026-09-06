# 真实数据操作指南

## 1. 保留两个相互独立的空间

`data/demo` / `results/demo` 仅存放合成测试；`data/processed` / `results/real`
用于真实研究。切换到 `config/production.yaml` 后，程序会拒绝带 SYNTHETIC
标记的数据。不要改掉标记冒充真实数据；完成的软件测试也不等于生物学验证。

## 2. 获取公开研究数据并冻结来源

```bash
python scripts/fetch_public.py --part metadata --extract
python scripts/fetch_public.py --part spatial --extract
python scripts/fetch_public.py --part singlecell --extract
```

空间文件来自 Zenodo 4739739，下载器校验发布方 MD5；同时保存 URL、字节数和
本地 SHA-256。单细胞合并包来自 GEO，当前脚本记录本地 SHA-256，但不能把它
写成已通过“发布方 SHA-256”。重新下载无发布方校验值的已有文件会停止，
需要先核对本地 provenance，或使用新的目录保存新版本。

原始 scRNA 读段受 EGA 控制。公开过滤矩阵可以开展下游分析，不能据此声称
自己重新完成了原始 FASTQ 比对、empty-droplet 判定或完整环境 RNA 校正。

## 3. 人工审核是有依据的数据整理步骤

下载后先查看文件列表、矩阵方向、条码和特征文件，再依照 `input_contract.md`
写入规范化输入，不凭文件名猜测分组。保留原始字段和映射表。

对于空间元数据，已实际核实的字段为：行索引 spot barcode、`nCount_RNA`、
`nFeature_RNA`、`subtype`、`patientid`、`Classification`。至少要做以下审核：

- patientid → patient_id；每个空间切片单独 section_id，spot ID 加前缀避免碰撞。
- subtype=TNBC 的来源定义保留；subtype=ER 的记录需回到临床表核查 HER2/PR。
  不能直接把 ER 改成 HRpos_HER2neg 让程序顺利运行。
- `Classification` 复制到 pathology_original。`Stroma` 可作为待核对间质标签；
  `Invasive cancer + stroma + lymphocytes` 必须保留 mixed 状态或用额外病理掩膜，
  不自动变成纯 tumor。Artefact、坏死、未确定区域按预设规则排除并记录。
- 从空间目录的 positions 文件读取 in_tissue 与全分辨率坐标。
  旧版无表头顺序为 barcode,in_tissue,array_row,array_col,pxl_row_in_fullres,
  pxl_col_in_fullres；映射 x=col，y=row。新版按字段名读取。
- 参考基因版本、特征 ID 和 gene_symbol 必须显式统一。若 GEO 文件只有 symbol，
  通过冻结的注释表保留唯一映射；输出映射成功率和丢弃原因。

本仓库提供通用 `prepare_matrix.py` 接入器，并提供空白元数据模板。
它不自动替研究者做不可验证的病理/临床判断。对于同一格式的多样本，可以先
为每个样本产生合规 h5ad，再用 `scripts/merge_assays.py` 合并，脚本严格检查基因
一致性，不做隐式多对多基因联接。

## 4. 获取 TCGA-BRCA bulk

```bash
python scripts/gdc_bulk.py discover
```

生成 `data/gdc/candidates.tsv`、查询内容、原始 API 快照和 GDC release 状态。
仅发现 open、Primary Tumor、STAR - Counts 文件，不会主动下载整个队列。
从正式临床来源审核 receptor、年龄、分期和批次，生成 `selected.tsv`：

```text
file_id    sample_id    patient_id    condition    batch
```

实际文件使用 tab 分隔；每行须填真实已审核标识，不复制示意标识。
每位患者只选一个肿瘤 aliquot，事先冻结重复样本的处理规则并保留排除表。
不能只凭 sample barcode 推断临床 HER2/ER，不使用表达结果挑选“最有利”样本。

```bash
python scripts/gdc_bulk.py fetch --selection data/gdc/selected.tsv
```

下载器验证选择的 file→patient→sample 对应，核对 GDC MD5，提取原始 unstranded
计数和基因注释。若真实分析需要调节其他协变量，在 selected.tsv 中加入完整
字段，并在观察结果之前改好 design。

## 5. 建立独立的空间参考与 cell2location 环境

参考患者应与需要验证的目标患者不重叠。公开同研究数据可用于探索性匹配，
但验证阶段要按患者留出；参考与目标共享患者时，当前适配器会停止。
可以逐个患者建立留出折，分别拟合后按 spot ID 合并**后验均值**。

在 Linux / GPU 主机单独创建 cell2location 官方推荐环境，按使用时的官方
安装说明记录精确版本与 CUDA/torch 驱动。本仓库 CPU 环境不包含 torch/
cell2location，避免把尚未解算和运行的 GPU 依赖称为已锁定环境。

```bash
# 以下程序在已建立并验证的 cell2location 环境中执行
python scripts/cell2location_fit.py \
  --reference data/curated/reference_donor_disjoint.h5ad \
  --spatial data/processed/spatial.h5ad \
  --cells-per-spot 10 \
  --output-dir results/cell2location \
  --accelerator gpu
```

上例的 10 是命令示例，必须用对应 H&E 的核计数估计替换。
检查训练历史、后验预测、病理一致性、未覆盖细胞类型、参考偏倚，并围绕
cells-per-spot 与 detection_alpha 做预设敏感性分析。默认长训练不等于已经收敛。
主流程的 `abundance_file` 指向 `abundance_means.csv`；保留 q05/q95 和模型文件。
适配器仅完成来源/API 审查，未在本次交付中执行 GPU 拟合。

## 6. 运行与审阅

```bash
bcflow run --config config/production.yaml
# 或分步运行 singlecell / pseudobulk / bulk / spatial / integrate / figures
```

程序先报错也可能是正确行为：当前真实数据达不到独立患者阈值、缺少受审核
标签、临床分组不足或批次混杂时，不应该生成显著性结论。修订研究设计或
补充数据后再运行，不把执行阈值降到 1 来获得期刊图。

审阅次序：样本纳排→单细胞/空间 QC→注释→伪 bulk 有效患者数→模型设计→
差异表达→参考映射诊断→空间病理→患者区间→跨模态一致性→独立实验验证。
所有真实数据图片仍标为探索性，只有完成这些检查才能用于正式稿件。

## 7. 资源规划（估计，不是实测基准）

CPU 演示一般只需数 GB 内存；真实约 10 万细胞及空间矩阵建议 Linux 64–128 GB
内存，单独准备 GPU 训练节点与数十 GB 数据空间。资源随基因数、细胞数、
空间样本数和训练参数变化，先对一个样本做资源试跑。GDC 全队列下载耗时
取决于文件清单大小和网络。不在 Git 仓库存储大规模表达矩阵或患者身份信息。

