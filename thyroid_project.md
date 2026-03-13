import os
import tarfile

# 设置路径
data_path = r"C:\Users\13918\Desktop\onek1k\bulk"
tar_file = os.path.join(data_path, "GSE248205_Processed_data.tar.gz")

# 解压文件
with tarfile.open(tar_file, "r:gz") as tar:
    tar.extractall(path=data_path)
    print("解压完成！")

# 查看解压后的文件结构
extract_dir = os.path.join(data_path)
for root, dirs, files in os.walk(extract_dir):
    level = root.replace(extract_dir, '').count(os.sep)
    indent = ' ' * 2 * level
    print(f'{indent}{os.path.basename(root)}/')
    subindent = ' ' * 2 * (level + 1)
    for file in files:
        file_path = os.path.join(root, file)
        file_size = os.path.getsize(file_path) / (1024*1024)  # MB
        print(f'{subindent}{file} ({file_size:.2f} MB)')
# 更精确地查看解压后新增的文件
import glob

# 查看所有解压出来的文件
all_files = []
for root, dirs, files in os.walk(data_path):
    for f in files:
        if f != "GSE248205_Processed_data.tar.gz":  # 排除原始压缩包
            full_path = os.path.join(root, f)
            all_files.append(full_path)
            print(full_path)

print(f"\n共找到 {len(all_files)} 个文件")
import os
import shutil
import json
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# ============================================
# 第一步：把每个样本的文件重组为标准 Visium 格式
# ============================================

base_dir = r"C:\Users\13918\Desktop\onek1k\bulk\Processed files"
output_base = r"C:\Users\13918\Desktop\onek1k\bulk\visium_formatted"

sample_names = ["C1", "C2", "GD1", "GD2", "GD3", "HT1", "HT2", "HT3"]

for sample in sample_names:
    src_dir = os.path.join(base_dir, sample)
    
    # 创建标准 Space Ranger 输出目录结构
    # outs/
    #   filtered_feature_bc_matrix/
    #     barcodes.tsv.gz
    #     features.tsv.gz
    #     matrix.mtx.gz
    #   spatial/
    #     scalefactors_json.json
    #     tissue_hires_image.png
    #     tissue_lowres_image.png
    #     tissue_positions_list.csv
    
    matrix_dir = os.path.join(output_base, sample, "filtered_feature_bc_matrix")
    spatial_dir = os.path.join(output_base, sample, "spatial")
    os.makedirs(matrix_dir, exist_ok=True)
    os.makedirs(spatial_dir, exist_ok=True)
    
    # 复制 matrix 文件（去掉样本名前缀）
    shutil.copy(os.path.join(src_dir, f"{sample}_barcodes.tsv.gz"),
                os.path.join(matrix_dir, "barcodes.tsv.gz"))
    shutil.copy(os.path.join(src_dir, f"{sample}_features.tsv.gz"),
                os.path.join(matrix_dir, "features.tsv.gz"))
    shutil.copy(os.path.join(src_dir, f"{sample}_matrix.mtx.gz"),
                os.path.join(matrix_dir, "matrix.mtx.gz"))
    
    # 复制 spatial 文件（去掉样本名前缀）
    shutil.copy(os.path.join(src_dir, f"{sample}_scalefactors_json.json"),
                os.path.join(spatial_dir, "scalefactors_json.json"))
    shutil.copy(os.path.join(src_dir, f"{sample}_tissue_hires_image.png"),
                os.path.join(spatial_dir, "tissue_hires_image.png"))
    shutil.copy(os.path.join(src_dir, f"{sample}_tissue_lowres_image.png"),
                os.path.join(spatial_dir, "tissue_lowres_image.png"))
    shutil.copy(os.path.join(src_dir, f"{sample}_tissue_positions_list.csv"),
                os.path.join(spatial_dir, "tissue_positions_list.csv"))
    
    print(f"✓ {sample} 文件重组完成")

print("\n所有样本重组完成！")
import scanpy as sc
import pandas as pd
import numpy as np
import json
from PIL import Image
import os

output_base = r"C:\Users\13918\Desktop\onek1k\bulk\visium_formatted"
sample_names = ["C1", "C2", "GD1", "GD2", "GD3", "HT1", "HT2", "HT3"]

condition_map = {
    "C1": "Control", "C2": "Control",
    "GD1": "GD", "GD2": "GD", "GD3": "GD",
    "HT1": "HT", "HT2": "HT", "HT3": "HT"
}

adata_dict = {}

for sample in sample_names:
    sample_dir = os.path.join(output_base, sample)
    matrix_dir = os.path.join(sample_dir, "filtered_feature_bc_matrix")
    spatial_dir = os.path.join(sample_dir, "spatial")
    
    print(f"\n{'='*50}")
    print(f"正在读取: {sample} ({condition_map[sample]})")
    print(f"{'='*50}")
    
    try:
        # 1. 读取表达矩阵
        adata = sc.read_10x_mtx(matrix_dir)
        adata.var_names_make_unique()
        
        # 2. 读取空间位置信息
        positions_file = os.path.join(spatial_dir, "tissue_positions_list.csv")
        
        # 先看一下文件头部，判断有没有列名
        with open(positions_file, 'r') as f:
            first_line = f.readline().strip()
        print(f"  positions文件首行: {first_line[:80]}")
        
        # 判断是否有header
        if 'barcode' in first_line.lower() or 'in_tissue' in first_line.lower():
            positions = pd.read_csv(positions_file, index_col=0)
        else:
            positions = pd.read_csv(positions_file, header=None, index_col=0,
                                     names=['in_tissue', 'array_row', 'array_col', 
                                            'pxl_row_in_fullres', 'pxl_col_in_fullres'])
        
        print(f"  positions形状: {positions.shape}")
        print(f"  positions列: {list(positions.columns)}")
        
        # 3. 只保留在表达矩阵中的barcode
        common_barcodes = adata.obs_names.intersection(positions.index)
        print(f"  表达矩阵barcodes: {adata.n_obs}")
        print(f"  positions barcodes: {positions.shape[0]}")
        print(f"  重叠barcodes: {len(common_barcodes)}")
        
        adata = adata[common_barcodes, :].copy()
        positions = positions.loc[common_barcodes]
        
        # 4. 添加空间坐标到 obsm
        adata.obsm['spatial'] = positions[['pxl_row_in_fullres', 'pxl_col_in_fullres']].values
        
        # 5. 读取 scale factors
        with open(os.path.join(spatial_dir, "scalefactors_json.json"), 'r') as f:
            scalefactors = json.load(f)
        print(f"  scalefactors: {scalefactors}")
        
        # 6. 读取组织图像
        hires_img = np.array(Image.open(os.path.join(spatial_dir, "tissue_hires_image.png")))
        lowres_img = np.array(Image.open(os.path.join(spatial_dir, "tissue_lowres_image.png")))
        
        # 7. 构建 uns['spatial'] 结构（scanpy标准格式）
        adata.uns['spatial'] = {
            sample: {
                'images': {
                    'hires': hires_img / 255.0 if hires_img.max() > 1 else hires_img,
                    'lowres': lowres_img / 255.0 if lowres_img.max() > 1 else lowres_img,
                },
                'scalefactors': scalefactors,
                'metadata': {}
            }
        }
        
        # 8. 添加 in_tissue 信息
        adata.obs['in_tissue'] = positions['in_tissue'].values
        adata.obs['array_row'] = positions['array_row'].values
        adata.obs['array_col'] = positions['array_col'].values
        
        # 9. 添加样本和条件信息
        adata.obs['sample'] = sample
        adata.obs['condition'] = condition_map[sample]
        
        # 设置 library_id
        adata.uns['spatial'][sample]['metadata']['library_id'] = sample
        
        print(f"  Spots数量: {adata.n_obs}")
        print(f"  基因数量: {adata.n_vars}")
        print(f"  ✓ 读取成功！")
        
        adata_dict[sample] = adata
        
    except Exception as e:
        import traceback
        print(f"  ✗ 读取失败: {e}")
        traceback.print_exc()

print(f"\n\n{'='*50}")
print(f"成功读取 {len(adata_dict)}/{len(sample_names)} 个样本")
for name, ad in adata_dict.items():
    print(f"  {name}: {ad.n_obs} spots × {ad.n_vars} genes [{condition_map[name]}]")
import matplotlib.pyplot as plt

# ============================================
# 3.1 先看一下每个样本的组织图像和空间分布
# ============================================

fig, axes = plt.subplots(2, 4, figsize=(24, 12))
axes = axes.flatten()

for idx, (sample, adata) in enumerate(adata_dict.items()):
    ax = axes[idx]
    
    # 获取低分辨率图像
    img = adata.uns['spatial'][sample]['images']['lowres']
    scale = adata.uns['spatial'][sample]['scalefactors']['tissue_lowres_scalef']
    
    # 显示组织图像
    ax.imshow(img)
    
    # 叠加spots（只显示in_tissue=1的）
    in_tissue_mask = adata.obs['in_tissue'] == 1
    coords = adata.obsm['spatial'][in_tissue_mask]
    ax.scatter(coords[:, 1] * scale, coords[:, 0] * scale, 
               s=3, alpha=0.5, c='red')
    
    in_tissue_count = in_tissue_mask.sum()
    ax.set_title(f"{sample} ({adata.obs['condition'].iloc[0]})\n"
                 f"Total: {adata.n_obs}, In tissue: {in_tissue_count}", fontsize=12)
    ax.axis('off')

plt.suptitle("GSE248205 - 所有样本空间分布", fontsize=16, fontweight='bold')
plt.tight_layout()
plt.savefig(r"C:\Users\13918\Desktop\onek1k\bulk\all_samples_spatial.png", dpi=150, bbox_inches='tight')
plt.show()

print("图片已保存！")
# ============================================
# 3.2 只保留 in_tissue=1 的 spots，然后合并
# ============================================

adata_filtered = {}

for sample, adata in adata_dict.items():
    # 过滤只保留组织上的spots
    mask = adata.obs['in_tissue'] == 1
    adata_f = adata[mask].copy()
    adata_filtered[sample] = adata_f
    print(f"{sample}: {adata.n_obs} -> {adata_f.n_obs} spots (过滤in_tissue)")

# 合并所有样本
import anndata as ad

# 合并前需要统一 uns 结构
adata_list = []
spatial_uns = {}

for sample, adata in adata_filtered.items():
    # 收集所有样本的空间信息
    spatial_uns[sample] = adata.uns['spatial'][sample]
    # 清理 uns 避免合并冲突
    adata.uns = {}
    adata_list.append(adata)

# 合并
adata_combined = ad.concat(adata_list, join='inner', merge='same')
adata_combined.obs_names_make_unique()

# 重新添加空间信息
adata_combined.uns['spatial'] = spatial_uns

print(f"\n{'='*50}")
print(f"合并后数据:")
print(f"  总spots: {adata_combined.n_obs}")
print(f"  总基因: {adata_combined.n_vars}")
print(f"\n各组样本数:")
print(adata_combined.obs['condition'].value_counts())
print(f"\n各样本spots数:")
print(adata_combined.obs['sample'].value_counts())


# ============================================
# 4.1 基本质控 QC
# ============================================

# 计算QC指标
adata_combined.var['mt'] = adata_combined.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata_combined, qc_vars=['mt'], inplace=True)

# 可视化QC
fig, axes = plt.subplots(1, 4, figsize=(20, 4))

sc.pl.violin(adata_combined, keys='total_counts', groupby='sample', 
             rotation=45, ax=axes[0], show=False)
axes[0].set_title('Total Counts per Spot')

sc.pl.violin(adata_combined, keys='n_genes_by_counts', groupby='sample', 
             rotation=45, ax=axes[1], show=False)
axes[1].set_title('Genes per Spot')

sc.pl.violin(adata_combined, keys='pct_counts_mt', groupby='sample', 
             rotation=45, ax=axes[2], show=False)
axes[2].set_title('MT% per Spot')

# 各条件的总counts分布
sc.pl.violin(adata_combined, keys='total_counts', groupby='condition', 
             rotation=45, ax=axes[3], show=False)
axes[3].set_title('Total Counts by Condition')

plt.tight_layout()
plt.savefig(r"C:\Users\13918\Desktop\onek1k\bulk\QC_violin.png", dpi=150, bbox_inches='tight')
plt.show()

# 打印QC统计
print("\n各样本QC统计:")
qc_stats = adata_combined.obs.groupby('sample').agg({
    'total_counts': ['median', 'mean'],
    'n_genes_by_counts': ['median', 'mean'],
    'pct_counts_mt': ['median', 'mean']
}).round(1)
print(qc_stats)


# ============================================
# 4.2 检查目标基因是否存在
# ============================================

target_genes = [
    'POU5F1', 'HCG11', 'HMGN4', 'ABHD16A', 'HIST1H3H',
    'ZKSCAN4', 'HLA-F', 'HLA-DQB2', 'ITPR3', 'HLA-DQB1', 'FAM134B'
]

print("目标基因检查:")
print(f"{'基因名':<15} {'是否存在':<10} {'检测到的spots数':<15} {'平均表达量'}")
print("-" * 60)

found_genes = []
missing_genes = []

for gene in target_genes:
    if gene in adata_combined.var_names:
        # 计算该基因在多少spots中有表达
        gene_expr = adata_combined[:, gene].X
        if hasattr(gene_expr, 'toarray'):
            gene_expr = gene_expr.toarray().flatten()
        else:
            gene_expr = np.array(gene_expr).flatten()
        
        n_expressed = (gene_expr > 0).sum()
        mean_expr = gene_expr.mean()
        pct = n_expressed / adata_combined.n_obs * 100
        
        print(f"{gene:<15} {'✓':<10} {n_expressed} ({pct:.1f}%)      {mean_expr:.3f}")
        found_genes.append(gene)
    else:
        # 尝试模糊匹配
        fuzzy = [v for v in adata_combined.var_names if gene in v]
        if fuzzy:
            print(f"{gene:<15} {'✗':<10} 可能匹配: {fuzzy[:5]}")
        else:
            print(f"{gene:<15} {'✗':<10} 未找到")
        missing_genes.append(gene)

print(f"\n找到 {len(found_genes)}/{len(target_genes)} 个目标基因")
print(f"找到的基因: {found_genes}")
if missing_genes:
    print(f"缺失的基因: {missing_genes}")

# ============================================
# 4.2.1 先查找 FAM134B 的新名称
# ============================================

# FAM134B 已被HUGO更名为 RETREG1
retreg_match = [v for v in adata_combined.var_names if 'RETREG' in v or 'FAM134' in v]
print(f"RETREG/FAM134相关基因: {retreg_match}")

# 更新目标基因列表
found_genes = ['POU5F1', 'HCG11', 'HMGN4', 'ABHD16A', 'HIST1H3H',
               'ZKSCAN4', 'HLA-F', 'HLA-DQB2', 'ITPR3', 'HLA-DQB1']

if 'RETREG1' in adata_combined.var_names:
    found_genes.append('RETREG1')
    print("✓ FAM134B 在数据中名为 RETREG1，已添加到目标基因列表")

# ============================================
# 4.3 预处理（修正版）
# ============================================

# 保存原始counts
adata_combined.layers['counts'] = adata_combined.X.copy()

# 过滤
sc.pp.filter_genes(adata_combined, min_cells=10)
print(f"过滤后基因数: {adata_combined.n_vars}")

# 标准化
sc.pp.normalize_total(adata_combined, target_sum=1e4)
sc.pp.log1p(adata_combined)

# 保存标准化后的数据
adata_combined.layers['normalized'] = adata_combined.X.copy()

# 高变基因（用 seurat 而非 seurat_v3，更稳定）
sc.pp.highly_variable_genes(adata_combined, flavor='seurat', 
                             min_mean=0.0125, max_mean=3, min_disp=0.5)

n_hvg = adata_combined.var['highly_variable'].sum()
print(f"高变基因数: {n_hvg}")

# 检查目标基因是否为高变基因
print("\n目标基因是否为高变基因:")
for gene in found_genes:
    if gene in adata_combined.var_names:
        is_hvg = adata_combined.var.loc[gene, 'highly_variable']
        print(f"  {gene}: {'HVG ✓' if is_hvg else '非HVG'}")

# PCA
sc.pp.scale(adata_combined, max_value=10)
sc.tl.pca(adata_combined, n_comps=50)

# 邻域图 + UMAP + 聚类
sc.pp.neighbors(adata_combined, n_neighbors=15, n_pcs=30)
sc.tl.umap(adata_combined)
sc.tl.leiden(adata_combined, resolution=0.5)

print(f"\nLeiden聚类数: {adata_combined.obs['leiden'].nunique()}")
print(adata_combined.obs['leiden'].value_counts().sort_index())
print(f"\n预处理完成！adata_combined: {adata_combined.shape}")


# ============================================
# 4.4 UMAP可视化：按样本、条件、聚类
# ============================================

fig, axes = plt.subplots(1, 3, figsize=(21, 6))

sc.pl.umap(adata_combined, color='condition', ax=axes[0], show=False, 
           title='Condition', frameon=False)
sc.pl.umap(adata_combined, color='sample', ax=axes[1], show=False, 
           title='Sample', frameon=False)
sc.pl.umap(adata_combined, color='leiden', ax=axes[2], show=False, 
           title='Leiden Clusters', frameon=False, legend_loc='on data')

plt.tight_layout()
plt.savefig(r"C:\Users\13918\Desktop\onek1k\bulk\UMAP_overview.png", dpi=150, bbox_inches='tight')
plt.show()

# ============================================
# 5.1 UMAP 总览
# ============================================

fig, axes = plt.subplots(1, 3, figsize=(21, 6))

sc.pl.umap(adata_combined, color='condition', ax=axes[0], show=False, 
           title='Condition', frameon=False)
sc.pl.umap(adata_combined, color='sample', ax=axes[1], show=False, 
           title='Sample', frameon=False)
sc.pl.umap(adata_combined, color='leiden', ax=axes[2], show=False, 
           title='Leiden Clusters', frameon=False, legend_loc='on data')

plt.tight_layout()
plt.savefig(r"C:\Users\13918\Desktop\onek1k\bulk\UMAP_overview.png", dpi=600, bbox_inches='tight')
plt.show()

# ============================================
# 5.2 目标基因在UMAP上的表达
# ============================================

# 用normalized层来展示基因表达
found_genes = ['POU5F1', 'HCG11', 'HMGN4', 'ABHD16A', 'HIST1H3H',
               'ZKSCAN4', 'HLA-F', 'HLA-DQB2', 'ITPR3', 'HLA-DQB1', 'RETREG1']

fig, axes = plt.subplots(3, 4, figsize=(24, 16))
axes = axes.flatten()

for i, gene in enumerate(found_genes):
    sc.pl.umap(adata_combined, color=gene, ax=axes[i], show=False,
               layer='normalized', cmap='Reds', frameon=False,
               title=f'{gene}', vmin=0)

# 隐藏多余的子图
axes[11].axis('off')

plt.suptitle('Target Genes Expression on UMAP', fontsize=16, fontweight='bold')
plt.tight_layout()
plt.savefig(r"C:\Users\13918\Desktop\onek1k\bulk\target_genes_UMAP.png", dpi=150, bbox_inches='tight')
plt.show()

# ============================================
# 5.3 目标基因按 Condition 的表达量比较 (Violin Plot)
# ============================================

fig, axes = plt.subplots(3, 4, figsize=(28, 18))
axes = axes.flatten()

for i, gene in enumerate(found_genes):
    sc.pl.violin(adata_combined, keys=gene, groupby='condition',
                 layer='normalized', ax=axes[i], show=False,
                 order=['Control', 'GD', 'HT'],
                 rotation=0)
    axes[i].set_title(gene, fontsize=14, fontweight='bold')

axes[11].axis('off')

plt.suptitle('Target Genes: Control vs GD vs HT', fontsize=16, fontweight='bold')
plt.tight_layout()
plt.savefig(r"C:\Users\13918\Desktop\onek1k\bulk\target_genes_violin.png", dpi=150, bbox_inches='tight')
plt.show()

# ============================================
# 5.4 目标基因 DotPlot：按条件
# ============================================

sc.pl.dotplot(adata_combined, var_names=found_genes, groupby='condition',
              layer='normalized', standard_scale='var',
              title='Target Genes Expression by Condition',
              save='_target_genes_condition.png')

# ============================================
# 6.1 带统计检验星号的 Violin Plot
# ============================================

from scipy import stats
import matplotlib.pyplot as plt
import numpy as np

found_genes = ['POU5F1', 'HCG11', 'HMGN4', 'ABHD16A', 'HIST1H3H',
               'ZKSCAN4', 'HLA-F', 'HLA-DQB2', 'ITPR3', 'HLA-DQB1', 'RETREG1']

# 获取normalized层数据
adata_plot = adata_combined.copy()
adata_plot.X = adata_plot.layers['normalized'].copy()

# 提取各组的index
ctrl_mask = adata_plot.obs['condition'] == 'Control'
gd_mask = adata_plot.obs['condition'] == 'GD'
ht_mask = adata_plot.obs['condition'] == 'HT'

def get_stars(pval):
    if pval < 0.0001:
        return '****'
    elif pval < 0.001:
        return '***'
    elif pval < 0.01:
        return '**'
    elif pval < 0.05:
        return '*'
    else:
        return 'ns'

fig, axes = plt.subplots(3, 4, figsize=(28, 20))
axes = axes.flatten()

# 存储统计结果
stat_results = []

for i, gene in enumerate(found_genes):
    ax = axes[i]
    
    # 提取表达值
    if hasattr(adata_plot[:, gene].X, 'toarray'):
        ctrl_expr = adata_plot[ctrl_mask, gene].X.toarray().flatten()
        gd_expr = adata_plot[gd_mask, gene].X.toarray().flatten()
        ht_expr = adata_plot[ht_mask, gene].X.toarray().flatten()
    else:
        ctrl_expr = np.array(adata_plot[ctrl_mask, gene].X).flatten()
        gd_expr = np.array(adata_plot[gd_mask, gene].X).flatten()
        ht_expr = np.array(adata_plot[ht_mask, gene].X).flatten()
    
    # Wilcoxon rank-sum test
    stat_gd, pval_gd = stats.mannwhitneyu(ctrl_expr, gd_expr, alternative='two-sided')
    stat_ht, pval_ht = stats.mannwhitneyu(ctrl_expr, ht_expr, alternative='two-sided')
    stat_gd_ht, pval_gd_ht = stats.mannwhitneyu(gd_expr, ht_expr, alternative='two-sided')
    
    stat_results.append({
        'gene': gene,
        'ctrl_mean': ctrl_expr.mean(),
        'gd_mean': gd_expr.mean(),
        'ht_mean': ht_expr.mean(),
        'pval_GD_vs_Ctrl': pval_gd,
        'pval_HT_vs_Ctrl': pval_ht,
        'pval_GD_vs_HT': pval_gd_ht,
        'stars_GD_vs_Ctrl': get_stars(pval_gd),
        'stars_HT_vs_Ctrl': get_stars(pval_ht),
        'stars_GD_vs_HT': get_stars(pval_gd_ht)
    })
    
    # 画violin
    data = [ctrl_expr, gd_expr, ht_expr]
    positions = [0, 1, 2]
    
    parts = ax.violinplot(data, positions=positions, showmeans=True, showmedians=False)
    colors = ['#2ecc71', '#e74c3c', '#3498db']
    for pc, color in zip(parts['bodies'], colors):
        pc.set_facecolor(color)
        pc.set_alpha(0.7)
    
    ax.set_xticks(positions)
    ax.set_xticklabels(['Control', 'GD', 'HT'], fontsize=11)
    ax.set_title(gene, fontsize=14, fontweight='bold')
    ax.set_ylabel('Normalized Expression')
    
    # 添加显著性标注
    y_max = max(ctrl_expr.max(), gd_expr.max(), ht_expr.max())
    h = y_max * 0.08  # 标注线的高度
    
    # Control vs GD
    y1 = y_max + h
    ax.plot([0, 0, 1, 1], [y1, y1+h*0.3, y1+h*0.3, y1], 'k-', lw=1)
    ax.text(0.5, y1+h*0.3, get_stars(pval_gd), ha='center', fontsize=11, color='red')
    
    # Control vs HT
    y2 = y1 + h * 1.5
    ax.plot([0, 0, 2, 2], [y2, y2+h*0.3, y2+h*0.3, y2], 'k-', lw=1)
    ax.text(1, y2+h*0.3, get_stars(pval_ht), ha='center', fontsize=11, color='blue')
    
    # GD vs HT
    y3 = y2 + h * 1.5
    ax.plot([1, 1, 2, 2], [y3, y3+h*0.3, y3+h*0.3, y3], 'k-', lw=1)
    ax.text(1.5, y3+h*0.3, get_stars(pval_gd_ht), ha='center', fontsize=11, color='purple')

axes[11].axis('off')

# 添加图例说明
axes[11].text(0.1, 0.7, 'Significance:', fontsize=12, fontweight='bold', transform=axes[11].transAxes)
axes[11].text(0.1, 0.55, '**** p < 0.0001', fontsize=11, transform=axes[11].transAxes)
axes[11].text(0.1, 0.42, '***  p < 0.001', fontsize=11, transform=axes[11].transAxes)
axes[11].text(0.1, 0.29, '**   p < 0.01', fontsize=11, transform=axes[11].transAxes)
axes[11].text(0.1, 0.16, '*    p < 0.05', fontsize=11, transform=axes[11].transAxes)
axes[11].text(0.1, 0.03, 'ns   not significant', fontsize=11, transform=axes[11].transAxes)

plt.suptitle('Target Genes Expression: Control vs GD vs HT\n(Mann-Whitney U test)', 
             fontsize=16, fontweight='bold')
plt.tight_layout()
plt.savefig(r"C:\Users\13918\Desktop\onek1k\bulk\target_genes_violin_significance.png", 
            dpi=200, bbox_inches='tight')
plt.show()

# 打印统计表格
stat_df = pd.DataFrame(stat_results)
print("\n统计检验结果:")
print(stat_df.to_string(index=False))


# ============================================
# 6.2 空间可视化（修正版）
# ============================================

key_genes = ['HLA-DQB1', 'HLA-F', 'ITPR3', 'RETREG1', 'HMGN4', 'ABHD16A']

for gene in key_genes:
    fig, axes = plt.subplots(2, 4, figsize=(28, 14))
    axes = axes.flatten()
    
    for idx, sample in enumerate(['C1','C2','GD1','GD2','GD3','HT1','HT2','HT3']):
        ax = axes[idx]
        
        # 用原始 adata_dict（保留了uns空间信息）
        adata_s = adata_dict[sample].copy()
        
        # 标准化
        sc.pp.normalize_total(adata_s, target_sum=1e4)
        sc.pp.log1p(adata_s)
        
        # 获取表达值
        if gene in adata_s.var_names:
            if hasattr(adata_s[:, gene].X, 'toarray'):
                expr = adata_s[:, gene].X.toarray().flatten()
            else:
                expr = np.array(adata_s[:, gene].X).flatten()
        else:
            expr = np.zeros(adata_s.n_obs)
        
        # 从 adata_combined.uns 获取空间信息
        img = adata_combined.uns['spatial'][sample]['images']['lowres']
        scale = adata_combined.uns['spatial'][sample]['scalefactors']['tissue_lowres_scalef']
        coords = adata_s.obsm['spatial']
        
        # 绘制
        ax.imshow(img)
        scatter = ax.scatter(
            coords[:, 1] * scale, 
            coords[:, 0] * scale,
            c=expr, cmap='Reds', s=12, alpha=0.8,
            vmin=0, vmax=np.percentile(expr, 99) if expr.max() > 0 else 1
        )
        plt.colorbar(scatter, ax=ax, shrink=0.6)
        
        cond = condition_map[sample]
        ax.set_title(f'{sample} ({cond})', fontsize=13, fontweight='bold')
        ax.axis('off')
    
    plt.suptitle(f'Spatial Expression: {gene}', fontsize=18, fontweight='bold')
    plt.tight_layout()
    plt.savefig(rf"C:\Users\13918\Desktop\onek1k\bulk\spatial_{gene}.png", 
                dpi=600, bbox_inches='tight')
    plt.show()
    print(f"✓ {gene} 空间图完成")

# ============================================
# 6.3 空间聚类可视化（修正版）
# ============================================

fig, axes = plt.subplots(2, 4, figsize=(28, 14))
axes = axes.flatten()

for idx, sample in enumerate(['C1','C2','GD1','GD2','GD3','HT1','HT2','HT3']):
    ax = axes[idx]
    
    # 获取该样本在合并数据中的leiden信息
    sample_mask = adata_combined.obs['sample'] == sample
    clusters = adata_combined.obs.loc[sample_mask, 'leiden']
    
    # 从 adata_combined.uns 获取空间信息
    img = adata_combined.uns['spatial'][sample]['images']['lowres']
    scale = adata_combined.uns['spatial'][sample]['scalefactors']['tissue_lowres_scalef']
    
    # 从原始 adata_dict 获取坐标
    coords = adata_dict[sample].obsm['spatial']
    
    ax.imshow(img)
    
    # 按leiden cluster着色
    unique_clusters = sorted(clusters.unique(), key=int)
    cmap = plt.cm.tab20(np.linspace(0, 1, 20))
    
    for cl in unique_clusters:
        cl_mask = (clusters == cl).values
        cl_coords = coords[cl_mask]
        color = cmap[int(cl) % 20]
        ax.scatter(cl_coords[:, 1] * scale, cl_coords[:, 0] * scale,
                   s=12, alpha=0.7, label=cl, color=color)
    
    cond = condition_map[sample]
    ax.set_title(f'{sample} ({cond})', fontsize=13, fontweight='bold')
    ax.axis('off')

handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, loc='center right', title='Leiden',
           fontsize=9, title_fontsize=11, bbox_to_anchor=(1.02, 0.5))

plt.suptitle('Spatial Distribution of Leiden Clusters', fontsize=18, fontweight='bold')
plt.tight_layout()
plt.savefig(r"C:\Users\13918\Desktop\onek1k\bulk\spatial_leiden_clusters.png",
            dpi=600, bbox_inches='tight')
plt.show()

# ============================================
# 7.1 定义你MR分析中涉及的细胞类型的marker基因
# ============================================

# 根据OneK1K数据集的细胞类型定义marker
cell_type_markers = {
    'Memory_B_Cell': ['CD27', 'CD19', 'MS4A1', 'IGHM', 'IGHG1', 'CD40', 'TNFRSF13B'],
    'Dendritic_Cell': ['ITGAX', 'HLA-DRA', 'HLA-DRB1', 'CLEC10A', 'CD1C', 'FCER1A', 'IRF8'],
    'NK_Cell': ['NKG7', 'GNLY', 'KLRD1', 'KLRB1', 'NCAM1', 'FCGR3A', 'PRF1'],
    'Naive_B_Cell': ['MS4A1', 'CD19', 'IGHD', 'IL4R', 'TCL1A', 'FCER2', 'CD22'],
    'CD8_T_Cell': ['CD8A', 'CD8B', 'CD3E', 'CD3D', 'GZMK', 'GZMA', 'CCL5'],
    'CD4_T_Cell': ['CD4', 'CD3E', 'CD3D', 'IL7R', 'LTB', 'LDHB', 'MAL'],
    'Plasma_Cell': ['JCHAIN', 'MZB1', 'SDC1', 'XBP1', 'IGHG1', 'IGHA1', 'PRDM1'],
    'CD8_S100B_T_Cell': ['CD8A', 'CD8B', 'S100B', 'CD3E'],
}

# 检查哪些marker在数据中存在
print("各细胞类型marker基因在数据中的存在情况:")
print("="*60)

valid_markers = {}
for ct, markers in cell_type_markers.items():
    found = [g for g in markers if g in adata_combined.var_names]
    missing = [g for g in markers if g not in adata_combined.var_names]
    valid_markers[ct] = found
    print(f"\n{ct}:")
    print(f"  找到 {len(found)}/{len(markers)}: {found}")
    if missing:
        print(f"  缺失: {missing}")

# ============================================
# 7.1 整理MR分析结果
# ============================================

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scanpy as sc

# MR分析结果：细胞类型 -> 基因 的因果关系
mr_results = pd.DataFrame({
    'cell_type': [
        'Memory B Cell',
        'Dendritic Cell', 
        'NK Recruiting Cell',
        'Naïve/Immature B Cell',
        'CD8 NC T Cell',
        'CD8 S100B T Cell',
        'CD4 NC T Cell',
        'CD8 NC T Cell',
        'Plasma Cell',
        'Plasma Cell',
        'CD4 NC T Cell',
        'CD8 NC T Cell'
    ],
    'gene': [
        'POU5F1', 'HCG11', 'HMGN4', 'ABHD16A',
        'HIST1H3H', 'ZKSCAN4', 'HLA-F', 'HLA-DQB2',
        'ITPR3', 'HLA-DQB1', 'RETREG1', 'RETREG1'
    ]
})

print("MR分析结果汇总:")
print(mr_results.to_string(index=False))
print(f"\n涉及细胞类型数: {mr_results['cell_type'].nunique()}")
print(f"涉及基因数: {mr_results['gene'].nunique()}")
print(f"\n各细胞类型涉及的基因:")
for ct in mr_results['cell_type'].unique():
    genes = mr_results[mr_results['cell_type']==ct]['gene'].tolist()
    print(f"  {ct}: {genes}")


# ============================================
# 7.2 定义各细胞类型的 marker 基因
# 用于空间转录组细胞类型打分/解卷积
# ============================================

# 标准免疫细胞 marker 基因（来自文献）
cell_type_markers = {
    # B细胞系
    'Memory B Cell':        ['CD19', 'CD27', 'CD38', 'MS4A1', 'IGHG1', 'IGHG3',
                             'AIM2', 'FCRL4', 'FCRL5'],
    'Naïve/Immature B Cell':['CD19', 'MS4A1', 'IGHD', 'IGHM', 'IL4R', 'TCL1A',
                             'FCER2', 'CD24'],
    'Plasma Cell':          ['MZB1', 'SDC1', 'CD38', 'IGHG1', 'JCHAIN', 'XBP1',
                             'PRDM1', 'IRF4', 'IGKC'],
    
    # T细胞系
    'CD4 NC T Cell':        ['CD3D', 'CD3E', 'CD4', 'TCF7', 'CCR7', 'LEF1',
                             'IL7R', 'SELL'],
    'CD8 NC T Cell':        ['CD3D', 'CD3E', 'CD8A', 'CD8B', 'TCF7', 'CCR7',
                             'SELL', 'LEF1'],
    'CD8 S100B T Cell':     ['CD3D', 'CD8A', 'S100B', 'GZMK', 'GZMH', 'NKG7',
                             'PRF1', 'GNLY'],
    
    # NK细胞
    'NK Recruiting Cell':   ['NCAM1', 'NKG7', 'GNLY', 'PRF1', 'GZMB', 'KLRD1',
                             'XCL1', 'XCL2', 'CCL3', 'CCL4'],
    
    # 树突状细胞
    'Dendritic Cell':       ['FCER1A', 'CD1C', 'CLEC10A', 'ITGAX', 'HLA-DRA',
                             'HLA-DRB1', 'LILRA4', 'IRF4', 'CLEC9A'],
}

# 检查这些marker基因在数据中是否存在
print("Marker基因可用性检查:")
available_markers = {}
for ct, markers in cell_type_markers.items():
    available = [g for g in markers if g in adata_combined.var_names]
    missing = [g for g in markers if g not in adata_combined.var_names]
    available_markers[ct] = available
    print(f"\n  {ct}:")
    print(f"    可用: {available}")
    if missing:
        print(f"    缺失: {missing}")

# ============================================
# 7.3 用 sc.tl.score_genes 对每个细胞类型打分
# 量化每个spot中各免疫细胞类型的富集程度
# ============================================

# 恢复normalized层
adata_combined.X = adata_combined.layers['normalized'].copy()

print("计算细胞类型得分...")
for ct, markers in available_markers.items():
    if len(markers) >= 3:  # 至少需要3个marker
        score_name = ct.replace('/', '_').replace(' ', '_')
        sc.tl.score_genes(
            adata_combined,
            gene_list=markers,
            score_name=score_name,
            random_state=42
        )
        print(f"  ✓ {ct} -> '{score_name}' 得分计算完成")
    else:
        print(f"  ✗ {ct}: marker基因不足 ({len(markers)}个)，跳过")

print(f"\n所有细胞类型得分已添加到 adata_combined.obs")
print("得分列:")
score_cols = [col for col in adata_combined.obs.columns if any(
    ct.replace('/', '_').replace(' ', '_') in col 
    for ct in cell_type_markers.keys()
)]
print(score_cols)


# ============================================
# 8.1 各细胞类型得分的空间分布
# ============================================

score_cols = ['Memory_B_Cell', 'Naïve_Immature_B_Cell', 'Plasma_Cell',
              'CD4_NC_T_Cell', 'CD8_NC_T_Cell', 'CD8_S100B_T_Cell',
              'NK_Recruiting_Cell', 'Dendritic_Cell']

# 每个细胞类型画一行（8个样本）
for score_name in score_cols:
    fig, axes = plt.subplots(2, 4, figsize=(28, 14))
    axes = axes.flatten()
    
    # 获取全局颜色范围
    all_scores = adata_combined.obs[score_name].values
    vmin_val = np.percentile(all_scores, 2)
    vmax_val = np.percentile(all_scores, 98)
    
    for idx, sample in enumerate(['C1','C2','GD1','GD2','GD3','HT1','HT2','HT3']):
        ax = axes[idx]
        
        # 获取该样本数据
        sample_mask = adata_combined.obs['sample'] == sample
        scores = adata_combined.obs.loc[sample_mask, score_name].values
        
        # 空间信息
        img = adata_combined.uns['spatial'][sample]['images']['lowres']
        scale = adata_combined.uns['spatial'][sample]['scalefactors']['tissue_lowres_scalef']
        coords = adata_dict[sample].obsm['spatial']
        
        ax.imshow(img)
        scatter = ax.scatter(
            coords[:, 1] * scale,
            coords[:, 0] * scale,
            c=scores, cmap='RdYlBu_r', s=12, alpha=0.8,
            vmin=vmin_val, vmax=vmax_val
        )
        plt.colorbar(scatter, ax=ax, shrink=0.6)
        
        cond = condition_map[sample]
        ax.set_title(f'{sample} ({cond})', fontsize=13, fontweight='bold')
        ax.axis('off')
    
    ct_display = score_name.replace('_', ' ')
    plt.suptitle(f'Spatial Distribution: {ct_display} Enrichment Score', 
                 fontsize=18, fontweight='bold')
    plt.tight_layout()
    plt.savefig(rf"C:\Users\13918\Desktop\onek1k\bulk\spatial_score_{score_name}.png",
                dpi=150, bbox_inches='tight')
    plt.show()
    print(f"✓ {ct_display} 空间图完成")

# ============================================
# 8.2 各细胞类型得分在三组间的比较（带星号）
# ============================================

from scipy import stats

def get_stars(pval):
    if pval < 0.0001: return '****'
    elif pval < 0.001: return '***'
    elif pval < 0.01: return '**'
    elif pval < 0.05: return '*'
    else: return 'ns'

fig, axes = plt.subplots(2, 4, figsize=(28, 14))
axes = axes.flatten()

ctrl_mask = adata_combined.obs['condition'] == 'Control'
gd_mask = adata_combined.obs['condition'] == 'GD'
ht_mask = adata_combined.obs['condition'] == 'HT'

for i, score_name in enumerate(score_cols):
    ax = axes[i]
    
    ctrl_vals = adata_combined.obs.loc[ctrl_mask, score_name].values
    gd_vals = adata_combined.obs.loc[gd_mask, score_name].values
    ht_vals = adata_combined.obs.loc[ht_mask, score_name].values
    
    # Violin plot
    data = [ctrl_vals, gd_vals, ht_vals]
    parts = ax.violinplot(data, positions=[0,1,2], showmeans=True)
    colors = ['#2ecc71', '#e74c3c', '#3498db']
    for pc, color in zip(parts['bodies'], colors):
        pc.set_facecolor(color)
        pc.set_alpha(0.7)
    
    ax.set_xticks([0,1,2])
    ax.set_xticklabels(['Control', 'GD', 'HT'], fontsize=11)
    ct_display = score_name.replace('_', ' ')
    ax.set_title(ct_display, fontsize=13, fontweight='bold')
    ax.set_ylabel('Enrichment Score')
    
    # 统计检验
    _, p_gd = stats.mannwhitneyu(ctrl_vals, gd_vals, alternative='two-sided')
    _, p_ht = stats.mannwhitneyu(ctrl_vals, ht_vals, alternative='two-sided')
    _, p_gh = stats.mannwhitneyu(gd_vals, ht_vals, alternative='two-sided')
    
    y_max = max(ctrl_vals.max(), gd_vals.max(), ht_vals.max())
    h = y_max * 0.08 if y_max > 0 else 0.05
    
    # Control vs GD
    y1 = y_max + h
    ax.plot([0,0,1,1], [y1, y1+h*0.3, y1+h*0.3, y1], 'k-', lw=1)
    ax.text(0.5, y1+h*0.3, get_stars(p_gd), ha='center', fontsize=10, color='red')
    
    # Control vs HT
    y2 = y1 + h*1.5
    ax.plot([0,0,2,2], [y2, y2+h*0.3, y2+h*0.3, y2], 'k-', lw=1)
    ax.text(1, y2+h*0.3, get_stars(p_ht), ha='center', fontsize=10, color='blue')
    
    # GD vs HT
    y3 = y2 + h*1.5
    ax.plot([1,1,2,2], [y3, y3+h*0.3, y3+h*0.3, y3], 'k-', lw=1)
    ax.text(1.5, y3+h*0.3, get_stars(p_gh), ha='center', fontsize=10, color='purple')

plt.suptitle('Immune Cell Type Enrichment Scores: Control vs GD vs HT',
             fontsize=16, fontweight='bold')
plt.tight_layout()
plt.savefig(r"C:\Users\13918\Desktop\onek1k\bulk\celltype_scores_violin.png",
            dpi=200, bbox_inches='tight')
plt.show()


# ============================================
# 8.3 修正版：在富集的细胞类型区域中验证MR靶基因
# ============================================

from scipy import stats

mr_pairs = [
    ('Memory_B_Cell',          'POU5F1'),
    ('Dendritic_Cell',         'HCG11'),
    ('NK_Recruiting_Cell',     'HMGN4'),
    ('Naïve_Immature_B_Cell',  'ABHD16A'),
    ('CD8_NC_T_Cell',          'HIST1H3H'),
    ('CD8_S100B_T_Cell',       'ZKSCAN4'),
    ('CD4_NC_T_Cell',          'HLA-F'),
    ('CD8_NC_T_Cell',          'HLA-DQB2'),
    ('Plasma_Cell',            'ITPR3'),
    ('Plasma_Cell',            'HLA-DQB1'),
    ('CD4_NC_T_Cell',          'RETREG1'),
    ('CD8_NC_T_Cell',          'RETREG1'),
]

adata_combined.X = adata_combined.layers['normalized'].copy()

fig, axes = plt.subplots(3, 4, figsize=(28, 20))
axes = axes.flatten()

stat_results_mr = []

for i, (ct_score, gene) in enumerate(mr_pairs):
    ax = axes[i]
    
    # 获取该细胞类型得分 top 25% 的spots
    threshold = adata_combined.obs[ct_score].quantile(0.75)
    enriched_mask = adata_combined.obs[ct_score] >= threshold
    adata_enriched = adata_combined[enriched_mask]
    
    # 提取各组表达值
    expr_dict = {}
    for cond in ['Control', 'GD', 'HT']:
        cond_mask = adata_enriched.obs['condition'] == cond
        if cond_mask.sum() > 0:
            if hasattr(adata_enriched[cond_mask, gene].X, 'toarray'):
                expr_dict[cond] = adata_enriched[cond_mask, gene].X.toarray().flatten()
            else:
                expr_dict[cond] = np.array(adata_enriched[cond_mask, gene].X).flatten()
        else:
            expr_dict[cond] = np.array([])
    
    ctrl_expr = expr_dict['Control']
    gd_expr = expr_dict['GD']
    ht_expr = expr_dict['HT']
    
    # 统计检验（只有两组都非空时才做）
    p_gd = stats.mannwhitneyu(ctrl_expr, gd_expr, alternative='two-sided')[1] if len(ctrl_expr)>5 and len(gd_expr)>5 else 1.0
    p_ht = stats.mannwhitneyu(ctrl_expr, ht_expr, alternative='two-sided')[1] if len(ctrl_expr)>5 and len(ht_expr)>5 else 1.0
    p_gh = stats.mannwhitneyu(gd_expr, ht_expr, alternative='two-sided')[1] if len(gd_expr)>5 and len(ht_expr)>5 else 1.0
    
    stat_results_mr.append({
        'cell_type': ct_score.replace('_', ' '),
        'gene': gene,
        'n_enriched_spots': int(enriched_mask.sum()),
        'n_ctrl': len(ctrl_expr), 'n_gd': len(gd_expr), 'n_ht': len(ht_expr),
        'ctrl_mean': float(ctrl_expr.mean()) if len(ctrl_expr)>0 else 0,
        'gd_mean': float(gd_expr.mean()) if len(gd_expr)>0 else 0,
        'ht_mean': float(ht_expr.mean()) if len(ht_expr)>0 else 0,
        'p_GD_vs_Ctrl': p_gd,
        'p_HT_vs_Ctrl': p_ht,
        'p_GD_vs_HT': p_gh,
    })
    
    # ---- 用 boxplot 代替 violinplot（更稳健，不怕少量数据）----
    plot_data = []
    plot_labels = []
    plot_colors = []
    
    color_map = {'Control': '#2ecc71', 'GD': '#e74c3c', 'HT': '#3498db'}
    
    for j, (cond, expr) in enumerate([('Control', ctrl_expr), ('GD', gd_expr), ('HT', ht_expr)]):
        if len(expr) > 0:
            plot_data.append(expr)
            plot_labels.append(f'{cond}\n(n={len(expr)})')
            plot_colors.append(color_map[cond])
        else:
            # 放一个空占位
            plot_data.append([0])
            plot_labels.append(f'{cond}\n(n=0)')
            plot_colors.append(color_map[cond])
    
    bp = ax.boxplot(plot_data, positions=[0,1,2], widths=0.5,
                     patch_artist=True, showfliers=False)
    for patch, color in zip(bp['boxes'], plot_colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.6)
    
    # 叠加散点（jitter）
    for j, expr in enumerate(plot_data):
        if len(expr) > 1:
            jitter = np.random.normal(0, 0.08, size=len(expr))
            ax.scatter(j + jitter, expr, s=3, alpha=0.3, color='black', zorder=3)
    
    ax.set_xticks([0,1,2])
    ax.set_xticklabels(plot_labels, fontsize=9)
    
    ct_display = ct_score.replace('_', ' ')
    ax.set_title(f'{gene} in {ct_display}', fontsize=11, fontweight='bold')
    ax.set_ylabel('Normalized Expr')
    
    # 星号标注
    all_vals = np.concatenate([e for e in [ctrl_expr, gd_expr, ht_expr] if len(e)>0])
    if len(all_vals) > 0:
        y_max = np.percentile(all_vals, 98)
        h = max(y_max * 0.1, 0.05)
        
        # Control vs GD
        y1 = y_max + h
        ax.plot([0,0,1,1], [y1, y1+h*0.3, y1+h*0.3, y1], 'k-', lw=1)
        ax.text(0.5, y1+h*0.3, get_stars(p_gd), ha='center', fontsize=9, color='red')
        
        # Control vs HT
        y2 = y1 + h*1.3
        ax.plot([0,0,2,2], [y2, y2+h*0.3, y2+h*0.3, y2], 'k-', lw=1)
        ax.text(1, y2+h*0.3, get_stars(p_ht), ha='center', fontsize=9, color='blue')
        
        # GD vs HT
        y3 = y2 + h*1.3
        ax.plot([1,1,2,2], [y3, y3+h*0.3, y3+h*0.3, y3], 'k-', lw=1)
        ax.text(1.5, y3+h*0.3, get_stars(p_gh), ha='center', fontsize=9, color='purple')

plt.suptitle('MR-identified Gene Expression in Cell-Type-Enriched Regions\n'
             '(Top 25% enriched spots per cell type)',
             fontsize=16, fontweight='bold')
plt.tight_layout()
plt.savefig(r"C:\Users\13918\Desktop\onek1k\bulk\MR_validation_celltype_gene.png",
            dpi=200, bbox_inches='tight')
plt.show()

# 打印汇总表
stat_mr_df = pd.DataFrame(stat_results_mr)
stat_mr_df['stars_GD'] = stat_mr_df['p_GD_vs_Ctrl'].apply(get_stars)
stat_mr_df['stars_HT'] = stat_mr_df['p_HT_vs_Ctrl'].apply(get_stars)
stat_mr_df['stars_GH'] = stat_mr_df['p_GD_vs_HT'].apply(get_stars)

print("\n" + "="*100)
print("MR靶基因在细胞类型富集区域的差异表达验证")
print("="*100)
print(stat_mr_df[['cell_type','gene','n_ctrl','n_gd','n_ht',
                   'ctrl_mean','gd_mean','ht_mean',
                   'stars_GD','stars_HT','stars_GH']].to_string(index=False))


# ============================================
# 9.1 汇总热图：Fold Change + 显著性
# ============================================

import matplotlib.patches as mpatches
from matplotlib.colors import TwoSlopeNorm

# 构建数据
stat_mr_df = pd.DataFrame(stat_results_mr)

# 计算 log2 fold change (加伪计数避免除零)
pseudo = 0.001
stat_mr_df['log2FC_GD_vs_Ctrl'] = np.log2((stat_mr_df['gd_mean'] + pseudo) / (stat_mr_df['ctrl_mean'] + pseudo))
stat_mr_df['log2FC_HT_vs_Ctrl'] = np.log2((stat_mr_df['ht_mean'] + pseudo) / (stat_mr_df['ctrl_mean'] + pseudo))

# 创建标签
stat_mr_df['label'] = stat_mr_df['gene'] + '\nin ' + stat_mr_df['cell_type']

# 构建热图矩阵
heatmap_data = stat_mr_df[['log2FC_GD_vs_Ctrl', 'log2FC_HT_vs_Ctrl']].values
labels_y = stat_mr_df['label'].values

fig, ax = plt.subplots(figsize=(8, 14))

# 对称色标
vmax = np.max(np.abs(heatmap_data[np.isfinite(heatmap_data)]))
norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)

im = ax.imshow(heatmap_data, cmap='RdBu_r', aspect='auto', norm=norm)

# 添加星号和数值
for row in range(heatmap_data.shape[0]):
    for col in range(heatmap_data.shape[1]):
        val = heatmap_data[row, col]
        
        if col == 0:
            star = get_stars(stat_mr_df.iloc[row]['p_GD_vs_Ctrl'])
        else:
            star = get_stars(stat_mr_df.iloc[row]['p_HT_vs_Ctrl'])
        
        # 显示fold change值和星号
        text_color = 'white' if abs(val) > vmax*0.6 else 'black'
        ax.text(col, row, f'{val:.2f}\n{star}', ha='center', va='center',
                fontsize=9, fontweight='bold', color=text_color)

ax.set_xticks([0, 1])
ax.set_xticklabels(['GD vs Control', 'HT vs Control'], fontsize=13, fontweight='bold')
ax.set_yticks(range(len(labels_y)))
ax.set_yticklabels(labels_y, fontsize=10)

plt.colorbar(im, ax=ax, label='log2(Fold Change)', shrink=0.6)
ax.set_title('MR-identified Genes: Fold Change in Cell-Type-Enriched Regions\n'
             '(Spatial Transcriptomics Validation)',
             fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig(r"C:\Users\13918\Desktop\onek1k\bulk\MR_heatmap_foldchange.png",
            dpi=200, bbox_inches='tight')
plt.show()

print("\n完整统计表:")
print(stat_mr_df[['cell_type','gene','n_ctrl','n_gd','n_ht',
                   'ctrl_mean','gd_mean','ht_mean',
                   'log2FC_GD_vs_Ctrl','log2FC_HT_vs_Ctrl',
                   'stars_GD','stars_HT']].to_string(index=False))


# ============================================
# 9.2 细胞类型富集得分与靶基因表达的相关性
#     → 证明基因表达确实随细胞类型富集增加
# ============================================

adata_combined.X = adata_combined.layers['normalized'].copy()

fig, axes = plt.subplots(3, 4, figsize=(28, 20))
axes = axes.flatten()

corr_results = []

for i, (ct_score, gene) in enumerate(mr_pairs):
    ax = axes[i]
    
    # 获取得分和表达值
    scores = adata_combined.obs[ct_score].values
    
    if hasattr(adata_combined[:, gene].X, 'toarray'):
        expr = adata_combined[:, gene].X.toarray().flatten()
    else:
        expr = np.array(adata_combined[:, gene].X).flatten()
    
    # 按condition分组绘制散点
    color_map = {'Control': '#2ecc71', 'GD': '#e74c3c', 'HT': '#3498db'}
    
    for cond in ['Control', 'GD', 'HT']:
        mask = adata_combined.obs['condition'].values == cond
        ax.scatter(scores[mask], expr[mask], s=2, alpha=0.2, 
                   color=color_map[cond], label=cond)
    
    # 整体Spearman相关
    from scipy.stats import spearmanr
    r, p = spearmanr(scores, expr)
    corr_results.append({
        'cell_type': ct_score.replace('_', ' '),
        'gene': gene,
        'spearman_r': r,
        'spearman_p': p
    })
    
    ct_display = ct_score.replace('_', ' ')
    ax.set_xlabel(f'{ct_display} Score', fontsize=9)
    ax.set_ylabel(f'{gene} Expression', fontsize=9)
    ax.set_title(f'{gene} vs {ct_display}\nr={r:.3f}, p={p:.1e}',
                 fontsize=10, fontweight='bold')
    ax.legend(fontsize=8, markerscale=3)

plt.suptitle('Correlation: Cell Type Enrichment Score vs Target Gene Expression',
             fontsize=0, fontweight='bold')
plt.tight_layout()
plt.savefig(r"C:\Users\13918\Desktop\onek1k\bulk\MR_correlation_score_gene.png",
            dpi=600, bbox_inches='tight')
plt.show()

# 打印相关性结果
corr_df = pd.DataFrame(corr_results)
corr_df['significance'] = corr_df['spearman_p'].apply(get_stars)
print("\n相关性分析结果:")
print(corr_df.to_string(index=False))


# ============================================
# 9.3 免疫细胞浸润程度比较（关键发现！）
#     → 证明AITD中免疫细胞浸润增加
# ============================================

fig, axes = plt.subplots(2, 4, figsize=(28, 14))
axes = axes.flatten()

# 计算每个样本中各细胞类型得分 top25% 的spots比例
score_cols = ['Memory_B_Cell', 'Naïve_Immature_B_Cell', 'Plasma_Cell',
              'CD4_NC_T_Cell', 'CD8_NC_T_Cell', 'CD8_S100B_T_Cell',
              'NK_Recruiting_Cell', 'Dendritic_Cell']

infiltration_data = []

for ct_score in score_cols:
    threshold = adata_combined.obs[ct_score].quantile(0.75)
    
    for sample in ['C1','C2','GD1','GD2','GD3','HT1','HT2','HT3']:
        sample_mask = adata_combined.obs['sample'] == sample
        sample_scores = adata_combined.obs.loc[sample_mask, ct_score]
        
        pct_enriched = (sample_scores >= threshold).sum() / len(sample_scores) * 100
        
        infiltration_data.append({
            'cell_type': ct_score.replace('_', ' '),
            'sample': sample,
            'condition': condition_map[sample],
            'pct_enriched': pct_enriched
        })

infil_df = pd.DataFrame(infiltration_data)

for i, ct_score in enumerate(score_cols):
    ax = axes[i]
    ct_display = ct_score.replace('_', ' ')
    
    ct_data = infil_df[infil_df['cell_type'] == ct_display]
    
    colors_bar = {'Control': '#2ecc71', 'GD': '#e74c3c', 'HT': '#3498db'}
    
    for j, (_, row) in enumerate(ct_data.iterrows()):
        ax.bar(j, row['pct_enriched'], 
               color=colors_bar[row['condition']], 
               alpha=0.8, edgecolor='black', linewidth=0.5)
    
    ax.set_xticks(range(len(ct_data)))
    ax.set_xticklabels(ct_data['sample'], rotation=45, fontsize=9)
    ax.set_ylabel('% Enriched Spots')
    ax.set_title(ct_display, fontsize=12, fontweight='bold')
    ax.axhline(y=25, color='gray', linestyle='--', alpha=0.5, label='Expected 25%')
    
    if i == 0:
        legend_patches = [mpatches.Patch(color=c, label=l) for l, c in colors_bar.items()]
        ax.legend(handles=legend_patches, fontsize=8)

plt.suptitle('Immune Cell Infiltration: % of Spots Enriched for Each Cell Type\n'
             '(Dashed line = expected 25% if uniform distribution)',
             fontsize=16, fontweight='bold')
plt.tight_layout()
plt.savefig(r"C:\Users\13918\Desktop\onek1k\bulk\immune_infiltration_barplot.png",
            dpi=200, bbox_inches='tight')
plt.show()

# 按condition汇总
print("\n各条件下免疫细胞浸润比例 (%):")
pivot = infil_df.groupby(['cell_type', 'condition'])['pct_enriched'].mean().unstack()
pivot = pivot[['Control', 'GD', 'HT']]
print(pivot.round(1).to_string())



# ============================================
# 10.3 共定位空间图：选3个最显著的MR配对
#     展示细胞类型得分和基因表达在同一组织上的叠加
# ============================================

# 选择最有意义的3对
top_pairs = [
    ('CD4_NC_T_Cell',    'HLA-F',    'CD4 NC T Cell'),
    ('Plasma_Cell',      'HLA-DQB1', 'Plasma Cell'),
    ('Naïve_Immature_B_Cell', 'ABHD16A', 'Naïve/Immature B Cell'),
]

# 选择代表性样本：C1(Control), GD1(GD), HT1(HT)
repr_samples = ['C1', 'GD1', 'HT1']

for ct_score, gene, ct_display in top_pairs:
    fig, axes = plt.subplots(2, 3, figsize=(21, 14))
    
    for col, sample in enumerate(repr_samples):
        sample_mask = adata_combined.obs['sample'] == sample
        
        # 空间信息
        img = adata_combined.uns['spatial'][sample]['images']['lowres']
        scale = adata_combined.uns['spatial'][sample]['scalefactors']['tissue_lowres_scalef']
        coords = adata_dict[sample].obsm['spatial']
        
        # 上排：细胞类型得分
        ax = axes[0, col]
        scores = adata_combined.obs.loc[sample_mask, ct_score].values
        ax.imshow(img)
        sc1 = ax.scatter(coords[:, 1] * scale, coords[:, 0] * scale,
                         c=scores, cmap='YlOrRd', s=15, alpha=0.8)
        plt.colorbar(sc1, ax=ax, shrink=0.5, label='Score')
        cond = condition_map[sample]
        ax.set_title(f'{sample} ({cond})\n{ct_display} Score', fontsize=12, fontweight='bold')
        ax.axis('off')
        
        # 下排：基因表达
        ax = axes[1, col]
        
        # 从原始数据获取标准化表达
        adata_s = adata_dict[sample].copy()
        sc.pp.normalize_total(adata_s, target_sum=1e4)
        sc.pp.log1p(adata_s)
        
        if hasattr(adata_s[:, gene].X, 'toarray'):
            expr = adata_s[:, gene].X.toarray().flatten()
        else:
            expr = np.array(adata_s[:, gene].X).flatten()
        
        ax.imshow(img)
        sc2 = ax.scatter(coords[:, 1] * scale, coords[:, 0] * scale,
                         c=expr, cmap='Purples', s=15, alpha=0.8,
                         vmin=0, vmax=np.percentile(expr, 99) if expr.max()>0 else 1)
        plt.colorbar(sc2, ax=ax, shrink=0.5, label='Expression')
        ax.set_title(f'{gene} Expression', fontsize=12, fontweight='bold')
        ax.axis('off')
    
    plt.suptitle(f'Spatial Co-localization: {ct_display} Enrichment & {gene} Expression\n'
                 f'(MR-identified causal pair)',
                 fontsize=16, fontweight='bold')
    plt.tight_layout()
    plt.savefig(rf"C:\Users\13918\Desktop\onek1k\bulk\colocalization_{ct_score}_{gene}.png",
                dpi=200, bbox_inches='tight')
    plt.show()
    print(f"✓ {ct_display} - {gene} 共定位图完成")


# ============================================
# 10.1 修正版：重新构建完整的汇总表格
# ============================================

# 从 stat_results_mr 重新构建
stat_mr_df = pd.DataFrame(stat_results_mr)

# 添加stars列
stat_mr_df['stars_GD'] = stat_mr_df['p_GD_vs_Ctrl'].apply(get_stars)
stat_mr_df['stars_HT'] = stat_mr_df['p_HT_vs_Ctrl'].apply(get_stars)
stat_mr_df['stars_GH'] = stat_mr_df['p_GD_vs_HT'].apply(get_stars)

# 计算 log2 fold change
pseudo = 0.001
stat_mr_df['log2FC_GD_vs_Ctrl'] = np.log2(
    (stat_mr_df['gd_mean'] + pseudo) / (stat_mr_df['ctrl_mean'] + pseudo))
stat_mr_df['log2FC_HT_vs_Ctrl'] = np.log2(
    (stat_mr_df['ht_mean'] + pseudo) / (stat_mr_df['ctrl_mean'] + pseudo))

print("stat_mr_df 列名:", list(stat_mr_df.columns))
print(stat_mr_df[['cell_type','gene','stars_GD','stars_HT',
                   'log2FC_GD_vs_Ctrl','log2FC_HT_vs_Ctrl']].to_string(index=False))

# 构建 summary_df
summary_df = stat_mr_df[['cell_type','gene','n_ctrl','n_gd','n_ht',
                          'ctrl_mean','gd_mean','ht_mean',
                          'log2FC_GD_vs_Ctrl','log2FC_HT_vs_Ctrl',
                          'stars_GD','stars_HT']].copy()

# 加入相关性
summary_df['spearman_r'] = corr_df['spearman_r'].values
summary_df['corr_sig'] = corr_df['significance'].values

# 加入浸润比例
infil_pivot = infil_df.groupby(['cell_type','condition'])['pct_enriched'].mean().unstack()

for idx, row in summary_df.iterrows():
    ct = row['cell_type']
    if ct in infil_pivot.index:
        summary_df.loc[idx, 'infil_Ctrl'] = infil_pivot.loc[ct, 'Control']
        summary_df.loc[idx, 'infil_GD'] = infil_pivot.loc[ct, 'GD']
        summary_df.loc[idx, 'infil_HT'] = infil_pivot.loc[ct, 'HT']

print("\n" + "="*120)
print("综合结果汇总表 (Table for manuscript)")
print("="*120)
print(summary_df.round(3).to_string(index=False))

# 保存CSV
summary_df.to_csv(r"C:\Users\13918\Desktop\onek1k\bulk\MR_ST_validation_summary.csv", index=False)
print("\n✓ 表格已保存为CSV！")


# ============================================
# 10.2 修正版：综合四面板图
# ============================================

import matplotlib.patches as mpatches
from matplotlib.colors import TwoSlopeNorm

fig, axes = plt.subplots(1, 4, figsize=(28, 12), 
                          gridspec_kw={'width_ratios': [3, 3, 1.5, 2]})

y_labels = [f"{row['gene']}\nin {row['cell_type']}" for _, row in summary_df.iterrows()]
n = len(summary_df)

# ---- 左1：GD vs Control log2FC ----
ax = axes[0]
fc_gd = summary_df['log2FC_GD_vs_Ctrl'].values
colors_gd = ['#e74c3c' if v > 0 else '#3498db' for v in fc_gd]

ax.barh(range(n), fc_gd, color=colors_gd, alpha=0.7, edgecolor='black')
for j, (val, star) in enumerate(zip(fc_gd, summary_df['stars_GD'].values)):
    offset = 0.5 if val >= 0 else -0.5
    ax.text(val + offset, j, star, ha='center', va='center', fontsize=9, fontweight='bold')

ax.set_yticks(range(n))
ax.set_yticklabels(y_labels, fontsize=9)
ax.set_xlabel('log2(Fold Change)', fontsize=11)
ax.set_title('GD vs Control', fontsize=14, fontweight='bold')
ax.axvline(x=0, color='black', linewidth=0.8)
ax.invert_yaxis()

# ---- 左2：HT vs Control log2FC ----
ax = axes[1]
fc_ht = summary_df['log2FC_HT_vs_Ctrl'].values
colors_ht = ['#e74c3c' if v > 0 else '#3498db' for v in fc_ht]

ax.barh(range(n), fc_ht, color=colors_ht, alpha=0.7, edgecolor='black')
for j, (val, star) in enumerate(zip(fc_ht, summary_df['stars_HT'].values)):
    offset = 0.5 if val >= 0 else -0.5
    ax.text(val + offset, j, star, ha='center', va='center', fontsize=9, fontweight='bold')

ax.set_yticks(range(n))
ax.set_yticklabels([])
ax.set_xlabel('log2(Fold Change)', fontsize=11)
ax.set_title('HT vs Control', fontsize=14, fontweight='bold')
ax.axvline(x=0, color='black', linewidth=0.8)
ax.invert_yaxis()

# ---- 中：Spearman r ----
ax = axes[2]
r_vals = summary_df['spearman_r'].values
colors_r = ['#e74c3c' if r > 0 else '#3498db' for r in r_vals]

ax.barh(range(n), r_vals, color=colors_r, alpha=0.7, edgecolor='black')
for j, (val, star) in enumerate(zip(r_vals, summary_df['corr_sig'].values)):
    offset = 0.03 if val >= 0 else -0.03
    ax.text(val + offset, j, star, ha='center', va='center', fontsize=8)

ax.set_yticks(range(n))
ax.set_yticklabels([])
ax.set_xlabel('Spearman r', fontsize=11)
ax.set_title('Score-Gene\nCorrelation', fontsize=12, fontweight='bold')
ax.axvline(x=0, color='black', linewidth=0.8)
ax.invert_yaxis()

# ---- 右：免疫浸润比例 ----
ax = axes[3]
x = np.arange(n)
width = 0.25

ctrl_infil = summary_df['infil_Ctrl'].values
gd_infil = summary_df['infil_GD'].values
ht_infil = summary_df['infil_HT'].values

ax.barh(x - width, ctrl_infil, width, label='Control', color='#2ecc71', alpha=0.7, edgecolor='black')
ax.barh(x, gd_infil, width, label='GD', color='#e74c3c', alpha=0.7, edgecolor='black')
ax.barh(x + width, ht_infil, width, label='HT', color='#3498db', alpha=0.7, edgecolor='black')

ax.axvline(x=25, color='gray', linestyle='--', alpha=0.5, label='Expected 25%')
ax.set_yticks(range(n))
ax.set_yticklabels([])
ax.set_xlabel('% Enriched Spots', fontsize=11)
ax.set_title('Immune Cell\nInfiltration', fontsize=12, fontweight='bold')
ax.legend(fontsize=8, loc='lower right')
ax.invert_yaxis()

plt.suptitle('Integrative Spatial Transcriptomics Validation of MR-identified\n'
             'Cell-Type-Specific Causal Genes in Autoimmune Thyroid Disease',
             fontsize=16, fontweight='bold')
plt.tight_layout()
plt.savefig(r"C:\Users\13918\Desktop\onek1k\bulk\MR_ST_integrated_summary.png",
            dpi=300, bbox_inches='tight')
plt.show()
print("✓ 综合图完成！")
