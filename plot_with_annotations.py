#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu

def extract_scores(tsv_file):
    scores = []
    smiles_list = []
    with open(tsv_file, 'r') as f:
        lines = f.readlines()
        if not lines:
            return [], []
        for line in lines[1:]:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                try:
                    scores.append(float(parts[1]))
                    smiles_list.append(parts[0])
                except:
                    scores.append(1.0)
                    smiles_list.append(parts[0])
    return scores, smiles_list

normal_scores, normal_smiles = extract_scores("normal_raw.tsv")
problematic_scores, _ = extract_scores("abnormal_raw.tsv")   # 同理，文件名按实际情况改

# 找出正常组中得分 > 0.01 的分子（略高于0）
outliers = [(smiles, score) for smiles, score in zip(normal_smiles, normal_scores) if score > 0.01]

print("正常组中得分 > 0.01 的分子：")
for smiles, score in outliers:
    print(f"  {smiles}: {score:.4f}")

df_normal = pd.DataFrame({'Score': normal_scores, 'Group': 'Normal'})
df_problematic = pd.DataFrame({'Score': problematic_scores, 'Group': 'Problematic'})
data = pd.concat([df_normal, df_problematic], ignore_index=True)

plt.figure(figsize=(8, 10))
ax = sns.boxplot(x='Group', y='Score', data=data, palette=['#1f77b4', '#ff7f0e'])
sns.stripplot(x='Group', y='Score', data=data, color='black', alpha=0.3, size=2)

# 在正常组箱线图上方添加注释，标注离群点
if outliers:
    x_pos = 0
    max_score = max(normal_scores)
    ax.annotate('Small rings (e.g., C1OC1, C1NCC1)\nwith ring strain',
                xy=(x_pos, max_score), xytext=(x_pos + 0.3, max_score + 0.1),
                arrowprops=dict(arrowstyle='->', color='red'),
                fontsize=9, color='red')

plt.title('Atomic-SCS Scores: Normal vs. Problematic Molecules\n(Normal group includes small rings with non-zero scores due to ring strain)')
plt.ylabel('Compliance Score (lower is better)')
plt.grid(axis='y', linestyle='--', alpha=0.7)
plt.tight_layout()
plt.savefig('boxplot_annotated.png', dpi=300)
plt.savefig('boxplot_annotated.pdf')
print("带注释的箱线图已保存")

# 统计检验
u_stat, p_value = mannwhitneyu(normal_scores, problematic_scores, alternative='less')
print(f"Mann-Whitney U test: U = {u_stat:.2f}, p = {p_value:.2e}")