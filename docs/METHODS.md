# 增强因果森林分析模块说明

## 概述

本项目在原有因果森林分析基础上，新增了四个专门的分析模块，旨在突出因果森林方法的创新性和临床实用性，与传统方法形成差异化。

## 新增模块结构

### 1. 个体化治疗策略模块 (`individualized_treatment.R`)

**核心功能：**
- 基于CATE值的患者聚类分析
- 个体化治疗建议生成
- 临床案例模拟
- 治疗效应分布可视化

**主要函数：**
- `perform_patient_clustering()`: 患者聚类分析
- `generate_treatment_recommendations()`: 生成个体化建议
- `plot_cate_distribution()`: 可视化CATE分布
- `simulate_clinical_cases()`: 模拟典型临床案例

### 2. 治疗异质性分析模块 (`treatment_heterogeneity.R`)

**核心功能：**
- 异质性统计检验
- 变量间交互效应分析
- 连续变量与治疗效应关系
- 异质性程度量化

**主要函数：**
- `test_treatment_heterogeneity()`: 异质性检验
- `analyze_interaction_effects()`: 交互效应分析
- `analyze_continuous_relationships()`: 连续变量关系分析
- `quantify_heterogeneity()`: 异质性量化

### 3. 临床决策支持模块 (`clinical_decision_support.R`)

**核心功能：**
- 简化的临床评分系统
- 决策树构建
- 决策曲线分析
- 不确定性量化

**主要函数：**
- `create_clinical_score()`: 创建评分系统
- `build_decision_tree()`: 构建决策树
- `perform_decision_curve_analysis()`: 决策曲线分析
- `quantify_prediction_uncertainty()`: 不确定性评估

### 4. 验证和稳健性分析模块 (`validation_analysis.R`)

**核心功能：**
- Bootstrap验证
- 交叉验证
- 敏感性分析
- 校准分析

**主要函数：**
- `bootstrap_validation()`: Bootstrap验证
- `cross_validation_subgroups()`: 交叉验证
- `hyperparameter_sensitivity()`: 超参数敏感性
- `calibration_analysis()`: 校准分析

## 使用方法

### 1. 基础准备

```r
# 在R中运行
source("enhanced_main.R")

# 执行完整的增强分析
main_results <- main_enhanced_analysis()
```

### 2. 单独运行特定分析

```r
# 加载必要的模块
source("individualized_treatment.R")

# 执行患者聚类
clustering_results <- perform_patient_clustering(matched_data, n_clusters = 4)

# 生成治疗建议
recommendations <- generate_treatment_recommendations(matched_data)
```

### 3. 结果输出位置

所有分析结果将保存在 `results/enhanced_analysis/` 目录下，按模块分类：

```
results/enhanced_analysis/
├── individualized_treatment/    # 个体化治疗分析结果
├── heterogeneity_analysis/      # 异质性分析结果
├── clinical_decision_support/   # 决策支持工具结果
├── validation_analysis/         # 验证分析结果
└── visualizations/              # 所有可视化图表
```

## 与传统方法的差异化优势

### 1. **从群体到个体**
- 传统方法：关注平均治疗效应
- 因果森林：提供个体化治疗效应估计

### 2. **从静态到动态**
- 传统方法：固定的亚组分析
- 因果森林：数据驱动的异质性发现

### 3. **从定性到定量**
- 传统方法：基于临床经验的决策
- 因果森林：量化的决策支持工具

### 4. **从假设到验证**
- 传统方法：预设的亚组假设
- 因果森林：稳健的验证框架

## 新文章的核心贡献点

1. **精准医疗视角**：从"是否停药"转向"谁应该停药"
2. **个体化决策工具**：开发实用的临床评分系统
3. **异质性深度挖掘**：识别影响治疗效应的关键因素
4. **临床实施路径**：提供具体的决策支持框架

## 技术要求

### R包依赖
```r
install.packages(c("grf", "dplyr", "ggplot2", "rpart", "rpart.plot", 
                   "cluster", "factoextra", "stringr", "viridis", 
                   "gridExtra", "pROC", "boot"))
```

### 数据要求
- 匹配后的数据集（包含CATE估计）
- 训练好的因果森林模型
- 标准化的变量命名

## 故障排除

### 常见问题
1. **内存不足**：减少Bootstrap迭代次数
2. **运行时间长**：使用较少的树数量进行测试
3. **缺少变量**：检查数据列名是否匹配

### 调试建议
- 使用小样本数据测试各模块
- 逐步运行而非一次性执行所有分析
- 检查中间结果的合理性

## 联系信息

如有技术问题，请检查：
1. R包版本兼容性
2. 数据格式要求
3. 内存和计算资源