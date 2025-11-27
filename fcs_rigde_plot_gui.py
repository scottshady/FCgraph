#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PyQt6-based GUI application for reading FCS flow cytometry files, automatically detecting fluorescence channels, calculating MFI (mean/median) statistics, and generating customizable ridge plots with color selection and file ordering.
"""

import os
import glob
import numpy as np
import pandas as pd
import warnings
import sys
from PyQt6.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
                             QLabel, QLineEdit, QPushButton, QTableWidget, QTableWidgetItem,
                             QFileDialog, QMessageBox, QDialog, QSpinBox, QCheckBox, QColorDialog,
                             QHeaderView, QProgressDialog, QAbstractItemView, QStyleOptionButton, QStyle)
from PyQt6.QtCore import Qt, QThread, pyqtSignal, QMimeData, QRect
from PyQt6.QtGui import QColor, QFont, QPalette, QDrag, QPainter, QPixmap, QPen
warnings.filterwarnings('ignore')

# 延迟导入matplotlib，避免NumPy版本冲突
try:
    import matplotlib
    matplotlib.use('Agg')  # 使用非交互式后端
    import matplotlib.pyplot as plt
    import seaborn as sns
    MATPLOTLIB_AVAILABLE = True
except ImportError as e:
    print(f"警告: matplotlib导入失败: {e}")
    MATPLOTLIB_AVAILABLE = False
except Exception as e:
    print(f"警告: matplotlib初始化失败: {e}")
    MATPLOTLIB_AVAILABLE = False

# 尝试导入不同的FCS读取库
USE_FLOWCAL = False
USE_FLOWIO = False
USE_FLOWUTILS = False
USE_FCSPARSER = False

# 优先使用fcsparser，因为它最简单且不需要复杂的依赖
try:
    import fcsparser
    USE_FCSPARSER = True
except ImportError:
    try:
        from flowcal import io as flowcal_io
        USE_FLOWCAL = True
    except ImportError:
        try:
            from flowio import FlowIO
            USE_FLOWIO = True
        except ImportError:
            try:
                import flowutils
                USE_FLOWUTILS = True
            except ImportError:
                print("错误: 需要安装以下库之一:")
                print("  pip install fcsparser  (推荐)")

# 设置中文字体（如果matplotlib可用）
if MATPLOTLIB_AVAILABLE:
    plt.rcParams['font.sans-serif'] = ['Arial Unicode MS', 'SimHei', 'DejaVu Sans']
    plt.rcParams['axes.unicode_minus'] = False

def detect_fitc_channel(channel_names):
    """
    自动识别FITC通道（优先）
    
    参数:
        channel_names: 通道名列表
        
    返回:
        通道名（字符串）或None
    """
    # 排除散射光和时间通道
    exclude_patterns = ['FSC', 'SSC', 'Time', 'Width']
    fluorescence_candidates = [
        ch for ch in channel_names 
        if not any(pattern in ch for pattern in exclude_patterns)
    ]
    
    # 优先识别FITC/488通道
    fitc_keywords = ['FITC', 'FL1', 'B530', '488']
    matched_channels = []
    
    for keyword in fitc_keywords:
        matched = [ch for ch in fluorescence_candidates if keyword.upper() in ch.upper()]
        if matched:
            matched_channels.extend(matched)
    
    if matched_channels:
        # 优先选择高度通道（H），如果没有则选择面积通道（A）
        height_channels = [ch for ch in matched_channels if '-H' in ch or 'H-' in ch]
        if height_channels:
            return height_channels[0]
        else:
            # 如果没有高度通道，使用面积通道
            area_channels = [ch for ch in matched_channels if '-A' in ch or 'A-' in ch]
            if area_channels:
                return area_channels[0]
            else:
                return matched_channels[0]
    
    # 如果找不到FITC，返回第一个非散射光通道
    if fluorescence_candidates:
        return fluorescence_candidates[0]
    
    return None

def read_fcs_file(file_path):
    """
    读取FCS文件（使用不同的库）
    
    返回:
        (data, channel_names) 或 None
    """
    try:
        if USE_FLOWCAL:
            fcs_data = flowcal_io.FCSData(file_path)
            channel_names = fcs_data.channels
            data = fcs_data
            return data, channel_names
        elif USE_FLOWIO:
            flowio = FlowIO(file_path)
            channel_names = [p['PnN'] for p in flowio.text['$P']]
            data = flowio.events
            return data, channel_names
        elif USE_FLOWUTILS:
            try:
                from flowutils import io as flowutils_io
                data = flowutils_io.load_fcs(file_path)
                channel_names = data.channels if hasattr(data, 'channels') else list(data.columns)
                return data, channel_names
            except:
                return None, None
        elif USE_FCSPARSER:
            # 使用fcsparser
            meta, data = fcsparser.parse(file_path, reformat_meta=True)
            # 从metadata中获取通道名
            channel_names = []
            for i in range(1, len(meta) + 1):
                pn_key = f'$P{i}N'
                if pn_key in meta:
                    channel_names.append(meta[pn_key])
            # 如果无法从metadata获取，使用DataFrame的列名
            if not channel_names and isinstance(data, pd.DataFrame):
                channel_names = list(data.columns)
            return data, channel_names
        else:
            return None, None
    except Exception as e:
        print(f"错误: 读取文件 {os.path.basename(file_path)} 时出错: {str(e)}")
        return None, None

def process_fcs_file_for_mfi(file_path, fluorescence_channel):
    """
    处理单个FCS文件，计算MFI的mean和median值
    
    参数:
        file_path: FCS文件路径
        fluorescence_channel: 荧光通道名
        
    返回:
        (MFI mean值, MFI median值) 或 (None, None)
    """
    try:
        # 读取FCS文件
        fcs_data, channel_names = read_fcs_file(file_path)
        
        if fcs_data is None:
            return None, None
        
        # 检查通道是否存在（处理可能的格式差异）
        channel_to_use = None
        channel_index = None
        possible_names = [
            fluorescence_channel,
            fluorescence_channel.replace('-', '.'),
            fluorescence_channel.replace('-', ''),
        ]
        
        for name in possible_names:
            if name in channel_names:
                channel_to_use = name
                channel_index = channel_names.index(name)
                break
        
        if channel_to_use is None:
            return None, None
        
        # 根据使用的库提取数据
        if USE_FCSPARSER:
            if isinstance(fcs_data, pd.DataFrame):
                # 应用门控
                fsc_a_col = None
                ssc_a_col = None
                
                for col in channel_names:
                    if 'FSC-A' in col or 'FSC.A' in col:
                        fsc_a_col = col
                    if 'SSC-A' in col or 'SSC.A' in col:
                        ssc_a_col = col
                
                if fsc_a_col and ssc_a_col:
                    gate_mask = (
                        (fcs_data[fsc_a_col] >= 0) & 
                        (fcs_data[fsc_a_col] <= 10000000) &
                        (fcs_data[ssc_a_col] >= 0) & 
                        (fcs_data[ssc_a_col] <= 1500000)
                    )
                    gated_data = fcs_data[gate_mask]
                else:
                    gated_data = fcs_data
                
                # 提取荧光通道数据
                fluorescence_values = gated_data[channel_to_use].values
            else:
                return None, None
        elif USE_FLOWCAL:
            fsc_a_col = None
            ssc_a_col = None
            
            for col in channel_names:
                if 'FSC-A' in col or 'FSC.A' in col:
                    fsc_a_col = col
                if 'SSC-A' in col or 'SSC.A' in col:
                    ssc_a_col = col
            
            if fsc_a_col and ssc_a_col:
                gate_mask = (
                    (fcs_data[:, fsc_a_col] >= 0) & 
                    (fcs_data[:, fsc_a_col] <= 10000000) &
                    (fcs_data[:, ssc_a_col] >= 0) & 
                    (fcs_data[:, ssc_a_col] <= 1500000)
                )
                gated_data = fcs_data[gate_mask]
            else:
                gated_data = fcs_data
            
            fluorescence_values = gated_data[:, channel_to_use]
        else:
            if isinstance(fcs_data, np.ndarray):
                data_array = fcs_data
            else:
                data_array = np.array(fcs_data)
            
            fsc_a_idx = None
            ssc_a_idx = None
            
            for i, col in enumerate(channel_names):
                if 'FSC-A' in col or 'FSC.A' in col:
                    fsc_a_idx = i
                if 'SSC-A' in col or 'SSC.A' in col:
                    ssc_a_idx = i
            
            if fsc_a_idx is not None and ssc_a_idx is not None:
                gate_mask = (
                    (data_array[:, fsc_a_idx] >= 0) & 
                    (data_array[:, fsc_a_idx] <= 10000000) &
                    (data_array[:, ssc_a_idx] >= 0) & 
                    (data_array[:, ssc_a_idx] <= 1500000)
                )
                gated_data = data_array[gate_mask]
            else:
                gated_data = data_array
            
            fluorescence_values = gated_data[:, channel_index]
        
        # 移除负值和NaN
        valid_mask = (fluorescence_values >= 0) & (~np.isnan(fluorescence_values))
        fluorescence_values = fluorescence_values[valid_mask]
        
        if len(fluorescence_values) == 0:
            return None, None
        
        # 计算MFI的mean和median值
        mfi_mean = np.mean(fluorescence_values)
        mfi_median = np.median(fluorescence_values)
        return mfi_mean, mfi_median
        
    except Exception as e:
        print(f"错误: 处理文件 {os.path.basename(file_path)} 时出错: {str(e)}")
        return None, None

def process_fcs_file(file_path, fluorescence_channel):
    """
    处理单个FCS文件
    
    参数:
        file_path: FCS文件路径
        fluorescence_channel: 荧光通道名
        
    返回:
        DataFrame包含荧光通道数据和样本名
    """
    try:
        # 读取FCS文件
        fcs_data, channel_names = read_fcs_file(file_path)
        
        if fcs_data is None:
            return None
        
        # 检查通道是否存在（处理可能的格式差异）
        channel_to_use = None
        channel_index = None
        possible_names = [
            fluorescence_channel,
            fluorescence_channel.replace('-', '.'),
            fluorescence_channel.replace('-', ''),
        ]
        
        for name in possible_names:
            if name in channel_names:
                channel_to_use = name
                channel_index = channel_names.index(name)
                break
        
        if channel_to_use is None:
            return None
        
        # 根据使用的库提取数据
        if USE_FCSPARSER:
            if isinstance(fcs_data, pd.DataFrame):
                # 应用门控
                fsc_a_col = None
                ssc_a_col = None
                
                for col in channel_names:
                    if 'FSC-A' in col or 'FSC.A' in col:
                        fsc_a_col = col
                    if 'SSC-A' in col or 'SSC.A' in col:
                        ssc_a_col = col
                
                if fsc_a_col and ssc_a_col:
                    gate_mask = (
                        (fcs_data[fsc_a_col] >= 0) & 
                        (fcs_data[fsc_a_col] <= 10000000) &
                        (fcs_data[ssc_a_col] >= 0) & 
                        (fcs_data[ssc_a_col] <= 1500000)
                    )
                    gated_data = fcs_data[gate_mask]
                else:
                    gated_data = fcs_data
                
                # 提取荧光通道数据
                fluorescence_values = gated_data[channel_to_use].values
            else:
                return None
        elif USE_FLOWCAL:
            fsc_a_col = None
            ssc_a_col = None
            
            for col in channel_names:
                if 'FSC-A' in col or 'FSC.A' in col:
                    fsc_a_col = col
                if 'SSC-A' in col or 'SSC.A' in col:
                    ssc_a_col = col
            
            if fsc_a_col and ssc_a_col:
                gate_mask = (
                    (fcs_data[:, fsc_a_col] >= 0) & 
                    (fcs_data[:, fsc_a_col] <= 10000000) &
                    (fcs_data[:, ssc_a_col] >= 0) & 
                    (fcs_data[:, ssc_a_col] <= 1500000)
                )
                gated_data = fcs_data[gate_mask]
            else:
                gated_data = fcs_data
            
            fluorescence_values = gated_data[:, channel_to_use]
        else:
            if isinstance(fcs_data, np.ndarray):
                data_array = fcs_data
            else:
                data_array = np.array(fcs_data)
            
            fsc_a_idx = None
            ssc_a_idx = None
            
            for i, col in enumerate(channel_names):
                if 'FSC-A' in col or 'FSC.A' in col:
                    fsc_a_idx = i
                if 'SSC-A' in col or 'SSC.A' in col:
                    ssc_a_idx = i
            
            if fsc_a_idx is not None and ssc_a_idx is not None:
                gate_mask = (
                    (data_array[:, fsc_a_idx] >= 0) & 
                    (data_array[:, fsc_a_idx] <= 10000000) &
                    (data_array[:, ssc_a_idx] >= 0) & 
                    (data_array[:, ssc_a_idx] <= 1500000)
                )
                gated_data = data_array[gate_mask]
            else:
                gated_data = data_array
            
            fluorescence_values = gated_data[:, channel_index]
        
        # 移除负值和NaN
        valid_mask = (fluorescence_values >= 0) & (~np.isnan(fluorescence_values))
        fluorescence_values = fluorescence_values[valid_mask]
        
        if len(fluorescence_values) == 0:
            return None
        
        # 创建DataFrame
        sample_name = os.path.splitext(os.path.basename(file_path))[0]
        df = pd.DataFrame({
            fluorescence_channel: fluorescence_values,
            'Sample': sample_name
        })
        
        return df
        
    except Exception as e:
        print(f"错误: 处理文件 {os.path.basename(file_path)} 时出错: {str(e)}")
        return None

def generate_plot(work_dir, file_configs, fluorescence_channel, output_file):
    """
    生成山脊图
    
    参数:
        work_dir: 工作目录
        file_configs: 文件配置列表，每个元素为 {'file_path': str, 'display_name': str, 'color': str, 'order': int}
        fluorescence_channel: 荧光通道名
        output_file: 输出文件路径
    """
    if not MATPLOTLIB_AVAILABLE:
        return False
    
    # 按order排序
    file_configs_sorted = sorted(file_configs, key=lambda x: x['order'])
    
    # 处理所有文件
    all_data_list = []
    for config in file_configs_sorted:
        df = process_fcs_file(config['file_path'], fluorescence_channel)
        if df is not None:
            # 更新样本名为显示名称
            df['Sample'] = config['display_name']
            all_data_list.append(df)
    
    if not all_data_list:
        return False
    
    # 合并所有数据
    combined_data = pd.concat(all_data_list, ignore_index=True)
    
    samples = combined_data['Sample'].unique()
    n_samples = len(samples)
    
    # 创建颜色映射
    sample_colors = {}
    for config in file_configs_sorted:
        if config['display_name'] in samples:
            sample_colors[config['display_name']] = config['color']
    
    # 计算合理的数据范围（使用对数坐标）
    all_values = combined_data[fluorescence_channel].values
    all_values = all_values[all_values > 0]
    
    if len(all_values) == 0:
        return False
    
    min_val = all_values.min()
    max_val = all_values.max()
    
    log_min = np.floor(np.log10(min_val))
    log_max = np.ceil(np.log10(max_val))
    
    if log_min < 1:
        log_min = 1
    x_min = 10 ** log_min
    x_max = 10 ** log_max
    
    if max_val < x_max * 0.8:
        x_max = max_val * 1.5
        log_max = np.ceil(np.log10(x_max))
    
    # 固定高度配置
    fixed_group_height = 1000
    height_scale = 0.5
    overlap_ratio = 0.75
    
    # 计算每个样本的基线
    baselines = [0]
    for i in range(1, len(samples)):
        prev_height_scaled = fixed_group_height * height_scale
        baseline = baselines[i-1] + prev_height_scaled * overlap_ratio
        baselines.append(baseline)
    
    last_height_scaled = fixed_group_height * height_scale
    total_height = baselines[-1] + last_height_scaled
    
    # 根据样本数量动态调整图片高度
    fig_height = max(4, n_samples * 0.8)
    fig_width = 4  # 横坐标长度（3英寸）
    
    # 创建图
    fig, ax = plt.subplots(1, 1, figsize=(fig_width, fig_height))
    
    # 反转绘制顺序
    for i, sample in enumerate(reversed(samples)):
        original_i = len(samples) - 1 - i
        sample_data = combined_data[combined_data['Sample'] == sample]
        sample_values = sample_data[sample_data[fluorescence_channel] > 0][fluorescence_channel].values
        
        if len(sample_values) == 0:
            continue
        
        # 使用直方图计算counts，然后轻微平滑
        bins = np.logspace(log_min, log_max, 100)
        hist, bin_edges = np.histogram(sample_values, bins=bins, density=False)
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
        hist_max = hist.max() if len(hist) > 0 else 0
        
        try:
            from scipy.interpolate import interp1d
            from scipy.ndimage import gaussian_filter1d
            
            non_zero_mask = hist > 0
            if np.sum(non_zero_mask) > 1:
                hist_smooth = gaussian_filter1d(hist, sigma=0.5)
                x_range = np.logspace(log_min, log_max, 300)
                interp_func = interp1d(bin_centers, hist_smooth, 
                                      kind='linear', 
                                      bounds_error=False, 
                                      fill_value=0)
                y_range = interp_func(x_range)
                y_range = np.maximum(y_range, 0)
            else:
                x_range = bin_centers
                y_range = hist
            
            y_max = y_range.max()
            
        except Exception as e:
            x_range = bin_centers
            y_range = hist
            y_max = hist_max
        
        # 获取颜色
        if sample in sample_colors:
            color = sample_colors[sample]
        else:
            color = plt.cm.viridis(original_i / max(1, n_samples - 1))
        
        # 归一化到固定高度
        y_max_current = y_range.max()
        if y_max_current > 0:
            y_range_normalized = (y_range / y_max_current) * fixed_group_height
            y_range = y_range_normalized * height_scale
        else:
            y_range = y_range * height_scale
        
        baseline = baselines[original_i]
        y_range_stacked = y_range + baseline
        
        # 填充曲线下方区域
        ax.fill_between(x_range, baseline, y_range_stacked, 
                        alpha=0.6, color=color, 
                        edgecolor='none',
                        where=y_range > 0)
        
        # 绘制平滑曲线
        ax.plot(x_range, y_range_stacked, 
                color=color, linewidth=2, 
                alpha=1.0)
    
    # 设置对数坐标
    ax.set_xscale('log')
    ax.set_xlim(x_min, x_max)
    
    # 设置x轴刻度
    major_ticks = [10**i for i in range(int(log_min), int(log_max) + 1)]
    ax.set_xticks(major_ticks)
    
    minor_ticks = []
    for j in range(int(log_min), int(log_max)):
        minor_ticks.extend([10**(j + 0.5)])
    if len(minor_ticks) > 0:
        ax.set_xticks(minor_ticks, minor=True)
    
    # 格式化x轴标签
    from matplotlib.ticker import FuncFormatter
    def log_formatter(x, pos):
        if x <= 0:
            return '0'
        exp = int(np.log10(x))
        if exp == 0:
            return '1'
        return f'10$^{{{exp}}}$'
    
    ax.xaxis.set_major_formatter(FuncFormatter(log_formatter))
    
    # 设置y轴范围
    if total_height > 0:
        ax.set_ylim(0, total_height * 1.25)
    
    # 设置字体
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['font.weight'] = 'bold'
    
    # 设置标签
    ax.set_ylabel('Count', fontsize=14, fontweight='bold', family='Arial')
    ax.set_xlabel(f'MFI ({fluorescence_channel})', 
                 fontsize=14, fontweight='bold', family='Arial')
    
    # 设置y轴刻度（间隔为300）
    if total_height > 0:
        y_ticks = np.arange(0, total_height * 1.25 + 300, 300)
        ax.set_yticks(y_ticks)
        ax.set_yticklabels([f'{int(t)}' for t in y_ticks], 
                          fontsize=11, fontweight='bold', family='Arial')
    
    # 加粗轴线
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(2)
    ax.spines['bottom'].set_linewidth(2)
    
    # 加粗刻度线
    ax.tick_params(axis='both', which='major', width=2, length=6)
    ax.tick_params(axis='both', which='minor', width=1.5, length=4)
    
    # 去掉网格线
    ax.grid(False)
    
    # 设置x轴刻度标签字体
    for label in ax.get_xticklabels():
        label.set_fontweight('bold')
        label.set_fontfamily('Arial')
        label.set_fontsize(11)
    
    # 创建图例
    legend_elements = []
    for i, sample in enumerate(samples):
        if sample in sample_colors:
            color = sample_colors[sample]
        else:
            color = plt.cm.viridis(i / max(1, n_samples - 1))
        from matplotlib.patches import Rectangle
        legend_patch = Rectangle((0, 0), 1, 1, 
                                 facecolor=color, 
                                 edgecolor=color, 
                                 linewidth=2,
                                 alpha=0.6)
        legend_elements.append((legend_patch, sample))
    
    handles = [elem[0] for elem in legend_elements]
    labels = [elem[1] for elem in legend_elements]
    legend = ax.legend(handles, labels, loc='upper right', 
                      fontsize=11,
                      frameon=False)
    for text in legend.get_texts():
        text.set_fontweight('bold')
        text.set_fontfamily('Arial')
    
    plt.tight_layout()
    
    # 保存图片
    plt.savefig(output_file, dpi=600, bbox_inches='tight')
    plt.close()
    
    return True

class CustomCheckBox(QCheckBox):
    """自定义复选框，显示对号"""
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setFixedSize(22, 22)
    
    def paintEvent(self, event):
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        
        # 绘制复选框边框
        rect = QRect(1, 1, 20, 20)
        if self.isChecked():
            painter.setBrush(QColor("#4CAF50"))
            painter.setPen(QPen(QColor("#4CAF50"), 2))
        else:
            painter.setBrush(QColor("white"))
            painter.setPen(QPen(QColor("#666"), 2))
        
        painter.drawRoundedRect(rect, 4, 4)
        
        # 如果选中，绘制对号
        if self.isChecked():
            pen = QPen(QColor("white"), 2.5)
            pen.setCapStyle(Qt.PenCapStyle.RoundCap)
            pen.setJoinStyle(Qt.PenJoinStyle.RoundJoin)
            painter.setPen(pen)
            # 绘制对号的两条线
            painter.drawLine(6, 11, 9, 14)
            painter.drawLine(9, 14, 16, 7)

# EditDialog类已移除，现在所有编辑都在主界面完成

class FCSRidgePlotGUI(QMainWindow):
    def __init__(self):
        super().__init__()
        self.work_dir = ""
        self.fcs_files = []
        self.file_configs = []  # 每个元素为 {'file_path': str, 'display_name': str, 'color': str, 'order': int, 'selected': bool}
        self.fluorescence_channel = None
        
        self.init_ui()
        self.apply_styles()
    
    def init_ui(self):
        self.setWindowTitle("FCS山脊图生成工具")
        self.setGeometry(100, 100, 1200, 800)  # 恢复原来的宽度
        
        # 中央部件
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        
        # 主布局
        main_layout = QVBoxLayout()
        main_layout.setSpacing(15)
        main_layout.setContentsMargins(15, 15, 15, 15)
        
        # 顶部：路径选择
        top_layout = QHBoxLayout()
        top_layout.addWidget(QLabel("输入路径:"))
        self.path_edit = QLineEdit()
        self.path_edit.setPlaceholderText("请选择包含FCS文件的目录...")
        top_layout.addWidget(self.path_edit)
        
        self.browse_btn = QPushButton("浏览")
        self.browse_btn.clicked.connect(self.browse_path)
        top_layout.addWidget(self.browse_btn)
        
        self.detect_btn = QPushButton("检测文件")
        self.detect_btn.clicked.connect(self.detect_files)
        top_layout.addWidget(self.detect_btn)
        
        main_layout.addLayout(top_layout)
        
        # 中间：文件列表表格
        self.table = QTableWidget()
        self.table.setColumnCount(6)  # 选中, 文件名, 显示名称, 颜色, MFI Mean, MFI Median
        self.table.setHorizontalHeaderLabels(['选中', '文件名', '显示名称', '颜色', 'MFI Mean', 'MFI Median'])
        self.table.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeMode.ResizeToContents)
        self.table.horizontalHeader().setSectionResizeMode(1, QHeaderView.ResizeMode.ResizeToContents)
        self.table.setColumnWidth(1, 200)  # 文件名列固定宽度（恢复原宽度）
        self.table.horizontalHeader().setSectionResizeMode(2, QHeaderView.ResizeMode.Stretch)
        self.table.setColumnWidth(2, 200)  # 显示名称列最小宽度（恢复原宽度）
        self.table.horizontalHeader().setSectionResizeMode(3, QHeaderView.ResizeMode.Fixed)
        self.table.setColumnWidth(3, 280)  # 颜色列固定宽度（恢复原宽度）
        self.table.horizontalHeader().setSectionResizeMode(4, QHeaderView.ResizeMode.ResizeToContents)
        self.table.horizontalHeader().setSectionResizeMode(5, QHeaderView.ResizeMode.ResizeToContents)
        self.table.setSelectionBehavior(QTableWidget.SelectionBehavior.SelectRows)
        self.table.setSelectionMode(QTableWidget.SelectionMode.SingleSelection)
        self.table.setAlternatingRowColors(True)
        # 设置行高为原来的5倍（约175px）
        self.table.verticalHeader().setDefaultSectionSize(100)
        self.table.verticalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Fixed)
        # 允许编辑显示名称列
        self.table.setEditTriggers(QTableWidget.EditTrigger.DoubleClicked | QTableWidget.EditTrigger.SelectedClicked)
        # 连接编辑完成信号
        self.table.itemChanged.connect(self.on_item_changed)
        
        main_layout.addWidget(self.table)
        
        # 底部：按钮
        bottom_layout = QHBoxLayout()
        bottom_layout.addStretch()
        
        self.csv_btn = QPushButton("生成CSV报告")
        self.csv_btn.clicked.connect(self.generate_csv_report)
        bottom_layout.addWidget(self.csv_btn)
        
        self.plot_btn = QPushButton("生成图像")
        self.plot_btn.clicked.connect(self.generate_plot)
        bottom_layout.addWidget(self.plot_btn)
        
        self.move_up_btn = QPushButton("上移")
        self.move_up_btn.clicked.connect(self.move_up)
        bottom_layout.addWidget(self.move_up_btn)
        
        self.move_down_btn = QPushButton("下移")
        self.move_down_btn.clicked.connect(self.move_down)
        bottom_layout.addWidget(self.move_down_btn)
        
        main_layout.addLayout(bottom_layout)
        
        central_widget.setLayout(main_layout)
    
    def apply_styles(self):
        """应用美化样式"""
        self.setStyleSheet("""
            QMainWindow {
                background: qlineargradient(x1:0, y1:0, x2:0, y2:1,
                    stop:0 #f0f4f8, stop:1 #e1e8ed);
            }
            QLabel {
                font-family: "Arial", "SimHei", "Microsoft YaHei";
                font-size: 14px;
                font-weight: bold;
                color: #2c3e50;
            }
            QLineEdit {
                padding: 10px;
                border: 2px solid #bdc3c7;
                border-radius: 6px;
                font-family: "Arial", "SimHei", "Microsoft YaHei";
                font-size: 14px;
                font-weight: bold;
                background-color: white;
            }
            QLineEdit:focus {
                border: 2px solid #3498db;
                background-color: #ecf0f1;
            }
            QPushButton {
                padding: 12px 24px;
                border: none;
                border-radius: 8px;
                font-family: "Arial", "SimHei", "Microsoft YaHei";
                font-size: 14px;
                font-weight: bold;
                background: qlineargradient(x1:0, y1:0, x2:0, y2:1,
                    stop:0 #3498db, stop:1 #2980b9);
                color: white;
                min-width: 120px;
            }
            QPushButton:hover {
                background: qlineargradient(x1:0, y1:0, x2:0, y2:1,
                    stop:0 #5dade2, stop:1 #3498db);
            }
            QPushButton:pressed {
                background: qlineargradient(x1:0, y1:0, x2:0, y2:1,
                    stop:0 #2980b9, stop:1 #1f618d);
            }
            QTableWidget {
                border: 3px solid #bdc3c7;
                border-radius: 10px;
                background-color: white;
                gridline-color: #ecf0f1;
                font-family: "Arial", "SimHei", "Microsoft YaHei";
                font-size: 16px;
                font-weight: bold;
                selection-background-color: #d5e8f4;
            }
            QTableWidget::item {
                padding: 15px;
                border: none;
                min-height: 35px;
                font-family: "Arial", "SimHei", "Microsoft YaHei";
                font-weight: bold;
            }
            QTableWidget::item:selected {
                background-color: #d5e8f4;
                color: #2c3e50;
            }
            QTableWidget::item:focus {
                background-color: #d5e8f4;
                outline: none;
            }
            QTableWidget::item:editable {
                font-size: 18px;
                font-weight: bold;
                min-height: 35px;
            }
            QHeaderView::section {
                background: qlineargradient(x1:0, y1:0, x2:0, y2:1,
                    stop:0 #34495e, stop:1 #2c3e50);
                color: white;
                padding: 18px;
                border: none;
                font-family: "Arial", "SimHei", "Microsoft YaHei";
                font-weight: bold;
                font-size: 18px;
                border-radius: 0px;
            }
            QHeaderView::section:first {
                border-top-left-radius: 8px;
            }
            QHeaderView::section:last {
                border-top-right-radius: 8px;
            }
            QCheckBox {
                spacing: 8px;
                font-family: "Arial", "SimHei", "Microsoft YaHei";
                font-weight: bold;
            }
            QCheckBox::indicator {
                width: 24px;
                height: 24px;
                border: 3px solid #7f8c8d;
                border-radius: 5px;
                background-color: white;
            }
            QCheckBox::indicator:checked {
                background-color: #27ae60;
                border: 3px solid #27ae60;
            }
            QCheckBox::indicator:hover {
                border: 3px solid #27ae60;
            }
            QDialog {
                background-color: white;
            }
            QSpinBox {
                padding: 8px;
                border: 2px solid #bdc3c7;
                border-radius: 5px;
                font-family: "Arial", "SimHei", "Microsoft YaHei";
                font-size: 14px;
                font-weight: bold;
            }
            QPushButton#colorBtn {
                min-width: 60px;
                max-width: 60px;
                min-height: 35px;
                max-height: 35px;
                border: 3px solid #bdc3c7;
                border-radius: 6px;
                padding: 0px;
            }
            QPushButton#colorBtn:hover {
                border: 3px solid #3498db;
            }
            QLineEdit#colorHex {
                padding: 10px;
                border: 2px solid #bdc3c7;
                border-radius: 6px;
                font-family: "Courier New", "Consolas", "Monaco", "Menlo", monospace;
                font-size: 16px;
                font-weight: bold;
                min-width: 220px;
                min-height: 30px;
            }
            QLineEdit#colorHex:focus {
                border: 2px solid #3498db;
                background-color: #ecf0f1;
            }
        """)
    
    def browse_path(self):
        path = QFileDialog.getExistingDirectory(self, "选择FCS文件目录")
        if path:
            self.path_edit.setText(path)
            self.work_dir = path
    
    def detect_files(self):
        self.work_dir = self.path_edit.text()
        if not self.work_dir or not os.path.isdir(self.work_dir):
            QMessageBox.critical(self, "错误", "请先选择有效的目录")
            return
        
        # 获取所有FCS文件
        self.fcs_files = glob.glob(os.path.join(self.work_dir, '*.fcs')) + \
                         glob.glob(os.path.join(self.work_dir, '*.FCS'))
        
        if not self.fcs_files:
            QMessageBox.warning(self, "警告", f"在目录 {self.work_dir} 下没有找到FCS文件")
            return
        
        # 从第一个文件识别通道
        try:
            first_fcs, channel_names = read_fcs_file(self.fcs_files[0])
            if first_fcs is None:
                QMessageBox.critical(self, "错误", "无法读取第一个文件")
                return
            self.fluorescence_channel = detect_fitc_channel(channel_names)
            
            if self.fluorescence_channel is None:
                QMessageBox.critical(self, "错误", "无法识别荧光通道")
                return
            
        except Exception as e:
            QMessageBox.critical(self, "错误", f"读取第一个文件时出错: {str(e)}")
            return
        
        # 显示进度对话框
        progress = QProgressDialog("正在处理文件...", "取消", 0, len(self.fcs_files), self)
        progress.setWindowModality(Qt.WindowModality.WindowModal)
        progress.setMinimumDuration(0)
        
        # 清空现有数据
        self.table.setRowCount(0)
        self.file_configs = []
        
        # 处理所有文件并计算MFI
        for i, fcs_file in enumerate(self.fcs_files):
            if progress.wasCanceled():
                break
            
            progress.setValue(i)
            progress.setLabelText(f"正在处理: {os.path.basename(fcs_file)}")
            QApplication.processEvents()
            
            file_name = os.path.basename(fcs_file)
            display_name = os.path.splitext(file_name)[0]
            
            # 计算MFI mean和median值
            mfi_result = process_fcs_file_for_mfi(fcs_file, self.fluorescence_channel)
            if mfi_result is not None:
                mfi_mean, mfi_median = mfi_result
            else:
                mfi_mean, mfi_median = None, None
            mfi_mean_str = f"{mfi_mean:.2f}" if mfi_mean is not None else "N/A"
            mfi_median_str = f"{mfi_median:.2f}" if mfi_median is not None else "N/A"
            
            # 创建配置
            config = {
                'file_path': fcs_file,
                'display_name': display_name,
                'color': '#%02x%02x%02x' % tuple(int(c * 255) for c in plt.cm.viridis(i / max(1, len(self.fcs_files) - 1))[:3]),
                'order': i,
                'selected': True,
                'mfi_mean': mfi_mean,
                'mfi_median': mfi_median
            }
            self.file_configs.append(config)
            
            # 添加到表格
            row = self.table.rowCount()
            self.table.insertRow(row)
            
            # 选中复选框（自定义样式）
            checkbox = CustomCheckBox()
            checkbox.setChecked(True)
            file_path = fcs_file
            checkbox.stateChanged.connect(lambda state, fp=file_path: self.update_selected_by_path(fp, state == Qt.CheckState.Checked.value))
            self.table.setCellWidget(row, 0, checkbox)
            
            # 文件名列（只读）
            file_name_item = QTableWidgetItem(file_name)
            file_name_item.setFont(QFont("Arial", 16, QFont.Weight.Bold))
            self.table.setItem(row, 1, file_name_item)
            self.table.item(row, 1).setFlags(self.table.item(row, 1).flags() & ~Qt.ItemFlag.ItemIsEditable)
            
            # 显示名称列（可编辑，增大字体）
            name_item = QTableWidgetItem(display_name)
            name_item.setFont(QFont("Arial", 18, QFont.Weight.Bold))
            self.table.setItem(row, 2, name_item)
            
            # 颜色列：颜色块按钮 + Hex编辑框
            color_widget = QWidget()
            color_layout = QHBoxLayout()
            color_layout.setContentsMargins(5, 5, 5, 5)
            color_layout.setSpacing(8)
            
            color_btn = QPushButton()
            color_btn.setObjectName("colorBtn")
            color_btn.setFixedSize(60, 35)  # 颜色按钮高度（恢复原大小）
            try:
                color = QColor(config['color'])
                if color.isValid():
                    color_btn.setStyleSheet(f"background-color: {color.name()};")
            except:
                pass
            color_btn.clicked.connect(lambda checked, fp=file_path: self.choose_color_for_file(fp))
            color_layout.addWidget(color_btn)
            
            color_hex = QLineEdit(config['color'])
            color_hex.setObjectName("colorHex")
            color_hex.setPlaceholderText("#RRGGBB")
            color_hex.setFont(QFont("Courier New", 16, QFont.Weight.Bold))
            color_hex.setMinimumWidth(220)  # Hex编辑框宽度
            color_hex.setMinimumHeight(30)  # Hex编辑框高度（恢复原大小）
            color_hex.textChanged.connect(lambda text, fp=file_path: self.update_color_by_hex(fp, text))
            color_layout.addWidget(color_hex)
            
            color_widget.setLayout(color_layout)
            self.table.setCellWidget(row, 3, color_widget)
            
            # MFI Mean列（只读）
            mfi_mean_item = QTableWidgetItem(mfi_mean_str)
            mfi_mean_item.setFont(QFont("Arial", 16, QFont.Weight.Bold))
            mfi_mean_item.setFlags(mfi_mean_item.flags() & ~Qt.ItemFlag.ItemIsEditable)
            self.table.setItem(row, 4, mfi_mean_item)
            
            # MFI Median列（只读）
            mfi_median_item = QTableWidgetItem(mfi_median_str)
            mfi_median_item.setFont(QFont("Arial", 16, QFont.Weight.Bold))
            mfi_median_item.setFlags(mfi_median_item.flags() & ~Qt.ItemFlag.ItemIsEditable)
            self.table.setItem(row, 5, mfi_median_item)
        
        progress.setValue(len(self.fcs_files))
        QMessageBox.information(self, "成功", f"检测到 {len(self.fcs_files)} 个FCS文件\n使用通道: {self.fluorescence_channel}")
    
    def update_selected_by_path(self, file_path, selected):
        """根据文件路径更新选中状态"""
        config = next((c for c in self.file_configs if c['file_path'] == file_path), None)
        if config:
            config['selected'] = selected
    
    def on_item_changed(self, item):
        """处理表格项编辑"""
        row = item.row()
        column = item.column()
        
        if column == 2:  # 显示名称列
            file_name_item = self.table.item(row, 1)
            if file_name_item:
                file_name = file_name_item.text()
                config = next((c for c in self.file_configs if os.path.basename(c['file_path']) == file_name), None)
                if config:
                    config['display_name'] = item.text()
    
    def move_up(self):
        """上移选中项"""
        current_row = self.table.currentRow()
        if current_row < 0:
            QMessageBox.warning(self, "警告", "请先选择要移动的项目")
            return
        
        if current_row == 0:
            return
        
        # 找到对应的配置
        file_name = self.table.item(current_row, 1).text()
        config = next((c for c in self.file_configs if os.path.basename(c['file_path']) == file_name), None)
        
        if config and config['order'] > 0:
            # 找到前一个项目
            prev_config = next((c for c in self.file_configs if c['order'] == config['order'] - 1), None)
            if prev_config:
                config['order'], prev_config['order'] = prev_config['order'], config['order']
                self.refresh_table()
                self.table.selectRow(current_row - 1)
    
    def move_down(self):
        """下移选中项"""
        current_row = self.table.currentRow()
        if current_row < 0:
            QMessageBox.warning(self, "警告", "请先选择要移动的项目")
            return
        
        if current_row >= self.table.rowCount() - 1:
            return
        
        # 找到对应的配置
        file_name = self.table.item(current_row, 1).text()
        config = next((c for c in self.file_configs if os.path.basename(c['file_path']) == file_name), None)
        
        max_order = max(c['order'] for c in self.file_configs) if self.file_configs else 0
        if config and config['order'] < max_order:
            # 找到后一个项目
            next_config = next((c for c in self.file_configs if c['order'] == config['order'] + 1), None)
            if next_config:
                config['order'], next_config['order'] = next_config['order'], config['order']
                self.refresh_table()
                self.table.selectRow(current_row + 1)
    
    def choose_color_for_file(self, file_path):
        """为文件选择颜色"""
        config = next((c for c in self.file_configs if c['file_path'] == file_path), None)
        if not config:
            return
        
        try:
            initial_color = QColor(config['color'])
        except:
            initial_color = QColor('#000000')
        
        color = QColorDialog.getColor(initial_color, self, "选择颜色")
        if color.isValid():
            config['color'] = color.name()
            # 更新颜色控件
            self.update_color_widget(file_path, color.name())
    
    def update_color_by_hex(self, file_path, hex_text):
        """通过Hex文本更新颜色"""
        config = next((c for c in self.file_configs if c['file_path'] == file_path), None)
        if not config:
            return
        
        # 验证Hex格式
        if hex_text.startswith('#'):
            try:
                color = QColor(hex_text)
                if color.isValid():
                    config['color'] = color.name()
                    # 更新颜色按钮
                    self.update_color_widget(file_path, color.name())
            except:
                pass
    
    def update_color_widget(self, file_path, color_hex):
        """更新颜色控件显示"""
        for row in range(self.table.rowCount()):
            file_name_item = self.table.item(row, 1)
            if not file_name_item:
                continue
            file_name = file_name_item.text()
            config = next((c for c in self.file_configs if os.path.basename(c['file_path']) == file_name), None)
            if config and config['file_path'] == file_path:
                color_widget = self.table.cellWidget(row, 3)
                if color_widget:
                    # 更新颜色按钮
                    color_btn = color_widget.findChild(QPushButton)
                    if color_btn:
                        try:
                            color = QColor(color_hex)
                            if color.isValid():
                                color_btn.setStyleSheet(f"background-color: {color.name()}; border: 2px solid #ddd; border-radius: 4px;")
                        except:
                            pass
                    # 更新Hex编辑框
                    color_hex_edit = color_widget.findChild(QLineEdit)
                    if color_hex_edit:
                        # 临时断开信号，避免循环更新
                        color_hex_edit.blockSignals(True)
                        color_hex_edit.setText(color_hex)
                        color_hex_edit.blockSignals(False)
                break
    
    def refresh_table(self):
        """刷新表格显示（按order排序）"""
        # 按order排序
        sorted_configs = sorted(self.file_configs, key=lambda x: x['order'])
        
        # 临时断开信号，避免刷新时触发
        self.table.itemChanged.disconnect()
        
        self.table.setRowCount(0)
        
        for i, config in enumerate(sorted_configs):
            row = self.table.rowCount()
            self.table.insertRow(row)
            
            file_name = os.path.basename(config['file_path'])
            mfi_mean_str = f"{config['mfi_mean']:.2f}" if config['mfi_mean'] is not None else "N/A"
            mfi_median_str = f"{config['mfi_median']:.2f}" if config['mfi_median'] is not None else "N/A"
            
            # 更新order为当前索引
            config['order'] = i
            
            # 选中复选框（自定义样式）
            checkbox = CustomCheckBox()
            checkbox.setChecked(config['selected'])
            file_path = config['file_path']
            checkbox.stateChanged.connect(lambda state, fp=file_path: self.update_selected_by_path(fp, state == Qt.CheckState.Checked.value))
            self.table.setCellWidget(row, 0, checkbox)
            
            # 文件名列（只读）
            file_name_item = QTableWidgetItem(file_name)
            file_name_item.setFont(QFont("Arial", 16, QFont.Weight.Bold))
            self.table.setItem(row, 1, file_name_item)
            self.table.item(row, 1).setFlags(self.table.item(row, 1).flags() & ~Qt.ItemFlag.ItemIsEditable)
            
            # 显示名称列（可编辑，增大字体）
            name_item = QTableWidgetItem(config['display_name'])
            name_item.setFont(QFont("Arial", 18, QFont.Weight.Bold))
            self.table.setItem(row, 2, name_item)
            
            # 颜色列：颜色块按钮 + Hex编辑框
            color_widget = QWidget()
            color_layout = QHBoxLayout()
            color_layout.setContentsMargins(5, 5, 5, 5)
            color_layout.setSpacing(8)
            
            color_btn = QPushButton()
            color_btn.setObjectName("colorBtn")
            color_btn.setFixedSize(60, 35)  # 颜色按钮高度（恢复原大小）
            try:
                color = QColor(config['color'])
                if color.isValid():
                    color_btn.setStyleSheet(f"background-color: {color.name()};")
            except:
                pass
            color_btn.clicked.connect(lambda checked, fp=file_path: self.choose_color_for_file(fp))
            color_layout.addWidget(color_btn)
            
            color_hex = QLineEdit(config['color'])
            color_hex.setObjectName("colorHex")
            color_hex.setPlaceholderText("#RRGGBB")
            color_hex.setFont(QFont("Courier New", 16, QFont.Weight.Bold))
            color_hex.setMinimumWidth(220)  # Hex编辑框宽度
            color_hex.setMinimumHeight(30)  # Hex编辑框高度（恢复原大小）
            color_hex.textChanged.connect(lambda text, fp=file_path: self.update_color_by_hex(fp, text))
            color_layout.addWidget(color_hex)
            
            color_widget.setLayout(color_layout)
            self.table.setCellWidget(row, 3, color_widget)
            
            # MFI Mean列（只读）
            mfi_mean_item = QTableWidgetItem(mfi_mean_str)
            mfi_mean_item.setFont(QFont("Arial", 16, QFont.Weight.Bold))
            mfi_mean_item.setFlags(mfi_mean_item.flags() & ~Qt.ItemFlag.ItemIsEditable)
            self.table.setItem(row, 4, mfi_mean_item)
            
            # MFI Median列（只读）
            mfi_median_item = QTableWidgetItem(mfi_median_str)
            mfi_median_item.setFont(QFont("Arial", 16, QFont.Weight.Bold))
            mfi_median_item.setFlags(mfi_median_item.flags() & ~Qt.ItemFlag.ItemIsEditable)
            self.table.setItem(row, 5, mfi_median_item)
        
        # 重新连接信号
        self.table.itemChanged.connect(self.on_item_changed)
    
    
    
    def generate_csv_report(self):
        """生成CSV报告"""
        if not self.file_configs:
            QMessageBox.warning(self, "警告", "没有文件数据，请先检测文件")
            return
        
        if not self.fluorescence_channel:
            QMessageBox.critical(self, "错误", "未识别到荧光通道")
            return
        
        # 选择保存位置
        output_file, _ = QFileDialog.getSaveFileName(
            self, "保存CSV报告", "", "CSV文件 (*.csv);;所有文件 (*.*)"
        )
        
        if not output_file:
            return
        
        # 生成报告
        report_data = []
        for config in sorted(self.file_configs, key=lambda x: x['order']):
            report_data.append({
                '文件名': os.path.basename(config['file_path']),
                '显示名称': config['display_name'],
                f'MFI Mean ({self.fluorescence_channel})': config['mfi_mean'] if config['mfi_mean'] is not None else 'N/A',
                f'MFI Median ({self.fluorescence_channel})': config['mfi_median'] if config['mfi_median'] is not None else 'N/A'
            })
        
        df = pd.DataFrame(report_data)
        df.to_csv(output_file, index=False, encoding='utf-8-sig')
        
        QMessageBox.information(self, "成功", f"CSV报告已保存到:\n{output_file}")
    
    def generate_plot(self):
        """生成图像"""
        if not self.file_configs:
            QMessageBox.warning(self, "警告", "没有文件数据，请先检测文件")
            return
        
        # 获取选中的文件
        selected_configs = [c for c in self.file_configs if c['selected']]
        
        if not selected_configs:
            QMessageBox.warning(self, "警告", "请至少选择一个文件")
            return
        
        if not self.fluorescence_channel:
            QMessageBox.critical(self, "错误", "未识别到荧光通道")
            return
        
        # 选择保存位置
        output_file, _ = QFileDialog.getSaveFileName(
            self, "保存图像", "", "PNG文件 (*.png);;所有文件 (*.*)"
        )
        
        if not output_file:
            return
        
        # 生成图像
        try:
            success = generate_plot(self.work_dir, selected_configs, self.fluorescence_channel, output_file)
            if success:
                QMessageBox.information(self, "成功", f"图像已保存到:\n{output_file}")
            else:
                QMessageBox.critical(self, "错误", "图像生成失败")
        except Exception as e:
            QMessageBox.critical(self, "错误", f"生成图像时出错:\n{str(e)}")

def main():
    app = QApplication(sys.argv)
    app.setStyle('Fusion')  # 使用Fusion样式，更现代
    
    window = FCSRidgePlotGUI()
    window.show()
    
    sys.exit(app.exec())

if __name__ == '__main__':
    main()