# FCS Ridge Plot GUI

PyQt6-based GUI application for reading FCS flow cytometry files, automatically detecting fluorescence channels, calculating MFI (mean/median) statistics, and generating customizable ridge plots with color selection and file ordering.

## Features

- **FCS File Reading**: Supports multiple FCS parsing libraries (fcsparser, flowcal, flowio, flowutils)
- **Automatic Channel Detection**: Automatically detects FITC/fluorescence channels from FCS files
- **MFI Statistics**: Calculates Mean Fluorescence Intensity (MFI) mean and median values
- **Customizable Ridge Plots**: Generate publication-quality ridge plots with:
  - Custom color selection for each sample
  - File ordering control (move up/down)
  - Display name editing
  - Log-scale x-axis
- **CSV Report Generation**: Export MFI statistics to CSV format
- **User-Friendly Interface**: Modern PyQt6 GUI with intuitive controls

## Requirements

### Required Dependencies
- Python 3.7+
- PyQt6
- numpy
- pandas
- matplotlib
- seaborn

### Optional FCS Parsing Libraries (at least one required)
- `fcsparser` (recommended)
- `flowcal`
- `flowio`
- `flowutils`

### Optional Dependencies
- `scipy` (for advanced smoothing in plots)

## Installation

1. Install required dependencies:
pip install PyQt6 numpy pandas matplotlib seaborn2. Install at least one FCS parsing library (recommended: fcsparser):
pip install fcsparser3. (Optional) Install scipy for enhanced plot smoothing:sh
pip install scipy## Usage

1. Run the application:
python fcs_rigde_plot_gui.py2. **Select Directory**: Click "浏览" (Browse) to select a directory containing FCS files

3. **Detect Files**: Click "检测文件" (Detect Files) to scan for FCS files and automatically detect fluorescence channels

4. **Configure Files**:
   - Check/uncheck files to include/exclude from plots
   - Edit display names by double-clicking the "显示名称" (Display Name) column
   - Customize colors by clicking the color button or editing the hex code
   - Reorder files using "上移" (Move Up) and "下移" (Move Down) buttons

5. **Generate Output**:
   - Click "生成图像" (Generate Plot) to create a ridge plot (PNG format, 600 DPI)
   - Click "生成CSV报告" (Generate CSV Report) to export MFI statistics

## Output

- **Ridge Plot**: High-resolution PNG image (600 DPI) with:
  - Log-scale x-axis for MFI values
  - Overlapping density distributions
  - Custom colors for each sample
  - Legend with sample names
  - Publication-ready formatting

- **CSV Report**: Contains:
  - File names
  - Display names
  - MFI Mean values
  - MFI Median values

## Channel Detection

The application automatically detects fluorescence channels by:
1. Excluding scatter chan
