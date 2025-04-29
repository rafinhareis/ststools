<p align="center">
  <img src="logo.png" alt="STSTools Logo" width="180"/>
</p>

<h1 align="center">STSTools</h1>
<p align="center"><strong>Data Analysis for STS Measurements from Nanosurf, Omicron, and Nanonis</strong></p>

---

## 🧪 What is STSTools?

**STSTools** is a Python library for analyzing and visualizing Scanning Tunneling Spectroscopy (STS) data from:

- 🧬 **Nanosurf**: `.nid` files
- ⚛️ **Omicron v6.4.0**: `.txt`, `.ibw` (IGOR grid maps)
- 🧲 **Nanonis**: `.dat` files

Built to help researchers process, inspect, and interpret I(V) spectra and grid maps through an intuitive Jupyter-based interface.

> ⚠️ *Important: Do not rename exported `.ibw` or `.dat` files — file structure is metadata-dependent.*

---

## 📦 Installation

### 🔹 From PyPI (recommended)

```bash
pip install ststools
```

### 🔹 From source (development mode)

```bash
git clone https://github.com/rafinhareis/STSTools.git
cd STSTools
pip install -e .
```

---

## 🚀 Getting Started

Simply import `ststools` in a Jupyter notebook — the main interface will launch automatically:

```python
import ststools
```

You’ll see options to load STS or grid files, define grid size, and start visualizing your data interactively — without needing to write custom scripts.

---

### 🧭 Interface Overview

| Element                  | Description                                                                 |
|---------------------------|-----------------------------------------------------------------------------|
| **X/Y Bias Range**        | Define the voltage window (min/max) for generating the current or dI/dV map |
| **Threshold Range**       | Set intensity limits for map visualization (removes outliers from colormap) |
| **Δ (Delta)**             | Controls voltage resolution for the map (interpolation step)               |
| **Map Type**              | Select the map to display: `I(V)`, `dI/dV`, or `log(I)`                     |
| **Colormap**              | Choose a color scheme (`viridis`, `inferno`, etc.)                          |
| **Row / Column sliders**  | Select the exact spatial point to extract STS curves from the grid         |
| **STS Map (middle)**      | Displays a bias map (bias vs position) based on current or dI/dV            |
| **Bottom Panel**          | Shows both I(V) and dI/dV(V) curves at the selected grid coordinate         |

---

### 📈 What You’ll See

- 🔲 **Left Panel**: Color map showing signal intensity at a chosen bias across the grid
- 📊 **Center Panel**: Interactive map of bias vs spatial position (usually voltage slices)
- 🧪 **Right Panel** *(optional)*: Can show alternative views like dI/dV or log(I)
- 📍 **Bottom Panel**:
  - I(V) curve at the selected position
  - Derivative curve with spatial coordinate and voltage range

**Example Output:**
```text
Current at Point (5.56, 3.33) nm
Gap = 2.15 eV, Tipo: p
```

---

### 💾 Output Options

- Export `.csv` with all extracted curves
- Save map image with selected colormap and bias
- Use `Display()` function for detailed fit on any selected point

---

## ✍️ Author

Developed and maintained by **Rafael dos Reis Barreto**  
📧 Contact: rafinhareis17@gmail.com  
📦 PyPI: [https://pypi.org/project/ststools](https://pypi.org/project/ststools/)  
💻 GitHub: [https://github.com/rafinhareis/STSTools](https://github.com/rafinhareis/STSTools)

> ❤️ If this tool helps you in your research, feel free to cite or mention me in your paper =D

---

## 📜 License

MIT License — feel free to use, fork and improve.
