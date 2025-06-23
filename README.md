# DNA Sequence Bloom Filter Web Application

A modern, interactive web app for DNA sequence analysis, comparison, and visualization using advanced Bloom filter-based tools.

---

## Features

### 🌱 **Bloom Filter Analysis**
- **Upload DNA sequence files** (FASTA or plain text).
- **Choose k-mer size** (5–20).
- **Three Bloom filter types:** Standard, Scalable, Partitioned.
- **Performance metrics:** Insertion time, memory usage, false positive rate.
- **Interactive visualizations:** Real-time performance graphs, k-mer analysis, and comparison charts.
- **Efficient pattern search:** Search for DNA motifs or sequences, with match context and statistics.

### 🧬 **Species Tools**
- **Compare two species' DNA files** side-by-side.
- **K-mer uniqueness heatmap:** Visualize unique and shared k-mers between species.
- **Venn diagram:** See overlap and uniqueness of k-mers.
- **Pairwise statistics:** Overlap and Jaccard similarity.
- **GC content, motif, and palindromic k-mer analysis** for each species.
- **Modern, responsive UI** with dark mode and clear feedback.

### 🔗 **Direct Sequence Alignment**
- **Upload wild-type and GMO DNA files** for global alignment (Needleman-Wunsch).
- **Interactive alignment visualization:** Color-coded matches, mismatches, insertions, deletions.
- **DNA strand SVG diagram:** Modern, zoomable, with tooltips.
- **Alignment summary pie chart** and additional analysis graphs.

### 📊 **Analysis Dashboard**
- **Compare performance** of all Bloom filter types on your data.
- **View and download results** for further study.

### ℹ️ **Know More**
- **Educational section** explaining Bloom filters, their types, and real-world use cases.
- **Comparison table** and diagrams for quick understanding.

### 🌓 **Modern UI/UX**
- **Floating dark mode toggle** (persistent across sessions).
- **Consistent teal accent color and modern design.**
- **Accessible navigation:** "Back to Main" on all feature pages.

---

## Installation

### Prerequisites
- Python 3.12 or higher
- pip

### Clone the Repository
```bash
git clone <your-repository-url>
cd BLoomsFilter
```

### Install Dependencies
```bash
pip install -r requirements.txt
```

---

## Usage

1. **Start the Flask app:**
   ```bash
   python app.py
   ```
2. **Open your browser:**  
   [http://localhost:5000](http://localhost:5000)

3. **Main Navigation:**
   - **Know More:** Learn about Bloom filters and DNA analysis.
   - **Go to Analysis:** Upload a DNA file, select k-mer size, and analyze with Bloom filters.
   - **Species Tools:** Compare two species' DNA, visualize k-mer uniqueness, and analyze motifs.
   - **Direct Sequence Alignment:** Align wild-type and GMO DNA, view interactive results.

---

## Input File Format

- **FASTA or plain text** (one sequence per line or FASTA headers).
- Example:
  ```
  >Sequence_1
  ATCGATCGATCGATCGATCG
  >Sequence_2
  GCTAGCTAGCTAGCTAGCTA
  ```

- **Sample files** are provided in the `dataset/` directory.

---

## Directory Structure

- `app.py` — Main Flask app
- `templates/` — All HTML pages (main, analysis, species tools, alignment, know more)
- `uploads/` — Temporary file storage
- `dataset/` — Example DNA files
- `requirements.txt` — Python dependencies

---

## Technical Details

- **Bloom filter implementations:** Standard, Scalable, Partitioned (pybloom-live and custom).
- **K-mer extraction, motif search, palindromic analysis, GC content calculation.**
- **Global sequence alignment:** Needleman-Wunsch algorithm for direct comparison.
- **Interactive visualizations:** Plotly, Chart.js, D3.js, and SVG.
- **Robust error handling** for file format, sequence validation, and upload issues.

---

## Contributing

- Issues and pull requests are welcome!
- Please ensure new features are well-documented and tested.

---

## License

MIT License (or specify your license here)
