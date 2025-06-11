# DNA Sequence Bloom Filter Application

This application implements three different types of Bloom filters (Standard, Scalable, and Partitioned) for efficient DNA sequence pattern matching and comparison. It provides a web interface for uploading DNA sequences and performing quick pattern searches.

## What is a Bloom Filter?

A Bloom filter is a space-efficient probabilistic data structure used to test whether an element is a member of a set. It can have false positives but never false negatives. This makes it perfect for DNA sequence matching where we want to quickly determine if a sequence pattern might exist in a larger dataset.

The application implements three types of Bloom filters:
1. **Standard Bloom Filter**: Fixed-size implementation with a predefined error rate
2. **Scalable Bloom Filter**: Automatically grows as more items are added while maintaining the desired false positive rate
3. **Partitioned Bloom Filter**: Divides the bit array into partitions to reduce false positives

## Features

- Web-based interface for easy interaction
- Support for FASTA format DNA sequence files
- Three different Bloom filter implementations
- Real-time performance metrics:
  - Insertion time
  - Memory usage
  - False positive rate
- Interactive visualizations using Plotly
- K-mer based sequence analysis
- Efficient DNA pattern searching
- Performance comparison between different Bloom filter types

## Installation

### Prerequisites
- Python 3.12 or higher
- pip (Python package installer)

### Clone the Repository
```bash
git clone <your-repository-url>
cd BLoomsFilter
```

### Install Dependencies
```bash
pip install -r requirements.txt
```

## Usage

1. Start the Flask application:
```bash
python app.py
```

2. Open your web browser and go to:
```
http://localhost:5000
```

3. Using the Web Interface:
   - Upload a DNA sequence file (FASTA format)
   - Set the k-mer size (between 5-20)
   - Use the search functionality to find patterns

### Input File Format
The application accepts FASTA format files containing DNA sequences. Example:
```
>Sequence_1
ATCGATCGATCGATCGATCG
>Sequence_2
GCTAGCTAGCTAGCTAGCTA
```

### Search Guidelines
- Search queries must be at least as long as the k-mer size
- Use only valid DNA characters (A, T, C, G)
- Longer sequences will provide more accurate results

## Technical Details

### Application Structure
- `app.py`: Main Flask application file
- `templates/`: HTML templates for the web interface
- `static/`: CSS, JavaScript, and other static files
- `uploads/`: Temporary storage for uploaded files

### Key Components
1. **File Processing**:
   - FASTA file parsing
   - K-mer extraction
   - Sequence validation

2. **Bloom Filter Implementation**:
   - Standard Bloom Filter using pybloom-live
   - Scalable Bloom Filter for growing datasets
   - Custom Partitioned Bloom Filter

3. **Performance Monitoring**:
   - Memory usage tracking
   - Insertion time measurement
   - False positive rate calculation

4. **Visualization**:
   - Real-time performance graphs
   - Comparison charts
   - Interactive plots

## Performance Considerations

- Memory Usage: The application optimizes memory usage through batch processing
- Search Speed: Utilizes efficient hashing and bit array operations
- Scalability: Handles large DNA sequences through chunked processing
- False Positive Rate: Maintained at approximately 0.1% (configurable)

## Error Handling

The application includes comprehensive error handling for:
- Invalid file formats
- Incorrect DNA sequences
- Memory limitations
- Invalid search queries

## Contributing

Feel free to submit issues, fork the repository, and create pull requests for any improvements.
