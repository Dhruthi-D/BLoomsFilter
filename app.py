import os
import time
import psutil
import hashlib
import numpy as np
from flask import Flask, render_template, request, jsonify, session
from werkzeug.utils import secure_filename
from pybloom_live import BloomFilter, ScalableBloomFilter
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import io
import base64
import re
import random
import plotly.graph_objects as go
import plotly.utils
import json
import logging
import traceback

# Configure logging
logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)

app = Flask(__name__)
app.secret_key = 'replace-this-with-a-very-secret-key-1234567890'  # Set a unique, secret key for session support

# Ensure upload directory exists with proper permissions
try:
    UPLOAD_FOLDER = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'uploads')
    if not os.path.exists(UPLOAD_FOLDER):
        os.makedirs(UPLOAD_FOLDER, exist_ok=True)
        logger.info(f"Created upload directory: {UPLOAD_FOLDER}")
    # Test write permissions
    test_file = os.path.join(UPLOAD_FOLDER, 'test.txt')
    with open(test_file, 'w') as f:
        f.write('test')
    os.remove(test_file)
    logger.info("Upload directory permissions verified")
except Exception as e:
    logger.error(f"Error setting up upload directory: {str(e)}")
    logger.error(traceback.format_exc())

app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
app.config['MAX_CONTENT_LENGTH'] = 16 * 1024 * 1024  # 16MB max file size

# Global variables
bloom_filters = {
    'standard': None,
    'scalable': None,
    'partitioned': None
}
performance_metrics = {
    'standard': {'insert_time': 0, 'search_time': 0, 'memory': 0, 'false_positives': 0, 'false_positive_rate': 0},
    'scalable': {'insert_time': 0, 'search_time': 0, 'memory': 0, 'false_positives': 0, 'false_positive_rate': 0},
    'partitioned': {'insert_time': 0, 'search_time': 0, 'memory': 0, 'false_positives': 0, 'false_positive_rate': 0}
}
k_value = 10  # Default k-value

# Store original kmers and metrics
original_kmers = set()

latest_results_data = {}

class PartitionedBloomFilter:
    def __init__(self, capacity, error_rate=0.001, num_partitions=4):
        self.num_partitions = num_partitions
        # Add 50% buffer to each partition
        partition_capacity = int((capacity * 1.5) / num_partitions)
        logger.debug(f"Creating partitioned Bloom filter with capacity {capacity}, partition capacity {partition_capacity}")
        self.partitions = [BloomFilter(capacity=partition_capacity, error_rate=error_rate) 
                         for _ in range(num_partitions)]
    
    def add(self, item):
        partition_index = hash(item) % self.num_partitions
        self.partitions[partition_index].add(item)
    
    def __contains__(self, item):
        partition_index = hash(item) % self.num_partitions
        return item in self.partitions[partition_index]

def clean_sequence(sequence):
    """Remove non-DNA characters and convert to uppercase."""
    return re.sub(r'[^ACGT]', '', sequence.upper())

def validate_dna_sequence(sequence):
    """Validate DNA sequence and ensure it's long enough for k-mers"""
    if not sequence:
        return False, "Sequence is empty"
    
    # Check for valid DNA characters
    if not all(base in 'ACGT' for base in sequence.upper()):
        return False, "Sequence contains invalid characters. Only A, C, G, T are allowed"
    
    # Check minimum length
    if len(sequence) < k_value:
        return False, f"Sequence must be at least {k_value} characters long for k={k_value}"
    
    return True, "Valid sequence"

def extract_kmers(sequence, k):
    """Extract k-mers from a DNA sequence, ensuring exact k length"""
    kmers = set()
    sequence = clean_sequence(sequence)
    if not sequence:
        return kmers
    
    # Only extract k-mers of exact length k
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i + k]
        if len(kmer) == k:  # Double check the length
            kmers.add(kmer)
    return kmers

def generate_random_kmers(num_kmers, k):
    """Generate random k-mers for false positive testing."""
    bases = 'ACGT'
    kmers = set()
    while len(kmers) < num_kmers:
        kmer = ''.join(random.choice(bases) for _ in range(k))
        kmers.add(kmer)
    return kmers

def measure_memory_usage(obj):
    """Measure memory usage of an object in MB."""
    try:
        return psutil.Process().memory_info().rss / 1024 / 1024
    except:
        return 0.0  # Return 0 if measurement fails

def create_performance_plots(metrics):
    """Create performance comparison plots using Plotly"""
    filter_types = ['standard', 'scalable', 'partitioned']
    
    # Time comparison plot
    time_fig = go.Figure(data=[
        go.Bar(
            x=filter_types,
            y=[metrics[ft]['insert_time'] for ft in filter_types],
            name='Insert Time',
            marker_color='rgb(55, 83, 109)'
        )
    ])
    
    # Add search time bars only if they exist
    if all('search_time' in metrics[ft] for ft in filter_types):
        time_fig.add_trace(go.Bar(
            x=filter_types,
            y=[metrics[ft]['search_time'] for ft in filter_types],
            name='Search Time',
            marker_color='rgb(26, 118, 255)'
        ))
    
    time_fig.update_layout(
        # title='Time Comparison by Bloom Filter Type',
        xaxis_title='Bloom Filter Type',
        yaxis_title='Time (seconds)',
        barmode='group',
        template='plotly_white'
    )
    
    # Memory usage comparison plot
    memory_fig = go.Figure(data=[
        go.Bar(
            x=filter_types,
            y=[metrics[ft]['memory'] for ft in filter_types],
            marker_color='rgb(26, 118, 255)'
        )
    ])
    memory_fig.update_layout(
        # title='Memory Usage by Bloom Filter Type',
        xaxis_title='Bloom Filter Type',
        yaxis_title='Memory Usage (MB)',
        template='plotly_white'
    )
    
    # False positive rate comparison plot
    fpr_fig = go.Figure(data=[
        go.Bar(
            x=filter_types,
            y=[metrics[ft]['false_positive_rate'] * 100 for ft in filter_types],
            marker_color='rgb(26, 118, 255)'
        )
    ])
    fpr_fig.update_layout(
        # title='False Positive Rate by Bloom Filter Type',
        xaxis_title='Bloom Filter Type',
        yaxis_title='False Positive Rate (%)',
        template='plotly_white'
    )
    
    return {
        'time_plot': json.dumps(time_fig, cls=plotly.utils.PlotlyJSONEncoder),
        'memory_plot': json.dumps(memory_fig, cls=plotly.utils.PlotlyJSONEncoder),
        'fpr_plot': json.dumps(fpr_fig, cls=plotly.utils.PlotlyJSONEncoder)
    }

def hash_function(data):
    """Custom hash function using hashlib"""
    return int(hashlib.md5(str(data).encode()).hexdigest(), 16)

def analyze_kmer_sizes(sequence, min_k=5, max_k=20):
    """Analyze different k-mer sizes and their characteristics."""
    results = []
    sequence = clean_sequence(sequence)
    
    for k in range(min_k, max_k + 1):
        kmers = set()
        for i in range(len(sequence) - k + 1):
            kmer = sequence[i:i + k]
            if len(kmer) == k:
                kmers.add(kmer)
        
        # Calculate metrics
        total_possible = len(sequence) - k + 1
        unique_kmers = len(kmers)
        if total_possible > 0:
            uniqueness_ratio = unique_kmers / total_possible
            coverage_ratio = unique_kmers / (4**k)  # Coverage compared to all possible k-mers
            # Calculate a score that balances uniqueness and practical usability
            practical_score = uniqueness_ratio * (1 - (k / max_k))  # Penalize very large k values
        else:
            uniqueness_ratio = 0
            coverage_ratio = 0
            practical_score = 0
        
        results.append({
            'k_size': k,
            'total_kmers': total_possible,
            'unique_kmers': unique_kmers,
            'uniqueness_ratio': uniqueness_ratio,
            'coverage_ratio': coverage_ratio,
            'practical_score': practical_score
        })
    
    # Find best k-mer sizes for different criteria
    best_uniqueness = max(results, key=lambda x: x['uniqueness_ratio'])
    best_coverage = max(results, key=lambda x: x['coverage_ratio'])
    best_practical = max(results, key=lambda x: x['practical_score'])
    
    return {
        'analysis': results,
        'recommendations': {
            'best_uniqueness': best_uniqueness['k_size'],
            'best_coverage': best_coverage['k_size'],
            'best_practical': best_practical['k_size'],
            'explanation': {
                'uniqueness': f"K={best_uniqueness['k_size']} gives highest uniqueness ratio of {best_uniqueness['uniqueness_ratio']:.2%}",
                'coverage': f"K={best_coverage['k_size']} gives best coverage of {best_coverage['coverage_ratio']:.2%}",
                'practical': f"K={best_practical['k_size']} gives best balance of uniqueness and practicality"
            }
        }
    }

def create_kmer_analysis_plot(kmer_data):
    """Create a plot showing k-mer size analysis."""
    analysis = kmer_data['analysis']
    recommendations = kmer_data['recommendations']
    
    k_sizes = [r['k_size'] for r in analysis]
    uniqueness_ratios = [r['uniqueness_ratio'] * 100 for r in analysis]
    coverage_ratios = [r['coverage_ratio'] * 100 for r in analysis]
    practical_scores = [r['practical_score'] * 100 for r in analysis]
    
    # Create the plot data
    data = [
        {
            'x': k_sizes,
            'y': uniqueness_ratios,
            'name': 'Uniqueness Ratio (%)',
            'type': 'scatter',
            'line': {'color': 'blue'},
            'hovertemplate': 'K-mer size: %{x}<br>Uniqueness: %{y:.2f}%'
        },
        {
            'x': k_sizes,
            'y': coverage_ratios,
            'name': 'Coverage Ratio (%)',
            'type': 'scatter',
            'line': {'color': 'red', 'dash': 'dash'},
            'hovertemplate': 'K-mer size: %{x}<br>Coverage: %{y:.2f}%'
        },
        {
            'x': k_sizes,
            'y': practical_scores,
            'name': 'Practical Score (%)',
            'type': 'scatter',
            'line': {'color': 'green', 'dash': 'dot'},
            'hovertemplate': 'K-mer size: %{x}<br>Practical Score: %{y:.2f}%'
        }
    ]
    
    # Create layout
    layout = {
        'title': 'K-mer Size Analysis and Recommendations',
        'xaxis': {
            'title': 'K-mer Size',
            'tickmode': 'linear',
            'dtick': 1
        },
        'yaxis': {
            'title': 'Percentage (%)',
            'range': [0, 100]
        },
        'hovermode': 'x unified',
        'showlegend': True,
        'legend': {
            'yanchor': 'top',
            'y': 0.99,
            'xanchor': 'left',
            'x': 0.01
        }
    }
    
    # Add vertical lines for recommended k-mer sizes
    shapes = []
    annotations = []
    colors = {'uniqueness': 'blue', 'coverage': 'red', 'practical': 'green'}
    
    y_positions = [95, 85, 75]  # Different y-positions for annotations
    i = 0
    for k_type, k_value in recommendations.items():
        if isinstance(k_value, int):  # Skip the 'explanation' key
            shapes.append({
                'type': 'line',
                'x0': k_value,
                'x1': k_value,
                'y0': 0,
                'y1': 100,
                'line': {
                    'color': colors[k_type.split('_')[1]],
                    'dash': 'dash'
                }
            })
            annotations.append({
                'x': k_value,
                'y': y_positions[i],
                'text': f"Best {k_type.split('_')[1]}: k={k_value}",
                'showarrow': False,
                'font': {'color': colors[k_type.split('_')[1]]}
            })
            i += 1
    
    layout['shapes'] = shapes
    layout['annotations'] = annotations
    
    return {'data': data, 'layout': layout}

@app.route('/')
def index():
    latest_results_data.clear()
    return render_template('index.html')

@app.route('/upload', methods=['POST'])
def upload_file():
    global k_value
    filepath = None
    try:
        logger.debug("Starting file upload process")
        
        # Verify upload directory exists and is writable
        if not os.path.exists(UPLOAD_FOLDER):
            os.makedirs(UPLOAD_FOLDER, exist_ok=True)
            logger.info(f"Created upload directory: {UPLOAD_FOLDER}")
        
        if 'file' not in request.files:
            logger.error("No file part in request")
            return jsonify({'error': 'No file provided'}), 400
        
        file = request.files['file']
        if file.filename == '':
            logger.error("No selected file")
            return jsonify({'error': 'No file selected'}), 400
        
        # Get k-value from request
        try:
            k_value = int(request.form.get('k_value', 10))
            if not 5 <= k_value <= 20:
                logger.error(f"Invalid k-value: {k_value}")
                return jsonify({'error': 'K-value must be between 5 and 20'}), 400
            logger.debug(f"Using k-value: {k_value}")
        except ValueError as e:
            logger.error(f"Error parsing k-value: {str(e)}")
            return jsonify({'error': 'Invalid k-value'}), 400

        # Save the file temporarily
        try:
            filename = secure_filename(file.filename)
            filepath = os.path.join(UPLOAD_FOLDER, filename)
            logger.debug(f"Attempting to save file to: {filepath}")
            
            # Ensure the file is not too large
            file.seek(0, os.SEEK_END)
            size = file.tell()
            file.seek(0)
            if size > app.config['MAX_CONTENT_LENGTH']:
                logger.error(f"File too large: {size} bytes")
                return jsonify({'error': 'File too large'}), 400
            
            file.save(filepath)
            logger.debug(f"File saved successfully to: {filepath}")
            
            # Verify file was saved correctly
            if not os.path.exists(filepath):
                raise Exception("File was not saved properly")
            
            # Verify file is readable
            with open(filepath, 'r', encoding='utf-8') as f:
                first_line = f.readline()
                logger.debug(f"First line of file: {first_line[:100]}...")
                
        except Exception as e:
            logger.error(f"Error saving file: {str(e)}")
            logger.error(traceback.format_exc())
            if filepath and os.path.exists(filepath):
                try:
                    os.remove(filepath)
                except:
                    pass
            return jsonify({'error': f'Error saving file: {str(e)}'}), 500
        
        # Process the file
        kmers = set()
        try:
            logger.debug(f"Reading file: {filepath}")
            with open(filepath, 'r', encoding='utf-8') as f:
                content = f.read()
                logger.debug(f"File content length: {len(content)}")
                
                if not content.strip():
                    raise Exception("File is empty")
                
                # Split content into lines and process each line
                lines = content.splitlines()
                logger.debug(f"Number of lines in file: {len(lines)}")
                
                for i, line in enumerate(lines, 1):
                    if not line.startswith('>'):  # Skip FASTA headers
                        sequence = line.strip().upper()
                        logger.debug(f"Processing line {i}, length: {len(sequence)}")
                        if len(sequence) >= k_value:
                            new_kmers = extract_kmers(sequence, k_value)
                            kmers.update(new_kmers)
                            logger.debug(f"Line {i}: Extracted {len(new_kmers)} k-mers")
            first_seq = next((line.strip() for line in lines if not line.startswith('>')), '')
        except Exception as e:
            logger.error(f"Error reading file: {str(e)}")
            logger.error(traceback.format_exc())
            if filepath and os.path.exists(filepath):
                try:
                    os.remove(filepath)
                except:
                    pass
            return jsonify({'error': f'Error reading file: {str(e)}'}), 500
        
        if not kmers:
            logger.error("No valid k-mers found in file")
            if filepath and os.path.exists(filepath):
                try:
                    os.remove(filepath)
                except:
                    pass
            return jsonify({'error': f'No valid {k_value}-mers found in the file'}), 400
        
        logger.debug(f"Total unique k-mers extracted: {len(kmers)}")
        
        # Initialize Bloom filters with increased capacity
        total_kmers = len(kmers)
        capacity = int(total_kmers * 3)  # Triple the capacity to be extra safe
        logger.debug(f"Initializing Bloom filters with capacity: {capacity}")
        
        try:
            # Clear existing filters
            bloom_filters['standard'] = None
            bloom_filters['scalable'] = None
            bloom_filters['partitioned'] = None
            
            # Initialize new filters
            bloom_filters['standard'] = BloomFilter(capacity=capacity, error_rate=0.001)
            bloom_filters['scalable'] = ScalableBloomFilter(initial_capacity=capacity, error_rate=0.001)
            bloom_filters['partitioned'] = PartitionedBloomFilter(capacity=capacity, error_rate=0.001)
            logger.debug("Bloom filters initialized successfully")
        except Exception as e:
            logger.error(f"Error initializing Bloom filters: {str(e)}")
            logger.error(traceback.format_exc())
            if filepath and os.path.exists(filepath):
                try:
                    os.remove(filepath)
                except:
                    pass
            return jsonify({'error': f'Error initializing Bloom filters: {str(e)}'}), 500
        
        # Measure insertion time and memory usage
        results = {}
        for filter_type, filter_instance in bloom_filters.items():
            try:
                logger.debug(f"Processing {filter_type} Bloom filter")
                start_time = time.time()
                
                # Add k-mers in smaller batches
                batch_size = 1000
                kmer_list = list(kmers)
                for i in range(0, len(kmer_list), batch_size):
                    batch = kmer_list[i:i + batch_size]
                    for kmer in batch:
                        if len(kmer) == k_value:
                            filter_instance.add(kmer)
                    logger.debug(f"Added batch {i//batch_size + 1} to {filter_type}")
                
                insert_time = time.time() - start_time
                logger.debug(f"Insertion time for {filter_type}: {insert_time:.4f}s")
                
                # Measure memory usage
                memory_usage = measure_memory_usage(filter_instance)
                logger.debug(f"Memory usage for {filter_type}: {memory_usage:.2f}MB")
                
                # Calculate false positive rate using decoy k-mers
                decoy_kmers = generate_random_kmers(1000, k_value)
                false_positives = sum(1 for kmer in decoy_kmers if kmer in filter_instance)
                false_positive_rate = false_positives / len(decoy_kmers)
                logger.debug(f"False positive rate for {filter_type}: {false_positive_rate:.4f}")
                
                results[filter_type] = {
                    'insert_time': insert_time,
                    'memory': memory_usage,
                    'false_positives': false_positives,
                    'false_positive_rate': false_positive_rate
                }
                
                # Update performance metrics
                performance_metrics[filter_type].update(results[filter_type])
            except Exception as e:
                logger.error(f"Error processing {filter_type} Bloom filter: {str(e)}")
                logger.error(traceback.format_exc())
                if filepath and os.path.exists(filepath):
                    try:
                        os.remove(filepath)
                    except:
                        pass
                return jsonify({'error': f'Error processing {filter_type} Bloom filter: {str(e)}'}), 500
        
        # Create performance plots
        plots = create_performance_plots(results)
        
        # After processing the file content
        kmer_analysis = analyze_kmer_sizes(content)
        kmer_plot = create_kmer_analysis_plot(kmer_analysis)
        
        # Update the results dictionary with the complete analysis
        results = {
            'message': f'Successfully processed {total_kmers} {k_value}-mers',
            'results': results,
            'plots': plots,
            'kmer_analysis': {
                'plot': kmer_plot,
                'recommendations': kmer_analysis['recommendations'],
                'analysis_data': kmer_analysis['analysis']
            }
        }
        
        # Store results in latest_results_data
        latest_results_data.update(results)
        
        return jsonify(results)
        
    except Exception as e:
        logger.error(f"Unexpected error in upload_file: {str(e)}")
        logger.error(traceback.format_exc())
        if filepath and os.path.exists(filepath):
            try:
                os.remove(filepath)
            except:
                pass
        return jsonify({'error': f'Error processing file: {str(e)}'}), 500

@app.route('/search', methods=['POST'])
def search():
    global k_value
    try:
        data = request.get_json()
        if not data or 'query' not in data:
            return jsonify({'error': 'No query provided'}), 400
        
        query = data['query'].strip().upper()
        is_valid, message = validate_dna_sequence(query)
        if not is_valid:
            return jsonify({'error': message}), 400
        
        # Check if filters are initialized
        if not all(bloom_filters.values()):
            return jsonify({'error': 'Please upload a file first'}), 400
        
        results = {}
        for filter_type, filter_instance in bloom_filters.items():
            start_time = time.time()
            kmers = extract_kmers(query, k_value)
            if not kmers:
                return jsonify({'error': f'No valid {k_value}-mers could be extracted from the query'}), 400
            
            # Find matching k-mers
            matching_kmers = [kmer for kmer in kmers if len(kmer) == k_value and kmer in filter_instance]
            matches = len(matching_kmers)
            search_time = time.time() - start_time
            
            results[filter_type] = {
                'search_time': search_time,
                'matches': matches,
                'total_kmers': len(kmers),
                'match_percentage': (matches / len(kmers) * 100) if kmers else 0,
                'found': matches > 0,
                'matching_kmers': matching_kmers,
                'query_length': len(query),
                'kmer_size': k_value
            }
            
            # Merge with previous metrics to ensure all columns are present
            merged = performance_metrics[filter_type].copy()
            merged.update(results[filter_type])
            results[filter_type] = merged
        
        # Create updated plots including search times
        plots = create_performance_plots(performance_metrics)
        
        # Store results, plots, and query_info in latest_results_data
        latest_results_data['results'] = results
        latest_results_data['plots'] = plots
        latest_results_data['query_info'] = {
            'sequence': query,
            'length': len(query),
            'kmer_size': k_value,
            'total_kmers': len(kmers)
        }
        
        return jsonify({
            'results': results,
            'plots': plots,
            'query_info': {
                'sequence': query,
                'length': len(query),
                'kmer_size': k_value,
                'total_kmers': len(kmers)
            }
        })
        
    except Exception as e:
        return jsonify({'error': str(e)}), 500

@app.route('/know-more')
def know_more():
    return render_template('know_more.html')

@app.route('/analysis')
def analysis():
    return render_template('analysis.html')

@app.route('/latest-results')
def latest_results():
    return jsonify(latest_results_data)

if __name__ == '__main__':
    app.run(debug=True) 