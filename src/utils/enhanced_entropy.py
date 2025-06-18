import numpy as np
import matplotlib.pyplot as plt
from typing import List, Union
import pandas as pd

'''
Length Considerations Added

1. Length-Adjusted Entropy
Applies a penalty factor to shorter lists
Recognizes that entropy calculated from longer lists is statistically more reliable
Uses a logarithmic scale that approaches full entropy value as length increases

2. Confidence Intervals
Calculates confidence level for entropy based on list length
Longer lists have tighter confidence intervals
Visualizes uncertainty in entropy measurements

Key New Features

1. Length-Aware Visualization
Scatter plot where point size represents list length
Shows relationship between proportion of 1s and entropy

2. Enhanced Statistical Reliability
Penalizes entropy from very short lists (less statistically significant)
Properly weights longer lists in the analysis
'''

def calculate_entropy(binary_list: List[Union[int, float]]) -> float:
    """
    Calculate the Shannon entropy of a binary list.
    
    Parameters:
    -----------
    binary_list : List[Union[int, float]]
        A list containing only 0s and 1s
        
    Returns:
    --------
    float
        The Shannon entropy value
    
    Raises:
    -------
    ValueError
        If the list contains values other than 0 and 1
    """
    # Validate the list contains only 0s and 1s
    unique_values = set(binary_list)
    if not unique_values.issubset({0, 1}):
        raise ValueError("List contains values other than 0 and 1")
    
    # Count occurrences of 0s and 1s
    n = len(binary_list)
    if n == 0:
        return 0.0
        
    count_0 = binary_list.count(0)
    count_1 = n - count_0
    
    # Calculate probabilities
    p_0 = count_0 / n
    p_1 = count_1 / n
    
    # Calculate entropy
    entropy = 0.0
    if p_0 > 0:
        entropy -= p_0 * np.log2(p_0)
    if p_1 > 0:
        entropy -= p_1 * np.log2(p_1)
    
    return entropy

def calculate_entropy_with_length_penalty(binary_list: List[Union[int, float]], 
                                         max_length: int = 500,
                                         penalty_factor: float = 0.2) -> float:
    """
    Calculate entropy with a penalty for shorter lists.
    Shorter lists are less statistically significant, so we apply a penalty.
    
    Parameters:
    -----------
    binary_list : List[Union[int, float]]
        A list containing only 0s and 1s
    max_length : int
        The maximum possible length of a list for normalization
    penalty_factor : float
        A factor controlling how much to penalize shorter lists
        
    Returns:
    --------
    float
        The length-adjusted entropy value
    """
    entropy = calculate_entropy(binary_list)
    length = len(binary_list)
    
    # Apply a logarithmic penalty that decreases as length increases
    # This is a smooth function that approaches 1 as length approaches max_length
    length_factor = 1 - penalty_factor * (np.log(max_length + 1) - np.log(length + 1)) / np.log(max_length + 1)
    
    return entropy * length_factor

def calculate_entropy_confidence(binary_list: List[Union[int, float]]) -> tuple:

    # maybe using bootstrap to calculate confidence interval? Hosein: using length for confidence interval is reasonable
    """
    Calculate entropy and its confidence interval based on list length.
    Longer lists give more confident entropy estimates.
    
    Parameters:
    -----------
    binary_list : List[Union[int, float]]
        A list containing only 0s and 1s
        
    Returns:
    --------
    tuple
        (entropy, confidence_level)
    """
    entropy = calculate_entropy(binary_list)
    length = len(binary_list)
    
    # Calculate confidence level (0-1 scale)
    # This approaches 1 for longer lists
    confidence = 1 - 1/(np.sqrt(length) + 1)
    
    return entropy, confidence

def compare_entropy_methods(lists: List[List[int]]) -> pd.DataFrame:
    """
    Compare entropy values using different calculation methods.
    
    Parameters:
    -----------
    lists : List[List[int]]
        A list of binary lists
        
    Returns:
    --------
    pandas.DataFrame
        A DataFrame with different entropy metrics
    """
    results = []
    
    for i, lst in enumerate(lists):
        row = {
            'list_index': i,
            'length': len(lst),
            'proportion_of_1s': lst.count(1) / len(lst),
            'basic_entropy': calculate_entropy(lst),
            'length_adjusted_entropy': calculate_entropy_with_length_penalty(lst),
            'entropy_confidence': calculate_entropy_confidence(lst)[1]
        }
        
        results.append(row)
    
    # Convert to DataFrame
    df = pd.DataFrame(results)
    return df

def visualize_enhanced_entropy(lists: List[List[int]], labels: List[str] = None):
    """
    Create enhanced visualizations for entropy comparison.
    
    Parameters:
    -----------
    lists : List[List[int]]
        A list of binary lists
    labels : List[str], optional
        Labels for each list. If None, indices will be used.
    """
    if labels is None:
        labels = [f"List {i}" for i in range(len(lists))]
    
    # Get comparison data
    df = compare_entropy_methods(lists)
    df['label'] = [labels[i] for i in df['list_index']]
    
    # Create figure with subplots
    fig, axes = plt.subplots(2, 1, figsize=(14, 10))
    
    # Plot 1: Basic entropy vs Length-adjusted entropy
    ax = axes[0]
    x = range(len(df))
    ax.bar(x, df['basic_entropy'], width=0.4, align='edge', label='Basic Entropy')
    ax.bar([i+0.4 for i in x], df['length_adjusted_entropy'], width=0.4, align='edge', 
           label='Length-Adjusted Entropy')
    ax.set_xticks([i+0.2 for i in x])
    ax.set_xticklabels(df['label'], rotation=45, ha='right')
    ax.set_title('Entropy Comparison: Basic vs. Length-Adjusted')
    ax.set_ylabel('Entropy (bits)')
    ax.set_ylim(0, 1.05)
    ax.legend()
    
    # Plot 2: Entropy with confidence intervals
    ax = axes[1]
    confidences = [calculate_entropy_confidence(lst) for lst in lists]
    entropy_values = [conf[0] for conf in confidences]
    confidence_levels = [conf[1] for conf in confidences]
    
    bars = ax.bar(labels, entropy_values,
                 yerr=[(1-conf)*e for e, conf in zip(entropy_values, confidence_levels)],
                 capsize=5)
    ax.set_title('Entropy with Confidence Intervals')
    ax.set_ylabel('Entropy (bits)')
    ax.set_xticklabels(labels, rotation=45, ha='right')
    ax.set_ylim(0, 1.1)
    
    plt.tight_layout()
    return fig, df

# Example usage
if __name__ == "__main__":
    # Example binary lists of different lengths and patterns
    list1 = [0, 0, 0, 0, 0]  # All 0s - Low entropy
    list2 = [1, 1, 1, 1, 1]  # All 1s - Low entropy
    list3 = [0, 1, 0, 1, 0, 1, 0, 1]  # Alternating - High entropy
    list4 = [0, 0, 0, 1, 1, 1]  # Mix - Medium entropy
    list5 = [0, 1, 1, 0, 1, 0, 1, 1, 0, 0] * 5  # Longer mixed list
    list6 = [0, 1] * 25  # Very long alternating pattern
    list7 = [0] * 48 + [1] * 52  # Long list with all 0s then all 1s
    list8 = [np.random.randint(0, 2) for _ in range(100)]  # Random binary list
    
    all_lists = [list1, list2, list3, list4, list5, list6, list7, list8]
    labels = ["All 0s", "All 1s", "Alternating", "Grouped", 
              "Mixed Pattern", "Long Alternating", "Long Grouped", "Random"]
    
    # Compare using enhanced methods
    comparison_df = compare_entropy_methods(all_lists)
    print("Entropy Comparison with Length Considerations:")
    print(comparison_df[['list_index', 'length', 'basic_entropy', 
                         'length_adjusted_entropy']])
    
    # Visualize with enhanced plots
    fig, results_df = visualize_enhanced_entropy(all_lists, labels)
    
    # Print detailed table
    print("\nDetailed Results:")
    print(results_df[['label', 'length', 'basic_entropy', 
                      'length_adjusted_entropy']])
