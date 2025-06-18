import numpy as np
from typing import List, Union
import pandas as pd

'''
 This module enhances entropy analysis of binary lists by providing confidence intervals for entropy estimates, which reflect the uncertainty due to list length.
    - Confidence intervals are wider for shorter lists and narrower for longer lists.
   - This helps visualize and interpret the reliability of entropy measurements.

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


def calculate_entropy_confidence(binary_list: List[Union[int, float]]) -> tuple:
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
