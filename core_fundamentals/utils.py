"""
utils.py
--------
Utility functions for Bioinformatics and Entrez API configurations.
Author: Surena Ghadam (surenaghadam-biotech)
"""

import os
from dotenv import load_dotenv

def get_ncbi_email() -> str:
    """
    Loads environment variables and retrieves the NCBI email securely.
    Raises ValueError if NCBI_EMAIL is missing to prevent silent failures.
    """
    
    load_dotenv()
    
    
    email = os.environ.get("NCBI_EMAIL")
    
    
    if not email:
        raise ValueError(
            "Configuration Error: 'NCBI_EMAIL' is missing from environment variables. "
            "Please ensure you have a .env file with your valid email address."
        )
        
    return email
