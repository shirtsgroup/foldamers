"""
Unit and regression test for the foldamers package.
"""

# Import package, test suite, and other packages as needed
import foldamers
import pytest
import sys

def test_foldamers_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "foldamers" in sys.modules
