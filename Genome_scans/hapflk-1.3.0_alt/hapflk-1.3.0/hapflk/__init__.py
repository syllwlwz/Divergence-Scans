import numpy as np

## Help pyinstaller find hidden imports
def dependencies_for_myprogram():
    from scipy.sparse.csgraph import _validation

missing=np.int16(-1)

def _complete_cases(x):
    return x!=missing

complete_cases=np.vectorize(_complete_cases)

def _is_na(x):
    return x==missing

is_na=np.vectorize(_is_na)

