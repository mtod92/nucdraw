import numpy as np
from typing import List, Tuple

def rotate(p: np.ndarray, origin: Tuple[float, float] =(0.0, 0.0), degrees: float =0) -> np.ndarray:
    # function from user: ImportanceOfBeingErnest retrieved from Stack Overflow
    angle = np.deg2rad(degrees)
    R = np.array([[np.cos(angle), -np.sin(angle)],
                  [np.sin(angle),  np.cos(angle)]])
    o = np.atleast_2d(origin)
    p = np.atleast_2d(p)
    return np.squeeze((R @ (p.T-o.T) + o.T).T)

def flatten(xss: List[List]) -> List:
    return [x for xs in xss for x in xs]
