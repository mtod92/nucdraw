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

def compute_repulsion_vectors(points, k=2, cutoff=None):
    N = len(points)
    arrows = np.zeros_like(points)

    for i in range(N):
        direction = np.zeros(2)
        for j in range(N):
            if i == j:
                continue
            diff = points[i] - points[j]
            dist = np.linalg.norm(diff)
            if cutoff is not None and dist > cutoff:
                continue
            if dist < 1e-5:
                continue  # avoid division by zero
            direction += diff / dist**k
        norm = np.linalg.norm(direction)
        if norm > 0:
            direction /= norm
        arrows[i] = direction
    return arrows