import RNA
import numpy as np
import matplotlib.pyplot as plt
from typing import List, Tuple

class nucleic_graph:
    def __init__(self, mfe, sz=5, degree=0, 
                 bckwargs={'lw':1, 'color':'k'}, 
                 bpkwargs={'lw':1, 'c':'red'}, 
                 scwargs={'s':10, 'c':'k'}):
        
        self.mfe = mfe
        self.sz = sz
        self.spacer = 1
        self.degree = degree
        
        self.generate_db_structure()
        self.generate_coordinates_and_pairs()

    def generate_db_structure(self):
        # Let's check whether it is a complex or single strand
        if '+' in self.mfe:
            complex = self.mfe.split('+')
            self.l = [len(entry) for entry in complex]
            self.dot_bracket_str = "".join(flatten([[entry, '.'*self.spacer] for entry in complex]))

        else:
            self.dot_bracket_str = self.mfe
            self.l = []

    def generate_coordinates_and_pairs(self):
        vrna_coords = RNA.get_xy_coordinates(self.dot_bracket_str)
        coords = []
        for i, _ in enumerate(self.dot_bracket_str):
            coord = (vrna_coords.get(i).X, vrna_coords.get(i).Y)
            coords.append(coord)
            
        coords = rotate(coords, np.mean(coords, axis=0), self.degree)

        self.coords = np.array(coords)
        self.pairs = self.parse_dot_bracket()

    def parse_dot_bracket(self):
        # Let's extract pairs of nucleotides as tuples
        stack = []
        pairs = []

        for i, char in enumerate(self.dot_bracket_str):
            if char == '(':
                stack.append(i)
            elif char == ')':
                if not stack:
                    raise ValueError(f"Unmatched closing bracket at position {i}")
                j = stack.pop()
                pairs.append((j, i))

        if stack:
            raise ValueError(f"Unmatched opening brackets at positions {stack}")
        
        return pairs

    def plotter(self, bpkwargs, bckwargs, scwargs):
        fig, ax = plt.subplots(figsize=(self.sz, self.sz))

        # Manually retrieve paired bases and relative coordinates to plot linkers
        for entry in self.pairs:
            ax.plot([self.coords[entry[0]][0], self.coords[entry[1]][0]], 
                    [self.coords[entry[0]][1], self.coords[entry[1]][1]], **bpkwargs, zorder=0)

        # Let's now plot the backbone
        if len(self.l) > 0: # it is a complex, it has to be treated differently
            # I need to remove the spacers I placed previously to generate the graph
            n = []
            for i in range(len(self.l)):
                n.append([sum(self.l[:i])+self.spacer*i, sum(self.l[:i+1])+self.spacer*i])

            revised_coordinates = []
            for entry in n:
                entry = list(range(entry[0], entry[1]))
                revised_coordinates.append(entry)
                ax.plot(self.coords[np.array(entry), 0], self.coords[np.array(entry), 1],
                        **bckwargs, zorder=1)
                ax.scatter(self.coords[np.array(entry), 0], self.coords[np.array(entry), 1],
                           **scwargs, zorder=2)
                
            # Let's update coordinates to get rid of spacers
            self.coords = self.coords[flatten(revised_coordinates)]

        else:
            ax.plot(self.coords[:,0], self.coords[:,1], **bckwargs)

        datalim = ((min(list(self.coords[:, 0]) + [ax.get_xlim()[0]]),
                    min(list(self.coords[:, 1]) + [ax.get_ylim()[0]])),
                   (max(list(self.coords[:, 0]) + [ax.get_xlim()[1]]),
                    max(list(self.coords[:, 1]) + [ax.get_ylim()[1]])))

        width = datalim[1][0] - datalim[0][0]
        height = datalim[1][1] - datalim[0][1]

        ax.set_aspect('equal', 'datalim')
        ax.update_datalim(datalim)
        ax.autoscale_view()
        ax.set_axis_off()
        
        self.ax = ax