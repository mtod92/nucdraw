import RNA
import numpy as np
import matplotlib.pyplot as plt
from typing import List, Tuple
from utils import rotate, flatten

class NucDraw:
    def __init__(self, dot_bracket_structure: str):
        self.dot_bracket_structure = dot_bracket_structure
    
    def generate(self, spacer: int=1, degree: float=0.0):
        self.generate_db_structure(spacer)
        self.generate_coordinates_and_pairs(degree)

    def generate_db_structure(self, spacer):
        self.spacer = spacer

        # Let's check whether it is a complex or single strand
        if '+' in self.dot_bracket_structure:
            complex = self.dot_bracket_structure.split('+')
            self.l = [len(entry) for entry in complex]
            self.dot_bracket_str = "".join(flatten([[entry, '.'*self.spacer] for entry in complex]))

        else:
            self.dot_bracket_str = self.dot_bracket_structure
            self.l = []

    def generate_coordinates_and_pairs(self, degree):
        vrna_coords = RNA.get_xy_coordinates(self.dot_bracket_str)
        coords = []
        for i, _ in enumerate(self.dot_bracket_str):
            coord = (vrna_coords.get(i).X, vrna_coords.get(i).Y)
            coords.append(coord)
            
        coords = rotate(coords, np.mean(coords, axis=0), degree)

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

    def plotter(self, sz: float=5.0, 
                 bckwargs: dict={'lw':1, 'color':'k'}, 
                 bpkwargs: dict={'lw':1, 'c':'red'}, 
                 scwargs: dict={'s':0, 'c':'k'}):
        fig, ax = plt.subplots(figsize=(sz, sz))

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

    def plot_circles(self, sequence: str = '', circle_size: float=5, circle_color = ['#1ea3eb', '#f5370c', '#22c716', '#e6d815']):

        if len(sequence) == 0: # no sequence is passed
            if isinstance(circle_color, list): # you are passing a list of colors, but they cannot be mapped because of no sequence
                color_table = circle_color[0]
            else:
                color_table = circle_color

            for i in range(len(self.coords)):
                self.ax.add_patch(plt.Circle((self.coords[i][0], self.coords[i][1]), circle_size, zorder=2, edgecolor='k', facecolor=color_table))

        else: # sequence is passed
            if len(sequence) != len(self.coords):
                raise ValueError("Sequence and structure must have the same length.")
        
            if isinstance(circle_color, list):
                color_table = {'A': circle_color[0], 'C': circle_color[1], 'G': circle_color[2], 'U': circle_color[3], 'T': circle_color[3]}
            else:
                color_table = {'A': circle_color, 'C': circle_color, 'G': circle_color, 'U': circle_color, 'T': circle_color}

            for i in range(len(self.coords)):
                self.ax.add_patch(plt.Circle((self.coords[i][0], self.coords[i][1]), circle_size, zorder=2, edgecolor='k', facecolor=color_table[sequence[i]]))

    def plot_sequence(self, sequence: str, kwargs: dict={'fontsize': 12, 'color': 'k'}):
        if len(sequence) != len(self.coords):
            raise ValueError("Sequence and structure must have the same length.")
    
        for i in range(len(sequence)):
            plt.text(self.coords[i][0], self.coords[i][1], sequence[i], ha='center', va='center', zorder=3, **kwargs)

    def multistrand_coloring(self, clr=[], bckwargs: dict={'lw':1}):

        if len(self.l) > 1: # check that it is indeed a complex
            for i in range(len(self.l)):
                idx = range(sum(self.l[:i]), sum(self.l[:i+1]))
                self.ax.plot(self.coords[idx, 0], self.coords[idx, 1], **bckwargs, color=clr[i], zorder=1)
            
        else:
            raise ValueError("This function should be used for multi-stranded complexes.")


