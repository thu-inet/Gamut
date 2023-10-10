from abc import ABC, abstractmethod

import matplotlib.pyplot as plt
import numpy as np

from Spectrum import Spectrum

from typing import Union, Callable


class Operator(ABC):

    def __init__(self, inp_num: int, label: str = None):
        self._inp_num = inp_num
        if label is None:
            label = 'Operator'
        self._label = label

    @property
    def label(self):
        return self._label

    @property
    def inp_num(self):
        return self._inp_num

    @label.setter
    def label(self, label: str):
        self._label = label

    def __str__(self):
        return f"Operator[{self.label}]"

    def __repr__(self):
        return f"Operator[{self.label}]\n \
|InpNum: {self.inp_num}\n"

    def __call__(self, spectra: list[Spectrum], check_type: bool = True, add_label: bool = True, *args, **kargs) -> Spectrum:

        if isinstance(spectra, Spectrum) and self.inp_num == 1:
            spectra = [spectra]

        if check_type and (not all([isinstance(spectrum, Spectrum) for spectrum in spectra])):
            raise TypeError('Input must be a list of Spectrum object.')
        if self.inp_num and len(spectra) != self.inp_num:
            raise TypeError(f'Input must be a list of length {self.inp_num}.')
        spectrum = self.__run__(spectra, *args, **kargs)
        if add_label:
            spectrum.label += ">>" + self._label
        return spectrum

    @abstractmethod
    def __run__(self, spectra: list[Spectrum], *args, **kargs) -> Spectrum:
        pass


class Pipe:

    def __init__(self, operators: list[Operator], label: str = None):
        for operator in operators[1:]:
            if operator.inp_num != 1:
                raise AttributeError(f'Subsequent operator of inp_num {operator.inp_num}')
        self._operators = operators
        self._inp_num = operators[0].inp_num

        if label is None:
            label = 'Pipe'
        self._label = label

    @property
    def operators(self):
        return self._operators

    @property
    def inp_num(self):
        return self._inp_num

    @property
    def label(self):
        return self._label

    @label.setter
    def label(self, label: str):
        self._label = label

    def __getattr__(self, attr):
        return getattr(self._operators, attr)

    def __getitem__(self, index):
        return self._operators[index]

    def __str__(self):
        return f"Pipe[{self.label}]"

    def __repr__(self):
        operators = '>>'.join([operator.label for operator in self._operators])
        return f"Pipe[{self.label}]\n \
|InpNum: {self.inp_num}\n \
|Operators: {operators}\n"

    def __call__(self, spectra: list[Spectrum]) -> Spectrum:

        if isinstance(spectra, Spectrum) and self.inp_num == 1: 
            spectra = [spectra]
        if not all([isinstance(spectrum, Spectrum) for spectrum in spectra]):
            raise TypeError('Input must be a list of Spectrum object.')
        self._spectra = []
        spectrum = spectra[0].copy()  # for pipe of zero operator
        for operator in self._operators:
            spectrum = operator(spectrum, check_type=False, add_label=True)
            self._spectra.append(spectrum)
        spectrum.label = spectra[0].label + ">>" + self._label
        return spectrum

    def append(self, operator):
        if isinstance(operator, Operator):
            if operator.inp_num != 1:
                raise AttributeError(f'Subsequent operator of inp_num {operator.inp_num}')
            self._operators.append(operator)
        elif isinstance(operator, Pipe):
            if Pipe.operators[0].inp_num != 1:
                raise AttributeError(f'Subsequent operator of inp_num {operator.inp_num}')
            self._operators += operator._operators
        else:
            raise TypeError('Input must be a Spectrum object.')

    @property
    def list_spectrum(self):
        return [spectrum for spectrum in self._spectra]

    def get_spectrum(self, index):
        return self._spectra[index]

class Node:

    def __init__(self, id: int, spectra: dict = None):
        self._id = id
        if spectra is None:
            spectra = {}
        self._spectra = spectra

    @property
    def id(self):
        return self._id

    @property
    def spectra(self):
        return self._spectra

    @spectra.setter
    def spectra(self, spectra: dict):
        self._spectra = spectra

    def __getitem__(self, key):
        return self._spectra[key]

    def __setitem__(self, key, value):
        self._spectra[key] = value

    def __str__(self):
        return f'Node[{self._id:>2d}]'

    def __repr__(self):
        return f'Node[{self._id:>2d}]\n \
|Spectra: {list(self.spectra.keys())}\n'


class PipeNet:

    def __init__(self, schemes: list[dict], label: str = None):

        if label is None:
            label = 'PipeNet'
        self._label = label

        self._schemes = schemes
        self._rank_nodes()

    @property
    def label(self):
        return self._label

    @label.setter
    def label(self, label: str):
        self._label = label

    def _rank_nodes(self):
        """
        Generate spectra nodes.
        """
        self._nodes = []
        ranked_nodes_id = []
        for scheme in self._schemes:
            nid = scheme["inp_id"]
            spectra = dict([(i, None) for i in scheme["inp_order"]])
            if nid not in ranked_nodes_id:
                self._nodes.append(Node(nid, spectra))
                ranked_nodes_id.append(scheme["inp_id"])
            else:
                node = [node for node in self._nodes if node.id == nid][0]
                node.spectra.update(spectra)

        for scheme in self._schemes:
            nid = scheme["outp_id"]
            spectra = {scheme['outp_order']: None}
            if nid not in ranked_nodes_id:
                self._nodes.append(Node(nid, spectra))
                ranked_nodes_id.append(scheme["inp_id"])
            else:
                node = [node for node in self._nodes if node.id == nid][0]
                if scheme['outp_order'] not in node.spectra.keys():
                    raise ValueError(f"Node {nid}'s output is not in node{node.id}'s input")

        if [node for node in self._nodes if node.id == 0][0] is None:
            raise ValueError("Node 0 does not exist.")

    def __call__(self, spectra: dict[Spectrum]):
        """
        Run the PipeNet.
        """
        node = [node for node in self._nodes if node.id == 0][0] # insert the input spectra into initial node
        if len(node.spectra) != len(spectra):
            raise ValueError(f"Input spectra' length {len(spectra)}  is not the same as Node 0's input length {len(node.spectra)}.")
        node.spectra = spectra

        print("================================")
        print(f"GAMUT: Start running {self.label}")
        print("================================")
        unexecuted_schemes = self._schemes.copy()
        while unexecuted_schemes:
            unexecuted_schemes_copy = unexecuted_schemes.copy()
            for scheme in unexecuted_schemes:
                if self._check_scheme(scheme):
                    print(f"|Execute: Node{scheme['inp_id']:>2d}|{scheme['inp_order']} >> {scheme['pipe']._label} >> Node{scheme['outp_id']:>2d}|[{scheme['outp_order']}]")
                    self._execute_scheme(scheme)
                    unexecuted_schemes_copy.remove(scheme)
                unexecuted_schemes = unexecuted_schemes_copy

    def _check_scheme(self, scheme: dict):
        """
        Check if the scheme is ready to execute.
        """
        inp_node = [node for node in self._nodes if node.id == scheme["inp_id"]][0]
        return all([inp_node[ind] is not None for ind in scheme['inp_order']])

    def _execute_scheme(self, scheme: dict):
        """
        Execute the scheme.
        """
        inp_node = [node for node in self._nodes if node.id == scheme["inp_id"]][0]
        outp_node = [node for node in self._nodes if node.id == scheme["outp_id"]][0]

        spectra = [inp_node[ind] for ind in scheme['inp_order']]
        spectrum = scheme['pipe'](spectra)
        outp_node[scheme['outp_order']] = spectrum

    def append(self, scheme: dict):
        """
        Append new schemes.
        """
        self._scheme.append(scheme)
        self._rank_nodes()

    def plot(self, verbose: bool = False):
        """
        Plot .mmd graph of current scheme.
        """
        lines = ['graph TD']
        for scheme in self._schemes:
            pipe = scheme['pipe']
            inp_node = [node for node in self._nodes if node.id == scheme["inp_id"]][0]
            outp_node = [node for node in self._nodes if node.id == scheme["outp_id"]][0]
            if verbose and isinstance(pipe, Pipe):
                scheme_repr = f"{scheme['inp_order']}>>" + ">>".join([operator.label for operator in pipe.operators]) + f">>[{scheme['outp_order']}]"
            else:
                scheme_repr = f"{scheme['inp_order']}>>{pipe._label}>>[{scheme['outp_order']}]"
            lines.append(f"{inp_node.id}[\"{inp_node}\"] -- \"{scheme_repr}\" --> {outp_node.id}[\"{outp_node}\"]")
        with open('pipenet.mmd', 'w') as fileopen:
            fileopen.write('\n'.join(lines))

    def get_node(self, nid):
        return [node for node in self._nodes if node.id == nid][0]

if __name__ == "__main__":

    from BasicOperator import FunctionalOperator, Multiplier

    # Single Input Operator Initialization
    linear = FunctionalOperator(func=lambda x: 3*x, label='linear')
    modulo = FunctionalOperator(func=lambda x: x % 2, label='modulo')
    divide = FunctionalOperator(func=lambda x: x / 2)

    # Multiple Input Operator Initialization
    multiply = Multiplier()

    # Initiate test spectrum
    spectrum = Spectrum(counts=np.arange(1, 5), label='original')

    # Linear Pipes
    pipe1 = Pipe([linear, modulo], label='linmod')
    pipe2 = Pipe([linear, divide], label='lindiv')

    # Nonlinear PipeNet
    schemes1 = [
                {'pipe': linear, 'inp_id': 0, 'inp_order': [0], 'outp_id': 1, 'outp_order': 0},
                {'pipe': modulo, 'inp_id': 0, 'inp_order': [0], 'outp_id': 2, 'outp_order': 0},
                {'pipe': pipe1, 'inp_id': 0, 'inp_order': [0], 'outp_id': 2, 'outp_order': 1},
                {'pipe': multiply, 'inp_id': 2, 'inp_order': [0, 1], 'outp_id': 3, 'outp_order': 0}
                ]
    net = PipeNet(schemes=schemes1)

    # Plot mmd graph of PipeNet
    net.plot(verbose=True)

    # Run spectrum on PipeNet
    net({0: spectrum})

    # Print spectra
    for node in net._nodes:
        print(node)
        for spectrum in node._spectra.values():
            print(repr(spectrum))