import numpy as np

from typing import Union, Callable
from collections import namedtuple
from abc import ABC, abstractmethod

from .Spectrum import Spectrum


class Operator(ABC):

    def __init__(self, inp_num: int | None, label: str = None):
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

    def __call__(self, spectra: list[Spectrum] | Spectrum, check_type: bool = True, add_label: bool = True, *args, **kargs) -> Spectrum:

        if isinstance(spectra, Spectrum) and self.inp_num == 1:
            spectra = [spectra]

        if check_type and (not all([isinstance(spectrum, Spectrum) for spectrum in spectra])):
            raise TypeError('Input must be a list of Spectrum object.')
        if self.inp_num and len(spectra) != self.inp_num:
            raise TypeError(f'Input must be a list of length {self.inp_num}.')

        spectrum = self.__run__(spectra, *args, **kargs)
        spectrum[:] = np.maximum(spectrum, 0)
        if add_label:
            spectrum.label += "|" + self._label
        return spectrum

    @abstractmethod
    def __run__(self, spectra: list[Spectrum], *args, **kargs) -> Spectrum:
        pass


class Pipe(Operator):

    def __init__(self, operators: list[Operator], label: str = None):
        for operator in operators[1:]:
            if operator.inp_num != 1:
                raise AttributeError(f'Subsequent operator of inp_num {operator.inp_num}')
        self._operators = operators

        if label is None:
            label = 'Pipe'
        super().__init__(operators[0].inp_num, label)

    @property
    def operators(self) -> list[Operator]:
        return self._operators

    @property
    def spectra(self) -> list[Spectrum]:
        return self._spectra

    def __getitem__(self, index):
        return self._operators[index]

    def __str__(self):
        return f"Pipe[{self.label}]"

    def __repr__(self):
        operators = '>>'.join([operator.label for operator in self._operators])
        return f"Pipe[{self.label}]\n \
|InpNum: {self.inp_num}\n \
|Operators: {operators}\n"

    def __run__(self, spectra: list[Spectrum], *args, **kargs) -> Spectrum:

        self._spectra = []
        label = spectra[0].label + ">>" + self._label

        spectrum = spectra[0].copy()
        for operator in self._operators:
            spectrum = operator(spectra, check_type=False, add_label=True)
            self._spectra.append(spectrum)
            spectra = [spectrum]

        spectrum.label = label
        return spectrum

    def append(self, operator):
        if isinstance(operator, Operator):
            if operator.inp_num != 1:
                raise AttributeError(f'Subsequent operator of inp_num {operator.inp_num}')
            self._operators.append(operator)
        elif isinstance(operator, Pipe):
            if operator[0].inp_num != 1:
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

    def __init__(self, id: int, spectra: dict[str, Spectrum]):
        self.id = id
        self.spectra = spectra
        
    def __getitem__(self, index):
        return self.spectra[index]

    def __setitem__(self, index, value):
        self.spectra[index] = value

class Flow:
    
    def __init__(self, operator: Operator, inp_id: int, inp_order: list[int] | int, outp_id: int, outp_order: int):
        self.operator = operator
        self.inp_id = inp_id
        self.inp_order = inp_order if isinstance(inp_order, list) else [inp_order]
        self.outp_id = outp_id
        self.outp_order = outp_order

class PipeNet:

    def __init__(self, flows: list[Flow], label: str = None):

        if label is None:
            label = 'PipeNet'
        self._label = label
        self._flows = flows
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
        for flow in self._flows:
            nid = flow.inp_id
            spectra = dict([(i, None) for i in flow.inp_order])
            if nid not in ranked_nodes_id:
                self._nodes.append(Node(nid, spectra))
                ranked_nodes_id.append(flow.inp_id)
            else:
                node = [node for node in self._nodes if node.id == nid][0]
                node.spectra.update(spectra)

        for flow in self._flows:
            nid = flow.outp_id
            spectra = {flow.outp_order: None}
            if nid not in ranked_nodes_id:
                self._nodes.append(Node(nid, spectra))
                ranked_nodes_id.append(flow.inp_id)
            else:
                node = [node for node in self._nodes if node.id == nid][0]
                if flow.outp_order not in node.spectra.keys():
                    raise ValueError(f"Node {nid}'s output is not in node{node.id}'s input")

        if [node for node in self._nodes if node.id == 0][0] is None:
            raise ValueError("Node 0 does not exist.")

    def __call__(self, spectra: list[Spectrum]):
        """
        Run the PipeNet.
        """
        node = [node for node in self._nodes if node.id == 0][0]  # insert the input spectra into initial node
        if len(node.spectra) != len(spectra):
            raise ValueError(f"Input spectra' length {len(spectra)}  is not the same as Node 0's input length {len(node.spectra)}.")
        node.spectra = spectra

        print("================================")
        print(f"GAMUT: Start running {self.label}")
        print("================================")
        unexecuted_flows = self._flows.copy()
        while unexecuted_flows:
            unexecuted_flows_copy = unexecuted_flows.copy()
            for flow in unexecuted_flows:
                if self._check_flow(flow):
                    print(f"|Execute: Node{flow.inp_id:>2d}|{flow.inp_order} >> {flow.operator.label} >> Node{flow.outp_id:>2d}|[{flow.outp_order}]")
                    self._execute_flow(flow)
                    unexecuted_flows_copy.remove(flow)
                unexecuted_flows = unexecuted_flows_copy

    def _check_flow(self, flow: Flow):
        """
        Check if the flow is ready to execute.
        """
        inp_node = [node for node in self._nodes if node.id == flow.inp_id][0]
        return all([inp_node[ind] is not None for ind in flow.inp_order])

    def _execute_flow(self, flow: Flow):
        """
        Execute the flow.
        """
        inp_node = [node for node in self._nodes if node.id == flow.inp_id][0]
        outp_node = [node for node in self._nodes if node.id == flow.outp_id][0]

        spectra = [inp_node[ind] for ind in flow.inp_order]
        spectrum = flow.operator(spectra)
        outp_node[flow.outp_order] = spectrum

    def append(self, flow: Flow):
        """
        Append new flows.
        """
        self._flows.append(flow)
        self._rank_nodes()

    def plot(self, verbose: bool = False):
        """
        Plot .mmd graph of current flow.
        """
        lines = ['graph TD']
        for flow in self._flows:
            pipe = flow.operator
            inp_node = [node for node in self._nodes if node.id == flow.inp_id][0]
            outp_node = [node for node in self._nodes if node.id == flow.outp_id][0]
            if verbose and isinstance(pipe, Pipe):
                flow_repr = f"{flow.inp_order}>>" + ">>".join([operator.label for operator in pipe.operators]) + f">>[{flow.outp_order}]"
            else:
                flow_repr = f"{flow.inp_order}>>{pipe._label}>>[{flow.outp_order}]"
            lines.append(f"{inp_node.id}[\"{inp_node}\"] -- \"{flow_repr}\" --> {outp_node.id}[\"{outp_node}\"]")
        with open('pipenet.mmd', 'w') as fileopen:
            fileopen.write('\n'.join(lines))

    def get_node(self, nid):
        return [node for node in self._nodes if node.id == nid][0]


