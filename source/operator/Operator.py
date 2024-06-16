import numpy as np

from collections import namedtuple
from abc import ABC, abstractmethod

from ..spectrum.Spectrum import Spectrum


class Operator(ABC):
    """
    Operator is the basic class in GAMUT to wrap all spectrum analysis algorithms.    
    """
    def __init__(self, inp_num: int | None, label: str | None = None):
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

        self.__dbg__(spectrum)
        spectrum[:] = np.maximum(spectrum, 0)
        if add_label:
            spectrum.label += "|" + self._label
        return spectrum

    @abstractmethod
    def __run__(self, spectra: list[Spectrum] | Spectrum, *args, **kargs) -> Spectrum:
        pass

    def __dbg__(self, spectrum: Spectrum) -> Spectrum:
        pass


class Pipe(Operator):

    def __init__(self, operators: list[Operator], label: str | None = None):
        for operator in operators[1:]:
            if operator.inp_num != 1:
                raise AttributeError(f'Subsequent operator of inp_num {operator.inp_num}')
        self._operators = operators

        if label is None:
            label = '>>'.join([operator.label for operator in operators])
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

    def __init__(self, id: int, spectra: dict[int, Spectrum | None]):
        self.id = id
        self.spectra = spectra
        
    def __getitem__(self, index: int) -> Spectrum | None:
        return self.spectra[index]

    def __setitem__(self, index, value):
        self.spectra[index] = value

    def __len__(self):
        return len(self.spectra)


class Flow:

    def __init__(self, operator: Operator | list[Operator] | Pipe, 
                 inp_id: int, inp_order: list[int] | int, 
                 outp_id: int, outp_order: int):
        if isinstance(operator, list):
            if all([isinstance(opr, Operator) for opr in operator]):
                operator = Pipe(operator)
            else:
                raise TypeError('Input must be a list of Operator object.')
        self.operator = operator
        self.inp_id = inp_id
        self.inp_order = inp_order if isinstance(inp_order, list) else [inp_order]
        self.outp_id = outp_id
        self.outp_order = outp_order


class PipeNet:

    def __init__(self, flows: list[Flow], label: str | None = None):

        if label is None:
            label = 'PipeNet'
        self._label = label
        self._flows = flows
        self._nodes = {}
        self._rank_nodes()

    @property
    def label(self):
        return self._label

    @label.setter
    def label(self, label: str):
        self._label = label

    @property
    def nodes(self) -> dict[int, Node]:
        return self._nodes

    def _rank_nodes(self):
        """
        Generate spectra nodes.
        """
        ranked_nodes_id = []
        for flow in self._flows:
            nid = flow.inp_id
            spectra = dict([(i, None) for i in flow.inp_order])
            if nid not in ranked_nodes_id:
                self._nodes[nid] = Node(nid, spectra)
                ranked_nodes_id.append(flow.inp_id)
            else:
                node = self._nodes[nid]
                node.spectra.update(spectra)

        ranked_nodes_id = []
        for flow in self._flows:
            nid = flow.outp_id
            spectra = {flow.outp_order: None}
            if nid not in ranked_nodes_id:
                self._nodes[nid] = Node(nid, spectra)
                ranked_nodes_id.append(flow.inp_id)
            else:
                node = self._nodes[nid]
                if flow.outp_order not in node.spectra.keys():
                    raise ValueError(f"Node {nid}'s output is not in node{node.id}'s input")

        if self.nodes[0] is None:
            raise ValueError("Node 0 does not exist.")

    def _check_flow(self, flow: Flow):
        """
        Check if the flow is ready to execute.
        """
        inp_node = self.nodes[flow.inp_id]
        return all([inp_node[ind] is not None for ind in flow.inp_order])

    def _execute_flow(self, flow: Flow):
        """
        Execute the flow.
        """
        inp_node = self.nodes[flow.inp_id]
        outp_node = self.nodes[flow.outp_id]
        spectra = [inp_node[ind] for ind in flow.inp_order]
        spectrum = flow.operator(spectra)
        outp_node[flow.outp_order] = spectrum

    def __call__(self, spectra: list[Spectrum]):
        """
        Run the PipeNet.
        """
        node = self.nodes[0]  # insert the input spectra into initial node
        if len(node.spectra) != len(spectra):
            raise ValueError(f"Input spectra' length {len(spectra)}  is not the same as Node 0's input length {len(node.spectra)}.")
        node.spectra = dict([(i, spectrum) for i, spectrum in enumerate(spectra)])

        print("================================")
        print(f"GAMUT: Start running {self.label}")
        print("================================")
        unexecuted_flows = self._flows.copy()
        while len(unexecuted_flows) > 0:
            unexecuted_flows_copy = unexecuted_flows.copy()
            for flow in unexecuted_flows:
                if self._check_flow(flow):
                    print(f"|Execute: Node{flow.inp_id:>2d}|{flow.inp_order} >> {flow.operator.label} >> Node{flow.outp_id:>2d}|[{flow.outp_order}]")
                    self._execute_flow(flow)
                    unexecuted_flows_copy.remove(flow)
                unexecuted_flows = unexecuted_flows_copy

        return self._nodes.copy()

    def plot(self, vertical: bool = False):
        """
        Plot .mmd graph of current flow.
        """
        rankdir = 'LR' if vertical else 'TD'
        maingraph = ["%%{ init: { 'flowchart': { 'curve': 'linear' } } }%%",
                     f'graph {rankdir}\n',
                     'classDef Node fill:#FBE5D5,stroke:#333,stroke-width:2px;',
                     'classDef Operator fill:#BDD7EE,stroke:#333,stroke-width:2px;',
                     f'subgraph {self.label}']
        subgraphs = []

        for nid, node in self._nodes.items():
            maingraph.append(f"\tNode{nid}[Node{nid}]:::Node")

        subgraph_id = 1
        for flow in self._flows:

            inp_node = self._nodes[flow.inp_id]
            outp_node = self._nodes[flow.outp_id]
            if not isinstance(flow.operator, Pipe):
                flow_repr = f"{flow.inp_order}>>{flow.operator.label}>>[{flow.outp_order}]"
            else:
                flow_repr = f"{flow.inp_order}>>Pipe{subgraph_id}>>[{flow.outp_order}]"
            
            if inp_node is not None and outp_node is not None:
                maingraph.append(f"\tNode{inp_node.id} -- \"{flow_repr}\" --> Node{outp_node.id}")
            else:
                raise ValueError(f"\tNode {flow.inp_id} or Node {flow.outp_id} does not exist.")

            if isinstance(flow.operator, Pipe):
                subgraph = [f"subgraph Pipe{subgraph_id}"]
                for i in range(len(flow.operator.operators)):
                    subgraph.append(f"\tOperator{i}[\"{flow.operator.operators[i].label}\"]:::Operator")
                for i in range(len(flow.operator.operators)-1):
                    subgraph.append(f"\tOperator{i} --> Operator{i+1}")
                subgraph.append('end\n')
                subgraphs.append(subgraph)

        maingraph.append('end\n')
        graphlines = maingraph + [line for subgraph in subgraphs for line in subgraph]
        with open('pipenet.mmd', 'w') as fileopen:
            fileopen.write('\n'.join(graphlines))


