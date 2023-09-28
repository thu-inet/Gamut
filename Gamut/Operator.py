from abc import ABC, abstractmethod
from Spectrum import Spectrum


class Operator(ABC):

    @abstractmethod
    def __init__(self, label: str = None):
        if label is None:
            label = 'operator'
        self._label = label
        pass

    @property
    def label(self):
        return self._label

    @abstractmethod
    def __run__(self) -> object:
        pass

    def __call__(self, spectrum: Spectrum):
        if not isinstance(spectrum, Spectrum):
            raise TypeError('Input must be a Spectrum object.')
        else:
            return self.__run__(spectrum)


class DualOperator(ABC):

    @abstractmethod
    def __init__(self, label='dualoperator'):
        self._label = label
        pass

    @property
    def label(self):
        return self._label

    @abstractmethod
    def __run__(self):
        pass

    def __call__(self, spectrum1, spectrum2):
        if not (isinstance(spectrum1, Spectrum)
                and (isinstance(spectrum2, Spectrum))):
            raise TypeError('Input must be a Spectrum object.')
        else:
            return self.__run__(spectrum1, spectrum2)


class OperatorPipe():

    def __init__(self, *args):
        self._operators = args
        self.__label__()

    @property
    def label(self):
        self.__label__()
        return self._label

    @property
    def operators(self):
        return self._operators

    def __label__(self):
        self._label = ''.join([operator.label for operator in self._operators])

    def __call__(self, spectrum):
        if isinstance(spectrum, Spectrum):
            for operator in self._operators:
                spectrum = operator(spectrum)
        else:
            raise TypeError('Input must be a Spectrum object.')
        return spectrum

    def add(self, operator):
        if isinstance(operator, Operator):
            self._operators.append(operator)
            self.__label__()
        elif isinstance(operator, OperatorPipe):
            self._operators += operator._operators
            self.__label__()
    
    # def plot(self):
    #     lines = ['graph TD']
    #     for i, scheme in enumerate(self._schemes):
    #         for j, operator in enumerate(
    #          self._operatorpipes[i].operators):
    #             start_node_id = f"{scheme[0]}-{scheme[1]}-{j}"
    #             end_node_id = f"{scheme[0]}-{scheme[1]}-{j+1}"
    #             if j == 0:
    #                 start_node_id = f"{scheme[0]}"
    #             if j == len(self._operatorpipes[i].operators)-1:
    #                 end_node_id = f"{scheme[1]}"
    #             lines.append(f'    {start_node_id}[{start_node_id}] -- "{operator._label}" --> {end_node_id}[{end_node_id}]')
    #     with open('pipenet.mmd', 'w') as fileopen:
    #         fileopen.write('\n'.join(lines))


class OperatorPipeNet():

    def __init__(self, *args, schemes=None):
        self._operatorpipes = args
        if schemes is None:
            schemes = [(i, i+1) for i in range(len(self._operatorpipes))]
        self._schemes = schemes
        self._nodes = self.__rank__()
        self._pipedspectrums = [[] for i in range(max(self._nodes)+1)]

    def __rank__(self):
        nodes = []
        for scheme in self._schemes:
            for node in scheme:
                if node not in nodes:
                    nodes.append(node)
        return sorted(nodes)

    def __call__(self, spectrum):

        nodes = self._nodes.copy()  # nodes han's been arrived
        schemes = self._schemes.copy()  # operators han's been used
        self._pipedspectrums[0].append(spectrum)

        print("================================")
        print("GAMUT: Start running pipenet")
        print("================================")
        while nodes:
            print("-------")
            rest_nodes = nodes.copy()
            for node in nodes:
                if node not in [scheme[1] for scheme in schemes]:
                    rest_nodes.remove(node)
                    print("|Target node: ", node)
            nodes = rest_nodes.copy()
            print("|Untargeted nodes: ", *nodes)

            rest_schemes = schemes.copy()
            for scheme in schemes:
                if scheme[0] not in nodes:
                    index_operator = self._schemes.index(scheme)
                    new_spectrum = self._operatorpipes[index_operator](
                        self._pipedspectrums[scheme[0]][0])
                    self._pipedspectrums[scheme[1]].append(new_spectrum)
                    rest_schemes.remove(scheme)
                    print("|Target scheme: ", scheme)
                    print("|Execute operator: ",
                          self._operatorpipes[index_operator].label)
            schemes = rest_schemes.copy()
            print("|Untargeted schemes: ", *schemes)

    def add(self, pipe, scheme=None):
        if isinstance(pipe, OperatorPipe):
            self._operatorpipes.append(pipe)
            self._scheme.append(scheme)

    def plot(self):
        lines = ['graph TD']
        for i, scheme in enumerate(self._schemes):
            for j, operator in enumerate(
             self._operatorpipes[i].operators):
                start_node_id = f"Pipe{i}{j}"
                start_node_name = f"Pipe[{i}]Spectrum[{j}]"
                end_node_id = f"Pipe{i}{j+1}"
                end_node_name = f"Pipe[{i}]Spectrum[{j+1}]"
                if j == 0:
                    start_node_id = f"Spectrum{scheme[0]}"
                    start_node_name = f"Spectrum[{scheme[0]}]"
                if j == len(self._operatorpipes[i].operators)-1:
                    end_node_id = f"Spectrum{scheme[1]}"
                    end_node_name = f"Spectrum[{scheme[1]}]"
                lines.append(f"    {start_node_id}[\"{start_node_name}\"] -- \"{operator._label[1:]}\" --> {end_node_id}[\"{end_node_name}\"]")
        with open('pipenet.mmd', 'w') as fileopen:
            fileopen.write('\n'.join(lines))