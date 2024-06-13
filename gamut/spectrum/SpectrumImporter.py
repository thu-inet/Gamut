import numpy as np
import pickle
import pandas as pd

from pathlib import Path
from xml.dom.minidom import parse as parse_xml
from re import match

from ..classes import Region, Calibration
from .Spectrum import Spectrum


class SpectrumImporter:
    
    @staticmethod
    def from_pickle(filename: str | Path) -> 'Spectrum':
        """
        Import spectrum from pickle files.
        
        :param filename: Path to the file where the spectrum is saved.
        
        :return: Spectrum object.
        """
        with open(filename, 'rb') as fileopen:
            spectrum = pickle.load(fileopen)
        return spectrum
    
    @staticmethod
    def from_excel(filename: str | Path, counts_column: str, energy_column: str | None = None) -> 'Spectrum':
        """
        Import spectrum from columns of excel files.
        
        :param filename: Path to the file where the spectrum is saved.
        :param counts_column: The name of the column containing the counts.
        :param energy_column: The name of the column containing the energies.
        
        :return: Spectrum object.
        """
        df = pd.read_excel(str(filename), index_col=0)
        if counts_column not in df.columns:
            raise ValueError(f"{counts_column} not in df.columns")
        else:
            spectrum = Spectrum(df[counts_column].to_numpy())
            if energy_column and energy_column not in df.columns:
                raise ValueError(f"{energy_column} not in df.columns")
            else:
                energies = df[energy_column].to_numpy()
                data = list(np.vstack((np.arange(energies.shape[0]), energies)))
                spectrum.ergcal = Calibration(method='linear', data=data)
        return spectrum

    @staticmethod
    def from_xml(filename: str | Path) -> 'Spectrum':
        """
        Import spectrum from xml files.
        
        :param filename: Path to the file where the spectrum is saved.
        
        :return: Spectrum object.
        """
        spec_xml = parse_xml(str(filename))
        if (spec_el := spec_xml.getElementsByTagName('counts')[0]) is None:
            raise ValueError(f"{filename} is not a valid xml file containing spectrum counts.")
        else:
            spectrum = Spectrum(counts=[float(val) for val in spec_el.firstChild.nodeValue.split()],
                    label=spec_el.getAttribute('name'))

        if (ergcal_el := spec_xml.getElementsByTagName('ergcal')[0]) is not None:
            spectrum.ergcal = Calibration(method = ergcal_el.getAttribute('method'),
                                      params=[float(val) for val in ergcal_el.getAttribute('params').split()])

        if (FWHMcal_el := spec_xml.getElementsByTagName('FWHMcal')[0]) is not None:
            spectrum.FWHMcal = Calibration(method=FWHMcal_el.getAttribute('method'),
                                       params=[float(val) for val in FWHMcal_el.getAttribute('params').split()])

        if (regions_el := spec_xml.getElementsByTagName('regions')[0]) is not None:
            for region_el in regions_el.getElementsByTagName('region'):
                region = Region(int(region_el.getAttribute('ind_left')), int(region_el.getAttribute('ind_right')))

                try:
                    region.slope = float(region_el.getAttribute('slope'))
                    region.offset = float(region_el.getAttribute('offset'))
                except:
                    pass

                for peak_el in region_el.getElementsByTagName('peak'):
                    peak = dict([(attr.name, float(attr.value)) for attr in peak_el.attributes.values()])
                    peak['location'] = int(peak['location'])
                    region.peaks.append(peak)
                    
                spectrum.regions.append(region)

        return spectrum

    @staticmethod
    def from_MCNP(filepath: str, tally_id: int = 8) -> 'Spectrum':
        """
        Import spectrum from MCNP output files.
        
        :param filepath: Path to the file where the spectrum is saved.
        
        :return: Spectrum object.
        """
        with open(filepath, 'r') as fileopen:
            filelines = fileopen.readlines()   

        index = 0
        while f'tally type {tally_id}' not in filelines[index]:
            index += 1
        index += 1

        energies, counts = [], []
        while 'energy' not in filelines[index]:
            index += 1
        index += 1
        while 'total' not in filelines[index]:
            erg, eff, err = filelines[index].split()
            energies.append(float(erg))
            counts.append(float(eff))
            index += 1

        while 'run terminated when' not in filelines[index]:
            index += 1
        match_results = match(r'\s*run terminated when \s* ([0-9]*) \s* particle histories were done.', filelines[index])
        if match_results:
            NPS = int(match_results.group(1))
        else:
            NPS = 1

        counts = np.array(counts) * NPS
        energies = np.array(energies) * 1000

        return Spectrum(counts, ergcal=Calibration(method='linear', data=np.stack([np.arange(len(counts)), energies], axis=1)))

    @staticmethod
    def from_GammaVision(filename: str | Path) -> 'Spectrum':
        """
        Import spectrum from Ortec .chn or .spe ASCII plain text file.
        
        :param filename: Path to the file where the spectrum is saved.
        
        :return: Spectrum object.
        """
        if Path(filename).suffix.lower() not in ['.spe', '.chn']:
            raise ValueError(f"{filename} is not a valid GammaVision file.")
        with open(filename, 'r') as fopen:
            filelines = fopen.readlines()
        data_index = next((i for i, line in enumerate(filelines) if line.startswith('$DATA')), None)
        if data_index is None:
            raise ValueError(f"{filename} is not a valid GammaVision file.")
        index = data_index + 2
        counts = []
        while filelines[index].strip().isdigit():
            counts.append(int(filelines[index].strip()))
            index += 1
            if index >= len(filelines):
                break
        self = 'Spectrum'(counts)

        region_index = next((i for i, line in enumerate(filelines) if line.startswith('$ROI')), None)
        if region_index is not None:
            index = region_index + 2
            self.regions = []
            while filelines[index].strip().isdigit():
                left, right = filelines[index].split()
                self.regions.append(Region(int(left), int(right)))
            index += 1

        ergcal_index = next((i for i, line in enumerate(filelines) if line.startswith('$ENER_FIT')), None)
        if ergcal_index is not None:
            incerp, slope = filelines[ergcal_index+1].split()
            self.ergcal = Calibration(method='linear', params=[float(incerp), float(slope)])
        return self

    